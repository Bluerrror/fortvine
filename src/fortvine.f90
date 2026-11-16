! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Shahab Aldin Shojaeezadeh. Part of fortvine (see LICENSE). If you use it, please cite:
! Shojaeezadeh, S. A. (2026). fortvine: a complete compendium to work on vine copula modeling.
! fortvine.f90 -- fortvine: vine copulas in Fortran (pair-copula eval, MLE fitting,
! vine density on prediction grids, Kendall's tau).
!
! Kernels. The typed Fortran API is in fortvine_api.f90; every routine here is
! also bind(C) so the library can be called from C/C++/Julia/R as well.
! Array convention of these kernels: observations are the LAST dimension,
! e.g. u(2, n) or u(d, n) (one observation's values contiguous).
!
! Family codes: 0 indep, 1 gaussian, 2 student, 3 clayton, 4 gumbel,
!               5 frank, 6 joe, 7 bb1, 8 bb6, 9 bb7, 10 bb8, 11 tawn,
!               12 plackett, 13 fgm, 14 amh, 15 galambos, 16 husler_reiss,
!               17 bb2, 18 bb3, 19 bb5, 20 bb9, 21 bb10, 22 bb4 (src/families_ext.inc),
!               23 tev (Student extreme value), 24 tll, 25 bernstein (nonparametric:
!               their estimate lives in a per-edge aux vector, see np_fit).
! Rotations follow vinecopulib / torchvine (counter-clockwise data rotation).
module fortvine
  use iso_c_binding
  use fortvine_ad
  use fortvine_ad1
  use fortvine_ad2
  use iso_fortran_env, only: int64
  !$ use omp_lib
  implicit none
  private
  integer, parameter :: dp = c_double
  real(dp), parameter :: EPSU = 1.0e-10_dp
  real(dp), parameter :: TINYV = tiny(1.0_dp)
  real(dp), parameter :: PI = 3.14159265358979323846_dp
  real(dp), parameter :: SQRT2 = 1.41421356237309504880_dp
  ! internal parameter vector: pp(1:3) parameters, pp(4:NPP) per-parameter constants (prep_par)
  integer, parameter :: NPP = 8
  real(dp), parameter :: LOGTINY = -708.3964185322641_dp   ! log(tiny(1d0))
  integer, parameter :: NSUB = 600        ! subsample size for the first fitting stage
  integer, parameter :: MS = 3            ! multi-start count for multi-parameter fits
  ! nonparametric pair copulas (24 tll, 25 bernstein): aux(1) family, aux(2) effective
  ! number of parameters, aux(3) grid size / order m, then for tll aux(4:33) the grid
  ! points and aux(34:933) the 30 x 30 density values; for bernstein aux(4:3+m*m)
  ! the cell weights
  integer, parameter, public :: NAUX = 3 + 30 + 900
  integer, parameter :: TLLM = 30
  integer, parameter :: NMEV = 60       ! Nelder-Mead budget of the subsample stage (value-only fits)
  ! likelihood passes sum NBK fixed blocks of observations in a fixed order: the result does not
  ! depend on how many threads a pass gets (bit-identical for any thread count)
  integer, parameter :: NBK = 128

  type :: ctx_t
    integer :: fam = 0, mode = 0, n = 0, nin = 1
    logical :: usew = .false.
    logical :: qn = .false.      ! newton_fit: BFGS-update the BHHH curvature (ridge-shaped likelihoods)
    logical :: lognu = .false.   ! t-EV: parameter 2 is log(nu) (the ridge is near-linear in log nu)
    real(dp) :: tol = 1.0e-8_dp  ! newton_fit: relative stopping tolerance
    ! newton_fit racing: give up once the bound on the attainable log-likelihood falls below
    ! tgt; cut returns that bound (-huge: not cut)
    real(dp) :: tgt = -huge(1.0_dp), cut = -huge(1.0_dp)
    ! racing pool of the selection (threads in use, fits running); nin is then the selection's
    ! thread budget and each pass takes a fair share of the idle threads (grab / release)
    integer, pointer :: pool(:) => null()
    real(dp) :: fixed = 0.0_dp, lbr = -1.0_dp, ubr = 1.0_dp, lastrho = 0.0_dp
    real(dp), pointer :: a(:) => null(), b(:) => null(), w(:) => null()
    real(dp), pointer :: x(:) => null(), y(:) => null()
    ! discrete margins: unrotated data a0/b0 (= F(x)) and am/bm (= F(x-)), flags da/db
    logical :: disc = .false., da = .false., db = .false.
    integer :: rot = 0
    real(dp), pointer :: a0(:) => null(), b0(:) => null(), am(:) => null(), bm(:) => null()
    ! disc_prep: row kinds and the parameter-free transforms of the two points of each row
    integer, allocatable :: dk(:)
    real(dp), allocatable :: xp(:, :), xq(:, :)
    ! t-EV (tev_prep): per observation t = log b / log(ab), log(ab), w = logit t; the 2n
    ! abscissae s = +-w sorted once (z = k (e^(s/nu) - rho) is increasing in s for every
    ! parameter value), sgn(j) = +i for s = w_i, -i for s = -w_i
    real(dp), allocatable :: tq(:), lq(:), wq(:), sq(:)
    integer, allocatable :: pq(:)
  end type ctx_t

  ! ---- pair-fit memo cache (fv_fitcache_*): exact reuse of (data, family, rotation,
  ! method, start, box) -> (parameters, loglik, nonparametric estimate) across repeated
  ! selections, e.g. tuning loops whose first-tree pairs are identical across trials
  type :: fc_entry
    logical :: used = .false.
    integer(int64) :: k1 = 0, k2 = 0
    integer :: n = 0
    real(dp) :: par(3) = 0.0_dp, ll = 0.0_dp, bnd = 0.0_dp    ! bnd > ll: a racing screen
    real(dp), allocatable :: ax(:)
  end type fc_entry
  type(fc_entry), allocatable, save :: fc_tab(:)
  integer, save :: fc_cap = 0
  integer(int64), save :: fc_hits = 0, fc_miss = 0
  ! host abort (another thread sets it): fits in flight stop at their next Newton step and
  ! new ones return at once; the results are meaningless and are not cached
  logical, public, volatile :: fv_abort = .false.

  public :: vine_rosenblatt_disc, bicop_cdf_x, ibeta_reg, pt, mle1, cand_spec
  public :: vine_pdf_x, vine_pdf_disc, vine_pdf_disc_x, vine_rosenblatt_x, vine_inverse_rosenblatt_x, vine_select_x, &
            bicop_eval_x, bicop_fit_x, bicop_select_x, criterion_w, par_bounds, vine_edge_uv, qnorm
  public :: fv_fitcache_enable, fv_fitcache_clear, fv_fitcache_disable, fv_fitcache_stats, fv_fitcache_capacity
  public :: fam_is_np, gen_is, gen_valid, gen_dup, gen_name, gen_codes, gen_npar, gen_bounds, gen_exch, gen_lt, gen_ut, gen_flip
  public :: fv_bicop_eval, fv_bicop_pdf_hfuncs, fv_fit_jobs, fv_vine_pdf, &
            fv_kendall_tau_c, fv_tau_matrix, fv_vine_select, fv_version, &
            fv_vine_rosenblatt, fv_vine_inverse_rosenblatt, fv_pseudo_obs, fv_bicop_fit, fv_bicop_select

contains

  ! ------------------------------------------------------------------ utils
  pure real(dp) function log1p_(x)
    real(dp), intent(in) :: x
    real(dp) :: u
    u = 1.0_dp + x
    if (u == 1.0_dp) then
      log1p_ = x
    else
      log1p_ = log(u) * x / (u - 1.0_dp)
    end if
  end function log1p_

  pure real(dp) function expm1_(x)
    real(dp), intent(in) :: x
    real(dp) :: u
    if (abs(x) > 0.5_dp) then
      expm1_ = exp(x) - 1.0_dp
      return
    end if
    u = exp(x)
    if (u == 1.0_dp) then
      expm1_ = x
    else
      expm1_ = (u - 1.0_dp) * x / log(u)
    end if
  end function expm1_

  pure real(dp) function clampu(u)
    real(dp), intent(in) :: u
    clampu = min(max(u, EPSU), 1.0_dp - EPSU)
  end function clampu

  pure logical function finite_(x)
    real(dp), intent(in) :: x
    finite_ = (x == x) .and. (abs(x) <= huge(x))
  end function finite_

  ! ------------------------------------------------------ normal distribution
  pure real(dp) function pnorm(x)
    real(dp), intent(in) :: x
    pnorm = 0.5_dp * erfc(-x / SQRT2)
  end function pnorm

  pure real(dp) function dnorm(x)
    real(dp), intent(in) :: x
    dnorm = exp(-0.5_dp * x * x) * 0.3989422804014327_dp
  end function dnorm

  ! Wichura (1988) AS241 PPND16, ~1e-16 relative accuracy.
  pure real(dp) function qnorm(p)
    real(dp), intent(in) :: p
    real(dp) :: q, r, v
    q = p - 0.5_dp
    if (abs(q) <= 0.425_dp) then
      r = 0.180625_dp - q * q
      v = q * (((((((r * 2509.0809287301226727_dp + 33430.575583588128105_dp) * r &
          + 67265.770927008700853_dp) * r + 45921.953931549871457_dp) * r &
          + 13731.693765509461125_dp) * r + 1971.5909503065514427_dp) * r &
          + 133.14166789178437745_dp) * r + 3.387132872796366608_dp) &
          / (((((((r * 5226.495278852545925_dp + 28729.085735721942674_dp) * r &
          + 39307.89580009271061_dp) * r + 21213.794301586595867_dp) * r &
          + 5394.1960214247511077_dp) * r + 687.1870074920579083_dp) * r &
          + 42.313330701600911252_dp) * r + 1.0_dp)
      qnorm = v
      return
    end if
    if (q < 0.0_dp) then
      r = p
    else
      r = 1.0_dp - p
    end if
    r = sqrt(-log(r))
    if (r <= 5.0_dp) then
      r = r - 1.6_dp
      v = (((((((r * 7.7454501427834140764e-4_dp + 0.0227238449892691845833_dp) * r &
          + 0.24178072517745061177_dp) * r + 1.27045825245236838258_dp) * r &
          + 3.64784832476320460504_dp) * r + 5.7694972214606914055_dp) * r &
          + 4.6303378461565452959_dp) * r + 1.42343711074968357734_dp) &
          / (((((((r * 1.05075007164441684324e-9_dp + 5.475938084995344946e-4_dp) * r &
          + 0.0151986665636164571966_dp) * r + 0.14810397642748007459_dp) * r &
          + 0.68976733498510000455_dp) * r + 1.6763848301838038494_dp) * r &
          + 2.05319162663775882187_dp) * r + 1.0_dp)
    else
      r = r - 5.0_dp
      v = (((((((r * 2.01033439929228813265e-7_dp + 2.71155556874348757815e-5_dp) * r &
          + 0.0012426609473880784386_dp) * r + 0.026532189526576123093_dp) * r &
          + 0.29656057182850489123_dp) * r + 1.7848265399172913358_dp) * r &
          + 5.4637849111641143699_dp) * r + 6.6579046435011037772_dp) &
          / (((((((r * 2.04426310338993978564e-15_dp + 1.4215117583164458887e-7_dp) * r &
          + 1.8463183175100546818e-5_dp) * r + 7.868691311456132591e-4_dp) * r &
          + 0.0148753612908506148525_dp) * r + 0.13692988092273580531_dp) * r &
          + 0.59983220655588793769_dp) * r + 1.0_dp)
    end if
    if (q < 0.0_dp) v = -v
    qnorm = v
  end function qnorm

  ! ---------------------------------------------------- Student-t distribution
  ! Regularised incomplete beta from its continued fraction (Abramowitz & Stegun 1964, eq. 26.5.8;
  ! DLMF 8.17.22):
  !   I_x(a, b) = x^a (1 - x)^b / (a B(a, b)) * K,   K = 1 / (1 + d_1 / (1 + d_2 / (1 + ...))),
  !   d_(2m+1) = -(a + m)(a + b + m) x / ((a + 2m)(a + 2m + 1)),   d_(2m) = m (b - m) x / ((a + 2m - 1)(a + 2m)).
  ! K is summed with the forward three-term (Wallis) recurrence of its convergents A_n / B_n
  ! (the fundamental recurrence formulas, e.g. Cuyt, Petersen, Verdonk, Waadeland & Jones 2008, Handbook
  ! of Continued Fractions for Special Functions, Springer, ch. 1): A_n = A_(n-1) + d_(n-1) A_(n-2),
  ! likewise B_n, every partial denominator being 1.  A and B are rescaled together when |B| leaves [1e-150, 1e150] (the convergent
  ! A/B is unchanged by a common factor), so nothing overflows.  Used where it converges fast,
  ! x < (a + 1) / (a + b + 2); the other side goes through I_x(a, b) = 1 - I_(1-x)(b, a).
  pure real(dp) function ibeta_cfrac(a, b, x) result(k)
    real(dp), intent(in) :: a, b, x
    real(dp) :: an_prev, an_cur, bn_prev, bn_cur, w, dn, kprev, sc
    integer :: n, m
    an_prev = 0.0_dp; bn_prev = 1.0_dp          ! A_0, B_0
    an_cur = 1.0_dp; bn_cur = 1.0_dp            ! A_1, B_1  (K_1 = 1)
    k = 1.0_dp
    do n = 1, 2000
      m = n / 2
      if (mod(n, 2) == 1) then
        dn = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1.0_dp))
      else
        dn = m * (b - m) * x / ((a + 2 * m - 1.0_dp) * (a + 2 * m))
      end if
      w = an_cur + dn * an_prev; an_prev = an_cur; an_cur = w
      w = bn_cur + dn * bn_prev; bn_prev = bn_cur; bn_cur = w
      if (abs(bn_cur) > 1.0e150_dp .or. (abs(bn_cur) < 1.0e-150_dp .and. bn_cur /= 0.0_dp)) then
        sc = 1.0_dp / abs(bn_cur)
        an_prev = an_prev * sc; an_cur = an_cur * sc; bn_prev = bn_prev * sc; bn_cur = bn_cur * sc
      end if
      if (bn_cur /= 0.0_dp) then
        kprev = k; k = an_cur / bn_cur
        if (abs(k - kprev) <= 1.0e-15_dp * abs(k)) exit
      end if
    end do
  end function ibeta_cfrac

  ! The same continued fraction and its directional derivative along (da, db, dx): the recurrence is
  ! differentiated term by term (forward mode), dA_n = dA_(n-1) + dd_(n-1) A_(n-2) + d_(n-1) dA_(n-2),
  ! and dK = (dA - K dB) / B.  Converged when value AND derivative are (an integer b ends the
  ! fraction's value exactly while the derivative in b still has a tail).
  pure subroutine ibeta_cfrac_d(a, b, x, da, db, dx, k, dk)
    real(dp), intent(in) :: a, b, x, da, db, dx
    real(dp), intent(out) :: k, dk
    real(dp) :: ap, ac, bp, bc, gap, gac, gbp, gbc, w, wg, num, dnum, den, dden, dn, ddn, kprev, dkprev, sc
    integer :: n, m
    ap = 0.0_dp; bp = 1.0_dp; ac = 1.0_dp; bc = 1.0_dp
    gap = 0.0_dp; gbp = 0.0_dp; gac = 0.0_dp; gbc = 0.0_dp
    k = 1.0_dp; dk = 0.0_dp
    do n = 1, 2000
      m = n / 2
      if (mod(n, 2) == 1) then
        num = -(a + m) * (a + b + m) * x
        dnum = -(da * (a + b + m) * x + (a + m) * (da + db) * x + (a + m) * (a + b + m) * dx)
        den = (a + 2 * m) * (a + 2 * m + 1.0_dp); dden = da * (2.0_dp * a + 4 * m + 1.0_dp)
      else
        num = m * (b - m) * x; dnum = m * (db * x + (b - m) * dx)
        den = (a + 2 * m - 1.0_dp) * (a + 2 * m); dden = da * (2.0_dp * a + 4 * m - 1.0_dp)
      end if
      dn = num / den; ddn = (dnum - dn * dden) / den
      w = ac + dn * ap; wg = gac + ddn * ap + dn * gap; ap = ac; gap = gac; ac = w; gac = wg
      w = bc + dn * bp; wg = gbc + ddn * bp + dn * gbp; bp = bc; gbp = gbc; bc = w; gbc = wg
      if (abs(bc) > 1.0e150_dp .or. (abs(bc) < 1.0e-150_dp .and. bc /= 0.0_dp)) then
        sc = 1.0_dp / abs(bc)
        ap = ap * sc; ac = ac * sc; bp = bp * sc; bc = bc * sc
        gap = gap * sc; gac = gac * sc; gbp = gbp * sc; gbc = gbc * sc
      end if
      if (bc /= 0.0_dp) then
        kprev = k; dkprev = dk
        k = ac / bc; dk = (gac - k * gbc) / bc
        if (abs(k - kprev) <= 1.0e-15_dp * abs(k) .and. abs(dk - dkprev) <= 1.0e-15_dp * (abs(dk) + abs(k))) exit
      end if
    end do
  end subroutine ibeta_cfrac_d

  ! d/dnu of the t_nu CDF at fixed x (exact): T = 1/2 -+ 1/2 I_y(a, b) with (a, b, y)
  ! depending on nu; differentiates log B (digamma), the prefactor y^a (1-y)^b and the
  ! continued fraction (ibeta_cfrac_d).  cs = [log B(nu/2, 1/2), psi(nu/2), psi((nu+1)/2)].
  pure real(dp) function dpt_dnu(x0, nu, cs) result(dT)
    real(dp), intent(in) :: x0, nu, cs(3)
    real(dp) :: x, t2, y, dy, a, b, da, db, lbt, dlbt, cf, dcf, sg, e
    x = -abs(x0); t2 = x * x
    if (nu > t2) then
      a = 0.5_dp; b = 0.5_dp * nu; da = 0.0_dp; db = 0.5_dp
      y = t2 / (nu + t2); dy = -t2 / ((nu + t2) * (nu + t2)); sg = -0.5_dp
    else
      a = 0.5_dp * nu; b = 0.5_dp; da = 0.5_dp; db = 0.0_dp
      y = nu / (nu + t2); dy = t2 / ((nu + t2) * (nu + t2)); sg = 0.5_dp
    end if
    if (y <= 0.0_dp .or. y >= 1.0_dp) then
      dT = 0.0_dp; return
    end if
    ! d log B(a, b)/dnu = 1/2 (psi(nu/2) - psi((nu+1)/2)) in both branches
    lbt = -cs(1) + a * log(y) + b * log1p_(-y)
    dlbt = -0.5_dp * (cs(2) - cs(3)) + da * log(y) + db * log1p_(-y) + (a / y - b / (1.0_dp - y)) * dy
    if (y < (a + 1.0_dp) / (a + b + 2.0_dp)) then
      call ibeta_cfrac_d(a, b, y, da, db, dy, cf, dcf)
      e = exp(lbt) * cf / a
      dT = sg * e * (dlbt + dcf / cf - da / a)
    else
      call ibeta_cfrac_d(b, a, 1.0_dp - y, db, da, -dy, cf, dcf)
      e = exp(lbt) * cf / b
      dT = -sg * e * (dlbt + dcf / cf - db / b)
    end if
    if (x0 > 0.0_dp) dT = -dT
  end function dpt_dnu

  ! Regularised incomplete beta I_x(a, b).
  pure real(dp) function ibeta_reg(a, b, x)
    real(dp), intent(in) :: a, b, x
    ibeta_reg = ibeta_reg_l(a, b, x, log_gamma(a) + log_gamma(b) - log_gamma(a + b))
  end function ibeta_reg

  ! I_x(a, b) given lbeta = log B(a, b) (precomputed once per parameter value)
  pure real(dp) function ibeta_reg_l(a, b, x, lbeta)
    real(dp), intent(in) :: a, b, x, lbeta
    real(dp) :: lpre
    if (x <= 0.0_dp) then
      ibeta_reg_l = 0.0_dp; return
    end if
    if (x >= 1.0_dp) then
      ibeta_reg_l = 1.0_dp; return
    end if
    lpre = a * log(x) + b * log1p_(-x) - lbeta          ! log[x^a (1-x)^b / B(a, b)]
    if (x * (a + b + 2.0_dp) < a + 1.0_dp) then
      ibeta_reg_l = exp(lpre) * ibeta_cfrac(a, b, x) / a
    else
      ibeta_reg_l = 1.0_dp - exp(lpre) * ibeta_cfrac(b, a, 1.0_dp - x) / b
    end if
  end function ibeta_reg_l

  ! Lower-tail CDF for x <= 0 (accurate in the far tail).
  pure real(dp) function pt_lower(x, nu)
    real(dp), intent(in) :: x, nu
    pt_lower = pt_lower_l(x, nu, lbeta_t(nu))
  end function pt_lower

  ! log B(nu/2, 1/2)
  pure real(dp) function lbeta_t(nu)
    real(dp), intent(in) :: nu
    lbeta_t = log_gamma(0.5_dp * nu) + log_gamma(0.5_dp) - log_gamma(0.5_dp * (nu + 1.0_dp))
  end function lbeta_t

  ! log normalising constant of the t density
  pure real(dp) function kdens_t(nu)
    real(dp), intent(in) :: nu
    kdens_t = log_gamma(0.5_dp * (nu + 1.0_dp)) - log_gamma(0.5_dp * nu) - 0.5_dp * log(nu * PI)
  end function kdens_t

  ! lower-tail t CDF for x <= 0 given lb = log B(nu/2, 1/2)
  pure real(dp) function pt_lower_l(x, nu, lb)
    real(dp), intent(in) :: x, nu, lb
    real(dp) :: t2
    t2 = x * x
    if (nu > t2) then
      ! 1/2 - 1/2 I_{t2/(nu+t2)}(1/2, nu/2): small t, avoid 1 - (nu/(nu+t2))
      pt_lower_l = 0.5_dp - 0.5_dp * ibeta_reg_l(0.5_dp, 0.5_dp * nu, t2 / (nu + t2), lb)
    else
      pt_lower_l = 0.5_dp * ibeta_reg_l(0.5_dp * nu, 0.5_dp, nu / (nu + t2), lb)
    end if
  end function pt_lower_l

  pure real(dp) function pt_l(x, nu, lb)
    real(dp), intent(in) :: x, nu, lb
    if (x <= 0.0_dp) then
      pt_l = pt_lower_l(x, nu, lb)
    else
      pt_l = 1.0_dp - pt_lower_l(-x, nu, lb)
    end if
  end function pt_l

  pure real(dp) function pt(x, nu)
    real(dp), intent(in) :: x, nu
    if (x <= 0.0_dp) then
      pt = pt_lower(x, nu)
    else
      pt = 1.0_dp - pt_lower(-x, nu)
    end if
  end function pt

  pure real(dp) function ldt(x, nu)
    real(dp), intent(in) :: x, nu
    ldt = log_gamma(0.5_dp * (nu + 1.0_dp)) - log_gamma(0.5_dp * nu) &
          - 0.5_dp * log(nu * PI) - 0.5_dp * (nu + 1.0_dp) * log1p_(x * x / nu)
  end function ldt

  ! Student-t quantile.  Start value: Hill, G. W. (1970) Algorithm 396: Student's t-quantiles,
  ! Commun. ACM 13(10), 619-620, doi:10.1145/355598.355600 (exact for nu = 1 and nu = 2; otherwise
  ! Hill's expansion about the normal deviate in the body and his tail formula), then a safeguarded
  ! Halley iteration on the upper-tail probability with the exact t CDF and density (qt_polish).
  pure real(dp) function qt(p, nu)
    real(dp), intent(in) :: p, nu
    qt = qt_l(p, nu, kdens_t(nu), lbeta_t(nu))
  end function qt

  ! kd = log density constant (kdens_t), lb = log B(nu/2, 1/2) (lbeta_t)
  pure real(dp) function qt_l(p, nu, kd, lb)
    real(dp), intent(in) :: p, nu, kd, lb
    real(dp) :: two_tail, t
    two_tail = 2.0_dp * min(p, 1.0_dp - p)            ! Hill's argument: P(|T| > t)
    if (two_tail >= 1.0_dp) then
      qt_l = 0.0_dp; return
    end if
    if (.not. (two_tail > 0.0_dp)) then                ! p = 0 or 1: the infinite quantile
      t = huge(1.0_dp)
    else
      t = qt_polish(hill396(two_tail, nu), 0.5_dp * two_tail, nu, kd, lb)
    end if
    if (p < 0.5_dp) t = -t
    qt_l = t
  end function qt_l

  ! Hill (1970), Algorithm 396, for the two-tail probability pp in (0, 1) and nu degrees of freedom
  ! (real nu is fine).  Hill's normal deviate normdev(pp/2) is the lower-tail one (negative); a
  ! numerical check against exact quantiles confirms that sign (start error ~1e-5 at nu = 3,
  ! ~1e-8 at nu = 50).  expm1_ replaces Hill's two-branch exp(y) - 1 / y + y^2/2 (same function).
  pure real(dp) function hill396(pp, nu) result(t)
    real(dp), intent(in) :: pp, nu
    real(dp) :: a, b, c, d, y, z
    if (abs(nu - 1.0_dp) < 1.0e-12_dp) then            ! Cauchy: t = cot(pi pp / 2)
      t = cos(0.5_dp * PI * pp) / sin(0.5_dp * PI * pp); return
    end if
    if (abs(nu - 2.0_dp) < 1.0e-12_dp) then            ! closed form for 2 df
      t = sqrt(2.0_dp / (pp * (2.0_dp - pp)) - 2.0_dp); return
    end if
    if (abs(nu - 0.5_dp) < 1.0e-3_dp) then             ! Hill's a = 1/(nu - 1/2) is singular here
      t = 1.0_dp; return                               ! (only a start: qt_polish brackets and converges)
    end if
    a = 1.0_dp / (nu - 0.5_dp)
    b = 48.0_dp / (a * a)
    c = ((20700.0_dp * a / b - 98.0_dp) * a - 16.0_dp) * a + 96.36_dp
    d = ((94.5_dp / (b + c) - 3.0_dp) / b + 1.0_dp) * sqrt(0.5_dp * PI * a) * nu
    y = (d * pp)**(2.0_dp / nu)
    if (y > 0.05_dp + a) then                          ! body: expansion about the normal deviate
      z = qnorm(0.5_dp * pp)
      y = z * z
      if (nu < 5.0_dp) c = c + 0.3_dp * (nu - 4.5_dp) * (z + 0.6_dp)
      c = (((0.05_dp * d * z - 5.0_dp) * z - 7.0_dp) * z - 2.0_dp) * z + b + c
      y = (((((0.4_dp * y + 6.3_dp) * y + 36.0_dp) * y + 94.5_dp) / c - y - 3.0_dp) / b + 1.0_dp) * z
      y = expm1_(a * y * y)
    else                                               ! far tail
      y = ((1.0_dp / (((nu + 6.0_dp) / (nu * y) - 0.089_dp * d - 0.822_dp) * (nu + 2.0_dp) * 3.0_dp) &
           + 0.5_dp / (nu + 4.0_dp)) * y - 1.0_dp) * (nu + 1.0_dp) / (nu + 2.0_dp) + 1.0_dp / y
    end if
    t = sqrt(nu * y)
  end function hill396

  ! Solve P(T > t) = tail (0 < tail < 1/2) for t > 0 from the start t0: Halley's method on
  ! g(t) = P(T > t) - tail, g' = -f(t), g'' = f(t) (nu + 1) t / (nu + t^2) (log-derivative of the t
  ! density), i.e. t <- t + h / (1 - h (nu + 1) t / (2 (nu + t^2))), h = g / f.  A bracket [lo, hi]
  ! is kept from the sign of g; a step that leaves it (or a non-finite start) bisects / doubles
  ! instead.  Halley converges cubically, so once an accepted Halley step is below 1e-6 t the error
  ! left is ~1e-18 t and the iteration stops; a step below rounding stops it as well.
  pure real(dp) function qt_polish(t0, tail, nu, kd, lb) result(t)
    real(dp), intent(in) :: t0, tail, nu, kd, lb
    real(dp) :: lo, hi, g, f, h, r, step
    integer :: it
    lo = 0.0_dp; hi = huge(1.0_dp)
    t = t0
    if (.not. (t > 0.0_dp .and. t < hi)) t = 1.0_dp
    do it = 1, 300
      g = pt_lower_l(-t, nu, lb) - tail
      if (g == 0.0_dp) exit
      if (g > 0.0_dp) then
        lo = t
      else
        hi = t
      end if
      f = exp(kd - 0.5_dp * (nu + 1.0_dp) * log1p_(t * t / nu))
      if (f > 0.0_dp) then
        h = g / f
        r = 1.0_dp - 0.5_dp * h * (nu + 1.0_dp) * t / (nu + t * t)
        if (r > 0.5_dp) then
          step = h / r                                 ! Halley
        else
          step = h                                     ! Newton when the Halley factor is unreliable
        end if
        if (abs(step) <= 4.0_dp * epsilon(1.0_dp) * t) exit       ! converged to rounding
        if (t + step > lo .and. t + step < hi) then
          t = t + step
          if (r > 0.5_dp .and. abs(step) <= 1.0e-6_dp * t) exit
          cycle
        end if
      end if
      if (hi < huge(1.0_dp)) then                      ! safeguard: bisect the bracket / expand it
        t = 0.5_dp * (lo + hi)
        if (hi - lo <= 4.0_dp * epsilon(1.0_dp) * t) exit
      else
        t = 2.0_dp * max(t, 1.0_dp)
      end if
    end do
  end function qt_polish

  ! Internal parameter vector with per-parameter constants (computed once per
  ! edge / likelihood evaluation, never per observation).
  pure function prep_par(fam, par) result(pp)
    integer, intent(in) :: fam
    real(dp), intent(in) :: par(3)
    real(dp) :: pp(NPP), nu
    pp = 0.0_dp
    pp(1:3) = par
    if (fam == 23) then
      nu = min(max(par(2), 1.0_dp), 60.0_dp)
      pp(7) = kdens_t(nu + 1.0_dp)
      pp(8) = lbeta_t(nu + 1.0_dp)
    end if
    if (fam == 2) then
      nu = min(max(par(2), 2.0_dp), 100.0_dp)
      pp(4) = log_gamma(0.5_dp * (nu + 2.0_dp)) + log_gamma(0.5_dp * nu) - 2.0_dp * log_gamma(0.5_dp * (nu + 1.0_dp))
      pp(5) = kdens_t(nu)
      pp(6) = lbeta_t(nu)
      pp(7) = kdens_t(nu + 1.0_dp)
      pp(8) = lbeta_t(nu + 1.0_dp)
    end if
  end function prep_par

  ! log(1 - exp(z)) for z < 0 (Maechler 2012)
  pure real(dp) function log1mexp(z)
    real(dp), intent(in) :: z
    if (z > -0.6931471805599453_dp) then
      log1mexp = log(-expm1_(z))
    else
      log1mexp = log1p_(-exp(z))
    end if
  end function log1mexp

  ! ----------------------------------------------- Archimedean generators (BB)
  pure subroutine arch(fam, th, de, x, phi, d1, d2)
    integer, intent(in) :: fam
    real(dp), intent(in) :: th, de, x
    real(dp), intent(out) :: phi, d1, d2
    real(dp) :: t1, tmp, l, ox, nrm
    select case (fam)
    case (7)  ! bb1
      t1 = x**(-th) - 1.0_dp
      phi = t1**de
      d1 = -de * th * x**(-(1.0_dp + th)) * t1**(de - 1.0_dp)
      d2 = de * th * t1**de * ((1.0_dp + de * th) - (1.0_dp + th) * x**th) &
           / max((x**th - 1.0_dp)**2 * x * x, TINYV)
    case (8)  ! bb6
      ox = max(1.0_dp - x, TINYV)
      tmp = ox**th
      l = log1p_(-tmp)
      phi = max(-l, TINYV)**de
      d1 = de * th * max(-l, TINYV)**(de - 1.0_dp) * ox**(th - 1.0_dp) / min(tmp - 1.0_dp, -TINYV)
      d2 = max(-l, TINYV)**(de - 2.0_dp) * de * th * ox**(th - 2.0_dp) &
           * ((de - 1.0_dp) * th * tmp - (tmp + th - 1.0_dp) * l) / (tmp - 1.0_dp)**2
    case (9)  ! bb7
      ox = max(1.0_dp - x, TINYV)
      tmp = ox**th
      phi = max(1.0_dp - tmp, TINYV)**(-de) - 1.0_dp
      d1 = -de * th * max(1.0_dp - tmp, TINYV)**(-1.0_dp - de) * ox**(th - 1.0_dp)
      d2 = de * th * max(1.0_dp - tmp, TINYV)**(-2.0_dp - de) * ox**(th - 2.0_dp) &
           * (th - 1.0_dp + (1.0_dp + de * th) * tmp)
    case (10) ! bb8
      nrm = max(1.0_dp - max(1.0_dp - de, TINYV)**th, TINYV)
      ox = max(1.0_dp - de * x, TINYV)
      tmp = ox**th
      phi = -log(max(1.0_dp - tmp, TINYV) / nrm)
      d1 = -de * th * ox**(th - 1.0_dp) / max(1.0_dp - tmp, TINYV)
      d2 = de * de * th * ox**(th - 2.0_dp) * (th - 1.0_dp + tmp) / (tmp - 1.0_dp)**2
    case default
      phi = 0.0_dp; d1 = -1.0_dp; d2 = 0.0_dp
    end select
  end subroutine arch

  pure real(dp) function arch_psi(fam, th, de, t)
    integer, intent(in) :: fam
    real(dp), intent(in) :: th, de, t
    real(dp) :: tt
    tt = max(t, TINYV)
    select case (fam)
    case (7)
      arch_psi = (tt**(1.0_dp / de) + 1.0_dp)**(-1.0_dp / th)
    case (8)
      arch_psi = 1.0_dp - max(-expm1_(-tt**(1.0_dp / de)), TINYV)**(1.0_dp / th)
    case (9)
      arch_psi = 1.0_dp - (1.0_dp - (1.0_dp + tt)**(-1.0_dp / de))**(1.0_dp / th)
    case (10)
      arch_psi = (1.0_dp - max(1.0_dp + exp(-tt) * (max(1.0_dp - de, TINYV)**th - 1.0_dp), TINYV) &
                 **(1.0_dp / th)) / de
    case default
      arch_psi = exp(-tt)
    end select
  end function arch_psi

  ! Tawn Pickands function A(t) and derivatives (vinecopulib tawn.ipp; A'' in the
  ! log form that stays finite for small psi).
  pure subroutine tawn_pickands(pp, t, aa, a1, a2)
    real(dp), intent(in) :: pp(3), t
    real(dp), intent(out) :: aa, a1, a2
    real(dp) :: p1, p2, th, tmp, lx, ly, hi
    p1 = pp(1); p2 = pp(2); th = pp(3)
    tmp = (p2 * t)**th + (p1 * (1.0_dp - t))**th
    aa = (1.0_dp - p1) * (1.0_dp - t) + (1.0_dp - p2) * t + tmp**(1.0_dp / th)
    a1 = p1 - p2 + tmp**(1.0_dp / th - 1.0_dp) * (p2 * (p2 * t)**(th - 1.0_dp) &
         - p1 * (p1 * (1.0_dp - t))**(th - 1.0_dp))
    if (p1 * (1.0_dp - t) <= 0.0_dp .or. p2 * t <= 0.0_dp) then
      a2 = 0.0_dp; return
    end if
    lx = th * log(p2 * t); ly = th * log(p1 * (1.0_dp - t))
    hi = max(lx, ly)
    a2 = exp(log(th - 1.0_dp) + (1.0_dp / th - 2.0_dp) * (hi + log1p_(exp(min(lx, ly) - hi))) &
             + lx + ly - 2.0_dp * log(t) - 2.0_dp * log1p_(-t))
  end subroutine tawn_pickands

  ! ------------------------------------------------- unrotated pair-copula kernels
  ! (a, b) already rotated and clamped. Returns log-density and both
  ! h-functions h1 = P(V<=b | U=a), h2 = P(U<=a | V=b).
  pure subroutine eval0(fam, pp, a, b, nh, lp, h1, h2)
    integer, intent(in) :: fam
    real(dp), intent(in) :: pp(NPP), a, b
    logical, intent(in) :: nh
    real(dp), intent(out) :: lp, h1, h2
    real(dp) :: r, s, x, y, nu, th, t1, t2, ls, st, lx, ly, ua, ub, eu, ev, ed, den
    real(dp) :: pa, pb, cc, da, db, dc, ddc, dum, p1, p2, la, lb2, tt, aa, a1, a2, lx2, ly2, hi2, g1, g2, gc
    p1 = pp(1); p2 = pp(2)
    h1 = b; h2 = a
    select case (fam)
    case (0)
      lp = 0.0_dp
    case (1)
      r = min(max(p1, -0.999999_dp), 0.999999_dp)
      x = qnorm(a); y = qnorm(b)
      s = 1.0_dp - r * r
      lp = -0.5_dp * log(s) - (r * r * (x * x + y * y) - 2.0_dp * r * x * y) / (2.0_dp * s)
      if (nh) then
        h1 = pnorm((y - r * x) / sqrt(s))
        h2 = pnorm((x - r * y) / sqrt(s))
      end if
    case (2)
      r = min(max(p1, -0.999999_dp), 0.999999_dp)
      nu = min(max(p2, 2.0_dp), 100.0_dp)
      x = qt_l(a, nu, pp(5), pp(6)); y = qt_l(b, nu, pp(5), pp(6))
      s = (1.0_dp - r) * (1.0_dp + r)
      lp = pp(4) - 0.5_dp * log(s) &
           + 0.5_dp * (nu + 1.0_dp) * (log1p_(x * x / nu) + log1p_(y * y / nu)) &
           - 0.5_dp * (nu + 2.0_dp) * log1p_((x * x - 2.0_dp * r * x * y + y * y) / (nu * s))
      if (nh) then
        h1 = pt_l((y - r * x) / sqrt((nu + x * x) * s / (nu + 1.0_dp)), nu + 1.0_dp, pp(8))
        h2 = pt_l((x - r * y) / sqrt((nu + y * y) * s / (nu + 1.0_dp)), nu + 1.0_dp, pp(8))
      end if
    case (3)
      th = p1
      if (th < 1.0e-10_dp) then
        lp = 0.0_dp; return
      end if
      t1 = -th * log(a); t2 = -th * log(b)
      ! log(a^-th + b^-th - 1), computed stably
      ls = max(t1, t2) + log(exp(t1 - max(t1, t2)) + exp(t2 - max(t1, t2)) - exp(-max(t1, t2)))
      lp = log1p_(th) + (-1.0_dp - th) * (log(a) + log(b)) + (-2.0_dp - 1.0_dp / th) * ls
      if (nh) then
        h1 = exp((-th - 1.0_dp) * log(a) + (-1.0_dp / th - 1.0_dp) * ls)
        h2 = exp((-th - 1.0_dp) * log(b) + (-1.0_dp / th - 1.0_dp) * ls)
      end if
    case (4)
      th = max(p1, 1.0_dp)
      x = -log(a); y = -log(b)
      lx = log(x); ly = log(y)
      t1 = th * lx; t2 = th * ly
      ls = max(t1, t2) + log1p_(exp(min(t1, t2) - max(t1, t2)))   ! log(x^th + y^th)
      st = exp(ls / th)
      lp = -st + (th - 1.0_dp) * (lx + ly) + x + y + (1.0_dp / th - 2.0_dp) * ls &
           + log(th - 1.0_dp + st)
      if (nh) then
        h1 = exp(-st + x + (th - 1.0_dp) * lx + (1.0_dp / th - 1.0_dp) * ls)
        h2 = exp(-st + y + (th - 1.0_dp) * ly + (1.0_dp / th - 1.0_dp) * ls)
      end if
    case (5)
      th = p1
      if (abs(th) < 1.0e-10_dp) then
        lp = 0.0_dp; return
      end if
      eu = expm1_(-th * a); ev = expm1_(-th * b); ed = expm1_(-th)
      den = ed + eu * ev
      lp = log(max(-th * ed * (eu + 1.0_dp) * (ev + 1.0_dp), TINYV)) - 2.0_dp * log(abs(den))
      if (nh) then
        h1 = (eu + 1.0_dp) * ev / den
        h2 = (ev + 1.0_dp) * eu / den
      end if
    case (6)
      th = max(p1, 1.0_dp)
      ua = (1.0_dp - a)**th; ub = (1.0_dp - b)**th
      s = ua + ub - ua * ub
      lp = (1.0_dp / th - 2.0_dp) * log(s) + (th - 1.0_dp) * (log(1.0_dp - a) + log(1.0_dp - b)) &
           + log(th - 1.0_dp + s)
      if (nh) then
        h1 = exp((th - 1.0_dp) * log(1.0_dp - a) + log1p_(-ub) + (1.0_dp / th - 1.0_dp) * log(s))
        h2 = exp((th - 1.0_dp) * log(1.0_dp - b) + log1p_(-ua) + (1.0_dp / th - 1.0_dp) * log(s))
      end if
    case (11)
      ! Tawn (extreme-value): vinecopulib ExtremeValueBicop::pdf_raw / hfunc*_raw
      la = log(a); lb2 = log(b)
      tt = lb2 / (la + lb2)
      call tawn_pickands(pp, tt, aa, a1, a2)
      lp = (la + lb2) * aa - la - lb2 + log(max(aa * aa + (1.0_dp - 2.0_dp * tt) * a1 * aa &
           - (1.0_dp - tt) * tt * (a1 * a1 + a2 / (la + lb2)), TINYV))
      if (nh) then
        h1 = exp((la + lb2) * aa) * (aa - tt * a1) / a
        h2 = exp((la + lb2) * aa) * (aa + (1.0_dp - tt) * a1) / b
      end if
    case (12:22)
      call ext_eval_r(fam, pp(1:3), a, b, log(a), log(b), nh, lp, h1, h2)
    case (23)
      call tev_eval(pp, a, b, nh, lp, h1, h2)
    case (7:10)
      ! Closed forms: log c = log psi''(s) + log|phi'(a)| + log|phi'(b)|,
      ! h1 = psi'(s) phi'(a), s = phi(a) + phi(b); C itself is never formed.
      call bb_closed(fam, p1, p2, a, b, log(a), log(b), log1p_(-a), log1p_(-b), nh, lp, h1, h2)
    case default
      lp = 0.0_dp
      if (gen_is(fam)) then                         ! generic families (families_gen.inc)
        call gen_eval_rv(fam, pp(1:3), a, b, lp, g1, g2, gc)
        if (nh) then
          h1 = exp(g1); h2 = exp(g2)
        end if
      end if
    end select
    h1 = min(max(h1, 0.0_dp), 1.0_dp)
    h2 = min(max(h2, 0.0_dp), 1.0_dp)
  end subroutine eval0

  ! BB1/BB6/BB7/BB8 log-density and h-functions in closed form (log space).
  ! la = log a, la1 = log(1-a) (precomputable per observation).
  pure subroutine bb_closed(fam, th, de, a, b, la, lb, la1, lb1, nh, lp, h1, h2)
    integer, intent(in) :: fam
    real(dp), intent(in) :: th, de, a, b, la, lb, la1, lb1
    logical, intent(in) :: nh
    real(dp), intent(out) :: lp, h1, h2
    real(dp) :: ta, tb, lta, ltb, ls, w, ua, ub, ra, rb, lra, lrb, z, lz, ps, lq, lpa, lpb, eta, ga, gb
    real(dp) :: l1ga, l1gb, l1da, l1db, emz, e1, e2
    select case (fam)
    case (7)   ! BB1: phi = (t^-th - 1)^de, psi = (1 + s^(1/de))^(-1/th)
      ta = expm1_(-th * la); tb = expm1_(-th * lb)
      lta = log(ta); ltb = log(tb)
      e1 = de * lta; e2 = de * ltb
      ls = max(e1, e2) + log1p_(exp(min(e1, e2) - max(e1, e2)))
      w = exp(ls / de)
      lp = log(th * (de - 1.0_dp) + (th * de + 1.0_dp) * w) + (-1.0_dp / th - 2.0_dp) * log1p_(w) &
           + (1.0_dp / de - 2.0_dp) * ls + (-th - 1.0_dp) * (la + lb) + (de - 1.0_dp) * (lta + ltb)
      if (nh) then
        ps = (-1.0_dp / th - 1.0_dp) * log1p_(w) + (1.0_dp / de - 1.0_dp) * ls
        h1 = exp(ps + (-th - 1.0_dp) * la + (de - 1.0_dp) * lta)
        h2 = exp(ps + (-th - 1.0_dp) * lb + (de - 1.0_dp) * ltb)
      end if
    case (8)   ! BB6: phi = (-log(1-(1-t)^th))^de, psi = 1 - (1 - exp(-s^(1/de)))^(1/th)
      ua = exp(th * la1); ub = exp(th * lb1)
      ra = -log1p_(-ua); rb = -log1p_(-ub)
      lra = log(ra); lrb = log(rb)
      e1 = de * lra; e2 = de * lrb
      ls = max(e1, e2) + log1p_(exp(min(e1, e2) - max(e1, e2)))
      w = exp(ls / de)
      z = -expm1_(-w); lz = log(z)
      lpa = (de - 1.0_dp) * lra + (th - 1.0_dp) * la1 + ra
      lpb = (de - 1.0_dp) * lrb + (th - 1.0_dp) * lb1 + rb
      lp = log(th) + (1.0_dp / th - 2.0_dp) * lz - w + log(w) - 2.0_dp * ls &
           + log(z * (w + de - 1.0_dp) + (1.0_dp - 1.0_dp / th) * exp(-w) * w) + lpa + lpb
      if (nh) then
        ps = (1.0_dp / th - 1.0_dp) * lz - w + log(w) - ls
        h1 = exp(ps + lpa); h2 = exp(ps + lpb)
      end if
    case (9)   ! BB7: phi = (1-(1-t)^th)^-de - 1, psi = 1 - (1 - (1+s)^(-1/de))^(1/th)
      ua = exp(th * la1); ub = exp(th * lb1)
      ra = log1p_(-ua); rb = log1p_(-ub)
      ls = log1p_(expm1_(-de * ra) + expm1_(-de * rb))     ! log(1 + s)
      lq = -ls / de
      z = -expm1_(lq); lz = log(z)
      lpa = (-de - 1.0_dp) * ra + (th - 1.0_dp) * la1
      lpb = (-de - 1.0_dp) * rb + (th - 1.0_dp) * lb1
      lp = log(th) + (1.0_dp / th - 2.0_dp) * lz + lq - 2.0_dp * ls &
           + log((1.0_dp - 1.0_dp / th) * exp(lq) + z * (1.0_dp + de)) + lpa + lpb
      if (nh) then
        ps = (1.0_dp / th - 1.0_dp) * lz + lq - ls
        h1 = exp(ps + lpa); h2 = exp(ps + lpb)
      end if
    case default   ! BB8: phi = -log((1-(1-de t)^th)/eta), psi = (1 - (1 - eta e^-s)^(1/th))/de
      eta = -expm1_(th * log1p_(-de))
      l1da = log1p_(-de * a); l1db = log1p_(-de * b)
      ga = exp(th * l1da); gb = exp(th * l1db)
      l1ga = log1p_(-ga); l1gb = log1p_(-gb)
      emz = exp(l1ga + l1gb) / eta              ! 1 - z = eta e^-s
      z = 1.0_dp - emz; lz = log(z)
      lpa = (th - 1.0_dp) * l1da - l1ga
      lpb = (th - 1.0_dp) * l1db - l1gb
      lp = log(de * th) + (1.0_dp / th - 2.0_dp) * lz + log(emz) + log(1.0_dp - emz / th) + lpa + lpb
      if (nh) then
        ps = (1.0_dp / th - 1.0_dp) * lz + log(emz)
        h1 = exp(ps + lpa); h2 = exp(ps + lpb)
      end if
    end select
  end subroutine bb_closed

  ! Student-t extreme-value copula (Demarta & McNeil 2005), rho in (-1,1), nu > 0:
  ! A(t) = t T(z1) + (1-t) T(z2), T = t_{nu+1} CDF, x = t/(1-t), k = sqrt((nu+1)/(1-rho^2)),
  ! z1 = k (x^(1/nu) - rho), z2 = k (x^(-1/nu) - rho); A', A'' in closed form with
  ! z' = +-k/nu x^(+-1/nu) L, L = 1/(t(1-t)).  T has no closed-form nu-derivative, so
  ! this family is fitted derivative-free (fit_one).
  pure subroutine tev_eval(pp, a, b, nh, lp, h1, h2)
    real(dp), intent(in) :: pp(NPP), a, b
    logical, intent(in) :: nh
    real(dp), intent(out) :: lp, h1, h2
    real(dp) :: r, nu, k, la, lb, lab, t, x1, x2, z1, z2, L, d1, d2, dd1, dd2, f1, f2, fp1, fp2
    real(dp) :: aa, a1, a2, c1, T1, T2
    r = min(max(pp(1), -0.999_dp), 0.999_dp)
    nu = min(max(pp(2), 1.0_dp), 60.0_dp)
    k = sqrt((nu + 1.0_dp) / ((1.0_dp - r) * (1.0_dp + r)))
    la = log(a); lb = log(b); lab = la + lb
    t = lb / lab
    t = min(max(t, 1.0e-300_dp), 1.0_dp - epsilon(1.0_dp))
    L = 1.0_dp / (t * (1.0_dp - t))
    x1 = exp((log(t) - log1p_(-t)) / nu); x2 = exp((log1p_(-t) - log(t)) / nu)   ! x^(1/nu), x^(-1/nu)
    z1 = k * (x1 - r); z2 = k * (x2 - r)
    d1 = k / nu * x1 * L; d2 = -k / nu * x2 * L
    dd1 = k / nu * x1 * L * L * (1.0_dp / nu - (1.0_dp - 2.0_dp * t))
    dd2 = k / nu * x2 * L * L * (1.0_dp / nu + (1.0_dp - 2.0_dp * t))
    c1 = nu + 1.0_dp
    f1 = exp(pp(7) - 0.5_dp * (c1 + 1.0_dp) * log1p_(z1 * z1 / c1))
    f2 = exp(pp(7) - 0.5_dp * (c1 + 1.0_dp) * log1p_(z2 * z2 / c1))
    fp1 = -f1 * (c1 + 1.0_dp) * z1 / (c1 + z1 * z1)
    fp2 = -f2 * (c1 + 1.0_dp) * z2 / (c1 + z2 * z2)
    T1 = pt_l(z1, c1, pp(8)); T2 = pt_l(z2, c1, pp(8))
    aa = t * T1 + (1.0_dp - t) * T2
    a1 = T1 - T2 + t * f1 * d1 + (1.0_dp - t) * f2 * d2
    a2 = 2.0_dp * (f1 * d1 - f2 * d2) + t * (fp1 * d1 * d1 + f1 * dd1) + (1.0_dp - t) * (fp2 * d2 * d2 + f2 * dd2)
    lp = lab * aa - lab + log(max(aa * aa + (1.0_dp - 2.0_dp * t) * a1 * aa &
         - (1.0_dp - t) * t * (a1 * a1 + a2 / lab), TINYV))
    h1 = b; h2 = a
    if (nh) then
      h1 = exp(lab * aa) * (aa - t * a1) / a
      h2 = exp(lab * aa) * (aa + (1.0_dp - t) * a1) / b
    end if
  end subroutine tev_eval

  ! t_{c1} CDF of a dual argument: d/dz T = t density (kd = log normaliser, lb = log B)
  pure function tcdf_hd(z, c1, kd, lb) result(r)
    type(hd), intent(in) :: z
    real(dp), intent(in) :: c1, kd, lb
    type(hd) :: r
    r = z
    r%v = pt_l(z%v, c1, lb)
    r%g = exp(kd - 0.5_dp * (c1 + 1.0_dp) * log1p_(z%v * z%v / c1)) * z%g
  end function tcdf_hd

  ! t-EV log-density over dual parameters (exact rho-score; nu is profiled by the caller,
  ! so P(2) carries no derivative and the df-dependence of T needs none either)
  pure function tev_lpdf_hd(P, a, b, kd, lbt) result(lp)
    type(hd), intent(in) :: P(3)
    real(dp), intent(in) :: a, b, kd, lbt
    type(hd) :: lp, k, z1, z2, d1, d2, dd1, dd2, f1, f2, fp1, fp2, t1, t2, aa, a1, a2
    real(dp) :: nu, c1, la, lb, lab, t, L, x1, x2
    nu = min(max(P(2)%v, 1.0_dp), 60.0_dp); c1 = nu + 1.0_dp
    la = log(a); lb = log(b); lab = la + lb
    t = min(max(lb / lab, 1.0e-300_dp), 1.0_dp - epsilon(1.0_dp))
    L = 1.0_dp / (t * (1.0_dp - t))
    x1 = exp((log(t) - log1p_(-t)) / nu); x2 = exp((log1p_(-t) - log(t)) / nu)
    k = hsqrt(c1 / ((1.0_dp - P(1)) * (1.0_dp + P(1))))
    z1 = k * (x1 - P(1)); z2 = k * (x2 - P(1))
    d1 = k * (x1 * L / nu); d2 = k * (-x2 * L / nu)
    dd1 = k * (x1 * L * L * (1.0_dp / nu - (1.0_dp - 2.0_dp * t)) / nu)
    dd2 = k * (x2 * L * L * (1.0_dp / nu + (1.0_dp - 2.0_dp * t)) / nu)
    f1 = hexp(kd - 0.5_dp * (c1 + 1.0_dp) * hlog1p(z1 * z1 / c1))
    f2 = hexp(kd - 0.5_dp * (c1 + 1.0_dp) * hlog1p(z2 * z2 / c1))
    fp1 = -(f1 * (c1 + 1.0_dp) * z1 / (c1 + z1 * z1))
    fp2 = -(f2 * (c1 + 1.0_dp) * z2 / (c1 + z2 * z2))
    t1 = tcdf_hd(z1, c1, kd, lbt); t2 = tcdf_hd(z2, c1, kd, lbt)
    aa = t * t1 + (1.0_dp - t) * t2
    a1 = t1 - t2 + t * f1 * d1 + (1.0_dp - t) * f2 * d2
    a2 = 2.0_dp * (f1 * d1 - f2 * d2) + t * (fp1 * d1 * d1 + f1 * dd1) + (1.0_dp - t) * (fp2 * d2 * d2 + f2 * dd2)
    lp = lab * aa - lab + hlog(aa * aa + (1.0_dp - 2.0_dp * t) * a1 * aa - (1.0_dp - t) * t * (a1 * a1 + a2 / lab))
  end function tev_lpdf_hd

  ! t_{c1} CDF of dual (z, c1): d/dz T = density f, d/dc1 T = dpt_dnu (exact)
  pure function tcdf_hd2(z, c1, f, cs) result(r)
    type(hd), intent(in) :: z, c1
    real(dp), intent(in) :: f, cs(3)
    type(hd) :: r
    r%v = pt_l(z%v, c1%v, cs(1))
    r%g = f * z%g + dpt_dnu(z%v, c1%v, cs) * c1%g
  end function tcdf_hd2

  ! t-EV log-density over dual (rho, nu): exact joint score, so (rho, nu) are fitted by
  ! one Fisher-scoring Newton instead of a Brent profile over nu.
  ! cs = [log B(c1/2, 1/2), psi(c1/2), psi((c1+1)/2), log normaliser of t_{c1}], c1 = nu+1
  pure function tev_lpdf_hd2(P, a, b, cs) result(lp)
    type(hd), intent(in) :: P(3)
    real(dp), intent(in) :: a, b, cs(4)
    type(hd) :: lp, nu, c1, kd, k, x1, x2, z1, z2, d1, d2, dd1, dd2, f1, f2, fp1, fp2, t1, t2, aa, a1, a2
    real(dp) :: la, lb, lab, t, L, lt
    nu = P(2); c1 = nu + 1.0_dp
    la = log(a); lb = log(b); lab = la + lb
    t = min(max(lb / lab, 1.0e-300_dp), 1.0_dp - epsilon(1.0_dp))
    L = 1.0_dp / (t * (1.0_dp - t))
    lt = log(t) - log1p_(-t)
    x1 = hexp(lt / nu); x2 = hexp(-lt / nu)
    k = hsqrt(c1 / ((1.0_dp - P(1)) * (1.0_dp + P(1))))
    z1 = k * (x1 - P(1)); z2 = k * (x2 - P(1))
    d1 = k * x1 * L / nu; d2 = -(k * x2 * L / nu)
    dd1 = k * x1 * (L * L) * (1.0_dp / nu - (1.0_dp - 2.0_dp * t)) / nu
    dd2 = k * x2 * (L * L) * (1.0_dp / nu + (1.0_dp - 2.0_dp * t)) / nu
    ! log normaliser of t_{c1}: d/dc = (psi((c+1)/2) - psi(c/2)) / 2 - 1/(2c)
    kd%v = cs(4); kd%g = (0.5_dp * (cs(3) - cs(2)) - 0.5_dp / c1%v) * c1%g
    f1 = hexp(kd - 0.5_dp * (c1 + 1.0_dp) * hlog1p(z1 * z1 / c1))
    f2 = hexp(kd - 0.5_dp * (c1 + 1.0_dp) * hlog1p(z2 * z2 / c1))
    fp1 = -(f1 * (c1 + 1.0_dp) * z1 / (c1 + z1 * z1))
    fp2 = -(f2 * (c1 + 1.0_dp) * z2 / (c1 + z2 * z2))
    t1 = tcdf_hd2(z1, c1, f1%v, cs(1:3)); t2 = tcdf_hd2(z2, c1, f2%v, cs(1:3))
    aa = t * t1 + (1.0_dp - t) * t2
    a1 = t1 - t2 + t * f1 * d1 + (1.0_dp - t) * f2 * d2
    a2 = 2.0_dp * (f1 * d1 - f2 * d2) + t * (fp1 * d1 * d1 + f1 * dd1) + (1.0_dp - t) * (fp2 * d2 * d2 + f2 * dd2)
    lp = lab * aa - lab + hlog(aa * aa + (1.0_dp - 2.0_dp * t) * a1 * aa - (1.0_dp - t) * t * (a1 * a1 + a2 / lab))
  end function tev_lpdf_hd2

  ! ---- t-EV likelihood by sorted cumulative quadrature ---------------------------
  ! A(t) needs the t_{nu+1} CDF at z = k (e^(+-w/nu) - rho) for every observation, which
  ! dominates the cost (2 incomplete-beta continued fractions per observation).  The 2n
  ! abscissae have a parameter-free order, so after one sort the CDF (and its exact
  ! df-derivative) follows from the previous value plus the integral of the density over
  ! the gap: Gauss-Legendre with 2-5 nodes, chosen by the gap so that each step's error is
  ! below 1e-17; gaps wider than 0.2 restart from the exact pt_l / dpt_dnu.  Exact to
  ! rounding like the continued fraction, ~3x cheaper; the df-derivative is ~10x cheaper.
  subroutine tev_prep(c)
    type(ctx_t), intent(inout) :: c
    integer :: i, n
    real(dp) :: la, lb, lab, t
    real(dp), allocatable :: key(:)
    integer, allocatable :: ix(:)
    n = c%n
    if (allocated(c%tq)) deallocate(c%tq, c%lq, c%wq, c%sq, c%pq)
    allocate(c%tq(n), c%lq(n), c%wq(n), c%sq(2 * n), c%pq(2 * n), key(2 * n), ix(2 * n))
    do i = 1, n
      la = log(c%a(i)); lb = log(c%b(i)); lab = la + lb
      t = min(max(lb / lab, 1.0e-300_dp), 1.0_dp - epsilon(1.0_dp))
      c%tq(i) = t; c%lq(i) = lab; c%wq(i) = log(t) - log1p_(-t)
      key(i) = c%wq(i); key(n + i) = -c%wq(i)
    end do
    do i = 1, 2 * n
      ix(i) = i
    end do
    call sort_idx(2 * n, key, ix)
    do i = 1, 2 * n
      c%sq(i) = key(ix(i))
      if (ix(i) <= n) then
        c%pq(i) = ix(i)
      else
        c%pq(i) = -(ix(i) - n)
      end if
    end do
  end subroutine tev_prep

  ! index sort (heap sort; ties keep any order)
  pure subroutine sort_idx(n, key, ix)
    integer, intent(in) :: n
    real(dp), intent(in) :: key(n)
    integer, intent(inout) :: ix(n)
    integer :: i, j, l, r, t
    if (n < 2) return
    l = n / 2 + 1; r = n
    do
      if (l > 1) then
        l = l - 1; t = ix(l)
      else
        t = ix(r); ix(r) = ix(1); r = r - 1
        if (r == 1) then
          ix(1) = t; return
        end if
      end if
      i = l; j = 2 * l
      do while (j <= r)
        if (j < r) then
          if (key(ix(j)) < key(ix(j + 1))) j = j + 1
        end if
        if (key(t) < key(ix(j))) then
          ix(i) = ix(j); i = j; j = 2 * j
        else
          j = r + 1
        end if
      end do
      ix(i) = t
    end do
  end subroutine sort_idx

  ! t_{c} CDF (and, if wd, d/dc at fixed z) at the sorted abscissae zs; chunks of 512 of
  ! the sorted range are independent (each restarts from the exact value), so nin threads
  ! share the pass.  cs = [log B(c/2,1/2), psi(c/2), psi((c+1)/2), log normaliser]
  subroutine tcdf_sorted(m, zs, c1, cs, wd, nin, ts, dts)
    integer, intent(in) :: m, nin
    real(dp), intent(in) :: zs(m), c1, cs(4)
    logical, intent(in) :: wd
    real(dp), intent(out) :: ts(m), dts(m)
    real(dp), parameter :: X2(1) = [0.5773502691896257645_dp], W2(1) = [1.0_dp]
    real(dp), parameter :: X3(2) = [0.0_dp, 0.7745966692414833770_dp], W3(2) = [8.0_dp / 9.0_dp, 5.0_dp / 9.0_dp]
    real(dp), parameter :: X4(2) = [0.3399810435848562648_dp, 0.8611363115940525752_dp], &
                           W4(2) = [0.6521451548625461427_dp, 0.3478548451374538574_dp]
    real(dp), parameter :: X5(3) = [0.0_dp, 0.5384693101056830910_dp, 0.9061798459386639928_dp], &
                           W5(3) = [0.5688888888888888889_dp, 0.4786286704993664680_dp, 0.2369268850561890875_dp]
    integer :: nch, ch, j0, j1, j
    real(dp) :: dk, h, mid, inc, dinc, csd(3)
    csd = cs(1:3)
    dk = 0.5_dp * (cs(3) - cs(2)) - 0.5_dp / c1
    nch = max(1, m / 512)                ! fixed chunks: the result does not depend on nin
    !$omp parallel do private(ch, j0, j1, j, h, mid, inc, dinc) num_threads(min(nin, nch)) if(min(nin, nch) > 1) &
    !$omp   schedule(static)
    do ch = 1, nch
      j0 = (ch - 1) * m / nch + 1; j1 = ch * m / nch
      ts(j0) = pt_l(zs(j0), c1, cs(1))
      dts(j0) = 0.0_dp
      if (wd) dts(j0) = dpt_dnu(zs(j0), c1, csd)
      do j = j0 + 1, j1
        h = zs(j) - zs(j - 1)
        if (h > 0.2_dp .or. .not. (h >= 0.0_dp)) then
          ts(j) = pt_l(zs(j), c1, cs(1))
          if (wd) dts(j) = dpt_dnu(zs(j), c1, csd)
          cycle
        end if
        mid = 0.5_dp * (zs(j) + zs(j - 1)); h = 0.5_dp * h
        if (h < 1.0e-3_dp) then
          call gl(X2, W2, .false., mid, h, inc, dinc)
        else if (h < 1.0e-2_dp) then
          call gl(X3, W3, .true., mid, h, inc, dinc)
        else if (h < 4.0e-2_dp) then
          call gl(X4, W4, .false., mid, h, inc, dinc)
        else
          call gl(X5, W5, .true., mid, h, inc, dinc)
        end if
        ts(j) = ts(j - 1) + inc
        if (wd) dts(j) = dts(j - 1) + dinc
      end do
    end do
    !$omp end parallel do
  contains
    ! symmetric Gauss-Legendre on [mid-h, mid+h]; zero0: the first node is the centre
    ! (everything passed explicitly: the caller's loop variables are OpenMP-private)
    subroutine gl(xn, wn, zero0, mid, h, inc, dinc)
      real(dp), intent(in) :: xn(:), wn(:), mid, h
      logical, intent(in) :: zero0
      real(dp), intent(out) :: inc, dinc
      integer :: q
      real(dp) :: u, l1, f, fac
      inc = 0.0_dp; dinc = 0.0_dp
      do q = 1, size(xn)
        fac = 1.0_dp
        if (.not. (zero0 .and. q == 1)) fac = 2.0_dp
        u = mid - h * xn(q)
        l1 = log1p_(u * u / c1); f = exp(cs(4) - 0.5_dp * (c1 + 1.0_dp) * l1)
        inc = inc + wn(q) * f
        if (wd) dinc = dinc + wn(q) * f * (dk - 0.5_dp * l1 + 0.5_dp * (c1 + 1.0_dp) * u * u / (c1 * (c1 + u * u)))
        if (fac > 1.5_dp) then
          u = mid + h * xn(q)
          l1 = log1p_(u * u / c1); f = exp(cs(4) - 0.5_dp * (c1 + 1.0_dp) * l1)
          inc = inc + wn(q) * f
          if (wd) dinc = dinc + wn(q) * f * (dk - 0.5_dp * l1 + 0.5_dp * (c1 + 1.0_dp) * u * u / (c1 * (c1 + u * u)))
        end if
      end do
      inc = h * inc; dinc = h * dinc
    end subroutine gl
  end subroutine tcdf_sorted

  ! t-EV log-likelihood (k = 0), with exact scores in rho (k = 1) or (rho, nu) (k = 2);
  ! c%lognu: p(2) is log nu.  BHHH curvature as in ll_hd.
  subroutine tev_ll(c, k, p, f, g, hm)
    type(ctx_t), intent(in) :: c
    integer, intent(in) :: k
    real(dp), intent(in) :: p(3)
    real(dp), intent(out) :: f, g(3), hm(3, 3)
    real(dp), allocatable :: zs(:), ts(:), dts(:), t1(:), t2(:), d1(:), d2(:)
    real(dp) :: rho, nu, c1, kk, cs(4), lp, tt, lab, w
    type(hd) :: P2(3), lpd
    integer :: m, j, i, n, jj, bk, nt
    real(dp) :: fb(NBK), gb(3, NBK), hb(3, 3, NBK)
    n = c%n; m = 2 * n
    rho = min(max(p(1), -0.999_dp), 0.999_dp)
    nu = p(2)
    if (c%lognu) nu = exp(p(2))
    nu = min(max(nu, 1.0_dp), 60.0_dp); c1 = nu + 1.0_dp
    cs = [lbeta_t(c1), digamma(0.5_dp * c1), digamma(0.5_dp * (c1 + 1.0_dp)), kdens_t(c1)]
    kk = sqrt(c1 / ((1.0_dp - rho) * (1.0_dp + rho)))
    allocate(zs(m), ts(m), dts(m), t1(n), t2(n), d1(n), d2(n))
    do j = 1, m
      zs(j) = kk * (exp(c%sq(j) / nu) - rho)
    end do
    nt = grab(c)
    call tcdf_sorted(m, zs, c1, cs, k >= 2, nt, ts, dts)
    call release(c, nt)
    do j = 1, m
      i = c%pq(j)
      if (i > 0) then
        t1(i) = ts(j); d1(i) = dts(j)
      else
        t2(-i) = ts(j); d2(-i) = dts(j)
      end if
    end do
    f = 0.0_dp; g = 0.0_dp; hm = 0.0_dp
    nt = grab(c)
    if (k == 0) then
      !$omp parallel do private(bk, i, lp) num_threads(nt) if(nt > 1) schedule(dynamic)
      do bk = 1, NBK
      fb(bk) = 0.0_dp
      do i = bk0(bk, n), bk0(bk + 1, n) - 1
        lp = tev_lp_T(rho, nu, c1, kk, cs(4), c%tq(i), c%lq(i), c%wq(i), t1(i), t2(i))
        if (.not. finite_(lp)) lp = LOGTINY
        if (c%usew) then
          fb(bk) = fb(bk) + c%w(i) * lp
        else
          fb(bk) = fb(bk) + lp
        end if
      end do
      end do
      !$omp end parallel do
      call release(c, nt)
      f = sum(fb)
      return
    end if
    do j = 1, 3
      if (j <= k) then
        P2(j) = hvar(p(j), j)
      else
        P2(j) = hconst(p(j))
      end if
    end do
    if (c%lognu) P2(2) = hexp(P2(2))
    !$omp parallel do private(bk, i, jj, lpd) num_threads(nt) if(nt > 1) schedule(dynamic)
    do bk = 1, NBK
    fb(bk) = 0.0_dp; gb(:, bk) = 0.0_dp; hb(:, :, bk) = 0.0_dp
    do i = bk0(bk, n), bk0(bk + 1, n) - 1
      lpd = tev_lp_Thd(P2, cs, c%tq(i), c%lq(i), c%wq(i), t1(i), d1(i), t2(i), d2(i))
      if (.not. (finite_(lpd%v) .and. all(abs(lpd%g) <= huge(1.0_dp)))) lpd = hconst(LOGTINY)
      if (c%usew) then
        fb(bk) = fb(bk) + c%w(i) * lpd%v; gb(:, bk) = gb(:, bk) + c%w(i) * lpd%g
        do jj = 1, k
          hb(1:k, jj, bk) = hb(1:k, jj, bk) - c%w(i) * lpd%g(1:k) * lpd%g(jj)
        end do
      else
        fb(bk) = fb(bk) + lpd%v; gb(:, bk) = gb(:, bk) + lpd%g
        do jj = 1, k
          hb(1:k, jj, bk) = hb(1:k, jj, bk) - lpd%g(1:k) * lpd%g(jj)
        end do
      end if
    end do
    end do
    !$omp end parallel do
    call release(c, nt)
    f = sum(fb); g = sum(gb, 2); hm = sum(hb, 3)
  end subroutine tev_ll

  ! t-EV log-density of one observation given T1 = T(z1), T2 = T(z2) (tev_eval)
  pure real(dp) function tev_lp_T(r, nu, c1, k, kd, t, lab, w, T1, T2) result(lp)
    real(dp), intent(in) :: r, nu, c1, k, kd, t, lab, w, T1, T2
    real(dp) :: L, x1, x2, z1, z2, d1, d2, dd1, dd2, f1, f2, fp1, fp2, aa, a1, a2
    L = 1.0_dp / (t * (1.0_dp - t))
    x1 = exp(w / nu); x2 = exp(-w / nu)
    z1 = k * (x1 - r); z2 = k * (x2 - r)
    d1 = k / nu * x1 * L; d2 = -k / nu * x2 * L
    dd1 = k / nu * x1 * L * L * (1.0_dp / nu - (1.0_dp - 2.0_dp * t))
    dd2 = k / nu * x2 * L * L * (1.0_dp / nu + (1.0_dp - 2.0_dp * t))
    f1 = exp(kd - 0.5_dp * (c1 + 1.0_dp) * log1p_(z1 * z1 / c1))
    f2 = exp(kd - 0.5_dp * (c1 + 1.0_dp) * log1p_(z2 * z2 / c1))
    fp1 = -f1 * (c1 + 1.0_dp) * z1 / (c1 + z1 * z1)
    fp2 = -f2 * (c1 + 1.0_dp) * z2 / (c1 + z2 * z2)
    aa = t * T1 + (1.0_dp - t) * T2
    a1 = T1 - T2 + t * f1 * d1 + (1.0_dp - t) * f2 * d2
    a2 = 2.0_dp * (f1 * d1 - f2 * d2) + t * (fp1 * d1 * d1 + f1 * dd1) + (1.0_dp - t) * (fp2 * d2 * d2 + f2 * dd2)
    lp = lab * aa - lab + log(max(aa * aa + (1.0_dp - 2.0_dp * t) * a1 * aa &
         - (1.0_dp - t) * t * (a1 * a1 + a2 / lab), TINYV))
  end function tev_lp_T

  ! dual counterpart: P(1) = rho, P(2) = nu (either may carry derivatives); dT1, dT2 are
  ! the exact df-derivatives of T at fixed z
  pure function tev_lp_Thd(P, cs, t, lab, w, T1v, dT1, T2v, dT2) result(lp)
    type(hd), intent(in) :: P(3)
    real(dp), intent(in) :: cs(4), t, lab, w, T1v, dT1, T2v, dT2
    type(hd) :: lp, nu, c1, kd, k, x1, x2, z1, z2, d1, d2, dd1, dd2, f1, f2, fp1, fp2, t1, t2, aa, a1, a2
    real(dp) :: L
    nu = P(2); c1 = nu + 1.0_dp
    L = 1.0_dp / (t * (1.0_dp - t))
    x1 = hexp(w / nu); x2 = hexp(-w / nu)
    k = hsqrt(c1 / ((1.0_dp - P(1)) * (1.0_dp + P(1))))
    z1 = k * (x1 - P(1)); z2 = k * (x2 - P(1))
    d1 = k * x1 * L / nu; d2 = -(k * x2 * L / nu)
    dd1 = k * x1 * (L * L) * (1.0_dp / nu - (1.0_dp - 2.0_dp * t)) / nu
    dd2 = k * x2 * (L * L) * (1.0_dp / nu + (1.0_dp - 2.0_dp * t)) / nu
    kd%v = cs(4); kd%g = (0.5_dp * (cs(3) - cs(2)) - 0.5_dp / c1%v) * c1%g
    f1 = hexp(kd - 0.5_dp * (c1 + 1.0_dp) * hlog1p(z1 * z1 / c1))
    f2 = hexp(kd - 0.5_dp * (c1 + 1.0_dp) * hlog1p(z2 * z2 / c1))
    fp1 = -(f1 * (c1 + 1.0_dp) * z1 / (c1 + z1 * z1))
    fp2 = -(f2 * (c1 + 1.0_dp) * z2 / (c1 + z2 * z2))
    t1%v = T1v; t1%g = f1%v * z1%g + dT1 * c1%g
    t2%v = T2v; t2%g = f2%v * z2%g + dT2 * c1%g
    aa = t * t1 + (1.0_dp - t) * t2
    a1 = t1 - t2 + t * f1 * d1 + (1.0_dp - t) * f2 * d2
    a2 = 2.0_dp * (f1 * d1 - f2 * d2) + t * (fp1 * d1 * d1 + f1 * dd1) + (1.0_dp - t) * (fp2 * d2 * d2 + f2 * dd2)
    lp = lab * aa - lab + hlog(aa * aa + (1.0_dp - 2.0_dp * t) * a1 * aa - (1.0_dp - t) * t * (a1 * a1 + a2 / lab))
  end function tev_lp_Thd

  ! inverse of h1 in the second argument: find v with h1(a, v) = w.
  pure real(dp) function hinv10(fam, pp, a, w0)
    integer, intent(in) :: fam
    real(dp), intent(in) :: pp(NPP), a, w0
    real(dp) :: p1, p2, w, r, nu, x, q, th, eu, ed, lo, hi, v, lp, h1, h2, s
    integer :: it
    p1 = pp(1); p2 = pp(2)
    w = clampu(w0)
    select case (fam)
    case (0)
      hinv10 = w
    case (1)
      r = min(max(p1, -0.999999_dp), 0.999999_dp)
      hinv10 = pnorm(qnorm(w) * sqrt(1.0_dp - r * r) + r * qnorm(a))
    case (2)
      r = min(max(p1, -0.999999_dp), 0.999999_dp)
      nu = min(max(p2, 2.0_dp), 100.0_dp)
      x = qt_l(a, nu, pp(5), pp(6)); q = qt_l(w, nu + 1.0_dp, pp(7), pp(8))
      hinv10 = pt_l(r * x + q * sqrt((nu + x * x) * (1.0_dp - r) * (1.0_dp + r) / (nu + 1.0_dp)), nu, pp(6))
    case (4)
      hinv10 = hinv_gumbel(max(p1, 1.0_dp), a, w)
    case (6)
      hinv10 = hinv_joe(max(p1, 1.0_dp), a, w)
    case (3)
      th = p1
      if (th < 1.0e-10_dp) then
        hinv10 = w; return
      end if
      ! v = (1 + a^-th * ((w^(-th/(th+1)) - 1)))^(-1/th), written without the
      ! cancellation of (w a^(th+1))^(-th/(th+1)) - a^-th for small a
      s = exp(-th * log(a)) * expm1_((-th / (th + 1.0_dp)) * log(w))
      hinv10 = exp(-log1p_(s) / th)
    case (5)
      th = p1
      if (abs(th) < 1.0e-10_dp) then
        hinv10 = w; return
      end if
      eu = expm1_(-th * a); ed = expm1_(-th)
      hinv10 = (-1.0_dp / th) * log1p_(w * ed / max(eu + 1.0_dp - w * eu, TINYV))
    case default
      hinv10 = newton_hinv(fam, pp, a, w, 1)
    end select
    hinv10 = min(max(hinv10, 0.0_dp), 1.0_dp)
  end function hinv10

  ! inverse of h2 in the first argument: find v with h2(v, b) = w.  Exchangeable
  ! families reuse hinv10 by symmetry; Tawn (asymmetric) is inverted directly.
  pure real(dp) function hinv20(fam, pp, b, w0)
    integer, intent(in) :: fam
    real(dp), intent(in) :: pp(NPP), b, w0
    real(dp) :: w, lo, hi, v, lp, h1, h2
    integer :: it
    if (fam /= 11 .and. .not. (gen_is(fam) .and. .not. gen_exch(fam))) then
      hinv20 = hinv10(fam, pp, b, w0); return
    end if
    w = clampu(w0)
    hinv20 = newton_hinv(fam, pp, b, w, 2)
  end function hinv20

  ! Gumbel hinv: with x = -log u, L = -log w and d = (x^th + y^th)^(1/th) - x,
  ! log h1 = -d - (th-1) log1p(d/x), so solve G(d) = d + (th-1) log1p(d/x) = L.
  ! G is increasing and concave: Newton from the lower bound L/G'(0) is monotone.
  pure real(dp) function hinv_gumbel(th, a, w)
    real(dp), intent(in) :: th, a, w
    real(dp) :: x, L, d, dn, g, gp, y
    integer :: it
    if (th <= 1.0_dp + 1.0e-12_dp) then
      hinv_gumbel = w; return
    end if
    x = -log(a); L = -log(w)
    d = L / (1.0_dp + (th - 1.0_dp) / x)
    do it = 1, 40
      g = d + (th - 1.0_dp) * log1p_(d / x) - L
      gp = 1.0_dp + (th - 1.0_dp) / (x + d)
      dn = d - g / gp
      if (abs(dn - d) <= 4.0e-16_dp * max(dn, 1.0e-300_dp)) then
        d = dn; exit
      end if
      d = dn
    end do
    y = x * exp(log(expm1_(th * log1p_(d / x))) / th)
    hinv_gumbel = exp(-y)
  end function hinv_gumbel

  ! Joe hinv: x = -log(1 - (1-u)^th), K = 1/expm1(x), L = -log w; solve
  ! G(d) = d + (1 - 1/th) log1p(K (1 - e^-d)) = L (concave), then
  ! (1-v)^th = 1 - e^-d.
  pure real(dp) function hinv_joe(th, a, w)
    real(dp), intent(in) :: th, a, w
    real(dp) :: x, k, L, c1, d, dn, g, gp, em, b
    integer :: it
    if (th <= 1.0_dp + 1.0e-12_dp) then
      hinv_joe = w; return
    end if
    x = -log1mexp(th * log1p_(-a))
    k = 1.0_dp / expm1_(x)
    L = -log(w)
    c1 = 1.0_dp - 1.0_dp / th
    d = L / (1.0_dp + c1 * k)
    do it = 1, 60
      em = -expm1_(-d)
      g = d + c1 * log1p_(k * em) - L
      gp = 1.0_dp + c1 * k * exp(-d) / (1.0_dp + k * em)
      dn = d - g / gp
      if (abs(dn - d) <= 4.0e-16_dp * max(dn, 1.0e-300_dp)) then
        d = dn; exit
      end if
      d = dn
    end do
    b = -expm1_(-d)
    hinv_joe = -expm1_(log(b) / th)
  end function hinv_joe

  ! Solve h(v) = w for v in (0,1): Newton with the copula density as the exact
  ! derivative (d h1(a,v)/dv = c(a,v), d h2(v,b)/dv = c(v,b)), bracketed and
  ! falling back to bisection; typically 3-6 evaluations instead of 52.
  pure real(dp) function newton_hinv(fam, pp, cond, w, which)
    integer, intent(in) :: fam, which
    real(dp), intent(in) :: pp(NPP), cond, w
    real(dp) :: lo, hi, v, vn, lp, h1, h2, f, dens
    integer :: it
    lo = 0.0_dp; hi = 1.0_dp
    v = w                                   ! independence solution
    do it = 1, 100
      if (which == 1) then
        call eval0(fam, pp, cond, clampu(v), .true., lp, h1, h2)
        f = h1 - w
      else
        call eval0(fam, pp, clampu(v), cond, .true., lp, h1, h2)
        f = h2 - w
      end if
      if (f < 0.0_dp) then
        lo = v
      else
        hi = v
      end if
      if (abs(f) <= 1.0e-15_dp .or. hi - lo <= 1.0e-15_dp) exit
      dens = exp(lp)
      vn = v - f / dens
      if (.not. (vn > lo .and. vn < hi)) vn = 0.5_dp * (lo + hi)
      if (abs(vn - v) <= 1.0e-15_dp * max(v, 1.0e-300_dp)) then
        v = vn; exit
      end if
      v = vn
    end do
    newton_hinv = v
  end function newton_hinv

  ! --------------------------------------------------------- rotation layer
  pure subroutine rotate(rot, u1, u2, a, b)
    integer, intent(in) :: rot
    real(dp), intent(in) :: u1, u2
    real(dp), intent(out) :: a, b
    select case (rot)
    case (90)
      a = u2; b = 1.0_dp - u1
    case (180)
      a = 1.0_dp - u1; b = 1.0_dp - u2
    case (270)
      a = 1.0_dp - u2; b = u1
    case default
      a = u1; b = u2
    end select
  end subroutine rotate

  pure subroutine eval_rot(fam, rot, pp, u1, u2, nh, lp, h1, h2)
    integer, intent(in) :: fam, rot
    real(dp), intent(in) :: pp(NPP), u1, u2
    logical, intent(in) :: nh
    real(dp), intent(out) :: lp, h1, h2
    real(dp) :: a, b, g1, g2
    call rotate(rot, clampu(u1), clampu(u2), a, b)
    call eval0(fam, pp, a, b, nh, lp, g1, g2)
    select case (rot)
    case (90)
      h1 = g2; h2 = 1.0_dp - g1
    case (180)
      h1 = 1.0_dp - g1; h2 = 1.0_dp - g2
    case (270)
      h1 = 1.0_dp - g2; h2 = g1
    case default
      h1 = g1; h2 = g2
    end select
  end subroutine eval_rot

  pure real(dp) function hinv_rot(fam, rot, pp, u1, u2, which)
    integer, intent(in) :: fam, rot, which
    real(dp), intent(in) :: pp(NPP), u1, u2
    real(dp) :: a, b, v
    call rotate(rot, clampu(u1), clampu(u2), a, b)
    ! hinv1_0(a,b) = hinv10(a,b);  hinv2_0(a,b) = hinv20(b,a)
    select case (rot*10 + which)
    case (1);    v = hinv10(fam, pp, a, b)
    case (2);    v = hinv20(fam, pp, b, a)
    case (901);  v = hinv20(fam, pp, b, a)
    case (902);  v = 1.0_dp - hinv10(fam, pp, a, b)
    case (1801); v = 1.0_dp - hinv10(fam, pp, a, b)
    case (1802); v = 1.0_dp - hinv20(fam, pp, b, a)
    case (2701); v = 1.0_dp - hinv20(fam, pp, b, a)
    case (2702); v = hinv10(fam, pp, a, b)
    case default; v = hinv10(fam, pp, a, b)
    end select
    hinv_rot = min(max(v, 0.0_dp), 1.0_dp)
  end function hinv_rot

  ! ======================================================================
  ! Analytic derivatives: log-density over dual-number parameters (exact
  ! scores), and a projected Fisher-scoring / Levenberg-Marquardt MLE.
  ! ======================================================================

  ! Parameter-independent per-observation transforms, computed once per fit.
  pure subroutine prep_xy(fam, a, b, x, y)
    integer, intent(in) :: fam
    real(dp), intent(in) :: a, b
    real(dp), intent(out) :: x, y
    select case (fam)
    case (1)
      x = qnorm(a); y = qnorm(b)
    case (4)
      x = -log(a); y = -log(b)
    case (7, 11, 12:22)
      x = log(a); y = log(b)
    case (8, 9)
      x = log1p_(-a); y = log1p_(-b)
    case default
      x = a; y = b
    end select
  end subroutine prep_xy

  pure subroutine arch_hd(fam, th, de, x, phi, d1, d2)
    integer, intent(in) :: fam
    type(hd), intent(in) :: th, de, x
    type(hd), intent(out) :: phi, d1, d2
    type(hd) :: t1, ox, tmp, l, nrm
    select case (fam)
    case (7)
      t1 = x**(-th) - 1.0_dp
      phi = t1**de
      d1 = -(de * th) * x**(-(th + 1.0_dp)) * t1**(de - 1.0_dp)
      d2 = de * th * t1**de * ((1.0_dp + de * th) - (1.0_dp + th) * x**th) &
           / ((x**th - 1.0_dp) * (x**th - 1.0_dp) * x * x)
    case (8)
      ox = 1.0_dp - x
      tmp = ox**th
      l = hlog1p(-tmp)
      phi = (-l)**de
      d1 = de * th * (-l)**(de - 1.0_dp) * ox**(th - 1.0_dp) / (tmp - 1.0_dp)
      d2 = (-l)**(de - 2.0_dp) * de * th * ox**(th - 2.0_dp) &
           * ((de - 1.0_dp) * th * tmp - (tmp + th - 1.0_dp) * l) / ((tmp - 1.0_dp) * (tmp - 1.0_dp))
    case (9)
      ox = 1.0_dp - x
      tmp = ox**th
      phi = (1.0_dp - tmp)**(-de) - 1.0_dp
      d1 = -(de * th) * (1.0_dp - tmp)**(-de - 1.0_dp) * ox**(th - 1.0_dp)
      d2 = de * th * (1.0_dp - tmp)**(-de - 2.0_dp) * ox**(th - 2.0_dp) &
           * (th - 1.0_dp + (1.0_dp + de * th) * tmp)
    case default  ! 10, bb8
      nrm = 1.0_dp - (1.0_dp - de)**th
      ox = 1.0_dp - de * x
      tmp = ox**th
      phi = -hlog((1.0_dp - tmp) / nrm)
      d1 = -(de * th) * ox**(th - 1.0_dp) / (1.0_dp - tmp)
      d2 = de * de * th * ox**(th - 2.0_dp) * (th - 1.0_dp + tmp) / ((tmp - 1.0_dp) * (tmp - 1.0_dp))
    end select
  end subroutine arch_hd

  pure function arch_psi_hd(fam, th, de, t) result(r)
    integer, intent(in) :: fam
    type(hd), intent(in) :: th, de, t
    type(hd) :: r
    select case (fam)
    case (7)
      r = (t**(1.0_dp / de) + 1.0_dp)**(-1.0_dp / th)
    case (8)
      r = 1.0_dp - (-hexpm1(-(t**(1.0_dp / de))))**(1.0_dp / th)
    case (9)
      r = 1.0_dp - (1.0_dp - (1.0_dp + t)**(-1.0_dp / de))**(1.0_dp / th)
    case default
      r = (1.0_dp - (1.0_dp + hexp(-t) * ((1.0_dp - de)**th - 1.0_dp))**(1.0_dp / th)) / de
    end select
  end function arch_psi_hd

  ! log(exp(p) + exp(q)) for dual numbers
  pure function lse(p, q) result(r)
    type(hd), intent(in) :: p, q
    type(hd) :: r
    if (p%v >= q%v) then
      r = p + hlog1p(hexp(q - p))
    else
      r = q + hlog1p(hexp(p - q))
    end if
  end function lse

  ! log c(a, b; P) with exact parameter derivatives, once per dual type
#define VT type(hd)
#define VCONST hconst
#define LPDF_NAME lpdf_hd
#define EXT_NAME ext_lpdf_hd
#define GEN_NAME gen_ext_lpdf_hd
#include "families_core.inc"
#include "families_ext.inc"
#include "families_gen.inc"
#undef VT
#undef VCONST
#undef LPDF_NAME
#undef EXT_NAME
#undef GEN_NAME
#define VT type(hd1)
#define VCONST h1const
#define LPDF_NAME lpdf_h1
#define EXT_NAME ext_lpdf_h1
#define GEN_NAME gen_ext_lpdf_h1
#include "families_core.inc"
#include "families_ext.inc"
#include "families_gen.inc"
#undef VT
#undef VCONST
#undef LPDF_NAME
#undef EXT_NAME
#undef GEN_NAME
#define VT type(hd2)
#define VCONST h2const
#define LPDF_NAME lpdf_h2
#define EXT_NAME ext_lpdf_h2
#define GEN_NAME gen_ext_lpdf_h2
#include "families_core.inc"
#include "families_ext.inc"
#include "families_gen.inc"
#undef VT
#undef VCONST
#undef LPDF_NAME
#undef EXT_NAME
#undef GEN_NAME


  ! first observation of block b of n (bk0(NBK + 1, n) = n + 1)
  pure integer function bk0(b, n)
    integer, intent(in) :: b, n
    bk0 = 1 + int((int(b - 1, int64) * n) / NBK)
  end function bk0

  ! threads for one likelihood pass: c%nin, or a fair share (budget / running fits, at most one
  ! per 2000 observations: smaller teams cost more than they save) of the selection's idle
  ! threads, reserved atomically so the budget is never exceeded
  integer function grab(c)
    type(ctx_t), intent(in) :: c
    integer :: nb, want, old
    grab = c%nin
    if (.not. associated(c%pool)) return
    !$omp atomic read
    nb = c%pool(2)
    want = max(0, min(c%nin / max(1, nb), c%n / 2000) - 1)
    !$omp atomic capture
    old = c%pool(1)
    c%pool(1) = c%pool(1) + want
    !$omp end atomic
    grab = max(0, min(want, c%nin - old))
    if (grab < want) then
      !$omp atomic
      c%pool(1) = c%pool(1) - (want - grab)
    end if
    grab = 1 + grab
  end function grab

  subroutine release(c, nt)
    type(ctx_t), intent(in) :: c
    integer, intent(in) :: nt
    if (associated(c%pool) .and. nt > 1) then
      !$omp atomic
      c%pool(1) = c%pool(1) - (nt - 1)
    end if
  end subroutine release

  ! Log-likelihood, exact scores and BHHH curvature of a one- or two-parameter family
  ! (k = npar) through the width-specialised dual types hd1 / hd2 (as ll_hd, minus the
  ! unused score slots of hd).
  subroutine ll_hk(c, k, p, f, g, hm)
    type(ctx_t), intent(in) :: c
    integer, intent(in) :: k
    real(dp), intent(in) :: p(3)
    real(dp), intent(out) :: f, g(3), hm(3, 3)
    type(hd1) :: P1(3), l1
    type(hd2) :: P2(3), l2
    real(dp) :: bp(3), wi, lv, lg(2), fb(NBK), gb(2, NBK), bb(3, NBK)
    integer :: i, bk, nt
    f = 0.0_dp; g = 0.0_dp; bp = 0.0_dp
    if (k == 1) then
      P1 = [h1var(p(1), 1), h1const(p(2)), h1const(p(3))]
    else
      P2 = [h2var(p(1), 1), h2var(p(2), 2), h2const(p(3))]
    end if
    nt = grab(c)
    !$omp parallel do private(bk, i, l1, l2, lv, lg, wi) num_threads(nt) if(nt > 1) schedule(dynamic)
    do bk = 1, NBK
    fb(bk) = 0.0_dp; gb(:, bk) = 0.0_dp; bb(:, bk) = 0.0_dp
    do i = bk0(bk, c%n), bk0(bk + 1, c%n) - 1
      lg = 0.0_dp
      if (k == 1) then
        l1 = lpdf_h1(c%fam, P1, c%a(i), c%b(i), c%x(i), c%y(i))
        lv = l1%v; lg(1) = l1%g(1)
      else
        l2 = lpdf_h2(c%fam, P2, c%a(i), c%b(i), c%x(i), c%y(i))
        lv = l2%v; lg = l2%g
      end if
      if (.not. (finite_(lv) .and. all(abs(lg) <= huge(1.0_dp)))) then
        lv = LOGTINY; lg = 0.0_dp
      end if
      wi = 1.0_dp
      if (c%usew) wi = c%w(i)
      fb(bk) = fb(bk) + wi * lv
      gb(:, bk) = gb(:, bk) + wi * lg
      bb(:, bk) = bb(:, bk) - wi * [lg(1) * lg(1), lg(2) * lg(1), lg(2) * lg(2)]
    end do
    end do
    !$omp end parallel do
    call release(c, nt)
    f = sum(fb); g(1:2) = sum(gb, 2); bp = sum(bb, 2)
    hm = 0.0_dp
    hm(1:k, 1:k) = reshape([bp(1), bp(2), bp(2), bp(3)], [k, k])
  end subroutine ll_hk

  ! total weighted log-likelihood, exact gradient, BHHH curvature (outer product of scores)
  subroutine ll_hd(c, k, p, f, g, hm)
    type(ctx_t), intent(in) :: c
    integer, intent(in) :: k
    real(dp), intent(in) :: p(3)
    real(dp), intent(out) :: f, g(3), hm(3, 3)
    type(hd) :: PP(3), lp
    integer :: i, j, bk, nt
    real(dp) :: kd, lbt, cs4(4), c1, p2, fb(NBK), gb(3, NBK), hb(3, 3, NBK)
    if (k == npar_of(c%fam) .and. k <= 2 .and. c%fam /= 2 .and. c%fam /= 23) then
      call ll_hk(c, k, p, f, g, hm)          ! one- and two-parameter families: narrower duals
      return
    end if
    if (c%fam == 23 .and. allocated(c%pq)) then
      call tev_ll(c, k, p, f, g, hm)
      return
    end if
    do j = 1, 3
      if (j <= k) then
        PP(j) = hvar(p(j), j)
      else
        PP(j) = hconst(p(j))
      end if
    end do
    p2 = p(2)
    if (c%lognu) then
      PP(2) = hexp(PP(2)); p2 = exp(p(2))
    end if
    f = 0.0_dp; g = 0.0_dp; hm = 0.0_dp
    kd = 0.0_dp; lbt = 0.0_dp
    if (c%fam == 23) then
      c1 = min(max(p2, 1.0_dp), 60.0_dp) + 1.0_dp
      kd = kdens_t(c1); lbt = lbeta_t(c1)
      cs4 = [lbt, digamma(0.5_dp * c1), digamma(0.5_dp * (c1 + 1.0_dp)), kd]
    end if
    nt = grab(c)
    !$omp parallel do private(bk, i, j, lp) num_threads(nt) if(nt > 1) schedule(dynamic)
    do bk = 1, NBK
    fb(bk) = 0.0_dp; gb(:, bk) = 0.0_dp; hb(:, :, bk) = 0.0_dp
    do i = bk0(bk, c%n), bk0(bk + 1, c%n) - 1
      if (c%fam == 23 .and. k >= 2) then
        lp = tev_lpdf_hd2(PP, c%a(i), c%b(i), cs4)
      else if (c%fam == 23) then
        lp = tev_lpdf_hd(PP, c%a(i), c%b(i), kd, lbt)
      else
        lp = lpdf_hd(c%fam, PP, c%a(i), c%b(i), c%x(i), c%y(i))
      end if
      if (.not. (finite_(lp%v) .and. all(abs(lp%g) <= huge(1.0_dp)))) lp = hconst(LOGTINY)
      ! BHHH / Fisher scoring: curvature = -sum of score outer products
      if (c%usew) then
        fb(bk) = fb(bk) + c%w(i) * lp%v; gb(:, bk) = gb(:, bk) + c%w(i) * lp%g
        do j = 1, k
          hb(1:k, j, bk) = hb(1:k, j, bk) - c%w(i) * lp%g(1:k) * lp%g(j)
        end do
      else
        fb(bk) = fb(bk) + lp%v; gb(:, bk) = gb(:, bk) + lp%g
        do j = 1, k
          hb(1:k, j, bk) = hb(1:k, j, bk) - lp%g(1:k) * lp%g(j)
        end do
      end if
    end do
    end do
    !$omp end parallel do
    call release(c, nt)
    f = sum(fb); g = sum(gb, 2); hm = sum(hb, 3)
  end subroutine ll_hd

  ! Maximise the log-likelihood in the box [lb, ub]: projected Fisher-scoring
  ! (BHHH) steps with exact analytic scores, Levenberg-Marquardt damping.  For BB8, Husler-Reiss,
  ! BB2, BB9 and BB10 a positive parameter whose box spans >= 2 decades is optimised in log scale
  ! (z = log x, score x dl/dx): the same MLE in fewer steps (BB9's delta took 40 steps of x1.3
  ! from 0.03 to 12).  Measured per family: Galambos, BB3, BB4 and BB5 converge faster in x.
  subroutine newton_fit(c, k, x0, lb, ub, xb, fb)
    type(ctx_t), intent(inout) :: c
    integer, intent(in) :: k
    real(dp), intent(in) :: x0(3), lb(3), ub(3)
    real(dp), intent(out) :: xb(3), fb
    real(dp) :: x(3), lo(3), hi(3), f, g(3), hm(3, 3), xn(3), fn, gn(3), hn(3, 3), a(3, 3), st(3), lam, sc
    real(dp) :: bm(3, 3), sv(3), yv(3), bs(3), sbs, ys, sp(3), r, xe(3), fe, ge(3), he(3, 3)
    logical :: act(3), ok, ls(3)
    integer :: it, tr, i
    lo = lb; hi = ub
    sp = 0.0_dp
    ls = .false.
    do i = 1, k
      sc = 1.0e-7_dp * max(ub(i) - lb(i), 1.0e-8_dp)
      lo(i) = lb(i) + sc; hi(i) = ub(i) - sc
      ls(i) = lb(i) > 0.0_dp .and. ub(i) >= 100.0_dp * lb(i) .and. any(c%fam == [10, 16, 17, 20, 21])
    end do
    x = x0
    x(1:k) = min(max(x0(1:k), lo(1:k)), hi(1:k))
    x = zof(x); lo = zof(lo); hi = zof(hi)
    call ev(x, f, g, hm)
    bm = -hm                    ! curvature model (negative Hessian of the log-likelihood)
    lam = 0.0_dp
    do it = 1, 40
      if (fv_abort) exit
      act = .false.
      do i = 1, k
        act(i) = (x(i) <= lo(i) .and. g(i) < 0.0_dp) .or. (x(i) >= hi(i) .and. g(i) > 0.0_dp)
      end do
      if (all(act(1:k) .or. abs(g(1:k)) <= 0.1_dp * c%tol * (1.0_dp + abs(f)))) exit
      ok = .false.
      do tr = 1, 8
        a = bm
        do i = 1, k
          a(i, i) = a(i, i) + lam * max(abs(bm(i, i)), 1.0e-8_dp)
        end do
        call solve_spd(k, a, g, act, st, ok)
        if (ok) then
          ! undamped Newton decrement already negligible -> converged, no trial needed
          if (lam == 0.0_dp .and. 0.5_dp * abs(sum(g(1:k) * st(1:k))) <= 0.1_dp * c%tol * (1.0_dp + abs(f))) then
            ok = .false.; exit
          end if
          if (lam == 0.0_dp .and. it > 1) then
            ! racing: remaining gain <= max(5 d, d + 0.8 |f| + b), d the predicted gain of this step,
            ! b = 5 (1 for one parameter).  2339 traced polishes, n = 1500-30000: never above
            ! max(5 d, d + 0.37 |f| + 5); the 653 one-parameter ones never above 5 d + 0.08
            sc = 0.5_dp * sum(g(1:k) * st(1:k))
            sc = f + max(5.0_dp * sc, sc + 0.8_dp * abs(f) + merge(1.0_dp, 5.0_dp, k == 1))
            if (sc < c%tgt) then
              c%cut = sc; ok = .false.; exit
            end if
          end if
          xn = x
          xn(1:k) = min(max(x(1:k) + st(1:k), lo(1:k)), hi(1:k))
          if (c%qn) then
            call ev(xn, fn, gn, hn)   ! scores are cheap relative to the value here
          else
            fn = ll_par(c, xof(xn))     ! value-only trial; scores only for accepted steps
          end if
          if (finite_(fn) .and. fn >= f - 1.0e-12_dp * (1.0_dp + abs(f))) exit
        end if
        ! rejected: damp to about the length t (in [0.1, 0.5]) at which the quadratic through f,
        ! the step's directional derivative and fn peaks (backtracking by interpolation; the
        ! old x10 from 1e-3 re-tried near-identical steps)
        sc = 0.1_dp
        if (ok .and. finite_(fn)) then
          sc = sum(g(1:k) * st(1:k))
          if (sc > 0.0_dp) then
            sc = min(0.5_dp, max(0.1_dp, sc / (2.0_dp * (sc - (fn - f)))))
          else
            sc = 0.1_dp
          end if
        end if
        lam = (1.0_dp + lam) / sc - 1.0_dp
        ok = .false.
      end do
      if (.not. ok) exit
      if (.not. c%qn) call ev(xn, fn, gn, hn)
      ! Newton decrement: predicted gain g's/2 is negligible -> converged
      if (abs(fn - f) <= c%tol * (1.0_dp + abs(f)) .and. &
          0.5_dp * abs(sum(g(1:k) * st(1:k))) <= c%tol * (1.0_dp + abs(f))) then
        x = xn; f = fn; exit
      end if
      if (maxval(abs(xn(1:k) - x(1:k)) / (1.0_dp + abs(x(1:k)))) < 1.0e-9_dp) then
        x = xn; f = fn; exit
      end if
      ! Aitken: steps of one direction shrinking by a constant ratio r (linear convergence, e.g.
      ! a crawl along a ridge or towards a bound) -> try their limit x + s r / (1 - r); kept
      ! only if it improves the log-likelihood
      sv(1:k) = xn(1:k) - x(1:k)
      r = sqrt(sum(sv(1:k)**2) / max(sum(sp(1:k)**2), tiny(1.0_dp)))
      if (r > 0.5_dp .and. r < 1.0_dp .and. sum(sv(1:k) * sp(1:k)) > 0.99_dp * r * sum(sp(1:k)**2)) then
        xe = xn
        xe(1:k) = min(max(xn(1:k) + sv(1:k) * (r / (1.0_dp - r)), lo(1:k)), hi(1:k))
        call ev(xe, fe, ge, he)
        if (finite_(fe) .and. fe > fn) then
          xn = xe; fn = fe; gn = ge; hn = he
        end if
      end if
      sp = sv
      if (c%qn) then
        ! BFGS secant update of the curvature model from the exact score change: BHHH alone
        ! zig-zags along strongly correlated ridges (t-EV rho / log nu); this restores
        ! superlinear convergence.  Reset to BHHH when the curvature condition fails.
        sv = 0.0_dp; yv = 0.0_dp
        sv(1:k) = xn(1:k) - x(1:k); yv(1:k) = g(1:k) - gn(1:k)
        ys = sum(yv(1:k) * sv(1:k))
        bs(1:k) = matmul(bm(1:k, 1:k), sv(1:k)); sbs = sum(sv(1:k) * bs(1:k))
        if (ys > 1.0e-12_dp * sqrt(sum(yv(1:k)**2) * sum(sv(1:k)**2)) .and. sbs > 0.0_dp) then
          do i = 1, k
            bm(1:k, i) = bm(1:k, i) - bs(1:k) * bs(i) / sbs + yv(1:k) * yv(i) / ys
          end do
        else
          bm = -hn
        end if
      else
        bm = -hn
      end if
      x = xn; f = fn; g = gn; hm = hn
      lam = 0.1_dp * lam
      if (lam < 1.0e-6_dp) lam = 0.0_dp
    end do
    xb = xof(x); fb = f
  contains
    function zof(y) result(z)               ! model parameters -> optimisation coordinates
      real(dp), intent(in) :: y(3)
      real(dp) :: z(3)
      z = y
      where (ls) z = log(y)
    end function zof
    function xof(z) result(y)
      real(dp), intent(in) :: z(3)
      real(dp) :: y(3)
      y = z
      where (ls) y = exp(z)
    end function xof
    subroutine ev(z, fz, gz, hz)            ! log-likelihood, scores and BHHH in z
      real(dp), intent(in) :: z(3)
      real(dp), intent(out) :: fz, gz(3), hz(3, 3)
      real(dp) :: y(3), dz(3)
      y = xof(z)
      call ll_hd(c, k, y, fz, gz, hz)
      dz = merge(y, 1.0_dp, ls)
      gz = gz * dz
      do i = 1, 3
        hz(:, i) = hz(:, i) * dz * dz(i)
      end do
    end subroutine ev
  end subroutine newton_fit

  ! All-n log-likelihood f at x and the Fisher-scoring (BHHH) predicted gain d of a full step from
  ! x (bounds active as in newton_fit); d = -1 if the curvature is not SPD.
  subroutine pred_gain(c, k, x, lb, ub, f, d)
    type(ctx_t), intent(in) :: c
    integer, intent(in) :: k
    real(dp), intent(in) :: x(3), lb(3), ub(3)
    real(dp), intent(out) :: f, d
    real(dp) :: g(3), hm(3, 3), st(3), sc
    logical :: act(3), ok
    integer :: i
    call ll_hd(c, k, x, f, g, hm)
    act = .false.
    do i = 1, k
      sc = 2.0e-7_dp * max(ub(i) - lb(i), 1.0e-8_dp)
      act(i) = (x(i) <= lb(i) + sc .and. g(i) < 0.0_dp) .or. (x(i) >= ub(i) - sc .and. g(i) > 0.0_dp)
    end do
    hm = -hm
    call solve_spd(k, hm, g, act, st, ok)
    d = -1.0_dp
    if (ok) d = 0.5_dp * sum(g(1:k) * st(1:k))
  end subroutine pred_gain

  ! Solve A s = g on the free coordinates (Cholesky); fails if A is not SPD.
  pure subroutine solve_spd(k, a, g, act, s, ok)
    integer, intent(in) :: k
    real(dp), intent(in) :: a(3, 3), g(3)
    logical, intent(in) :: act(3)
    real(dp), intent(out) :: s(3)
    logical, intent(out) :: ok
    real(dp) :: l(3, 3), b(3), z(3)
    integer :: idx(3), m, i, j, q
    m = 0
    do i = 1, k
      if (.not. act(i)) then
        m = m + 1; idx(m) = i
      end if
    end do
    s = 0.0_dp; ok = .true.
    if (m == 0) return
    l = 0.0_dp
    do j = 1, m
      l(j, j) = a(idx(j), idx(j)) - sum(l(j, 1:j - 1)**2)
      if (.not. (l(j, j) > 0.0_dp)) then
        ok = .false.; return
      end if
      l(j, j) = sqrt(l(j, j))
      do i = j + 1, m
        l(i, j) = (a(idx(i), idx(j)) - sum(l(i, 1:j - 1) * l(j, 1:j - 1))) / l(j, j)
      end do
    end do
    do i = 1, m
      b(i) = g(idx(i))
    end do
    do i = 1, m
      z(i) = (b(i) - sum(l(i, 1:i - 1) * z(1:i - 1))) / l(i, i)
    end do
    do i = m, 1, -1
      b(i) = (z(i) - sum(l(i + 1:m, i) * b(i + 1:m))) / l(i, i)
    end do
    do q = 1, m
      s(idx(q)) = b(q)
    end do
  end subroutine solve_spd

  ! ---------------------------------------------------- extended families (generic)
#define VT real(dp)
#define REALMODE
#define EXT_NAME ext_eval_r
#include "families_ext.inc"
#undef EXT_NAME
#undef REALMODE
#undef VT
#define VT real(dp)
#define REALMODE
#define VCONST
#define GEN_NAME gen_eval_rv
#include "families_gen.inc"
#undef GEN_NAME
#undef VCONST
#undef REALMODE
#undef VT

  ! ---------------------------------------------------------------- C API: eval
  integer(c_int) function fv_version() bind(C, name="fv_version")
    fv_version = 1
  end function fv_version

  ! what: 0 logpdf, 1 pdf, 2 hfunc1, 3 hfunc2, 4 hinv1, 5 hinv2.  u is torch (n,2).
  subroutine fv_bicop_eval(fam, rot, par, n, u, what, out, nthreads) bind(C, name="fv_bicop_eval")
    integer(c_int), value :: fam, rot, n, what, nthreads
    real(c_double), intent(in) :: par(3), u(2, n)
    real(c_double), intent(out) :: out(n)
    real(dp) :: axd(1)
    call bicop_eval_x(fam, rot, par, axd, n, u, what, out, nthreads)
  end subroutine fv_bicop_eval

  ! what: 0 logpdf, 1 pdf, 2 hfunc1, 3 hfunc2, 4 hinv1, 5 hinv2; ax = aux vector (fam 24, 25)
  subroutine bicop_eval_x(fam, rot, par, ax, n, u, what, out, nthreads)
    integer, intent(in) :: fam, rot, n, what, nthreads
    real(dp), intent(in) :: par(3), u(2, n), ax(:)
    real(dp), intent(out) :: out(n)
    integer :: i
    real(dp) :: lp, h1, h2, pp(NPP)
    pp = prep_par(fam, par)
    !$omp parallel do private(i, lp, h1, h2) num_threads(nthreads) if(n > 4096) schedule(static)
    do i = 1, n
      select case (what)
      case (0, 1)
        call eval_x(fam, rot, pp, ax, u(1, i), u(2, i), .false., lp, h1, h2)
        if (what == 0) then
          out(i) = lp
        else
          out(i) = max(exp(lp), TINYV)
        end if
      case (2, 3)
        call eval_x(fam, rot, pp, ax, u(1, i), u(2, i), .true., lp, h1, h2)
        if (what == 2) then
          out(i) = h1
        else
          out(i) = h2
        end if
      case default
        out(i) = hinv_x(fam, rot, pp, ax, u(1, i), u(2, i), what - 3)
      end select
    end do
    !$omp end parallel do
  end subroutine bicop_eval_x

  subroutine fv_bicop_pdf_hfuncs(fam, rot, par, n, u, pdf, h1, h2, nthreads) &
      bind(C, name="fv_bicop_pdf_hfuncs")
    integer(c_int), value :: fam, rot, n, nthreads
    real(c_double), intent(in) :: par(3), u(2, n)
    real(c_double), intent(out) :: pdf(n), h1(n), h2(n)
    integer :: i
    real(dp) :: lp, pp(NPP)
    pp = prep_par(fam, par)
    !$omp parallel do private(i, lp) num_threads(nthreads) if(n > 4096) schedule(static)
    do i = 1, n
      call eval_rot(fam, rot, pp, u(1, i), u(2, i), .true., lp, h1(i), h2(i))
      pdf(i) = max(exp(lp), TINYV)
    end do
    !$omp end parallel do
  end subroutine fv_bicop_pdf_hfuncs

  ! ----------------------------------------------------------- log-likelihood
  real(dp) function ll_par(c, pp)
    type(ctx_t), intent(in) :: c
    real(dp), intent(in) :: pp(3)
    integer :: i
    real(dp) :: lp, h1, h2, s, ppx(NPP), q1, q2, axd(1), gd(3), hd3(3, 3), sb(NBK)
    integer :: bk, nt
    if (c%fam == 23 .and. allocated(c%pq) .and. .not. c%disc) then
      call tev_ll(c, 0, pp, s, gd, hd3)
      ll_par = s
      return
    end if
    if (c%lognu) then
      ppx = prep_par(c%fam, [pp(1), exp(pp(2)), pp(3)])
    else
      ppx = prep_par(c%fam, pp)
    end if
    nt = grab(c)
    !$omp parallel do private(bk, i, lp, h1, h2, q1, q2) num_threads(nt) if(nt > 1) schedule(dynamic)
    do bk = 1, NBK
      sb(bk) = 0.0_dp
      do i = bk0(bk, c%n), bk0(bk + 1, c%n) - 1
        if (c%disc) then
          call eval_disc(c%fam, c%rot, ppx, axd, c%a0(i), c%b0(i), c%am(i), c%bm(i), c%da, c%db, lp, h1, h2, &
                         q1, q2, need_h=.false.)
        else
          call eval0(c%fam, ppx, c%a(i), c%b(i), .false., lp, h1, h2)
        end if
        ! an invalid/zero density is scored as log(tiny) -- the same floor the vine
        ! evaluator clamps to -- so invalid parameter regions are never rewarded
        if (.not. finite_(lp)) lp = LOGTINY
        if (c%usew) then
          sb(bk) = sb(bk) + c%w(i) * lp
        else
          sb(bk) = sb(bk) + lp
        end if
      end do
    end do
    !$omp end parallel do
    call release(c, nt)
    ll_par = sum(sb)
  end function ll_par

  ! Student loglik for cached t-scores c%x, c%y at (rho, nu).
  real(dp) function ll_student_xy(c, rho, nu)
    type(ctx_t), intent(in) :: c
    real(dp), intent(in) :: rho, nu
    real(dp) :: k0, s, lp, r, tb(NBK)
    integer :: i, bk, nt
    r = min(max(rho, -0.999999_dp), 0.999999_dp)
    s = 1.0_dp - r * r
    k0 = log_gamma(0.5_dp * (nu + 2.0_dp)) + log_gamma(0.5_dp * nu) &
         - 2.0_dp * log_gamma(0.5_dp * (nu + 1.0_dp)) - 0.5_dp * log(s)
    nt = grab(c)
    !$omp parallel do private(bk, i, lp) num_threads(nt) if(nt > 1) schedule(dynamic)
    do bk = 1, NBK
    tb(bk) = 0.0_dp
    do i = bk0(bk, c%n), bk0(bk + 1, c%n) - 1
      lp = k0 + 0.5_dp * (nu + 1.0_dp) * (log1p_(c%x(i)**2 / nu) + log1p_(c%y(i)**2 / nu)) &
           - 0.5_dp * (nu + 2.0_dp) * log1p_((c%x(i)**2 - 2.0_dp * r * c%x(i) * c%y(i) + c%y(i)**2) / (nu * s))
      if (.not. finite_(lp)) lp = LOGTINY
      if (c%usew) then
        tb(bk) = tb(bk) + c%w(i) * lp
      else
        tb(bk) = tb(bk) + lp
      end if
    end do
    end do
    !$omp end parallel do
    call release(c, nt)
    ll_student_xy = sum(tb)
  end function ll_student_xy

  subroutine student_scores(c, nu)
    type(ctx_t), intent(inout) :: c
    real(dp), intent(in) :: nu
    integer :: i, nt
    real(dp) :: kd, lb
    kd = kdens_t(nu); lb = lbeta_t(nu)
    nt = grab(c)
    !$omp parallel do private(i) num_threads(nt) if(nt > 1) schedule(static)
    do i = 1, c%n
      c%x(i) = qt_l(c%a(i), nu, kd, lb)
      c%y(i) = qt_l(c%b(i), nu, kd, lb)
    end do
    !$omp end parallel do
    call release(c, nt)
  end subroutine student_scores

  ! 1-D objective (minimised).  mode 1: one-parameter family; 2: student rho
  ! with cached scores (nu = c%fixed); 3: student nu, rho profiled out;
  ! 4: student nu with rho fixed (= c%fixed).
  recursive real(dp) function f1(c, x) result(f)
    type(ctx_t), intent(inout) :: c
    real(dp), intent(in) :: x
    type(ctx_t) :: ci
    real(dp) :: rb, fb, rp(3)
    select case (c%mode)
    case (1)
      f = -ll_par(c, [x, 0.0_dp, 0.0_dp])
    case (2)
      f = -ll_student_xy(c, x, c%fixed)
    case (3)
      call student_scores(c, x)
      ! inner rho: Newton with the analytic score, warm-started at the last rho
      rb = student_rho_newton(c, x, c%lastrho, c%lbr, c%ubr)
      c%lastrho = rb
      f = -ll_student_xy(c, rb, x)
    case (4)
      call student_scores(c, x)
      f = -ll_student_xy(c, c%fixed, x)
    case (5)
      ! t-EV, x = log(nu): inner rho by projected Newton on the exact (dual-number) score
      call newton_fit(c, 1, [c%lastrho, exp(x), 0.0_dp], [c%lbr, exp(x), 0.0_dp], [c%ubr, exp(x), 0.0_dp], rp, fb)
      c%lastrho = rp(1)
      f = -fb
    case default
      f = 0.0_dp
    end select
  end function f1

  ! MLE of rho for fixed nu given cached t-scores (c%x, c%y): Newton on the
  ! analytic score dl/drho = rho/s - (nu+2)/2 Q'/(nu+Q), BHHH curvature,
  ! bisection-safe within [lo, hi].  Typically 2-4 cheap passes.
  real(dp) function student_rho_newton(c, nu, r0, lo, hi) result(r)
    type(ctx_t), intent(in) :: c
    real(dp), intent(in) :: nu, r0, lo, hi
    real(dp) :: g, h, s, nq, q, dq, gi, a, b, rn, x, y, gb(NBK), hb(NBK)
    integer :: it, i, bk, nt
    a = lo; b = hi
    r = min(max(r0, a + 1.0e-9_dp), b - 1.0e-9_dp)
    do it = 1, 30
      s = (1.0_dp - r) * (1.0_dp + r)
      nt = grab(c)
      !$omp parallel do private(bk, i, x, y, nq, q, dq, gi) num_threads(nt) if(nt > 1) schedule(dynamic)
      do bk = 1, NBK
      gb(bk) = 0.0_dp; hb(bk) = 0.0_dp
      do i = bk0(bk, c%n), bk0(bk + 1, c%n) - 1
        x = c%x(i); y = c%y(i)
        nq = x * x - 2.0_dp * r * x * y + y * y
        q = nq / s
        dq = (-2.0_dp * x * y * s + 2.0_dp * r * nq) / (s * s)
        gi = r / s - 0.5_dp * (nu + 2.0_dp) * dq / (nu + q)
        if (c%usew) gi = c%w(i) * gi
        gb(bk) = gb(bk) + gi
        hb(bk) = hb(bk) + gi * gi
      end do
      end do
      !$omp end parallel do
      call release(c, nt)
      g = sum(gb); h = sum(hb)
      if (g > 0.0_dp) then
        a = r
      else
        b = r
      end if
      if (h <= 0.0_dp) exit
      rn = r + g / h
      if (.not. (rn > a .and. rn < b)) rn = 0.5_dp * (a + b)
      if (abs(rn - r) <= 1.0e-10_dp) then
        r = rn; exit
      end if
      r = rn
    end do
  end function student_rho_newton

  ! Brent (1973) localmin on [ax, bx].
  recursive subroutine brent(c, ax, bx, tol, xmin, fmin)
    type(ctx_t), intent(inout) :: c
    real(dp), intent(in) :: ax, bx, tol
    real(dp), intent(out) :: xmin, fmin
    real(dp), parameter :: CG = 0.3819660112501051_dp
    real(dp) :: a, b, d, e, p, q, r, u, v, w, x, fu, fv, fw, fx, xm, tol1, tol2
    integer :: it
    a = min(ax, bx); b = max(ax, bx)
    v = a + CG * (b - a); w = v; x = v
    e = 0.0_dp; d = 0.0_dp
    fx = f1(c, x); fv = fx; fw = fx
    do it = 1, 100
      xm = 0.5_dp * (a + b)
      tol1 = tol * abs(x) + 1.0e-10_dp
      tol2 = 2.0_dp * tol1
      if (abs(x - xm) <= tol2 - 0.5_dp * (b - a)) exit
      if (abs(e) > tol1) then
        r = (x - w) * (fx - fv)
        q = (x - v) * (fx - fw)
        p = (x - v) * q - (x - w) * r
        q = 2.0_dp * (q - r)
        if (q > 0.0_dp) p = -p
        q = abs(q)
        r = e; e = d
        if (abs(p) >= abs(0.5_dp * q * r) .or. p <= q * (a - x) .or. p >= q * (b - x)) then
          if (x >= xm) then
            e = a - x
          else
            e = b - x
          end if
          d = CG * e
        else
          d = p / q
          u = x + d
          if (u - a < tol2 .or. b - u < tol2) d = sign(tol1, xm - x)
        end if
      else
        if (x >= xm) then
          e = a - x
        else
          e = b - x
        end if
        d = CG * e
      end if
      if (abs(d) >= tol1) then
        u = x + d
      else
        u = x + sign(tol1, d)
      end if
      fu = f1(c, u)
      if (fu <= fx) then
        if (u >= x) then
          a = x
        else
          b = x
        end if
        v = w; fv = fw; w = x; fw = fx; x = u; fx = fu
      else
        if (u < x) then
          a = u
        else
          b = u
        end if
        if (fu <= fw .or. w == x) then
          v = w; fv = fw; w = u; fw = fu
        else if (fu <= fv .or. v == x .or. v == w) then
          v = u; fv = fu
        end if
      end if
    end do
    xmin = x; fmin = fx
  end subroutine brent

  ! Nelder-Mead (k = 2 or 3) in a logistic box transform, coarse multistart first.
  subroutine nmk(c, k, x0, lb, ub, xb, fb, polish, maxev)
    type(ctx_t), intent(inout) :: c
    integer, intent(in) :: k
    real(dp), intent(in) :: x0(3), lb(3), ub(3)
    real(dp), intent(out) :: xb(3), fb
    logical, intent(in), optional :: polish     ! local polish: no grid, small simplex, few evaluations
    integer, intent(in), optional :: maxev      ! evaluation budget per restart
    real(dp) :: s(3, 4), f(4), z(3), zc(3), zr(3), ze(3), zk(3), fr, fe, fk, ft, fz, step
    real(dp), parameter :: grid(3) = [-2.0_dp, 0.0_dp, 2.0_dp]
    integer :: i, j, m, q, nev, rs, ilo, ihi, inx, g(3), ng, mev
    logical :: pol
    z = 0.0_dp
    do i = 1, k
      z(i) = logit_box(x0(i), lb(i), ub(i))
    end do
    nev = 0
    fz = fbox(z)
    ng = 3**k
    pol = .false.
    if (present(polish)) pol = polish
    if (pol) ng = 0
    do q = 0, ng - 1
      m = q
      do i = 1, k
        g(i) = mod(m, 3) + 1; m = m / 3
      end do
      zc = 0.0_dp
      do i = 1, k
        zc(i) = grid(g(i))
      end do
      ft = fbox(zc)
      if (ft < fz) then
        fz = ft; z = zc
      end if
    end do
    step = merge(0.1_dp, 0.8_dp, pol)
    do rs = 1, merge(1, 2, pol)
      s(:, 1) = z; f(1) = fz
      do j = 2, k + 1
        s(:, j) = z; s(j - 1, j) = z(j - 1) + step
        f(j) = fbox(s(:, j))
      end do
      mev = merge(5 * k + 10, 250 * rs, pol)
      if (present(maxev)) mev = min(mev, maxev * rs)
      do while (nev < mev)
        ilo = minloc(f(1:k + 1), 1); ihi = maxloc(f(1:k + 1), 1)
        inx = ilo
        do j = 1, k + 1
          if (j /= ihi .and. f(j) >= f(inx)) inx = j
        end do
        if (f(ihi) - f(ilo) <= 1.0e-8_dp * (1.0_dp + abs(f(ilo)))) exit
        zc = 0.0_dp
        do j = 1, k + 1
          if (j /= ihi) zc = zc + s(:, j)
        end do
        zc = zc / k
        zr = zc + (zc - s(:, ihi)); fr = fbox(zr)
        if (fr < f(ilo)) then
          ze = zc + 2.0_dp * (zc - s(:, ihi)); fe = fbox(ze)
          if (fe < fr) then
            s(:, ihi) = ze; f(ihi) = fe
          else
            s(:, ihi) = zr; f(ihi) = fr
          end if
        else if (fr < f(inx)) then
          s(:, ihi) = zr; f(ihi) = fr
        else
          if (fr < f(ihi)) then
            zk = zc + 0.5_dp * (zr - zc)
          else
            zk = zc + 0.5_dp * (s(:, ihi) - zc)
          end if
          fk = fbox(zk)
          if (fk < min(fr, f(ihi))) then
            s(:, ihi) = zk; f(ihi) = fk
          else
            do j = 1, k + 1
              if (j /= ilo) then
                s(:, j) = s(:, ilo) + 0.5_dp * (s(:, j) - s(:, ilo))
                f(j) = fbox(s(:, j))
              end if
            end do
          end if
        end if
      end do
      ilo = minloc(f(1:k + 1), 1)
      z = s(:, ilo); fz = f(ilo)
      step = 0.15_dp     ! restart with a small simplex to confirm convergence
    end do
    fb = fz
    xb = 0.0_dp
    do i = 1, k
      xb(i) = lb(i) + (ub(i) - lb(i)) / (1.0_dp + exp(-z(i)))
    end do
  contains
    real(dp) function fbox(zz)
      real(dp), intent(in) :: zz(3)
      real(dp) :: p(3)
      integer :: mm
      p = 0.0_dp
      do mm = 1, k
        p(mm) = lb(mm) + (ub(mm) - lb(mm)) / (1.0_dp + exp(-zz(mm)))
      end do
      fbox = -ll_par(c, p)
      if (.not. finite_(fbox)) fbox = huge(1.0_dp)
      nev = nev + 1
    end function fbox
  end subroutine nmk

  ! Start values for the extended families: best of a per-family grid (log-spaced
  ! for scale parameters), value-only likelihood evaluations.
  subroutine ext_grid_start(c, fam, lb, ub, xs, fs)
    type(ctx_t), intent(inout) :: c
    integer, intent(in) :: fam
    real(dp), intent(in) :: lb(3), ub(3)
    real(dp), intent(out) :: xs(3, MS), fs(MS)
    real(dp) :: g1(8), g2(4), p(3)
    integer :: n1, n2, i, j
    n2 = 1; g2(1) = 0.0_dp
    select case (fam)
    case (12); n1 = 8; g1 = [0.05_dp, 0.2_dp, 0.5_dp, 0.8_dp, 1.5_dp, 4.0_dp, 15.0_dp, 80.0_dp]
    case (13); n1 = 5; g1(1:5) = [-0.8_dp, -0.4_dp, 0.0_dp, 0.4_dp, 0.8_dp]
    case (14); n1 = 5; g1(1:5) = [-0.8_dp, -0.3_dp, 0.1_dp, 0.5_dp, 0.85_dp]
    case (15); n1 = 6; g1(1:6) = [0.1_dp, 0.3_dp, 0.7_dp, 1.3_dp, 2.5_dp, 6.0_dp]
    case (16); n1 = 6; g1(1:6) = [0.2_dp, 0.5_dp, 1.0_dp, 1.8_dp, 3.5_dp, 8.0_dp]
    case (17); n1 = 4; g1(1:4) = [0.1_dp, 0.5_dp, 1.5_dp, 5.0_dp]; n2 = 4; g2 = [0.1_dp, 0.5_dp, 1.5_dp, 5.0_dp]
    case (18, 19, 20); n1 = 4; g1(1:4) = [1.1_dp, 1.5_dp, 2.5_dp, 5.0_dp]; n2 = 4; g2 = [0.1_dp, 0.5_dp, 1.5_dp, 5.0_dp]
    case (21); n1 = 4; g1(1:4) = [0.2_dp, 0.8_dp, 2.0_dp, 6.0_dp]; n2 = 4; g2 = [0.1_dp, 0.4_dp, 0.75_dp, 0.98_dp]
    case (22); n1 = 4; g1(1:4) = [0.1_dp, 0.5_dp, 1.5_dp, 4.0_dp]; n2 = 4; g2 = [0.1_dp, 0.5_dp, 1.5_dp, 5.0_dp]
    case (23); n1 = 4; g1(1:4) = [-0.3_dp, 0.2_dp, 0.6_dp, 0.9_dp]; n2 = 4; g2 = [1.5_dp, 4.0_dp, 10.0_dp, 25.0_dp]
    case default; n1 = 1; g1(1) = 0.5_dp * (lb(1) + ub(1))
    end select
    fs = -huge(1.0_dp); xs = 0.0_dp
    do i = 1, n1
      do j = 1, n2
        p = 0.0_dp
        p(1) = min(max(g1(i), lb(1)), ub(1))
        if (npar_of(fam) >= 2) p(2) = min(max(g2(j), lb(2)), ub(2))
        call keep_top(p, ll_par(c, p), xs, fs)
      end do
    end do
  end subroutine ext_grid_start

  ! Start values for the generic families: 4 values per parameter (4^k grid, k <= 3), log-spaced
  ! in wide positive boxes, linear in unit boxes (where the upper end, e.g. psi = 1, is included).
  ! The grid is ranked by Blomqvist beta (one closed-form C(1/2, 1/2) per point, no likelihood
  ! pass); only the MS points closest to the data beta are scored by the likelihood.
  subroutine gen_grid_start(c, fam, lb, ub, xs, fs)
    type(ctx_t), intent(inout) :: c
    integer, intent(in) :: fam
    real(dp), intent(in) :: lb(3), ub(3)
    real(dp), intent(out) :: xs(3, MS), fs(MS)
    real(dp), parameter :: RW(4) = [0.1_dp, 0.35_dp, 0.6_dp, 0.85_dp], RU(4) = [0.2_dp, 0.5_dp, 0.8_dp, 1.0_dp]
    real(dp) :: p(3), bd, bm, xt(3, MS), ft(MS), lp, l1, l2, lc, sw
    integer :: k, q, m, i, j
    k = gen_npar(fam)
    fs = -huge(1.0_dp); xs = 0.0_dp; ft = -huge(1.0_dp); xt = 0.0_dp
    if (c%usew) then
      sw = sum(c%w(1:c%n))
      bd = 4.0_dp * sum(c%w(1:c%n), mask=c%a(1:c%n) <= 0.5_dp .and. c%b(1:c%n) <= 0.5_dp) / max(sw, TINYV) - 1.0_dp
    else
      bd = 4.0_dp * count(c%a(1:c%n) <= 0.5_dp .and. c%b(1:c%n) <= 0.5_dp) / real(max(c%n, 1), dp) - 1.0_dp
    end if
    do q = 0, 4**k - 1
      m = q; p = 0.0_dp
      do i = 1, k
        j = mod(m, 4) + 1; m = m / 4
        if (ub(i) <= 1.0_dp) then
          p(i) = lb(i) + (ub(i) - lb(i)) * RU(j)
        else if (lb(i) > 0.0_dp .and. ub(i) >= 20.0_dp * lb(i)) then
          p(i) = lb(i) * (ub(i) / lb(i))**RW(j)
        else
          p(i) = lb(i) + (ub(i) - lb(i)) * RW(j)
        end if
      end do
      call gen_eval_rv(fam, p, 0.5_dp, 0.5_dp, lp, l1, l2, lc)
      bm = 4.0_dp * exp(lc) - 1.0_dp
      if (bm == bm) call keep_top(p, -abs(bm - bd), xt, ft)
    end do
    do i = 1, MS
      if (ft(i) > -huge(1.0_dp)) call keep_top(xt(:, i), ll_par(c, xt(:, i)), xs, fs)
    end do
  end subroutine gen_grid_start

  ! Pre-screen of a generic candidate in vine selection: Blomqvist-beta start and one Newton
  ! fit on a stride subsample (~NSUB rows); returns the subsample log-likelihood scaled to n.
  ! Only the best few candidates per edge go on to the full-data racing and polish.
  subroutine gen_prescreen(n, u1, u2, fam, rot, lb, ub, usew, w, llsc)
    integer, intent(in) :: n, fam, rot
    real(dp), intent(in) :: u1(n), u2(n), lb(3), ub(3), w(:)
    logical, intent(in) :: usew
    real(dp), intent(out) :: llsc
    real(dp), allocatable, target :: a(:), b(:), ww(:)
    type(ctx_t) :: cs
    real(dp) :: xs0(3, MS), fs0(MS), pj(3), fj
    integer :: st, m, i
    st = max(1, n / NSUB); m = (n - 1) / st + 1     ! ~NSUB rows (NSUB/2 ranks worse: DE110020 CRPS +0.4 %)
    allocate(a(m), b(m), ww(m))
    do i = 1, m
      call rotate(rot, clampu(u1(1 + (i - 1) * st)), clampu(u2(1 + (i - 1) * st)), a(i), b(i))
    end do
    ww = 1.0_dp
    if (usew) ww = w(1:n:st)
    cs%fam = fam; cs%n = m; cs%usew = usew; cs%nin = 1
    cs%a => a; cs%b => b; cs%w => ww; cs%x => a; cs%y => b
    cs%tol = 1.0e-3_dp; cs%qn = .true.          ! a ranking only: the kept candidates are fitted exactly
    call gen_grid_start(cs, fam, lb, ub, xs0, fs0)
    llsc = -huge(1.0_dp)
    if (fs0(1) <= -huge(1.0_dp)) return
    call newton_fit(cs, gen_npar(fam), xs0(:, 1), lb, ub, pj, fj)
    if (finite_(fj)) llsc = fj * real(n, dp) / real(m, dp)
  end subroutine gen_prescreen

  ! Coarse logistic-grid start for multi-parameter families (value-only evals).
  subroutine grid_start(c, k, x0, lb, ub, xs, fs)
    type(ctx_t), intent(inout) :: c
    integer, intent(in) :: k
    real(dp), intent(in) :: x0(3), lb(3), ub(3)
    real(dp), intent(out) :: xs(3, MS), fs(MS)
    real(dp), parameter :: grid(3) = [-2.0_dp, 0.0_dp, 2.0_dp]
    real(dp) :: p(3)
    integer :: q, m, i
    fs = -huge(1.0_dp); xs = 0.0_dp
    call keep_top(x0, ll_par(c, x0), xs, fs)
    do q = 0, 3**k - 1
      m = q
      p = 0.0_dp
      do i = 1, k
        p(i) = lb(i) + (ub(i) - lb(i)) / (1.0_dp + exp(-grid(mod(m, 3) + 1)))
        m = m / 3
      end do
      call keep_top(p, ll_par(c, p), xs, fs)
    end do
  end subroutine grid_start

  ! insert (p, f) into the descending top-MS list
  pure subroutine keep_top(p, f, xs, fs)
    real(dp), intent(in) :: p(3), f
    real(dp), intent(inout) :: xs(3, MS), fs(MS)
    integer :: j, k
    if (.not. finite_(f)) return
    do j = 1, MS
      if (f > fs(j)) then
        do k = MS, j + 1, -1
          fs(k) = fs(k - 1); xs(:, k) = xs(:, k - 1)
        end do
        fs(j) = f; xs(:, j) = p
        return
      end if
    end do
  end subroutine keep_top

  pure real(dp) function logit_box(x, lb, ub)
    real(dp), intent(in) :: x, lb, ub
    real(dp) :: fr
    fr = (x - lb) / max(ub - lb, TINYV)
    fr = min(max(fr, 1.0e-6_dp), 1.0_dp - 1.0e-6_dp)
    logit_box = log(fr / (1.0_dp - fr))
  end function logit_box

  ! ----------------------------------------------------------- C API: fitting
  ! njob independent (dataset, family, rotation) fits.  U is (nds, n, 2) in
  ! torch order -> Fortran U(2, n, nds).  method: 0 = evaluate loglik at x0,
  ! 1 = MLE within [lb, ub], 2 = student with rho fixed at x0(1), nu by MLE.
  subroutine fv_fit_jobs(njob, n, nds, u, ds, fam, rot, method, x0, lb, ub, usew, w, &
                          par, ll, nthreads) bind(C, name="fv_fit_jobs")
    integer(c_int), value :: njob, n, nds, usew, nthreads
    real(c_double), intent(in) :: u(2, n, nds), w(n), x0(3, njob), lb(3, njob), ub(3, njob)
    integer(c_int), intent(in) :: ds(njob), fam(njob), rot(njob), method(njob)
    real(c_double), intent(out) :: par(3, njob), ll(njob)
    integer :: j
    !$ call omp_set_max_active_levels(2)
    !$omp parallel do private(j) num_threads(max(1, min(nthreads, njob))) schedule(dynamic, 1)
    do j = 1, njob
      call fit_one(n, u(:, :, ds(j) + 1), fam(j), rot(j), method(j), x0(:, j), lb(:, j), ub(:, j), &
                   usew /= 0, w, par(:, j), ll(j), max(1, nthreads / max(1, min(nthreads, njob))))
    end do
    !$omp end parallel do
  end subroutine fv_fit_jobs

  ! ---------------------------------------------------------- pair-fit cache
  ! fv_fitcache_enable(capacity): turn the cache on with room for capacity fits (a
  ! hash table; when full, colliding entries are replaced).  Keys are two independent
  ! 32-bit FNV-style hashes of every bit of the pair data (+ weights, F(x-) of discrete
  ! margins) together with family, rotation, method, start and box: a hit returns the
  ! very fit that was computed for identical inputs.  Thread-safe (OpenMP critical).
  subroutine fv_fitcache_enable(capacity)
    integer, intent(in) :: capacity
    !$omp critical (fv_fitcache)
    if (allocated(fc_tab)) deallocate(fc_tab)
    fc_cap = max(16, capacity)
    allocate(fc_tab(fc_cap))
    fc_hits = 0; fc_miss = 0
    !$omp end critical (fv_fitcache)
  end subroutine fv_fitcache_enable

  integer function fv_fitcache_capacity()
    fv_fitcache_capacity = fc_cap
  end function fv_fitcache_capacity

  subroutine fv_fitcache_clear()
    integer :: i
    !$omp critical (fv_fitcache)
    if (allocated(fc_tab)) then
      do i = 1, fc_cap
        fc_tab(i)%used = .false.
        if (allocated(fc_tab(i)%ax)) deallocate(fc_tab(i)%ax)
      end do
    end if
    fc_hits = 0; fc_miss = 0
    !$omp end critical (fv_fitcache)
  end subroutine fv_fitcache_clear

  subroutine fv_fitcache_disable()
    !$omp critical (fv_fitcache)
    if (allocated(fc_tab)) deallocate(fc_tab)
    fc_cap = 0
    !$omp end critical (fv_fitcache)
  end subroutine fv_fitcache_disable

  subroutine fv_fitcache_stats(hits, misses, used)
    integer(int64), intent(out) :: hits, misses
    integer, intent(out), optional :: used
    hits = fc_hits; misses = fc_miss
    if (present(used)) then
      used = 0
      if (allocated(fc_tab)) used = count(fc_tab(:)%used)
    end if
  end subroutine fv_fitcache_stats

  ! fold the bits of x(:) into the two hash lanes (each < 2^32; no integer overflow)
  pure subroutine fc_fold(m, x, h1, h2)
    integer, intent(in) :: m
    real(dp), intent(in) :: x(m)
    integer(int64), intent(inout) :: h1, h2
    integer(int64), parameter :: M32 = 4294967295_int64, P1 = 16777619_int64, P2 = 1099511_int64
    integer(int64) :: b, lo, hi
    integer :: i
    do i = 1, m
      b = transfer(x(i), 0_int64)
      lo = iand(b, M32); hi = iand(shiftr(b, 32), M32)
      h1 = iand(ieor(h1, lo) * P1, M32); h1 = iand(ieor(h1, hi) * P1, M32)
      h2 = iand(ieor(h2, hi) * P2 + 7_int64, M32); h2 = iand(ieor(h2, lo) * P2 + 7_int64, M32)
    end do
  end subroutine fc_fold

  ! ================================================================== analytical one-parameter fits
  ! Gaussian (1), Clayton (3), Gumbel (4), Frank (5) and Joe (6).  Every quantity is in closed form:
  ! the log-density, both h-functions and the copula CDF together with their first two derivatives in
  ! the parameter (kern_lp, kern_h, kern_c), so one pass over the data gives the log-likelihood, its
  ! score and its curvature -- for continuous pairs and for pairs with discrete margins alike.
  ! Estimator: Le Cam's one-step from the Kendall's-tau inversion (newton1); the Gaussian on a
  ! continuous pair is solved exactly (Cardano).  Kernel arguments t(4) per family:
  ! 1 (qnorm a, qnorm b) | 3 (log a, log b, a, b) | 4 (-log a, -log b, log(-log a), log(-log b)) |
  ! 5 (a, b) | 6 (log(1-a), log(1-b), a, b), all of the rotated, clamped pair.

  ! transforms of a pair shared by every family and rotation (continuous edges): 1 log a, 2 log b,
  ! 3 log(1-a), 4 log(1-b), 5 qnorm a, 6 qnorm b, 7 log(-log a), 8 log(-log b), 9 log(-log(1-a)),
  ! 10 log(-log(1-b))
  subroutine pc_transforms(n, u1, u2, t)
    integer, intent(in) :: n
    real(dp), intent(in) :: u1(n), u2(n)
    real(dp), intent(out) :: t(n, 10)
    real(dp) :: a, b
    integer :: i
    do i = 1, n
      a = clampu(u1(i)); b = clampu(u2(i))
      t(i, 1) = log(a); t(i, 2) = log(b); t(i, 3) = log1p_(-a); t(i, 4) = log1p_(-b)
      t(i, 5) = qnorm(a); t(i, 6) = qnorm(b)
      t(i, 7) = log(-t(i, 1)); t(i, 8) = log(-t(i, 2)); t(i, 9) = log(-t(i, 3)); t(i, 10) = log(-t(i, 4))
    end do
  end subroutine pc_transforms

  ! kernel arguments t(4) of family fam at a rotated, clamped point (a, b)
  pure subroutine kern_t(fam, a, b, t)
    integer, intent(in) :: fam
    real(dp), intent(in) :: a, b
    real(dp), intent(out) :: t(4)
    t = 0.0_dp
    select case (fam)
    case (1); t(1) = qnorm(a); t(2) = qnorm(b)
    case (3); t = [log(a), log(b), a, b]
    case (4); t(1) = -log(a); t(2) = -log(b); t(3) = log(t(1)); t(4) = log(t(2))
    case (5); t(1) = a; t(2) = b
    case (6); t = [log1p_(-a), log1p_(-b), a, b]
    end select
  end subroutine kern_t

  ! log-density lp and its theta-derivatives d1, d2 (eval0's density; flat in the clamped regions)
  pure subroutine kern_lp(fam, th0, t, lp, d1, d2)
    integer, intent(in) :: fam
    real(dp), intent(in) :: th0, t(4)
    real(dp), intent(out) :: lp, d1, d2
    real(dp) :: r, s, aa, bb, nn, n1, q, th, e, ls, r1, c1, it, it2, it3, t1, t2, m, mn, lm, ln, &
                l1, l2, gg, g1, g2, sst, s1, s2, den, eu, ev, ed, v, eu1, ev1, eu2, ev2, ed1, ed2, dn1, dn2, ua, ub, &
                ua1, ub1, ua2, ub2
    lp = 0.0_dp; d1 = 0.0_dp; d2 = 0.0_dp
    select case (fam)
    case (1)
      r = min(max(th0, -0.999999_dp), 0.999999_dp); s = 1.0_dp - r * r
      aa = t(1) * t(1) + t(2) * t(2); bb = t(1) * t(2)
      nn = r * r * aa - 2.0_dp * r * bb; n1 = 2.0_dp * r * aa - 2.0_dp * bb; q = n1 * s + 2.0_dp * r * nn
      lp = -0.5_dp * log(s) - nn / (2.0_dp * s)
      if (abs(th0) <= 0.999999_dp) then
        d1 = r / s - q / (2.0_dp * s * s)
        d2 = (s + 2.0_dp * r * r) / (s * s) - ((2.0_dp * aa * s + 2.0_dp * nn) * s * s + 4.0_dp * r * s * q) / (2.0_dp * s**4)
      end if
    case (3)
      th = max(th0, 0.0_dp); call clayton_l(th, -t(1), -t(2), l1, l2, r1)
      c1 = 1.0_dp / (1.0_dp + th)
      lp = log1p_(th) - (1.0_dp + th) * (t(1) + t(2)) - (1.0_dp + 2.0_dp * th) * l1
      d1 = c1 - (t(1) + t(2)) - 2.0_dp * l1 - (1.0_dp + 2.0_dp * th) * l2
      d2 = -c1 * c1 - 4.0_dp * l2 - (1.0_dp + 2.0_dp * th) * r1
    case (4)
      th = max(th0, 1.0_dp); it = 1.0_dp / th; it2 = it * it; it3 = it2 * it
      t1 = th * t(3); t2 = th * t(4); m = max(t1, t2); mn = min(t1, t2)
      e = exp(mn - m); ls = m + log1p_(e)
      lm = merge(t(3), t(4), t1 >= t2); ln = merge(t(4), t(3), t1 >= t2); q = 1.0_dp / (1.0_dp + e)
      l1 = (lm + ln * e) * q; l2 = (lm * lm + ln * ln * e) * q - l1 * l1
      gg = ls * it; g1 = l1 * it - ls * it2; g2 = l2 * it - 2.0_dp * l1 * it2 + 2.0_dp * ls * it3
      sst = exp(gg); s1 = sst * g1; s2 = sst * (g2 + g1 * g1)
      den = th - 1.0_dp + sst; q = 1.0_dp / den
      lp = -sst + (th - 1.0_dp) * (t(3) + t(4)) + t(1) + t(2) + (it - 2.0_dp) * ls + log(den)
      if (th0 >= 1.0_dp) then
        d1 = -s1 + (t(3) + t(4)) - ls * it2 + (it - 2.0_dp) * l1 + (1.0_dp + s1) * q
        d2 = -s2 + 2.0_dp * ls * it3 - 2.0_dp * l1 * it2 + (it - 2.0_dp) * l2 + s2 * q - ((1.0_dp + s1) * q)**2
      end if
    case (5)
      th = th0
      if (abs(th) < 1.0e-10_dp) return
      call frank_terms(th, t(1), t(2), eu, ev, ed, eu1, ev1, eu2, ev2, ed1, ed2, den, dn1, dn2)
      v = -th * ed * (eu + 1.0_dp) * (ev + 1.0_dp); q = 1.0_dp / den
      lp = log(max(v, TINYV) * q * q)
      d1 = -2.0_dp * dn1 * q; d2 = -2.0_dp * (dn2 * q - (dn1 * q)**2)
      if (v > TINYV) then
        d1 = d1 + 1.0_dp / th + ed1 / ed - (t(1) + t(2))
        d2 = d2 - 1.0_dp / (th * th) + ed2 / ed - (ed1 / ed)**2
      end if
    case (6)
      th = max(th0, 1.0_dp); it = 1.0_dp / th; it2 = it * it; it3 = it2 * it
      call joe_terms(th, t(1), t(2), ua, ub, ua1, ub1, ua2, ub2, s, s1, s2)
      ls = log(s); den = th - 1.0_dp + s; q = 1.0_dp / s; r = 1.0_dp / den
      lp = (it - 2.0_dp) * ls + (th - 1.0_dp) * (t(1) + t(2)) + log(den)
      if (th0 >= 1.0_dp) then
        d1 = -ls * it2 + (it - 2.0_dp) * s1 * q + (t(1) + t(2)) + (1.0_dp + s1) * r
        d2 = 2.0_dp * ls * it3 - 2.0_dp * s1 * it2 * q + (it - 2.0_dp) * (s2 * q - (s1 * q)**2) + s2 * r - ((1.0_dp + s1) * r)**2
      end if
    end select
  end subroutine kern_lp

  ! unrotated h-functions g(1) = h1 = P(V <= b | U = a), g(2) = h2 = P(U <= a | V = b) (eval0's) and
  ! their theta-derivatives gd, gdd
  pure subroutine kern_h(fam, th0, t, g, gd, gdd)
    integer, intent(in) :: fam
    real(dp), intent(in) :: th0, t(4)
    real(dp), intent(out) :: g(2), gd(2), gdd(2)
    real(dp) :: r, s, z(2), zd(2), zdd(2), th, e, ls, r1, k1, k2, lt(2), it, it2, it3, t1, t2, m, mn, lm, &
                ln, q, l1, l2, gg, g1, g2, sst, s1, s2, eu, ev, ed, eu1, ev1, eu2, ev2, ed1, ed2, den, dn1, dn2, nn(2), &
                n1(2), n2(2), ua, ub, ua1, ub1, ua2, ub2, uu(2), ud(2), udd(2)
    integer :: j
    g = 0.0_dp; gd = 0.0_dp; gdd = 0.0_dp
    select case (fam)
    case (1)
      r = min(max(th0, -0.999999_dp), 0.999999_dp); s = 1.0_dp - r * r
      z(1) = (t(2) - r * t(1)) / sqrt(s); z(2) = (t(1) - r * t(2)) / sqrt(s)
      zd(1) = (r * t(2) - t(1)) / s**1.5_dp; zd(2) = (r * t(1) - t(2)) / s**1.5_dp
      zdd(1) = t(2) / s**1.5_dp + 3.0_dp * r * (r * t(2) - t(1)) / s**2.5_dp
      zdd(2) = t(1) / s**1.5_dp + 3.0_dp * r * (r * t(1) - t(2)) / s**2.5_dp
      do j = 1, 2
        g(j) = pnorm(z(j))
        if (abs(th0) <= 0.999999_dp) then
          gd(j) = dnorm(z(j)) * zd(j); gdd(j) = dnorm(z(j)) * (zdd(j) - z(j) * zd(j)**2)
        end if
      end do
    case (3)
      th = max(th0, 0.0_dp); call clayton_l(th, -t(1), -t(2), l1, l2, r1)
      lt = [t(1), t(2)]
      do j = 1, 2
        k1 = -lt(j) - l1 - (1.0_dp + th) * l2
        k2 = -2.0_dp * l2 - (1.0_dp + th) * r1
        g(j) = exp(-(1.0_dp + th) * (lt(j) + l1)); gd(j) = g(j) * k1; gdd(j) = g(j) * (k2 + k1 * k1)
      end do
    case (4)
      th = max(th0, 1.0_dp); it = 1.0_dp / th; it2 = it * it; it3 = it2 * it
      t1 = th * t(3); t2 = th * t(4); m = max(t1, t2); mn = min(t1, t2)
      e = exp(mn - m); ls = m + log1p_(e)
      lm = merge(t(3), t(4), t1 >= t2); ln = merge(t(4), t(3), t1 >= t2); q = 1.0_dp / (1.0_dp + e)
      l1 = (lm + ln * e) * q; l2 = (lm * lm + ln * ln * e) * q - l1 * l1
      gg = ls * it; g1 = l1 * it - ls * it2; g2 = l2 * it - 2.0_dp * l1 * it2 + 2.0_dp * ls * it3
      sst = exp(gg); s1 = sst * g1; s2 = sst * (g2 + g1 * g1)
      lt = [t(3), t(4)]
      do j = 1, 2
        g(j) = exp(-sst + t(j) + (th - 1.0_dp) * lt(j) + (it - 1.0_dp) * ls)
        if (th0 >= 1.0_dp) then
          k1 = -s1 + lt(j) - ls * it2 + (it - 1.0_dp) * l1
          k2 = -s2 + 2.0_dp * ls * it3 - 2.0_dp * l1 * it2 + (it - 1.0_dp) * l2
          gd(j) = g(j) * k1; gdd(j) = g(j) * (k2 + k1 * k1)
        end if
      end do
    case (5)
      th = th0
      if (abs(th) < 1.0e-10_dp) then
        g = [t(2), t(1)]; return
      end if
      call frank_terms(th, t(1), t(2), eu, ev, ed, eu1, ev1, eu2, ev2, ed1, ed2, den, dn1, dn2)
      nn = [(eu + 1.0_dp) * ev, (ev + 1.0_dp) * eu]
      n1 = [eu1 * ev + (eu + 1.0_dp) * ev1, ev1 * eu + (ev + 1.0_dp) * eu1]
      n2 = [eu2 * ev + 2.0_dp * eu1 * ev1 + (eu + 1.0_dp) * ev2, ev2 * eu + 2.0_dp * ev1 * eu1 + (ev + 1.0_dp) * eu2]
      do j = 1, 2
        g(j) = nn(j) / den; gd(j) = (n1(j) - g(j) * dn1) / den; gdd(j) = (n2(j) - 2.0_dp * gd(j) * dn1 - g(j) * dn2) / den
      end do
    case (6)
      th = max(th0, 1.0_dp); it = 1.0_dp / th; it2 = it * it; it3 = it2 * it
      call joe_terms(th, t(1), t(2), ua, ub, ua1, ub1, ua2, ub2, s, s1, s2)
      ls = log(s)
      uu = [ub, ua]; ud = [ub1, ua1]; udd = [ub2, ua2]; lt = [t(1), t(2)]   ! h1 uses (1-b)^th, h2 (1-a)^th
      do j = 1, 2
        g(j) = exp((th - 1.0_dp) * lt(j) + log1p_(-uu(j)) + (it - 1.0_dp) * ls)
        if (th0 >= 1.0_dp) then
          k1 = lt(j) - ud(j) / (1.0_dp - uu(j)) - ls * it2 + (it - 1.0_dp) * s1 / s
          k2 = -udd(j) / (1.0_dp - uu(j)) - (ud(j) / (1.0_dp - uu(j)))**2 + 2.0_dp * ls * it3 - 2.0_dp * s1 * it2 / s &
               + (it - 1.0_dp) * (s2 / s - (s1 / s)**2)
          gd(j) = g(j) * k1; gdd(j) = g(j) * (k2 + k1 * k1)
        end if
      end do
    end select
    g = min(max(g, 0.0_dp), 1.0_dp)
  end subroutine kern_h

  ! copula CDF c (cdf0's) and its theta-derivatives cd, cdd at a rotated, clamped point (a, b)
  pure subroutine kern_c(fam, th0, a, b, c, cd, cdd)
    integer, intent(in) :: fam
    real(dp), intent(in) :: th0, a, b
    real(dp), intent(out) :: c, cd, cdd
    real(dp) :: r, s, x, y, qq, ph, th, e, ls, r1, lx, ly, it, it2, it3, t1, t2, m, mn, &
                lm, ln, q, l1, l2, gg, g1, g2, sst, s1, s2, eu, ev, ed, eu1, ev1, eu2, ev2, ed1, ed2, den, dn1, dn2, xx, &
                x1, x2, ll, ld, ldd, ua, ub, ua1, ub1, ua2, ub2, pp, kk, kd
    c = a * b; cd = 0.0_dp; cdd = 0.0_dp
    select case (fam)
    case (1)
      r = min(max(th0, -0.999999_dp), 0.999999_dp); s = 1.0_dp - r * r
      x = qnorm(a); y = qnorm(b); c = bvn_lower(x, y, r)
      if (abs(th0) <= 0.999999_dp) then                 ! Plackett: dPhi2/drho = phi2
        qq = x * x - 2.0_dp * r * x * y + y * y
        ph = exp(-qq / (2.0_dp * s)) / (2.0_dp * PI * sqrt(s))
        cd = ph; cdd = ph * (r / s + x * y / s - r * qq / (s * s))
      end if
    case (3)
      th = max(th0, 0.0_dp); call clayton_l(th, -log(a), -log(b), l1, l2, r1)
      c = exp(-l1); cd = -l2 * c; cdd = (l2 * l2 - r1) * c
    case (4)
      th = max(th0, 1.0_dp); it = 1.0_dp / th; it2 = it * it; it3 = it2 * it
      lx = log(-log(a)); ly = log(-log(b))
      t1 = th * lx; t2 = th * ly; m = max(t1, t2); mn = min(t1, t2)
      e = exp(mn - m); ls = m + log1p_(e)
      lm = merge(lx, ly, t1 >= t2); ln = merge(ly, lx, t1 >= t2); q = 1.0_dp / (1.0_dp + e)
      l1 = (lm + ln * e) * q; l2 = (lm * lm + ln * ln * e) * q - l1 * l1
      gg = ls * it; g1 = l1 * it - ls * it2; g2 = l2 * it - 2.0_dp * l1 * it2 + 2.0_dp * ls * it3
      sst = exp(gg); s1 = sst * g1; s2 = sst * (g2 + g1 * g1)
      c = exp(-sst)
      if (th0 >= 1.0_dp) then
        cd = -s1 * c; cdd = (s1 * s1 - s2) * c
      end if
    case (5)
      th = th0
      if (abs(th) < 1.0e-10_dp) return
      call frank_terms(th, a, b, eu, ev, ed, eu1, ev1, eu2, ev2, ed1, ed2, den, dn1, dn2)
      xx = eu * ev / ed
      x1 = (eu1 * ev + eu * ev1) / ed - eu * ev * ed1 / (ed * ed)
      x2 = (eu2 * ev + 2.0_dp * eu1 * ev1 + eu * ev2) / ed - 2.0_dp * (eu1 * ev + eu * ev1) * ed1 / (ed * ed) &
           - eu * ev * (ed2 / (ed * ed) - 2.0_dp * ed1 * ed1 / ed**3)
      ll = log1p_(xx); ld = x1 / (1.0_dp + xx); ldd = x2 / (1.0_dp + xx) - ld * ld
      c = -ll / th; cd = ll / (th * th) - ld / th; cdd = -2.0_dp * ll / th**3 + 2.0_dp * ld / (th * th) - ldd / th
    case (6)
      th = max(th0, 1.0_dp); it = 1.0_dp / th; it2 = it * it; it3 = it2 * it
      call joe_terms(th, log1p_(-a), log1p_(-b), ua, ub, ua1, ub1, ua2, ub2, s, s1, s2)
      ls = log(s); pp = exp(ls * it)                    ! C = 1 - s^(1/th)
      kk = s1 * it / s - ls * it2
      kd = (s2 / s - (s1 / s)**2) * it - 2.0_dp * s1 * it2 / s + 2.0_dp * ls * it3
      c = 1.0_dp - pp
      if (th0 >= 1.0_dp) then
        cd = -pp * kk; cdd = -pp * (kd + kk * kk)
      end if
    end select
    c = min(max(c, 0.0_dp), min(a, b))
  end subroutine kern_c

  ! Clayton's generator term L = log(a^-th + b^-th - 1) / th (x = -log a, y = -log b) and its
  ! theta-derivatives L1, L2: log1p/expm1 directly, and its Taylor series in th below 1e-4, where
  ! the direct derivatives cancel (L = x + y - th xy + th^2 xy(x+y)/2 - ...; independence at th = 0)
  pure subroutine clayton_l(th, x, y, l, l1, l2)
    real(dp), intent(in) :: th, x, y
    real(dp), intent(out) :: l, l1, l2
    real(dp) :: p, s1, c1, c2, c3, c4, e, e0, ls, r1, r2
    if (th < 1.0e-4_dp) then
      p = x * y; s1 = x + y
      c1 = -p; c2 = 0.5_dp * p * s1; c3 = -p * (2.0_dp * (x * x + y * y) + 9.0_dp * p) / 12.0_dp
      c4 = p * s1 * (x * x + y * y + 13.0_dp * p) / 24.0_dp
      l = s1 + th * (c1 + th * (c2 + th * (c3 + th * c4)))
      l1 = c1 + th * (2.0_dp * c2 + th * (3.0_dp * c3 + 4.0_dp * th * c4))
      l2 = 2.0_dp * c2 + th * (6.0_dp * c3 + 12.0_dp * th * c4)
    else
      e = exp(th * x); e0 = exp(th * y); ls = log1p_(expm1_(th * x) + expm1_(th * y))
      r1 = (x * e + y * e0) / (e + e0 - 1.0_dp); r2 = (x * x * e + y * y * e0) / (e + e0 - 1.0_dp) - r1 * r1
      l = ls / th; l1 = (r1 - l) / th; l2 = (r2 - 2.0_dp * l1) / th
    end if
  end subroutine clayton_l

  ! Frank: eu = expm1(-th a), ev, ed = expm1(-th), den = ed + eu ev and theta-derivatives
  pure subroutine frank_terms(th, a, b, eu, ev, ed, eu1, ev1, eu2, ev2, ed1, ed2, den, dn1, dn2)
    real(dp), intent(in) :: th, a, b
    real(dp), intent(out) :: eu, ev, ed, eu1, ev1, eu2, ev2, ed1, ed2, den, dn1, dn2
    eu = expm1_(-th * a); ev = expm1_(-th * b); ed = expm1_(-th)
    eu1 = -a * (eu + 1.0_dp); eu2 = a * a * (eu + 1.0_dp)
    ev1 = -b * (ev + 1.0_dp); ev2 = b * b * (ev + 1.0_dp)
    ed1 = -(ed + 1.0_dp); ed2 = ed + 1.0_dp
    den = ed + eu * ev
    dn1 = ed1 + eu1 * ev + eu * ev1
    dn2 = ed2 + eu2 * ev + 2.0_dp * eu1 * ev1 + eu * ev2
  end subroutine frank_terms

  ! Joe: ua = (1-a)^th = exp(th l1), ub, s = ua + ub - ua ub and theta-derivatives
  pure subroutine joe_terms(th, l1, l2, ua, ub, ua1, ub1, ua2, ub2, s, s1, s2)
    real(dp), intent(in) :: th, l1, l2
    real(dp), intent(out) :: ua, ub, ua1, ub1, ua2, ub2, s, s1, s2
    ua = exp(th * l1); ub = exp(th * l2)
    ua1 = l1 * ua; ub1 = l2 * ub; ua2 = l1 * ua1; ub2 = l2 * ub1
    s = ua + ub - ua * ub
    s1 = ua1 + ub1 - ua1 * ub - ua * ub1
    s2 = ua2 + ub2 - ua2 * ub - 2.0_dp * ua1 * ub1 - ua * ub2
  end subroutine joe_terms

  ! Newton's method on [lb, ub] from the closed-form start x0; each pass gives the log-likelihood f,
  ! score g and curvature h.  The iteration stops when the quadratic model predicts a remaining gain
  ! g^2 / 2|h| below 1e-6 log-likelihood units (the exact MLE to that accuracy), so a start inside
  ! the Le Cam one-step neighbourhood costs two passes.  Non-concave points step uphill; any step that
  ! lowers the likelihood is halved; the box is respected.
  subroutine newton1(pass, x0, lb, ub, th, f, ok)
    interface
      subroutine pass(t, f, g, h)
        import :: dp
        real(dp), intent(in) :: t
        real(dp), intent(out) :: f, g, h
      end subroutine pass
    end interface
    real(dp), intent(in) :: x0, lb, ub
    real(dp), intent(out) :: th, f
    logical, intent(out) :: ok
    real(dp) :: lo, hi, g, h, st, tn, fn, gn, hn
    integer :: it, tr
    logical :: acc
    ok = .false.
    lo = lb + 1.0e-10_dp * max(ub - lb, 1.0e-8_dp); hi = ub - 1.0e-10_dp * max(ub - lb, 1.0e-8_dp)
    th = min(max(x0, lo), hi); f = -huge(1.0_dp)
    if (.not. (hi > lo)) return
    call pass(th, f, g, h)
    if (.not. (finite_(f) .and. finite_(g) .and. finite_(h))) return
    do it = 1, 60
      if ((th <= lo .and. g < 0.0_dp) .or. (th >= hi .and. g > 0.0_dp)) exit
      if (abs(g) <= 1.0e-9_dp * (1.0_dp + abs(f))) exit
      if (h < 0.0_dp) then
        if (0.5_dp * g * g / (-h) <= 1.0e-6_dp) exit
        st = -g / h
      else
        st = sign(0.1_dp * max(abs(th), 1.0_dp), g)
      end if
      acc = .false.
      do tr = 1, 50
        tn = min(max(th + st, lo), hi)
        if (tn == th) exit
        call pass(tn, fn, gn, hn)
        if (finite_(fn) .and. finite_(gn) .and. finite_(hn) .and. fn >= f - 1.0e-12_dp * (1.0_dp + abs(f))) then
          acc = .true.; exit
        end if
        st = 0.5_dp * st
      end do
      if (.not. acc) exit
      th = tn; f = fn; g = gn; h = hn
    end do
    ok = .true.
  end subroutine newton1

  ! continuous pair (u1, u2) with transforms t (pc_transforms), rotation rot, start x0, box [lb, ub]
  subroutine mle1(fam, rot, n, u1, u2, t, usew, w, x0, lb, ub, par, ll, ok)
    integer, intent(in) :: fam, rot, n
    real(dp), intent(in) :: u1(n), u2(n), t(n, 10), w(:), x0(3), lb(3), ub(3)
    logical, intent(in) :: usew
    real(dp), intent(out) :: par(3), ll
    logical, intent(out) :: ok
    integer, parameter :: MAP(4, 0:3) = reshape([1, 2, 3, 4, 2, 3, 4, 1, 3, 4, 1, 2, 4, 1, 2, 3], [4, 4])
    integer, parameter :: LMAP(2, 0:3) = reshape([7, 8, 8, 9, 9, 10, 10, 7], [2, 4])
    real(dp), allocatable :: ww(:), ra(:), rb(:)
    real(dp) :: ss(3), rr(3), th, f, fr, lo, hi
    integer :: i, ir, k(4), kl(2), nr
    par = 0.0_dp; ll = -huge(1.0_dp)
    allocate(ww(n))
    if (usew) then
      ww = w(1:n)
    else
      ww = 1.0_dp
    end if
    ir = mod(rot / 90, 4); k = MAP(:, ir); kl = LMAP(:, ir)
    if (fam == 1) then
      ! exact MLE: the score vanishes at a root of W r^3 - B r^2 + (A - W) r - B in the sums
      ! W = sum w, A = sum w (x^2 + y^2), B = sum w x y of the normal scores (or at a box end)
      ss = [sum(ww), sum(ww * (t(:, 5)**2 + t(:, 6)**2)), sum(ww * t(:, 5) * t(:, 6))]
      lo = lb(1) + 1.0e-7_dp * (ub(1) - lb(1)); hi = ub(1) - 1.0e-7_dp * (ub(1) - lb(1))
      call cubic_roots(ss(1), -ss(3), ss(2) - ss(1), -ss(3), rr, nr)
      th = lo; f = gauss_ll(lo)
      fr = gauss_ll(hi)
      if (fr > f) then
        th = hi; f = fr
      end if
      do i = 1, nr
        if (rr(i) > lo .and. rr(i) < hi) then
          fr = gauss_ll(rr(i))
          if (fr > f) then
            th = rr(i); f = fr
          end if
        end if
      end do
      par(1) = th; ll = f; ok = finite_(f)
      return
    end if
    if (fam == 5) then                                    ! Frank (rotationless) uses the pair itself
      allocate(ra(n), rb(n))
      do i = 1, n
        call rotate(rot, clampu(u1(i)), clampu(u2(i)), ra(i), rb(i))
      end do
    end if
    call newton1(pass, x0(1), lb(1), ub(1), th, f, ok)
    par(1) = th; ll = f
  contains
    real(dp) function gauss_ll(r0)
      real(dp), intent(in) :: r0
      real(dp) :: r, s
      r = min(max(r0, -0.999999_dp), 0.999999_dp); s = 1.0_dp - r * r
      gauss_ll = -0.5_dp * ss(1) * log(s) - (r * r * ss(2) - 2.0_dp * r * ss(3)) / (2.0_dp * s)
    end function gauss_ll

    subroutine pass(th0, f, g, h)
      real(dp), intent(in) :: th0
      real(dp), intent(out) :: f, g, h
      real(dp) :: tk(4), lp, d1, d2
      integer :: i
      f = 0.0_dp; g = 0.0_dp; h = 0.0_dp
      do i = 1, n
        select case (fam)
        case (3); tk = [t(i, k(1)), t(i, k(2)), 0.0_dp, 0.0_dp]
        case (4); tk = [-t(i, k(1)), -t(i, k(2)), t(i, kl(1)), t(i, kl(2))]
        case (5); tk = [ra(i), rb(i), 0.0_dp, 0.0_dp]
        case default; tk = [t(i, k(3)), t(i, k(4)), 0.0_dp, 0.0_dp]
        end select
        call kern_lp(fam, th0, tk, lp, d1, d2)
        call acc_ll(ww(i), lp, d1, d2, f, g, h)
      end do
    end subroutine pass
  end subroutine mle1

  ! pair with discrete margins (ctx c after disc_prep): exact discrete likelihood of eval_disc --
  ! densities of untied rows, h-function differences of rows tied in one margin, rectangle
  ! probabilities of rows tied in both -- with closed-form score and curvature
  subroutine mle1d(c, x0, lb, ub, par, ll, ok)
    type(ctx_t), intent(in) :: c
    real(dp), intent(in) :: x0, lb, ub
    real(dp), intent(out) :: par, ll
    logical, intent(out) :: ok
    call newton1(pass, x0, lb, ub, par, ll, ok)
  contains
    subroutine pass(th0, f, g, h)
      real(dp), intent(in) :: th0
      real(dp), intent(out) :: f, g, h
      real(dp) :: lp, d1, d2, wi
      integer :: i
      f = 0.0_dp; g = 0.0_dp; h = 0.0_dp
      do i = 1, c%n
        wi = 1.0_dp
        if (c%usew) wi = c%w(i)
        call row(i, th0, lp, d1, d2)
        call acc_ll(wi, lp, d1, d2, f, g, h)
      end do
    end subroutine pass
    subroutine row(i, th0, lp, d1, d2)
      integer, intent(in) :: i
      real(dp), intent(in) :: th0
      real(dp), intent(out) :: lp, d1, d2
      real(dp) :: hp(2), hpd(2), hpdd(2), hq(2), hqd(2), hqdd(2), dd, ddd, dddd, cc(4), cd(4), cdd(4), sg
      integer :: j, kk
      select case (c%dk(i))
      case (0)
        call kern_lp(c%fam, th0, c%xp(:, i), lp, d1, d2)
      case (1, 2)
        ! h-function k of the original orientation = s0 + sg * g(j) of the rotated pair
        kk = c%dk(i)
        call hrot_map(c%rot, 3 - kk, j, sg)
        call kern_h(c%fam, th0, c%xp(:, i), hp, hpd, hpdd)
        call kern_h(c%fam, th0, c%xq(:, i), hq, hqd, hqdd)
        dd = sg * (hp(j) - hq(j)); ddd = sg * (hpd(j) - hqd(j)); dddd = sg * (hpdd(j) - hqdd(j))
        call log_mass(dd, ddd, dddd, merge(c%a0(i) - c%am(i), c%b0(i) - c%bm(i), kk == 1), lp, d1, d2)
      case default
        ! rectangle of the rotated pair: corners (lo, lo), (hi, lo), (lo, hi), (hi, hi) in xp
        call kern_c(c%fam, th0, c%xp(2, i), c%xp(4, i), cc(1), cd(1), cdd(1))
        call kern_c(c%fam, th0, c%xp(1, i), c%xp(4, i), cc(2), cd(2), cdd(2))
        call kern_c(c%fam, th0, c%xp(2, i), c%xp(3, i), cc(3), cd(3), cdd(3))
        call kern_c(c%fam, th0, c%xp(1, i), c%xp(3, i), cc(4), cd(4), cdd(4))
        dd = cc(1) - cc(2) - cc(3) + cc(4); ddd = cd(1) - cd(2) - cd(3) + cd(4); dddd = cdd(1) - cdd(2) - cdd(3) + cdd(4)
        call log_mass(dd, ddd, dddd, (c%a0(i) - c%am(i)) * (c%b0(i) - c%bm(i)), lp, d1, d2)
      end select
    end subroutine row
  end subroutine mle1d

  ! log(mass / width) and its derivatives (eval_disc's floor: tiny mass scores log(tiny), flat)
  pure subroutine log_mass(m, md, mdd, wd, lp, d1, d2)
    real(dp), intent(in) :: m, md, mdd, wd
    real(dp), intent(out) :: lp, d1, d2
    if (m / wd > TINYV) then
      lp = log(m / wd); d1 = md / m; d2 = mdd / m - d1 * d1
    else
      lp = log(TINYV); d1 = 0.0_dp; d2 = 0.0_dp
    end if
  end subroutine log_mass

  ! eval_rot's h-function k (1, 2) of the original orientation = const + sg * g(j) of the rotated pair
  pure subroutine hrot_map(rot, k, j, sg)
    integer, intent(in) :: rot, k
    integer, intent(out) :: j
    real(dp), intent(out) :: sg
    select case (rot * 10 + k)
    case (901); j = 2; sg = 1.0_dp
    case (902); j = 1; sg = -1.0_dp
    case (1801); j = 1; sg = -1.0_dp
    case (1802); j = 2; sg = -1.0_dp
    case (2701); j = 2; sg = -1.0_dp
    case (2702); j = 1; sg = 1.0_dp
    case default; j = k; sg = 1.0_dp
    end select
  end subroutine hrot_map

  ! weighted log-likelihood, score and curvature accumulation; a non-finite log-density scores
  ! LOGTINY (eval0's floor) and contributes no slope
  pure subroutine acc_ll(wi, lp, d1, d2, f, g, h)
    real(dp), intent(in) :: wi, lp, d1, d2
    real(dp), intent(inout) :: f, g, h
    if (abs(lp) < 1.0e300_dp .and. abs(d1) < 1.0e300_dp .and. abs(d2) < 1.0e300_dp) then
      f = f + wi * lp; g = g + wi * d1; h = h + wi * d2
    else
      f = f + wi * LOGTINY
    end if
  end subroutine acc_ll

  ! row kinds of a discrete pair and the kernel arguments of its points (mle1d): 0 untied (density at
  ! the rotated pair); 1 / 2 tied in the first / second margin (h-function at the pair and at
  ! (a-, b) / (a, b-)); 3 tied in both (rotated rectangle corners alo, ahi, blo, bhi in xp)
  subroutine disc_prep(c)
    type(ctx_t), intent(inout) :: c
    integer :: i
    real(dp) :: ra, rb, sa, sb
    logical :: xa, xb
    allocate(c%dk(c%n), c%xp(4, c%n), c%xq(4, c%n))
    c%xp = 0.0_dp; c%xq = 0.0_dp
    do i = 1, c%n
      xa = c%da .and. c%a0(i) - c%am(i) > 1.0e-12_dp
      xb = c%db .and. c%b0(i) - c%bm(i) > 1.0e-12_dp
      call rotate(c%rot, clampu(c%a0(i)), clampu(c%b0(i)), ra, rb)
      if (xa .and. xb) then
        c%dk(i) = 3
        call rotate(c%rot, clampu(c%am(i)), clampu(c%bm(i)), sa, sb)
        c%xp(:, i) = [min(ra, sa), max(ra, sa), min(rb, sb), max(rb, sb)]
        cycle
      end if
      call kern_t(c%fam, ra, rb, c%xp(:, i))
      if (xa) then
        c%dk(i) = 1
        call rotate(c%rot, clampu(c%am(i)), clampu(c%b0(i)), sa, sb)
      else if (xb) then
        c%dk(i) = 2
        call rotate(c%rot, clampu(c%a0(i)), clampu(c%bm(i)), sa, sb)
      else
        c%dk(i) = 0; cycle
      end if
      call kern_t(c%fam, sa, sb, c%xq(:, i))
    end do
  end subroutine disc_prep

  ! real roots of c3 x^3 + c2 x^2 + c1 x + c0 = 0 (c3 /= 0), Cardano / trigonometric form
  pure subroutine cubic_roots(c3, c2, c1, c0, x, nr)
    real(dp), intent(in) :: c3, c2, c1, c0
    real(dp), intent(out) :: x(3)
    integer, intent(out) :: nr
    real(dp) :: a, b, c, p, q, dsc, m, th, sh, u, v
    a = c2 / c3; b = c1 / c3; c = c0 / c3
    p = b - a * a / 3.0_dp; q = 2.0_dp * a**3 / 27.0_dp - a * b / 3.0_dp + c; sh = -a / 3.0_dp
    dsc = (q / 2.0_dp)**2 + (p / 3.0_dp)**3
    x = 0.0_dp
    if (dsc > 0.0_dp .or. p >= 0.0_dp) then                ! one real root
      u = -q / 2.0_dp + sqrt(max(dsc, 0.0_dp)); v = -q / 2.0_dp - sqrt(max(dsc, 0.0_dp))
      x(1) = sign(abs(u)**(1.0_dp / 3.0_dp), u) + sign(abs(v)**(1.0_dp / 3.0_dp), v) + sh; nr = 1
    else                                                   ! three real roots (p < 0)
      m = 2.0_dp * sqrt(-p / 3.0_dp)
      th = acos(min(max(3.0_dp * q / (p * m), -1.0_dp), 1.0_dp)) / 3.0_dp
      x(1) = m * cos(th) + sh; x(2) = m * cos(th - 2.0_dp * PI / 3.0_dp) + sh; x(3) = m * cos(th - 4.0_dp * PI / 3.0_dp) + sh
      nr = 3
    end if
  end subroutine cubic_roots



  ! bnd present: racing screen (see fit_one_raw); bnd = ll marks a final fit.  xr: polish from
  ! this screen point (the cache key stays the one of the plain fit).  Screens are cached under
  ! their own key, so a repeated selection is all cache hits.
  subroutine fit_one(n, u, fam, rot, method, x0, lb, ub, usew, w, par, ll, nin, ax, um, dsc, bnd, xr, tgt, pool)
    integer, intent(in) :: n, fam, rot, method, nin
    real(dp), intent(in) :: u(2, n), x0(3), lb(3), ub(3), w(:)
    logical, intent(in) :: usew
    real(dp), intent(out) :: par(3), ll
    real(dp), intent(inout), optional :: ax(:)
    real(dp), intent(in), optional :: um(:, :)
    logical, intent(in), optional :: dsc(2)
    real(dp), intent(out), optional :: bnd
    real(dp), intent(in), optional :: xr(3), tgt
    integer, intent(inout), target, optional :: pool(2)
    integer(int64) :: h1, h2
    logical :: dd, hit
    if (fv_abort) then
      par = x0; ll = -huge(1.0_dp)
      if (present(bnd)) bnd = ll
      return
    end if
    if (fc_cap == 0 .or. fam == 0 .or. method == 0) then
      call fit_one_raw(n, u, fam, rot, method, x0, lb, ub, usew, w, par, ll, nin, ax, um, dsc, bnd, xr, tgt, pool)
      if (present(bnd)) bnd = max(bnd, ll)
      return
    end if
    dd = .false.
    if (present(dsc)) dd = any(dsc)
    h1 = 2166136261_int64; h2 = 40503_int64
    call fc_fold(2 * n, u, h1, h2)
    if (usew) call fc_fold(n, w, h1, h2)
    if (dd) then
      call fc_fold(2 * n, um, h1, h2)
      call fc_fold(2, [merge(1.0_dp, 0.0_dp, dsc(1)), merge(1.0_dp, 0.0_dp, dsc(2))], h1, h2)
    end if
    call fc_fold(14, [real(fam, dp), real(rot, dp), real(method, dp), real(n, dp), x0, lb, ub, &
                  merge(1.0_dp, 0.0_dp, usew)], h1, h2)
    h1 = ior(shiftl(h1, 32), h2)
    h2 = ieor(h1, 6148914691236517205_int64)          ! key of the screen
    !$omp critical (fv_fitcache)
    hit = look(h1)
    if (.not. hit .and. present(bnd) .and. .not. present(xr)) hit = look(h2)
    if (hit) then
      fc_hits = fc_hits + 1
    else
      fc_miss = fc_miss + 1
    end if
    !$omp end critical (fv_fitcache)
    if (hit) return
    call fit_one_raw(n, u, fam, rot, method, x0, lb, ub, usew, w, par, ll, nin, ax, um, dsc, bnd, xr, tgt, pool)
    if (fv_abort) return
    if (present(bnd)) then
      if (bnd > ll) then                                 ! a screen, or a polish given up
        call store(h2, bnd); return
      end if
      bnd = ll
    end if
    call store(h1, ll)
  contains
    logical function look(key)                          ! (inside the critical section)
      integer(int64), intent(in) :: key
      integer :: i, j
      look = .false.
      if (fc_cap == 0) return
      do j = 0, 15
        i = modulo(int(modulo(key, int(fc_cap, int64))) + j, fc_cap) + 1
        if (.not. fc_tab(i)%used) return
        if (fc_tab(i)%k1 == key .and. fc_tab(i)%n == n) then
          par = fc_tab(i)%par; ll = fc_tab(i)%ll
          if (present(bnd)) bnd = fc_tab(i)%bnd
          if (present(ax) .and. allocated(fc_tab(i)%ax)) ax(1:size(fc_tab(i)%ax)) = fc_tab(i)%ax
          look = .true.; return
        end if
      end do
    end function look
    subroutine store(key, b)
      integer(int64), intent(in) :: key
      real(dp), intent(in) :: b
      integer :: i, j, slot
      !$omp critical (fv_fitcache)
      if (fc_cap > 0) then
        slot = int(modulo(key, int(fc_cap, int64))) + 1
        i = slot
        do j = 0, 15
          i = modulo(slot - 1 + j, fc_cap) + 1
          if (.not. fc_tab(i)%used .or. (fc_tab(i)%k1 == key .and. fc_tab(i)%n == n)) exit
          if (j == 15) i = slot              ! probe run full: replace the home slot
        end do
        fc_tab(i)%used = .true.; fc_tab(i)%k1 = key; fc_tab(i)%n = n
        fc_tab(i)%par = par; fc_tab(i)%ll = ll; fc_tab(i)%bnd = b
        if (allocated(fc_tab(i)%ax)) deallocate(fc_tab(i)%ax)
        if (present(ax) .and. fam_is_np(fam)) fc_tab(i)%ax = ax
      end if
      !$omp end critical (fv_fitcache)
    end subroutine store
  end subroutine fit_one

  subroutine fit_one_raw(n, u, fam, rot, method, x0, lb, ub, usew, w, par, ll, nin, ax, um, dsc, bnd, xr, tgt, pool)
    integer, intent(in) :: n, fam, rot, method, nin
    logical :: dsc_any, okf
    real(dp), intent(in) :: u(2, n), x0(3), lb(3), ub(3), w(:)
    logical, intent(in) :: usew
    real(dp), intent(out) :: par(3), ll
    real(dp), intent(inout), optional :: ax(:)        ! nonparametric estimate (fam 24, 25)
    real(dp), intent(in), optional :: um(:, :)        ! (2, n) F(x-) of discrete margins
    logical, intent(in), optional :: dsc(2)           ! discrete flags of the two margins
    ! Racing screen (bnd present): a two-stage fit stops at its subsample optimum and returns
    ! the exact all-n log-likelihood there (ll, a lower bound of the MLE's) and an upper bound
    ! bnd = ll + max(10 d, d + 0.08 |ll| + 5), d the predicted Fisher-scoring gain.  Over 23
    ! truth families x 25 candidates x both rotations (n = 5000, 1216 polishes) the polish never
    ! gained more than this minus 11 (0.05 |ll| would still hold; the miss of d grows with the
    ! pair's information), so a candidate whose bound cannot beat the best achieved score is
    ! never the selected one.  Every other path is a full fit and sets bnd = ll.
    real(dp), intent(out), optional :: bnd
    real(dp), intent(in), optional :: xr(3)            ! resume: polish from this screen point
    real(dp), intent(in), optional :: tgt              ! resume: give up below this log-likelihood (bnd > ll)
    integer, intent(inout), target, optional :: pool(2) ! the selection's thread pool (nin = its budget)
    real(dp), allocatable, target :: a0(:), b0(:), am(:), bm(:)
    real(dp) :: lbx(3), ubx(3)
    real(dp), allocatable, target :: a(:), b(:), xs(:), ys(:), ww(:)
    type(ctx_t) :: c, cs
    real(dp) :: xm, fm, lo2, hi2, xs0(3, MS), fs0(MS), pj(3), fj, opt(3, MS + 1), fopt(MS + 1)
    real(dp) :: fj2, gd
    integer :: nopt, j
    integer :: i, st
    logical :: two, scr
    if (present(bnd)) bnd = -huge(1.0_dp)
    scr = present(bnd) .and. .not. present(xr)          ! a racing screen (not a polish)
    allocate(a(n), b(n), ww(n), xs(n), ys(n))
    do i = 1, n
      call rotate(rot, clampu(u(1, i)), clampu(u(2, i)), a(i), b(i))
    end do
    if (usew) then
      ww = w
    else
      ww = 1.0_dp
    end if
    c%fam = fam; c%n = n; c%usew = usew; c%nin = max(1, nin)
    if (present(pool)) c%pool => pool
    c%a => a; c%b => b; c%w => ww
    c%x => xs; c%y => ys
    ! one-parameter families on continuous margins: analytical fit (mle1)
    dsc_any = .false.
    if (present(dsc)) dsc_any = any(dsc)
    if (method == 1 .and. .not. dsc_any .and. any(fam == [1, 3, 4, 5, 6])) then
      block
        real(dp), allocatable :: tt(:, :)
        allocate(tt(n, 10))
        call pc_transforms(n, u(1, :), u(2, :), tt)
        call mle1(fam, rot, n, u(1, :), u(2, :), tt, usew, w, x0, lb, ub, par, ll, okf)
      end block
      if (okf) return
    end if
    ! two-stage fitting: search on a deterministic stride subsample (~NSUB obs),
    ! then polish on all n.  The full-data stopping rule is unchanged, so the
    ! result is still the full-data MLE; it only needs 1-3 full-data iterations.
    two = n > 3 * NSUB .and. method /= 0
    if (two) then
      st = n / NSUB
      cs = c
      cs%a => a(1:n:st); cs%b => b(1:n:st); cs%w => ww(1:n:st)
      cs%x => xs(1:n:st); cs%y => ys(1:n:st)
      cs%n = size(cs%a); cs%nin = 1       ! ~NSUB observations: threads cost more than they save
      nullify(cs%pool)
      cs%tol = 1.0e-4_dp; cs%qn = .true.  ! only a start for the all-n polish: ~1/sqrt(NSUB) precise anyway
    end if
    par = x0
    if (fam == 0) then
      par = 0.0_dp; ll = 0.0_dp; return
    end if
    if (fam_is_np(fam)) then
      par = 0.0_dp
      if (present(dsc)) then
        if (any(dsc)) then
          ! discrete margins (vinecopulib): the estimate is fitted to the mid-ranks
          ! (F(x) + F(x-)) / 2 of the discrete margins; the log-likelihood is the
          ! discrete one (finite differences of h-functions / the exact copula CDF)
          do i = 1, n
            if (dsc(1)) a(i) = clampu(0.5_dp * (u(1, i) + um(1, i)))
            if (dsc(2)) b(i) = clampu(0.5_dp * (u(2, i) + um(2, i)))
          end do
          call np_fit(fam, n, a, b, usew, ww, ax, ll)
          ll = 0.0_dp
          do i = 1, n
            call eval_disc(fam, 0, prep_par(fam, par), ax, clampu(u(1, i)), clampu(u(2, i)), clampu(um(1, i)), &
                           clampu(um(2, i)), dsc(1), dsc(2), xm, fm, lo2, hi2, fj, need_h=.false.)
            if (.not. finite_(xm)) xm = LOGTINY
            ll = ll + ww(i) * xm
          end do
          return
        end if
      end if
      call np_fit(fam, n, a, b, usew, ww, ax, ll)
      return
    end if
    if (present(dsc)) then
      if (any(dsc)) then
        ! discrete margins: value-only likelihood (eval_disc), grid start + Nelder-Mead
        allocate(a0(n), b0(n), am(n), bm(n))
        do i = 1, n
          a0(i) = clampu(u(1, i)); b0(i) = clampu(u(2, i))
          am(i) = clampu(um(1, i)); bm(i) = clampu(um(2, i))
        end do
        c%disc = .true.; c%da = dsc(1); c%db = dsc(2); c%rot = rot
        c%a0 => a0; c%b0 => b0; c%am => am; c%bm => bm
        call par_bounds(fam, lbx, ubx)
        if (two) then
          cs%disc = .true.; cs%da = dsc(1); cs%db = dsc(2); cs%rot = rot
          cs%a0 => a0(1:n:st); cs%b0 => b0(1:n:st); cs%am => am(1:n:st); cs%bm => bm(1:n:st)
        end if
        if (method == 1 .and. any(fam == [1, 3, 4, 5, 6])) then      ! analytical (mle1d)
          call disc_prep(c)
          call mle1d(c, min(max(x0(1), lbx(1)), ubx(1)), lbx(1), ubx(1), xm, fm, okf)
          if (okf) then
            par = [xm, 0.0_dp, 0.0_dp]; ll = fm
            return
          end if
        end if
        if (npar_of(fam) == 1) then
          ! one parameter: Brent (superlinear, derivative-free) on the exact discrete likelihood
          c%mode = 1; lo2 = lbx(1); hi2 = ubx(1)
          if (two) then
            cs%mode = 1
            call brent(cs, lbx(1), ubx(1), 1.0e-3_dp, xm, fm)
            lo2 = max(lbx(1), xm - 0.05_dp * (ubx(1) - lbx(1))); hi2 = min(ubx(1), xm + 0.05_dp * (ubx(1) - lbx(1)))
          end if
          call brent(c, lo2, hi2, 1.0e-7_dp, xm, fm)
          par = [xm, 0.0_dp, 0.0_dp]; ll = -fm
          return
        end if
        if (two) then
          call ext_grid_start(cs, fam, lbx, ubx, xs0, fs0)
          call keep_top(min(max(x0, lbx), ubx), ll_par(cs, min(max(x0, lbx), ubx)), xs0, fs0)
          call nmk(cs, npar_of(fam), xs0(:, 1), lbx, ubx, pj, fm)
          call nmk(c, npar_of(fam), pj, lbx, ubx, par, fm, polish=.true.)
        else
          call ext_grid_start(c, fam, lbx, ubx, xs0, fs0)
          call keep_top(min(max(x0, lbx), ubx), ll_par(c, min(max(x0, lbx), ubx)), xs0, fs0)
          call nmk(c, npar_of(fam), xs0(:, 1), lbx, ubx, par, fm)
        end if
        ll = -fm
        return
      end if
    end if
    if (fam == 2) then
      if (method == 0) then
        call student_scores(c, x0(2))
        ll = ll_student_xy(c, x0(1), x0(2))
      else if (method == 2) then
        c%mode = 4; c%fixed = x0(1)
        lo2 = lb(2); hi2 = ub(2)
        if (two) then
          cs%mode = 4; cs%fixed = x0(1)
          call brent(cs, lb(2), ub(2), 1.0e-3_dp, xm, fm)
          lo2 = max(lb(2), xm / 2.5_dp); hi2 = min(ub(2), xm * 2.5_dp)
        end if
        call brent(c, lo2, hi2, 1.0e-4_dp, xm, fm)
        par = [x0(1), xm, 0.0_dp]; ll = -fm
      else
        c%mode = 3; c%lbr = lb(1); c%ubr = ub(1)
        lo2 = lb(2); hi2 = ub(2)
        if (present(xr)) then                   ! resume: the full-data bracket around the screen's nu
          lo2 = max(lb(2), xr(2) / 2.5_dp); hi2 = min(ub(2), xr(2) * 2.5_dp)
        else if (two) then
          cs%mode = 3; cs%lbr = lb(1); cs%ubr = ub(1); cs%lastrho = 0.5_dp * (lb(1) + ub(1))
          call brent(cs, lb(2), ub(2), 1.0e-3_dp, xm, fm)
          lo2 = max(lb(2), xm / 2.5_dp); hi2 = min(ub(2), xm * 2.5_dp)
          if (scr) then                         ! screen: exact all-n value at the subsample nu (rho solved)
            call student_scores(c, xm)
            par = [student_rho_newton(c, xm, cs%lastrho, lb(1), ub(1)), xm, 0.0_dp]
            ll = ll_student_xy(c, par(1), xm)
            bnd = ll + 0.15_dp * abs(ll) + 5.0_dp ! nu gains reached 0.055 |ll| (23 truths, n = 2000-20000)
            return
          end if
        end if
        c%lastrho = 0.5_dp * (lb(1) + ub(1))
        xm = -1.0_dp
        if (hi2 >= ub(2)) then                  ! optimum on the upper bound (the near-Gaussian end):
          fm = f1(c, 0.999_dp * hi2)            ! Brent would take ~17 evaluations to crawl there
          if (f1(c, hi2) < fm) xm = hi2
        end if
        if (xm < 0.0_dp) call brent(c, lo2, hi2, 1.0e-4_dp, xm, fm)
        par = [c%lastrho, xm, 0.0_dp]
        ! lastrho belongs to the last evaluated nu; refit rho at the optimum
        call student_scores(c, xm)
        par(1) = student_rho_newton(c, xm, par(1), lb(1), ub(1))
        ll = ll_student_xy(c, par(1), xm)
      end if
      return
    end if
    if (method == 0) then
      ll = ll_par(c, x0)
      return
    end if
    if (fam == 23) then
      ! t-EV: nu has no closed-form score (the t CDF's df-derivative), so nu is profiled
      ! by Brent while rho is solved by Newton on its exact score (as for the Student)
      ! (profile in log nu; subsample first, then all n around the subsample optimum).
      ! The subsample profile locates the basin; all n are then fitted jointly in
      ! (rho, nu) by Fisher scoring with the exact nu-score (dpt_dnu), which replaces the
      ! full-data Brent profile (~10 inner rho fits) by a few joint Newton steps.
      c%mode = 5; c%lbr = lb(1); c%ubr = ub(1); c%lastrho = 0.5_dp; c%qn = .true.
      call tev_prep(c)
      if (present(xr)) then
        c%lastrho = xr(1); xm = log(xr(2))
      else if (two) then
        cs%mode = 5; cs%lbr = lb(1); cs%ubr = ub(1); cs%lastrho = 0.5_dp; cs%qn = .true.
        call tev_prep(cs)
        call brent(cs, log(lb(2)), log(ub(2)), 5.0e-2_dp, xm, fm)
        fm = f1(cs, xm)                    ! lastrho of the optimum, not of the last trial
        c%lastrho = cs%lastrho
      else
        call brent(c, log(lb(2)), log(ub(2)), 1.0e-3_dp, xm, fm)
        fm = f1(c, xm)
      end if
      ! joint Newton in (rho, log nu): the likelihood ridge is near-linear in log nu
      c%lognu = .true.
      lbx = [lb(1), log(lb(2)), 0.0_dp]; ubx = [ub(1), log(ub(2)), 0.0_dp]
      if (two .and. scr) then
        par = [c%lastrho, xm, 0.0_dp]
        call pred_gain(c, 2, par, lbx, ubx, ll, gd)
        if (gd >= 0.0_dp) then
          bnd = ll + max(10.0_dp * gd, gd + 0.08_dp * abs(ll) + 5.0_dp)
          c%lognu = .false.; par(2) = exp(xm)
          return
        end if
      end if
      ! Near-independent pairs: the likelihood has an independence plateau (rho -> -1:
      ! T(z) -> 1, A -> 1) whose supremum in the box is the corner (rho_min, nu_max), and
      ! the subsample cannot tell it from a weak-dependence basin.  If the subsample
      ! optimum is on the plateau (rho at its bound), the corner replaces the slow walk
      ! along it and Newton starts from the dependence side; if the full-data fit is
      ! weak (ll < 5), both alternatives are tried as well.
      if (c%lastrho <= lb(1) + 1.0e-6_dp) then
        par = [lb(1), log(ub(2)), 0.0_dp]; ll = ll_par(c, par)
        call newton_fit(c, 2, [0.5_dp, xm, 0.0_dp], lbx, ubx, pj, fj)
        if (fj > ll) then
          par = pj; ll = fj
        end if
      else
        ! most fits end on the nu = nu_max (Husler-Reiss) face, which the subsample rarely
        ! resolves: a 1-D rho fit on that face is a few cheap steps, and the joint Newton
        ! starts from whichever start is better on all n
        c%lognu = .false.
        ! (from the dependence side: on the face the rho profile also has the plateau)
        call newton_fit(c, 1, [max(c%lastrho, 0.3_dp), ub(2), 0.0_dp], [lb(1), ub(2), 0.0_dp], &
                        [ub(1), ub(2), 0.0_dp], pj, fj)
        c%lognu = .true.
        fm = ll_par(c, [c%lastrho, xm, 0.0_dp])
        if (fj > fm) then
          call newton_fit(c, 2, [pj(1), log(ub(2)), 0.0_dp], lbx, ubx, par, ll)
        else
          call newton_fit(c, 2, [c%lastrho, xm, 0.0_dp], lbx, ubx, par, ll)
        end if
        if (ll < 5.0_dp) then
          fj = ll_par(c, [lb(1), log(ub(2)), 0.0_dp])
          if (fj > ll) then
            par = [lb(1), log(ub(2)), 0.0_dp]; ll = fj
          end if
        end if
      end if
      c%lognu = .false.
      par(2) = exp(par(2))
      par(3) = 0.0_dp
      return
    end if
    do i = 1, n
      call prep_xy(fam, a(i), b(i), xs(i), ys(i))
    end do
    if (two .and. present(xr)) then
      nopt = 1; opt(:, 1) = xr; fopt(1) = 0.0_dp
    else
      ! multi-start: Newton from the MS best grid points (on the subsample when
      ! two-stage), keep the best, then polish on all n.  Guards against ridge /
      ! boundary local optima (e.g. BB5 drifting to its delta -> 0 Gumbel limit).
      xs0(:, 1) = x0; fs0 = -huge(1.0_dp); fs0(1) = 0.0_dp
      if (two) then
        if (fam >= 7 .and. fam <= 10) call grid_start(cs, npar_of(fam), x0, lb, ub, xs0, fs0)
        if (fam >= 12 .and. .not. gen_is(fam)) call ext_grid_start(cs, fam, lb, ub, xs0, fs0)
        if (gen_is(fam)) call gen_grid_start(cs, fam, lb, ub, xs0, fs0)
      else
        if (fam >= 7 .and. fam <= 10) call grid_start(c, npar_of(fam), x0, lb, ub, xs0, fs0)
        if (fam >= 12 .and. .not. gen_is(fam)) call ext_grid_start(c, fam, lb, ub, xs0, fs0)
        if (gen_is(fam)) call gen_grid_start(c, fam, lb, ub, xs0, fs0)
      end if
      if (fam >= 12) then
        if (two) then
          call keep_top(x0, ll_par(cs, x0), xs0, fs0)
        else
          call keep_top(x0, ll_par(c, x0), xs0, fs0)
        end if
      end if
      ! stage 1: Newton from each start (subsample when two-stage)
      nopt = 0
      ! multi-start only where the likelihood is known to be multimodal / ridge-shaped
      ! (BB5, BB9); single-start Newton reproduces vinecopulib's optima elsewhere
      do i = 1, merge(MS, 1, fam == 19 .or. fam == 20 .or. gen_is(fam))
        if (fs0(i) <= -huge(1.0_dp)) exit
        if (two) then
          call newton_fit(cs, npar_of(fam), xs0(:, i), lb, ub, pj, fj)
        else
          call newton_fit(c, npar_of(fam), xs0(:, i), lb, ub, pj, fj)
        end if
        if (.not. finite_(fj)) cycle
        ! keep distinct optima only (the likelihood can be multimodal, e.g. BB5)
        if (any([(maxval(abs(pj(1:3) - opt(1:3, j)) / (1.0_dp + abs(opt(1:3, j)))) < 5.0e-2_dp, j = 1, nopt)])) then
          do j = 1, nopt      ! same mode: keep the better representative
            if (maxval(abs(pj(1:3) - opt(1:3, j)) / (1.0_dp + abs(opt(1:3, j)))) < 5.0e-2_dp .and. fj > fopt(j)) then
              opt(:, j) = pj; fopt(j) = fj
            end if
          end do
          cycle
        end if
        nopt = nopt + 1; opt(:, nopt) = pj; fopt(nopt) = fj
      end do
      if (nopt == 0) then
        nopt = 1; opt(:, 1) = xs0(:, 1)
      end if
    end if
    ! stage 2: polish every distinct optimum on all n, keep the best
    fm = -huge(1.0_dp)
    fj = maxval(fopt(1:nopt))
    if (two .and. scr) then                      ! racing screen: bound every competitive mode
      do i = 1, nopt
        if (fopt(i) < fj - 2.0_dp) cycle
        call pred_gain(c, npar_of(fam), opt(:, i), lb, ub, fj2, gd)
        if (gd < 0.0_dp .or. .not. finite_(fj2)) then
          bnd = -huge(1.0_dp); exit
        end if
        bnd = max(bnd, fj2 + max(10.0_dp * gd, gd + 0.08_dp * abs(fj2) + 5.0_dp))
        if (fj2 + gd > fm) then                   ! the resume polishes the best predicted mode
          fm = fj2 + gd; par = opt(:, i); ll = fj2
        end if
      end do
      if (bnd > -huge(1.0_dp)) return
      fm = -huge(1.0_dp)
    end if
    c%qn = .true.                                 ! BFGS polish: BHHH alone converges linearly off-model
    if (present(tgt)) c%tgt = tgt
    do i = 1, nopt
      ! only modes competitive on the subsample (within 2 log-lik units) are polished
      if (two .and. fopt(i) < fj - 2.0_dp) cycle
      if (two) then
        call newton_fit(c, npar_of(fam), opt(:, i), lb, ub, pj, fj2)
      else
        pj = opt(:, i); fj2 = fopt(i)
      end if
      if (finite_(fj2) .and. fj2 > fm) then
        fm = fj2; par = pj
      end if
    end do
    ll = fm
    if (present(bnd) .and. c%cut > -huge(1.0_dp)) bnd = c%cut
    if (.not. finite_(ll)) then
      ! ponytail: derivative-free fallback; only reached if Newton cannot start
      call nmk(c, npar_of(fam), x0, lb, ub, par, fm)
      ll = -fm
    end if
  end subroutine fit_one_raw

  ! --------------------------------------------------------- structure arrays
  ! From the natural-order struct array (edge, tree) derive min array and the
  ! needed-h-function flags (vinecopulib RVineStructure::compute_*).
  pure subroutine struct_arrays(d, trunc, st, mn, n1, n2)
    integer, intent(in) :: d, trunc, st(d, d)
    integer, intent(out) :: mn(d, d), n1(d, d), n2(d, d)
    integer :: i, j
    mn = st; n1 = 0; n2 = 0
    do j = 1, d - 1
      do i = 2, min(d - j, trunc)
        mn(j, i) = min(st(j, i), mn(j, i - 1))
      end do
    end do
    do i = 1, min(d - 2, max(trunc - 1, 0))
      do j = 1, d - 1 - i
        if (st(j, i + 1) /= mn(j, i + 1)) n1(mn(j, i + 1), i) = 1
        n2(j, i) = 1
        if (st(j, i + 1) == mn(j, i + 1)) n2(mn(j, i + 1), i) = 1
      end do
    end do
  end subroutine struct_arrays

  pure subroutine prep_all(d, trunc, fam, par, ppv)
    integer, intent(in) :: d, trunc, fam(d, d)
    real(dp), intent(in) :: par(3, d, d)
    real(dp), intent(out) :: ppv(NPP, d, d)
    integer :: t, e
    ppv = 0.0_dp
    do t = 1, trunc
      do e = 1, d - t
        ppv(:, e, t) = prep_par(fam(e, t), par(:, e, t))
      end do
    end do
  end subroutine prep_all

  ! --------------------------------------------------------- C API: vine pdf
  ! Vine density on (nsamp, bins, d) data in ORIGINAL variable order.  Every
  ! column except `tgt` (0-based original index, -1 = none) is constant across
  ! bins, so edges whose inputs are all bin-constant are evaluated once per
  ! sample.  order(d): 1-based structure order.  st/fam/rot: torch (d, d)
  ! [tree, edge] -> Fortran (edge, tree); par: (3, edge, tree).
  subroutine fv_vine_pdf(nsamp, bins, d, tgt, u, order, trunc, st, fam, rot, par, out, nthreads) &
      bind(C, name="fv_vine_pdf")
    integer(c_int), value :: nsamp, bins, d, tgt, trunc, nthreads
    real(c_double), intent(in) :: u(d, bins, nsamp), par(3, d, d)
    integer(c_int), intent(in) :: order(d), st(d, d), fam(d, d), rot(d, d)
    real(c_double), intent(out) :: out(bins, nsamp)
    real(dp) :: axd(1, d, d)
    call vine_pdf_x(nsamp, bins, d, tgt, u, order, trunc, st, fam, rot, par, 1, axd, out, nthreads)
  end subroutine fv_vine_pdf

  ! fv_vine_pdf with the aux vectors of nonparametric pair copulas, aux(na, d, d)
  ! erel(2, nrel) / prel(3, nrel, nsamp): non-simplified edges (edge, tree) and their parameters
  ! per sample (they must not depend on `tgt`: the bins share them)
  subroutine vine_pdf_x(nsamp, bins, d, tgt, u, order, trunc, st, fam, rot, par, na, aux, out, nthreads, erel, prel)
    integer, intent(in) :: nsamp, bins, d, tgt, trunc, nthreads, na
    real(dp), intent(in) :: u(d, bins, nsamp), par(3, d, d), aux(na, d, d)
    integer, intent(in) :: order(d), st(d, d), fam(d, d), rot(d, d)
    real(dp), intent(out) :: out(bins, nsamp)
    integer, intent(in), optional :: erel(:, :)
    real(dp), intent(in), optional :: prel(:, :, :)
    real(dp), allocatable :: hf1(:, :), hf2(:, :), lpf(:), t1(:), t2(:), ppl(:, :, :)
    logical, allocatable :: fu1(:), fu2(:)
    integer :: mn(d, d), nh1(d, d), nh2(d, d)
    real(dp) :: ppv(NPP, d, d)
    integer :: i, j, t, e, m, nb, k, kb1, kb2, ntgt, nrel, r
    logical :: fa, fb, full, useh2
    real(dp) :: lpc, lp, h1, h2, v1, v2
    call struct_arrays(d, trunc, st, mn, nh1, nh2)
    call prep_all(d, trunc, fam, par, ppv)
    nrel = 0
    if (present(erel) .and. present(prel)) nrel = size(erel, 2)
    ntgt = 0
    do j = 1, d
      if (order(j) - 1 == tgt) ntgt = j
    end do
    !$omp parallel private(hf1, hf2, lpf, t1, t2, fu1, fu2, i, j, t, e, m, nb, k, kb1, kb2, &
    !$omp&                 fa, fb, full, useh2, lpc, lp, h1, h2, v1, v2, ppl, r) num_threads(nthreads)
    allocate(hf1(bins, d), hf2(bins, d), lpf(bins), t1(bins), t2(bins), fu1(d), fu2(d))
    ppl = ppv
    !$omp do schedule(dynamic, 4)        ! hybrid P/E cores: no static split
    do i = 1, nsamp
      do r = 1, nrel
        ppl(:, erel(1, r), erel(2, r)) = prep_par(fam(erel(1, r), erel(2, r)), prel(:, r, i))
      end do
      do j = 1, d
        fu2(j) = (j == ntgt) .and. bins > 1
        if (fu2(j)) then
          do k = 1, bins
            hf2(k, j) = clampu(u(order(j), k, i))
          end do
        else
          hf2(1, j) = clampu(u(order(j), 1, i))
        end if
      end do
      hf1 = hf2; fu1 = fu2
      lpc = 0.0_dp; lpf = 0.0_dp
      do t = 1, trunc
        do e = 1, d - t
          m = mn(e, t)
          useh2 = (m == st(e, t))
          fa = fu2(e)
          if (useh2) then
            fb = fu2(m)
          else
            fb = fu1(m)
          end if
          full = fa .or. fb
          nb = merge(bins, 1, full)
          do k = 1, nb
            kb1 = merge(k, 1, fa)
            kb2 = merge(k, 1, fb)
            v1 = hf2(kb1, e)
            if (useh2) then
              v2 = hf2(kb2, m)
            else
              v2 = hf1(kb2, m)
            end if
            if (fam(e, t) == 0) then
              lp = 0.0_dp; h1 = v2; h2 = v1
            else
              call eval_x(fam(e, t), rot(e, t), ppl(:, e, t), aux(:, e, t), v1, v2, &
                          nh1(e, t) /= 0 .or. nh2(e, t) /= 0, lp, h1, h2)
            end if
            if (full) then
              lpf(k) = lpf(k) + lp
            else
              lpc = lpc + lp
            end if
            t1(k) = h1; t2(k) = h2
          end do
          if (nh1(e, t) /= 0) then
            hf1(1:nb, e) = t1(1:nb); fu1(e) = full
          end if
          if (nh2(e, t) /= 0) then
            hf2(1:nb, e) = t2(1:nb); fu2(e) = full
          end if
        end do
      end do
      do k = 1, bins
        out(k, i) = max(exp(lpc + lpf(k)), TINYV)
      end do
    end do
    !$omp end do
    deallocate(hf1, hf2, lpf, t1, t2, fu1, fu2)
    !$omp end parallel
  end subroutine vine_pdf_x

  ! Inputs (v1, v2) of every edge of tree tq at the samples u(d, n) (original variable order), as
  ! vine_pdf_x feeds them to the pair copulas; erel / prel as there.  v1o / v2o(e, i), e <= d - tq.
  subroutine vine_edge_uv(n, d, u, order, trunc, st, fam, rot, par, na, aux, tq, v1o, v2o, nthreads, erel, prel)
    integer, intent(in) :: n, d, trunc, na, tq, nthreads
    real(dp), intent(in) :: u(d, n), par(3, d, d), aux(na, d, d)
    integer, intent(in) :: order(d), st(d, d), fam(d, d), rot(d, d)
    real(dp), intent(out) :: v1o(d, n), v2o(d, n)
    integer, intent(in), optional :: erel(:, :)
    real(dp), intent(in), optional :: prel(:, :, :)
    integer :: mn(d, d), nh1(d, d), nh2(d, d), i, j, t, e, m, r, nrel
    real(dp) :: ppv(NPP, d, d), hf1(d), hf2(d), lp, h1, h2, v1, v2
    real(dp), allocatable :: ppl(:, :, :)
    call struct_arrays(d, trunc, st, mn, nh1, nh2)
    call prep_all(d, trunc, fam, par, ppv)
    nrel = 0
    if (present(erel) .and. present(prel)) nrel = size(erel, 2)
    v1o = 0.5_dp; v2o = 0.5_dp
    !$omp parallel private(i, j, t, e, m, r, hf1, hf2, lp, h1, h2, v1, v2, ppl) num_threads(nthreads)
    ppl = ppv
    !$omp do schedule(static)
    do i = 1, n
      do r = 1, nrel
        ppl(:, erel(1, r), erel(2, r)) = prep_par(fam(erel(1, r), erel(2, r)), prel(:, r, i))
      end do
      do j = 1, d
        hf2(j) = clampu(u(order(j), i))
      end do
      hf1 = hf2
      do t = 1, min(tq, trunc)
        do e = 1, d - t
          m = mn(e, t)
          v1 = hf2(e)
          if (m == st(e, t)) then
            v2 = hf2(m)
          else
            v2 = hf1(m)
          end if
          if (t == tq) then
            v1o(e, i) = v1; v2o(e, i) = v2
            cycle
          end if
          if (fam(e, t) == 0) then
            h1 = v2; h2 = v1
          else
            call eval_x(fam(e, t), rot(e, t), ppl(:, e, t), aux(:, e, t), v1, v2, .true., lp, h1, h2)
          end if
          if (nh1(e, t) /= 0) hf1(e) = h1
          if (nh2(e, t) /= 0) hf2(e) = h2
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine vine_edge_uv


  ! Vine density with discrete variables (vinecopulib var_types): u(d, n) = F(x),
  ! um(d, n) = F(x-) (used where disc(j)), original variable order.  Each stored
  ! pseudo-observation carries its minus version and a discreteness flag.
  subroutine vine_pdf_disc(n, d, u, um, disc, order, trunc, st, fam, rot, par, na, aux, out, nthreads)
    integer, intent(in) :: n, d, trunc, nthreads, na
    real(dp), intent(in) :: u(d, n), um(d, n), par(3, d, d), aux(na, d, d)
    logical, intent(in) :: disc(d)
    integer, intent(in) :: order(d), st(d, d), fam(d, d), rot(d, d)
    real(dp), intent(out) :: out(n)
    integer :: mn(d, d), nh1(d, d), nh2(d, d), i, j, t, e, m
    real(dp) :: hf1(d), hf2(d), hm1(d), hm2(d), ppv(NPP, d, d), v2, v2m, lp, h1, h2, q1, q2, s
    logical :: df1(d), df2(d), d2
    call struct_arrays(d, trunc, st, mn, nh1, nh2)
    call prep_all(d, trunc, fam, par, ppv)
    !$omp parallel do private(i, j, t, e, m, hf1, hf2, hm1, hm2, v2, v2m, lp, h1, h2, q1, q2, s, df1, df2, d2) &
    !$omp& num_threads(nthreads) schedule(dynamic, 256)
    do i = 1, n
      do j = 1, d
        hf2(j) = clampu(u(order(j), i)); df2(j) = disc(order(j))
        hm2(j) = hf2(j)
        if (df2(j)) hm2(j) = clampu(um(order(j), i))
      end do
      hf1 = hf2; hm1 = hm2; df1 = df2
      s = 0.0_dp
      do t = 1, trunc
        do e = 1, d - t
          m = mn(e, t)
          if (m == st(e, t)) then
            v2 = hf2(m); v2m = hm2(m); d2 = df2(m)
          else
            v2 = hf1(m); v2m = hm1(m); d2 = df1(m)
          end if
          call eval_disc(fam(e, t), rot(e, t), ppv(:, e, t), aux(:, e, t), hf2(e), v2, hm2(e), v2m, df2(e), d2, &
                         lp, h1, h2, q1, q2)
          s = s + lp
          if (nh1(e, t) /= 0) then
            hf1(e) = h1; hm1(e) = q1; df1(e) = d2
          end if
          hf2(e) = h2; hm2(e) = q2
        end do
      end do
      out(i) = max(exp(s), TINYV)
    end do
    !$omp end parallel do
  end subroutine vine_pdf_disc

  ! vine_pdf_disc on a prediction grid: u, um (d, bins, nsamp) in ORIGINAL variable order where
  ! only variable `tgt` (0-based) varies across bins.  As in vine_pdf_x, edges whose inputs are
  ! all bin-constant are evaluated once per sample, and h-functions only where a later tree
  ! reads them; every pair term is the same eval_disc call as in vine_pdf_disc.
  subroutine vine_pdf_disc_x(nsamp, bins, d, tgt, u, um, disc, order, trunc, st, fam, rot, par, na, aux, out, &
                             nthreads)
    integer, intent(in) :: nsamp, bins, d, tgt, trunc, nthreads, na
    real(dp), intent(in) :: u(d, bins, nsamp), um(d, bins, nsamp), par(3, d, d), aux(na, d, d)
    logical, intent(in) :: disc(d)
    integer, intent(in) :: order(d), st(d, d), fam(d, d), rot(d, d)
    real(dp), intent(out) :: out(bins, nsamp)
    real(dp), allocatable :: hf1(:, :), hf2(:, :), hm1(:, :), hm2(:, :), lpf(:), t1(:), t2(:), q1(:), q2(:)
    logical, allocatable :: fu1(:), fu2(:)
    logical :: df1(d), df2(d), d2, fa, fb, full, useh2, wh
    integer :: mn(d, d), nh1(d, d), nh2(d, d)
    real(dp) :: ppv(NPP, d, d), lpc, lp, v1, v1m, v2, v2m
    integer :: i, j, t, e, m, nb, k, kb1, kb2, ntgt
    call struct_arrays(d, trunc, st, mn, nh1, nh2)
    call prep_all(d, trunc, fam, par, ppv)
    ntgt = 0
    do j = 1, d
      if (order(j) - 1 == tgt) ntgt = j
    end do
    !$omp parallel private(hf1, hf2, hm1, hm2, lpf, t1, t2, q1, q2, fu1, fu2, df1, df2, d2, fa, fb, full, useh2, wh, &
    !$omp&                 lpc, lp, v1, v1m, v2, v2m, i, j, t, e, m, nb, k, kb1, kb2) num_threads(nthreads)
    allocate(hf1(bins, d), hf2(bins, d), hm1(bins, d), hm2(bins, d), lpf(bins), t1(bins), t2(bins), q1(bins), &
             q2(bins), fu1(d), fu2(d))
    !$omp do schedule(dynamic, 4)
    do i = 1, nsamp
      do j = 1, d
        fu2(j) = (j == ntgt) .and. bins > 1
        df2(j) = disc(order(j))
        do k = 1, merge(bins, 1, fu2(j))
          hf2(k, j) = clampu(u(order(j), k, i))
          hm2(k, j) = hf2(k, j)
          if (df2(j)) hm2(k, j) = clampu(um(order(j), k, i))
        end do
      end do
      hf1 = hf2; hm1 = hm2; df1 = df2; fu1 = fu2
      lpc = 0.0_dp; lpf = 0.0_dp
      do t = 1, trunc
        do e = 1, d - t
          m = mn(e, t)
          useh2 = (m == st(e, t))
          fa = fu2(e)
          if (useh2) then
            fb = fu2(m); d2 = df2(m)
          else
            fb = fu1(m); d2 = df1(m)
          end if
          full = fa .or. fb
          nb = merge(bins, 1, full)
          wh = nh1(e, t) /= 0 .or. nh2(e, t) /= 0
          do k = 1, nb
            kb1 = merge(k, 1, fa); kb2 = merge(k, 1, fb)
            v1 = hf2(kb1, e); v1m = hm2(kb1, e)
            if (useh2) then
              v2 = hf2(kb2, m); v2m = hm2(kb2, m)
            else
              v2 = hf1(kb2, m); v2m = hm1(kb2, m)
            end if
            call eval_disc(fam(e, t), rot(e, t), ppv(:, e, t), aux(:, e, t), v1, v2, v1m, v2m, df2(e), d2, &
                           lp, t1(k), t2(k), q1(k), q2(k), need_h=wh)
            if (full) then
              lpf(k) = lpf(k) + lp
            else
              lpc = lpc + lp
            end if
          end do
          if (nh1(e, t) /= 0) then
            hf1(1:nb, e) = t1(1:nb); hm1(1:nb, e) = q1(1:nb); df1(e) = d2; fu1(e) = full
          end if
          if (nh2(e, t) /= 0) then
            hf2(1:nb, e) = t2(1:nb); hm2(1:nb, e) = q2(1:nb); fu2(e) = full
          end if
        end do
      end do
      do k = 1, bins
        out(k, i) = max(exp(lpc + lpf(k)), TINYV)
      end do
    end do
    !$omp end do
    deallocate(hf1, hf2, hm1, hm2, lpf, t1, t2, q1, q2, fu1, fu2)
    !$omp end parallel
  end subroutine vine_pdf_disc_x

  ! Rosenblatt transform with discrete variables (vinecopulib Vinecop::rosenblatt):
  ! out = F(x_k | x_1..x_{k-1}) and outm = its left limit F(x_k- | ...), propagated with
  ! the h-functions at F(x) and at F(x-) of the discrete arguments (eval_disc, as in
  ! vine_pdf_disc).  Continuous variables: outm = out.  Brockwell's randomisation
  ! W out + (1 - W) outm is applied by the caller.
  subroutine vine_rosenblatt_disc(n, d, u, um, disc, order, trunc, st, fam, rot, par, na, aux, out, outm, nthreads)
    integer, intent(in) :: n, d, trunc, nthreads, na
    real(dp), intent(in) :: u(d, n), um(d, n), par(3, d, d), aux(na, d, d)
    logical, intent(in) :: disc(d)
    integer, intent(in) :: order(d), st(d, d), fam(d, d), rot(d, d)
    real(dp), intent(out) :: out(d, n), outm(d, n)
    integer :: mn(d, d), nh1(d, d), nh2(d, d), i, j, t, e, m
    real(dp) :: hf1(d), hf2(d), hm1(d), hm2(d), ppv(NPP, d, d), v2, v2m, lp, h1, h2, q1, q2
    logical :: df1(d), df2(d), d2
    call struct_arrays(d, trunc, st, mn, nh1, nh2)
    call prep_all(d, trunc, fam, par, ppv)
    !$omp parallel do private(i, j, t, e, m, hf1, hf2, hm1, hm2, v2, v2m, lp, h1, h2, q1, q2, df1, df2, d2) &
    !$omp& num_threads(nthreads) schedule(static)
    do i = 1, n
      do j = 1, d
        hf2(j) = clampu(u(order(j), i)); df2(j) = disc(order(j))
        hm2(j) = hf2(j)
        if (df2(j)) hm2(j) = clampu(um(order(j), i))
      end do
      hf1 = hf2; hm1 = hm2; df1 = df2
      do t = 1, trunc
        do e = 1, d - t
          m = mn(e, t)
          if (m == st(e, t)) then
            v2 = hf2(m); v2m = hm2(m); d2 = df2(m)
          else
            v2 = hf1(m); v2m = hm1(m); d2 = df1(m)
          end if
          call eval_disc(fam(e, t), rot(e, t), ppv(:, e, t), aux(:, e, t), hf2(e), v2, hm2(e), v2m, df2(e), d2, &
                         lp, h1, h2, q1, q2)
          if (nh1(e, t) /= 0) then
            hf1(e) = h1; hm1(e) = q1; df1(e) = d2
          end if
          hf2(e) = h2; hm2(e) = q2
        end do
      end do
      do j = 1, d
        out(order(j), i) = min(max(hf2(j), 1.0e-10_dp), 1.0_dp - 1.0e-10_dp)
        outm(order(j), i) = min(max(hm2(j), 1.0e-10_dp), 1.0_dp - 1.0e-10_dp)
      end do
    end do
    !$omp end parallel do
  end subroutine vine_rosenblatt_disc

  ! Rosenblatt transform (vinecopulib Vinecop::rosenblatt), u/out: torch (n, d).
  subroutine fv_vine_rosenblatt(n, d, u, order, trunc, st, fam, rot, par, out, nthreads) &
      bind(C, name="fv_vine_rosenblatt")
    integer(c_int), value :: n, d, trunc, nthreads
    real(c_double), intent(in) :: u(d, n), par(3, d, d)
    integer(c_int), intent(in) :: order(d), st(d, d), fam(d, d), rot(d, d)
    real(c_double), intent(out) :: out(d, n)
    real(dp) :: axd(1, d, d)
    call vine_rosenblatt_x(n, d, u, order, trunc, st, fam, rot, par, 1, axd, out, nthreads)
  end subroutine fv_vine_rosenblatt

  subroutine vine_rosenblatt_x(n, d, u, order, trunc, st, fam, rot, par, na, aux, out, nthreads)
    integer, intent(in) :: n, d, trunc, nthreads, na
    real(dp), intent(in) :: u(d, n), par(3, d, d), aux(na, d, d)
    integer, intent(in) :: order(d), st(d, d), fam(d, d), rot(d, d)
    real(dp), intent(out) :: out(d, n)
    integer :: mn(d, d), nh1(d, d), nh2(d, d), i, j, t, e, m
    real(dp) :: hf1(d), hf2(d), u2, lp, h1, h2, ppv(NPP, d, d)
    call struct_arrays(d, trunc, st, mn, nh1, nh2)
    call prep_all(d, trunc, fam, par, ppv)
    !$omp parallel do private(i, j, t, e, m, hf1, hf2, u2, lp, h1, h2) num_threads(nthreads) schedule(static)
    do i = 1, n
      do j = 1, d
        hf2(j) = clampu(u(order(j), i))
      end do
      hf1 = hf2
      do t = 1, trunc
        do e = 1, d - t
          m = mn(e, t)
          if (m == st(e, t)) then
            u2 = hf2(m)
          else
            u2 = hf1(m)
          end if
          if (fam(e, t) == 0) then
            h1 = u2; h2 = hf2(e)
          else
            call eval_x(fam(e, t), rot(e, t), ppv(:, e, t), aux(:, e, t), hf2(e), u2, .true., lp, h1, h2)
          end if
          if (nh1(e, t) /= 0) hf1(e) = h1
          hf2(e) = h2
        end do
      end do
      do j = 1, d
        out(order(j), i) = clampu(hf2(j))
      end do
    end do
    !$omp end parallel do
  end subroutine vine_rosenblatt_x

  ! Inverse Rosenblatt transform (vinecopulib Vinecop::inverse_rosenblatt).
  subroutine fv_vine_inverse_rosenblatt(n, d, u, order, trunc, st, fam, rot, par, out, nthreads) &
      bind(C, name="fv_vine_inverse_rosenblatt")
    integer(c_int), value :: n, d, trunc, nthreads
    real(c_double), intent(in) :: u(d, n), par(3, d, d)
    integer(c_int), intent(in) :: order(d), st(d, d), fam(d, d), rot(d, d)
    real(c_double), intent(out) :: out(d, n)
    real(dp) :: axd(1, d, d)
    call vine_inverse_rosenblatt_x(n, d, u, order, trunc, st, fam, rot, par, 1, axd, out, nthreads)
  end subroutine fv_vine_inverse_rosenblatt

  subroutine vine_inverse_rosenblatt_x(n, d, u, order, trunc, st, fam, rot, par, na, aux, out, nthreads)
    integer, intent(in) :: n, d, trunc, nthreads, na
    real(dp), intent(in) :: u(d, n), par(3, d, d), aux(na, d, d)
    integer, intent(in) :: order(d), st(d, d), fam(d, d), rot(d, d)
    real(dp), intent(out) :: out(d, n)
    integer :: mn(d, d), nh1(d, d), nh2(d, d), i, j, t0, v, m, ts
    real(dp) :: hv(0:d, d), hf(0:d, d), u0, u1, lp, h1, h2, ppv(NPP, d, d)
    call struct_arrays(d, trunc, st, mn, nh1, nh2)
    call prep_all(d, trunc, fam, par, ppv)
    !$omp parallel do private(i, j, t0, v, m, ts, hv, hf, u0, u1, lp, h1, h2) num_threads(nthreads) &
    !$omp& schedule(static)
    do i = 1, n
      do j = 1, d
        hv(min(trunc, d - j), j) = clampu(u(order(j), i))
      end do
      hf(0, d) = hv(0, d)
      do v = d - 1, 1, -1
        ts = min(trunc - 1, d - v - 1)
        do t0 = ts, 0, -1
          m = mn(v, t0 + 1)
          u0 = hv(t0 + 1, v)
          if (m == st(v, t0 + 1)) then
            u1 = hv(t0, m)
          else
            u1 = hf(t0, m)
          end if
          if (fam(v, t0 + 1) == 0) then
            hv(t0, v) = u0
          else
            hv(t0, v) = hinv_x(fam(v, t0 + 1), rot(v, t0 + 1), ppv(:, v, t0 + 1), aux(:, v, t0 + 1), u0, u1, 2)
          end if
          if (nh1(v, t0 + 1) /= 0) then
            if (fam(v, t0 + 1) == 0) then
              hf(t0 + 1, v) = u1
            else
              call eval_x(fam(v, t0 + 1), rot(v, t0 + 1), ppv(:, v, t0 + 1), aux(:, v, t0 + 1), hv(t0, v), u1, &
                          .true., lp, h1, h2)
              hf(t0 + 1, v) = h1
            end if
          end if
        end do
      end do
      do j = 1, d
        out(order(j), i) = hv(0, j)
      end do
    end do
    !$omp end parallel do
  end subroutine vine_inverse_rosenblatt_x

  ! Pseudo-observations: average ranks / (n + 1), per column. x/out torch (n, d).
  subroutine fv_pseudo_obs(n, d, x, out, nthreads) bind(C, name="fv_pseudo_obs")
    integer(c_int), value :: n, d, nthreads
    real(c_double), intent(in) :: x(d, n)
    real(c_double), intent(out) :: out(d, n)
    real(dp), allocatable :: col(:), r(:)
    integer, allocatable :: idx(:), buf(:)
    integer :: j, i, k, l
    !$omp parallel do private(j, i, k, l, col, r, idx, buf) num_threads(nthreads) schedule(dynamic, 1)
    do j = 1, d
      allocate(col(n), r(n), idx(n), buf(n))
      col = x(j, :)
      do i = 1, n
        idx(i) = i
      end do
      call msort_by(col, idx, buf, 1, n)
      i = 1
      do while (i <= n)
        k = i
        do while (k < n)
          if (col(idx(k + 1)) /= col(idx(i))) exit
          k = k + 1
        end do
        do l = i, k
          r(idx(l)) = 0.5_dp * real(i + k, dp)     ! average rank of the tie block
        end do
        i = k + 1
      end do
      out(j, :) = r / real(n + 1, dp)
      deallocate(col, r, idx, buf)
    end do
    !$omp end parallel do
  end subroutine fv_pseudo_obs

  recursive subroutine msort_by(key, ix, bf, lo, hi)
    real(dp), intent(in) :: key(:)
    integer, intent(inout) :: ix(:), bf(:)
    integer, intent(in) :: lo, hi
    integer :: mid, p, q, k
    if (hi <= lo) return
    mid = (lo + hi) / 2
    call msort_by(key, ix, bf, lo, mid)
    call msort_by(key, ix, bf, mid + 1, hi)
    p = lo; q = mid + 1; k = lo
    do while (p <= mid .and. q <= hi)
      if (key(ix(q)) < key(ix(p))) then
        bf(k) = ix(q); q = q + 1
      else
        bf(k) = ix(p); p = p + 1
      end if
      k = k + 1
    end do
    do while (p <= mid)
      bf(k) = ix(p); p = p + 1; k = k + 1
    end do
    do while (q <= hi)
      bf(k) = ix(q); q = q + 1; k = k + 1
    end do
    ix(lo:hi) = bf(lo:hi)
  end subroutine msort_by

  ! ------------------------------------------------------ C API: Kendall tau
  ! tau-b via Knight (1966): O(n log n), tie-corrected.
  real(dp) function tau_b(n, x, y)
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n), y(n)
    integer, allocatable :: idx(:), buf(:)
    integer(8) :: n0, n1, n2, n3, sw, tcount
    integer :: i, j
    real(dp) :: den
    if (n < 2) then
      tau_b = 0.0_dp; return
    end if
    allocate(idx(n), buf(n))
    do i = 1, n
      idx(i) = i
    end do
    sw = 0
    call msort(idx, buf, 1, n, 1, sw)   ! lexicographic (x, y)
    n0 = int(n, 8) * int(n - 1, 8) / 2
    n1 = 0; n3 = 0
    i = 1
    do while (i <= n)
      j = i
      do while (j < n)
        if (x(idx(j + 1)) /= x(idx(i))) exit
        j = j + 1
      end do
      tcount = j - i + 1
      n1 = n1 + tcount * (tcount - 1) / 2
      i = j + 1
    end do
    i = 1
    do while (i <= n)
      j = i
      do while (j < n)
        if (x(idx(j + 1)) /= x(idx(i)) .or. y(idx(j + 1)) /= y(idx(i))) exit
        j = j + 1
      end do
      tcount = j - i + 1
      n3 = n3 + tcount * (tcount - 1) / 2
      i = j + 1
    end do
    sw = 0
    call msort(idx, buf, 1, n, 2, sw)   ! by y, counting discordant swaps
    n2 = 0
    i = 1
    do while (i <= n)
      j = i
      do while (j < n)
        if (y(idx(j + 1)) /= y(idx(i))) exit
        j = j + 1
      end do
      tcount = j - i + 1
      n2 = n2 + tcount * (tcount - 1) / 2
      i = j + 1
    end do
    den = sqrt(real(n0 - n1, dp)) * sqrt(real(n0 - n2, dp))
    if (den <= 0.0_dp) then
      tau_b = 0.0_dp
    else
      tau_b = real(n0 - n1 - n2 + n3 - 2_8 * sw, dp) / den
    end if
    tau_b = min(max(tau_b, -1.0_dp), 1.0_dp)
  contains
    logical function lt(p, q, key)
      integer, intent(in) :: p, q, key
      if (key == 1) then
        lt = x(p) < x(q) .or. (x(p) == x(q) .and. y(p) < y(q))
      else
        lt = y(p) < y(q)
      end if
    end function lt
    recursive subroutine msort(ix, bf, lo, hi, key, cnt)
      integer, intent(inout) :: ix(:), bf(:)
      integer, intent(in) :: lo, hi, key
      integer(8), intent(inout) :: cnt
      integer :: mid, p, q, k
      if (hi <= lo) return
      mid = (lo + hi) / 2
      call msort(ix, bf, lo, mid, key, cnt)
      call msort(ix, bf, mid + 1, hi, key, cnt)
      p = lo; q = mid + 1; k = lo
      do while (p <= mid .and. q <= hi)
        if (lt(ix(q), ix(p), key)) then
          bf(k) = ix(q); q = q + 1
          cnt = cnt + (mid - p + 1)
        else
          bf(k) = ix(p); p = p + 1
        end if
        k = k + 1
      end do
      do while (p <= mid)
        bf(k) = ix(p); p = p + 1; k = k + 1
      end do
      do while (q <= hi)
        bf(k) = ix(q); q = q + 1; k = k + 1
      end do
      ix(lo:hi) = bf(lo:hi)
    end subroutine msort
  end function tau_b

  real(c_double) function fv_kendall_tau_c(n, x, y) bind(C, name="fv_kendall_tau")
    integer(c_int), value :: n
    real(c_double), intent(in) :: x(n), y(n)
    fv_kendall_tau_c = tau_b(n, x, y)
  end function fv_kendall_tau_c

  ! data is torch (n, d) -> Fortran X(d, n); T is (d, d) symmetric.
  subroutine fv_tau_matrix(n, d, xx, tm, nthreads) bind(C, name="fv_tau_matrix")
    integer(c_int), value :: n, d, nthreads
    real(c_double), intent(in) :: xx(d, n)
    real(c_double), intent(out) :: tm(d, d)
    real(dp), allocatable :: cols(:, :)
    integer :: p, j, k, npair
    allocate(cols(n, d))
    cols = transpose(xx)
    npair = d * (d - 1) / 2
    tm = 0.0_dp
    do j = 1, d
      tm(j, j) = 1.0_dp
    end do
    !$omp parallel do private(p, j, k) num_threads(nthreads) schedule(dynamic, 1)
    do p = 0, npair - 1
      call pair_index(p, d, j, k)
      tm(j, k) = tau_b(n, cols(:, j), cols(:, k))
      tm(k, j) = tm(j, k)
    end do
    !$omp end parallel do
  end subroutine fv_tau_matrix

  pure subroutine pair_index(p, d, j, k)
    integer, intent(in) :: p, d
    integer, intent(out) :: j, k
    integer :: r, rowlen
    r = p; j = 1
    rowlen = d - 1
    do while (r >= rowlen)
      r = r - rowlen
      j = j + 1
      rowlen = rowlen - 1
    end do
    k = j + 1 + r
  end subroutine pair_index

  ! ======================================================================
  ! Vine structure + pair-copula selection (Dissmann et al. 2013), a port of
  ! torchvine's VinecopSelectorTorch / vinecopulib's VinecopSelector.
  ! ======================================================================

  pure real(dp) function digamma(x0)
    real(dp), intent(in) :: x0
    real(dp) :: x, r, f
    x = x0; r = 0.0_dp
    do while (x < 12.0_dp)
      r = r - 1.0_dp / x
      x = x + 1.0_dp
    end do
    f = 1.0_dp / (x * x)
    digamma = r + log(x) - 0.5_dp / x - f * (1.0_dp / 12.0_dp - f * (1.0_dp / 120.0_dp &
              - f * (1.0_dp / 252.0_dp - f * (1.0_dp / 240.0_dp - f / 132.0_dp))))
  end function digamma

  ! integral_0^a t / (e^t - 1) dt, exact series (machine precision):
  !  a < 1:  a - a^2/4 + sum_k B_2k a^(2k+1) / ((2k+1) (2k)!)   (Bernoulli numbers)
  !  a >= 1: pi^2/6 - sum_k e^(-k a) (a/k + 1/k^2)
  pure real(dp) function debye_int(a)
    real(dp), intent(in) :: a
    real(dp), parameter :: bc(10) = [1.0_dp / 36.0_dp, -1.0_dp / 3600.0_dp, 1.0_dp / 211680.0_dp, &
      -1.0_dp / 10886400.0_dp, 1.0_dp / 526901760.0_dp, -691.0_dp / (2730.0_dp * 6227020800.0_dp), &
      7.0_dp / (6.0_dp * 1307674368000.0_dp), -3617.0_dp / (510.0_dp * 355687428096000.0_dp), &
      43867.0_dp / (798.0_dp * 121645100408832000.0_dp), -174611.0_dp / (330.0_dp * 51090942171709440000.0_dp)]
    integer :: k
    real(dp) :: s, a2, p, e, ek, term
    if (a < 1.0_dp) then
      a2 = a * a; p = a; s = 0.0_dp
      do k = 1, 10
        p = p * a2
        s = s + bc(k) * p
      end do
      debye_int = a - 0.25_dp * a2 + s
    else
      e = exp(-a); ek = 1.0_dp; s = 0.0_dp
      do k = 1, 200
        ek = ek * e
        term = ek * (a / k + 1.0_dp / (real(k, dp) * k))
        s = s + term
        if (term < 1.0e-17_dp * s) exit
      end do
      debye_int = PI * PI / 6.0_dp - s
    end if
  end function debye_int

  ! trigamma psi'(x), x > 0 (recurrence + asymptotic series)
  pure real(dp) function trigamma(x0)
    real(dp), intent(in) :: x0
    real(dp) :: x, r, f
    x = x0; r = 0.0_dp
    do while (x < 12.0_dp)
      r = r + 1.0_dp / (x * x)
      x = x + 1.0_dp
    end do
    f = 1.0_dp / (x * x)
    trigamma = r + 1.0_dp / x + 0.5_dp * f + f / x * (1.0_dp / 6.0_dp - f * (1.0_dp / 30.0_dp - f * (1.0_dp / 42.0_dp &
               - f * (1.0_dp / 30.0_dp - f * 5.0_dp / 66.0_dp))))
  end function trigamma

  pure real(dp) function frank_tau(th)
    real(dp), intent(in) :: th
    real(dp) :: a
    a = abs(th)
    if (a < 1.0e-5_dp) then
      frank_tau = 0.0_dp; return
    end if
    frank_tau = 1.0_dp - 4.0_dp / a + 4.0_dp * debye_int(a) / (a * a)
    if (th < 0.0_dp) frank_tau = -frank_tau
  end function frank_tau

  pure real(dp) function joe_tau(th)
    real(dp), intent(in) :: th
    joe_tau = min(max(1.0_dp + 2.0_dp * joe_r(th), 0.0_dp), 0.999999_dp)
  end function joe_tau

  ! R(th) = (psi(2) - psi(2/th + 1)) / (2 - th); near th = 2 the removable singularity is
  ! expanded: psi(2) - psi(2 + e) = -psi'(2) e - psi''(2) e^2/2 - psi'''(2) e^3/6, e = (2-th)/th
  pure real(dp) function joe_r(th)
    real(dp), intent(in) :: th
    real(dp), parameter :: t1 = 0.6449340668482264_dp, t2 = -0.4041138063191886_dp, t3 = 0.4939394022668291_dp
    real(dp) :: e
    if (abs(th - 2.0_dp) < 1.0e-3_dp) then
      e = (2.0_dp - th) / th
      joe_r = -(t1 + e * (0.5_dp * t2 + e * t3 / 6.0_dp)) / th
    else
      joe_r = (digamma(2.0_dp) - digamma(2.0_dp / th + 1.0_dp)) / (2.0_dp - th)
    end if
  end function joe_r

  ! d tau / d theta in closed form (Frank: Debye; Joe: trigamma), theta > 0
  pure real(dp) function dtau(fam, th)
    integer, intent(in) :: fam
    real(dp), intent(in) :: th
    real(dp) :: x, h
    if (fam == 5) then
      dtau = 4.0_dp / (th * th) + 4.0_dp / (th * expm1_(th)) - 8.0_dp * debye_int(th) / th**3
    else
      h = th
      if (abs(h - 2.0_dp) < 1.0e-3_dp) h = merge(2.001_dp, 1.999_dp, h >= 2.0_dp)
      x = 2.0_dp / h + 1.0_dp
      dtau = 2.0_dp * (trigamma(x) * 2.0_dp / (h * h) * (2.0_dp - h) &
             + (digamma(2.0_dp) - digamma(x))) / (2.0_dp - h)**2
    end if
  end function dtau

  ! Invert a monotone-increasing tau(theta) on [lo, hi]: Newton with the analytic
  ! derivative, bracketed (bisection fallback); 4-6 steps instead of 60 bisections.
  ! Frank is odd in theta, so it is solved for |tau| on theta > 0.
  pure real(dp) function invert_tau(fam, tau0, lo, hi)
    integer, intent(in) :: fam
    real(dp), intent(in) :: tau0, lo, hi
    real(dp) :: a, b, x, xn, f, tau, sg
    integer :: it
    sg = 1.0_dp; tau = tau0; a = lo; b = hi
    if (fam == 5) then
      sg = sign(1.0_dp, tau0); tau = abs(tau0); a = 1.0e-6_dp; b = max(abs(lo), abs(hi))
    end if
    if (fam == 5) then
      x = min(max(2.0_dp * tau / max(1.0_dp - tau, 1.0e-3_dp) * 3.0_dp, a), b)     ! ~ 6 tau for small tau
    else
      x = min(max(1.0_dp / max(1.0_dp - tau, 1.0e-3_dp), a), b)
    end if
    do it = 1, 100
      if (fam == 5) then
        f = frank_tau(x) - tau
      else
        f = 1.0_dp + 2.0_dp * joe_r(x) - tau
      end if
      if (f < 0.0_dp) then
        a = x
      else
        b = x
      end if
      xn = x - f / dtau(fam, x)
      if (.not. (xn > a .and. xn < b)) xn = 0.5_dp * (a + b)
      if (abs(xn - x) <= 1.0e-14_dp * x .or. b - a <= 1.0e-14_dp * x) then
        x = xn; exit
      end if
      x = xn
    end do
    invert_tau = sg * x
  end function invert_tau

  ! Bicop::tau_to_parameters (one-parameter families; student gets nu = 4).
  pure subroutine tau2par(fam, tau, p)
    integer, intent(in) :: fam
    real(dp), intent(in) :: tau
    real(dp), intent(out) :: p(3)
    real(dp) :: at
    p = 0.0_dp
    at = abs(tau)
    select case (fam)
    case (1)
      p(1) = sin(tau * PI / 2.0_dp)
    case (2)
      p(1) = sin(tau * PI / 2.0_dp); p(2) = 4.0_dp
    case (3)
      p(1) = min(max(2.0_dp * at / max(1.0e-12_dp, 1.0_dp - at), 1.0e-10_dp), 28.0_dp)
    case (4)
      p(1) = min(max(1.0_dp / max(1.0e-12_dp, 1.0_dp - at), 1.0_dp), 50.0_dp)
    case (5)
      if (abs(tau) >= 1.0e-12_dp) p(1) = invert_tau(5, tau, -35.0_dp + 1.0e-6_dp, 35.0_dp - 1.0e-5_dp)
    case (6)
      p(1) = invert_tau(6, at, 1.0_dp, 30.0_dp)
    end select
  end subroutine tau2par

  pure subroutine par_bounds(fam, lb, ub)
    integer, intent(in) :: fam
    real(dp), intent(out) :: lb(3), ub(3)
    lb = 0.0_dp; ub = 0.0_dp
    if (gen_is(fam)) then
      call gen_bounds(fam, lb, ub); return
    end if
    select case (fam)
    case (1); lb(1) = -1.0_dp; ub(1) = 1.0_dp
    case (2); lb(1:2) = [-1.0_dp, 2.0_dp]; ub(1:2) = [1.0_dp, 30.0_dp]
    case (3); lb(1) = 1.0e-10_dp; ub(1) = 28.0_dp
    case (4); lb(1) = 1.0_dp; ub(1) = 50.0_dp
    case (5); lb(1) = -35.0_dp; ub(1) = 35.0_dp
    case (6); lb(1) = 1.0_dp; ub(1) = 30.0_dp
    case (7); lb(1:2) = [0.0_dp, 1.0_dp]; ub(1:2) = [7.0_dp, 7.0_dp]
    case (8); lb(1:2) = [1.0_dp, 1.0_dp]; ub(1:2) = [6.0_dp, 8.0_dp]
    case (9); lb(1:2) = [1.0_dp, 0.01_dp]; ub(1:2) = [6.0_dp, 25.0_dp]
    case (10); lb(1:2) = [1.0_dp, 1.0e-4_dp]; ub(1:2) = [8.0_dp, 1.0_dp]
    case (11); lb = [0.0_dp, 0.0_dp, 1.0_dp]; ub = [1.0_dp, 1.0_dp, 60.0_dp]
    case (12); lb(1) = 1.0e-3_dp; ub(1) = 1.0e3_dp
    case (13); lb(1) = -1.0_dp; ub(1) = 1.0_dp
    case (14); lb(1) = -1.0_dp; ub(1) = 0.9999_dp
    case (15, 16); lb(1) = 0.01_dp; ub(1) = 30.0_dp
    case (17); lb(1:2) = [0.01_dp, 0.01_dp]; ub(1:2) = [30.0_dp, 30.0_dp]
    case (18, 19, 20); lb(1:2) = [1.0_dp, 0.01_dp]; ub(1:2) = [10.0_dp, 30.0_dp]
    case (21); lb(1:2) = [0.01_dp, 0.0_dp]; ub(1:2) = [30.0_dp, 1.0_dp]
    case (22); lb(1:2) = [0.01_dp, 0.01_dp]; ub(1:2) = [15.0_dp, 15.0_dp]
    case (23); lb(1:2) = [-0.95_dp, 1.0_dp]; ub(1:2) = [0.99_dp, 30.0_dp]
    end select
  end subroutine par_bounds

  ! Bicop::_adjust_bounds_like_vinecopulib: narrow 1-parameter bounds to tau +- 0.1.
  pure subroutine adjust_bounds(fam, tau, lb, ub)
    integer, intent(in) :: fam
    real(dp), intent(in) :: tau
    real(dp), intent(inout) :: lb(3), ub(3)
    real(dp) :: lo, hi, pl(3), ph(3)
    select case (fam)
    case (1, 3, 4, 5, 6)
      if (fam == 1 .or. fam == 5) then
        lo = max(tau - 0.1_dp, -0.99_dp); hi = min(tau + 0.1_dp, 0.99_dp)
      else
        lo = max(abs(tau) - 0.1_dp, 1.0e-10_dp); hi = min(abs(tau) + 0.1_dp, 0.95_dp)
      end if
      call tau2par(fam, lo, pl); call tau2par(fam, hi, ph)
      lb(1) = max(lb(1), pl(1)); ub(1) = min(ub(1), ph(1))
    case (2)
      lo = max(tau - 0.1_dp, -0.99_dp); hi = min(tau + 0.1_dp, 0.99_dp)
      lb(1:2) = [max(sin(lo * PI / 2.0_dp), -0.999_dp), 2.01_dp]
      ub(1:2) = [min(sin(hi * PI / 2.0_dp), 0.999_dp), 30.0_dp]
    case (11)
      ! vinecopulib ParBicop::adjust_parameters_bounds for tawn
      lb = [0.3_dp, 0.3_dp, 1.5_dp]; ub = [1.0_dp, 1.0_dp, 7.0_dp]
    end select
  end subroutine adjust_bounds

  pure real(dp) function winsor(tau)
    real(dp), intent(in) :: tau
    winsor = sign(min(max(abs(tau), 0.01_dp), 0.9_dp), tau)
  end function winsor

  pure integer function npar_of(fam)
    integer, intent(in) :: fam
    select case (fam)
    case (0, 24, 25); npar_of = 0      ! nonparametric: effective df is aux(2)
    case (2, 7, 8, 9, 10, 17:23); npar_of = 2
    case (11); npar_of = 3
    case (10000:); npar_of = gen_npar(fam)
    case default; npar_of = 1
    end select
  end function npar_of

  pure logical function can_rotate(fam)
    integer, intent(in) :: fam
    ! rotationless families cover both signs of dependence themselves
    can_rotate = .not. (fam == 0 .or. fam == 1 .or. fam == 2 .or. fam == 5 .or. (fam >= 12 .and. fam <= 14) &
                     .or. fam_is_np(fam))
  end function can_rotate

  ! nonparametric pair copulas (24 tll, 25 bernstein): their estimate lives in an aux vector.
  ! Family codes are open-ended (generic families start at 10000), so never test fam >= 24.
  pure elemental logical function fam_is_np(fam)
    integer, intent(in) :: fam
    fam_is_np = fam == 24 .or. fam == 25
  end function fam_is_np

  ! ======================================================================
  ! Generic families (src/families_gen.inc): fam = 10000 k + 100 og + ie
  !   k  = 1 plain, 2 Khoudraji asymmetrisation (two more parameters a, b in [0, 1])
  !   og = outer generator: 0 none (exp(-s)); 1 positive stable, 2 gamma, 3 Sibuya,
  !        4 log series, 5 tilted positive stable, 6 generalised Sibuya, 7 negative binomial;
  !        11..26 composition po.pin = psi_po(-log zeta_pin), og = 10 + 4 (po-1) + pin with
  !        po, pin in {1 positive stable, 2 gamma, 3 Sibuya, 4 log series}
  !   ie = inner: 0 independence (Archimedean); 1 logistic, 2 negative logistic, 3 Husler-Reiss,
  !        4 symmetric mixed, 5 asymmetric logistic, 6/7 its types 1/2, 8 asymmetric negative
  !        logistic, 9/10 its types 1/2, 11 asymmetric mixed, 12 Khoudraji-Husler-Reiss,
  !        13 Khoudraji-mixed (Pickands functions: extreme-value / Archimax)
  ! Parameters: outer first, then inner, then (k = 2) a, b; at most 3 in all.
  ! ======================================================================
  pure elemental logical function gen_is(fam)
    integer, intent(in) :: fam
    gen_is = fam >= 10000
  end function gen_is

  pure elemental integer function gen_og(fam)
    integer, intent(in) :: fam
    gen_og = mod(fam, 10000) / 100
  end function gen_og

  pure elemental integer function gen_ie(fam)
    integer, intent(in) :: fam
    gen_ie = mod(fam, 100)
  end function gen_ie

  pure elemental logical function gen_kh(fam)
    integer, intent(in) :: fam
    gen_kh = fam / 10000 == 2
  end function gen_kh

  pure elemental integer function gen_npo(og)
    integer, intent(in) :: og
    select case (og)
    case (0); gen_npo = 0
    case (1:4); gen_npo = 1
    case default; gen_npo = 2
    end select
  end function gen_npo

  pure elemental integer function gen_npi(ie)
    integer, intent(in) :: ie
    select case (ie)
    case (0); gen_npi = 0
    case (1:4); gen_npi = 1
    case (5, 8, 12, 13); gen_npi = 3
    case default; gen_npi = 2
    end select
  end function gen_npi

  pure elemental integer function gen_npar(fam)
    integer, intent(in) :: fam
    gen_npar = gen_npo(gen_og(fam)) + gen_npi(gen_ie(fam)) + merge(2, 0, gen_kh(fam))
  end function gen_npar

  ! a well-defined, identifiable code (degenerate products excluded)
  pure elemental logical function gen_valid(fam)
    integer, intent(in) :: fam
    integer :: og, ie, k
    gen_valid = .false.
    if (fam < 10000 .or. fam >= 30000) return
    k = fam / 10000; og = gen_og(fam); ie = gen_ie(fam)
    if (.not. (og <= 7 .or. (og >= 11 .and. og <= 26))) return
    if (ie > 13) return
    if (og == 11 .or. og == 21) return                 ! PS.PS = PS, Sibuya.Sibuya = Sibuya
    if (og == 14 .or. og == 22 .or. og == 26) return   ! log-series inner under PS / Sibuya / log-series:
                                                       ! not accurate in double precision (Joe 2014 A.1.2)
    if (ie > 0 .and. og >= 11 .and. mod(og - 11, 4) + 1 == 2) return   ! gamma-inner generator (BB2 / BB3 type)
                                                       ! with a Pickands inner: fails the consistency tests
    if (og == 0 .and. ie == 0) return                  ! independence
    if (og == 1 .and. ie == 1) return                  ! PS outer with logistic inner = logistic
    if (k == 2 .and. ie >= 5) return                   ! already asymmetric
    gen_valid = gen_npar(fam) >= 1 .and. gen_npar(fam) <= 3
  end function gen_valid

  ! the same copula is available under another code (a hand-coded family or another generic
  ! code); kept for checks, left out of the generic family sets
  pure elemental logical function gen_dup(fam)
    integer, intent(in) :: fam
    integer :: og, ie
    og = gen_og(fam); ie = gen_ie(fam)
    if (gen_kh(fam)) then
      gen_dup = og == 0 .or. (og == 1 .and. ie == 0)    ! Khoudraji-EV = ie 5 / 8 / 12 / 13; Khoudraji-Gumbel = Tawn
    else if (ie == 0) then                              ! Archimedean: Gumbel ... BB10, BB1-3, BB6-7
      gen_dup = og <= 7 .or. any(og == [12, 15, 16, 19, 20])
    else if (og == 0) then                              ! Gumbel, Galambos, Husler-Reiss, Tawn
      gen_dup = any(ie == [1, 2, 3, 5])
    else                                                ! BB5, BB4; psi with logistic = psi.PS
      gen_dup = (ie == 2 .and. (og == 1 .or. og == 2)) .or. (ie == 1 .and. (og <= 4 .or. any(og == [15, 19, 23])))
    end if
  end function gen_dup

  ! parameter box (outer, inner, Khoudraji in this order)
  pure subroutine gen_bounds(fam, lb, ub)
    integer, intent(in) :: fam
    real(dp), intent(out) :: lb(3), ub(3)
    integer :: og, ie, j, no, ni, po, pin
    real(dp) :: b1(2, 4), b2(2, 4)
    ! one-parameter LTs alone, and as the outer / inner part of a composition
    b1 = reshape([1.0_dp, 20.0_dp, 1.0e-3_dp, 20.0_dp, 1.0_dp, 20.0_dp, 1.0e-3_dp, 35.0_dp], [2, 4])
    b2 = reshape([1.0_dp, 10.0_dp, 1.0e-3_dp, 15.0_dp, 1.0_dp, 10.0_dp, 1.0e-3_dp, 30.0_dp], [2, 4])
    lb = 0.0_dp; ub = 0.0_dp
    og = gen_og(fam); ie = gen_ie(fam); no = gen_npo(og); ni = gen_npi(ie)
    select case (og)
    case (1:4); lb(1) = b1(1, og); ub(1) = b1(2, og)
    case (5); lb(1:2) = [1.0_dp, 0.01_dp]; ub(1:2) = [10.0_dp, 25.0_dp]
    case (6); lb(1:2) = [1.0_dp, 1.0e-3_dp]; ub(1:2) = [8.0_dp, 1.0_dp]
    case (7); lb(1:2) = [0.01_dp, 0.0_dp]; ub(1:2) = [30.0_dp, 0.99_dp]
    case (11:26)
      po = (og - 11) / 4 + 1; pin = mod(og - 11, 4) + 1
      lb(1:2) = [b2(1, po), b2(1, pin)]; ub(1:2) = [b2(2, po), b2(2, pin)]
    end select
    j = no
    select case (ie)
    case (1); lb(j + 1) = 1.0_dp; ub(j + 1) = 20.0_dp
    case (2, 3); lb(j + 1) = 0.01_dp; ub(j + 1) = 20.0_dp
    case (4); lb(j + 1) = 0.0_dp; ub(j + 1) = 1.0_dp
    case (5); lb(j + 1:j + 3) = [0.01_dp, 0.01_dp, 1.0_dp]; ub(j + 1:j + 3) = [1.0_dp, 1.0_dp, 20.0_dp]
    case (6, 7); lb(j + 1:j + 2) = [1.0_dp, 0.01_dp]; ub(j + 1:j + 2) = [20.0_dp, 1.0_dp]
    case (8); lb(j + 1:j + 3) = [0.01_dp, 0.01_dp, 0.01_dp]; ub(j + 1:j + 3) = [1.0_dp, 1.0_dp, 20.0_dp]
    case (9, 10); lb(j + 1:j + 2) = [0.01_dp, 0.01_dp]; ub(j + 1:j + 2) = [20.0_dp, 1.0_dp]
    case (11); lb(j + 1:j + 2) = [0.0_dp, 0.0_dp]; ub(j + 1:j + 2) = [1.0_dp, 1.0_dp]
    case (12); lb(j + 1:j + 3) = [0.01_dp, 0.01_dp, 0.01_dp]; ub(j + 1:j + 3) = [20.0_dp, 1.0_dp, 1.0_dp]
    case (13); lb(j + 1:j + 3) = [0.0_dp, 0.01_dp, 0.01_dp]; ub(j + 1:j + 3) = [1.0_dp, 1.0_dp, 1.0_dp]
    end select
    if (gen_kh(fam)) then
      lb(no + ni + 1:no + ni + 2) = 0.01_dp; ub(no + ni + 1:no + ni + 2) = 1.0_dp
    end if
  end subroutine gen_bounds

  ! exchangeable: h2 is h1 with the arguments swapped
  pure elemental logical function gen_exch(fam)
    integer, intent(in) :: fam
    gen_exch = .not. gen_kh(fam) .and. (gen_ie(fam) <= 4)
  end function gen_exch

  ! transpose a generic family in place, C(u, v) -> C(v, u) (vine edge flip). The outer generator
  ! is exchangeable, so only the inner Pickands function A(w) -> A(1-w) or the Khoudraji
  ! exponents change: asymmetric (negative) logistic swap psi1 / psi2, their types 1 / 2 swap
  ! codes, the asymmetric mixed model maps (th, phi) -> (th + 3 phi, -phi), Khoudraji swaps a / b
  pure subroutine gen_flip(fam, par)
    integer, intent(inout) :: fam
    real(dp), intent(inout) :: par(3)
    integer :: og, ie, no, ni, j
    real(dp) :: th, ph, th2, ph2
    og = gen_og(fam); ie = gen_ie(fam); no = gen_npo(og); ni = gen_npi(ie)
    if (gen_kh(fam)) then
      j = no + ni + 1
      par(j:j + 1) = par([j + 1, j])
      return
    end if
    select case (ie)
    case (5, 8, 12, 13)                      ! (psi1, psi2, th) / (lambda, a, b): swap the pair
      j = merge(no + 1, no + 2, ie == 5 .or. ie == 8)
      par(j:j + 1) = par([j + 1, j])
    case (6, 9); fam = fam + 1               ! type 1 <-> type 2, same (th, psi)
    case (7, 10); fam = fam - 1
    case (11)                                ! (th, kappa): phi = -th/3 + kappa ((1-th)/2 + th/3)
      th = par(no + 1); ph = -th / 3.0_dp + par(no + 2) * ((1.0_dp - th) / 2.0_dp + th / 3.0_dp)
      th2 = th + 3.0_dp * ph; ph2 = -ph
      par(no + 1) = th2
      par(no + 2) = (ph2 + th2 / 3.0_dp) / max((1.0_dp - th2) / 2.0_dp + th2 / 3.0_dp, 1.0e-300_dp)
    end select
  end subroutine gen_flip

  ! tail dependence (for family preselection): lower from a gamma outer or inner part (LT with
  ! polynomial decay), upper from a positive stable / Sibuya part (LT with 1 - psi ~ s^(1/th))
  ! or from any extreme-value inner
  pure elemental logical function gen_lt(fam)
    integer, intent(in) :: fam
    integer :: og
    og = gen_og(fam)
    gen_lt = og == 2 .or. (og >= 11 .and. ((og - 11) / 4 + 1 == 2 .or. mod(og - 11, 4) + 1 == 2))
  end function gen_lt

  pure elemental logical function gen_ut(fam)
    integer, intent(in) :: fam
    integer :: og, po, pin
    og = gen_og(fam)
    gen_ut = gen_ie(fam) > 0 .or. og == 1 .or. og == 3
    if (og >= 11) then
      po = (og - 11) / 4 + 1; pin = mod(og - 11, 4) + 1
      gen_ut = gen_ut .or. po == 1 .or. po == 3 .or. pin == 1 .or. pin == 3
    end if
  end function gen_ut

  ! every valid generic code, optionally only those without a duplicate elsewhere
  pure function gen_codes(new_only) result(cs)
    logical, intent(in) :: new_only
    integer, allocatable :: cs(:)
    integer :: k, og, ie, f, n
    integer :: buf(2 * 27 * 14)
    n = 0
    do k = 1, 2
      do og = 0, 26
        do ie = 0, 13
          f = 10000 * k + 100 * og + ie
          if (.not. gen_valid(f)) cycle
          if (new_only .and. gen_dup(f)) cycle
          n = n + 1; buf(n) = f
        end do
      end do
    end do
    cs = buf(1:n)
  end function gen_codes

  ! name: the literature name where the copula has one, otherwise built from its parts
  pure function gen_name(fam) result(nm)
    integer, intent(in) :: fam
    character(len=48) :: nm
    character(len=10), parameter :: PN(4) = [character(len=10) :: 'PS', 'gamma', 'Sibuya', 'logseries']
    character(len=10), parameter :: GN(7) = [character(len=10) :: 'PS', 'gamma', 'Sibuya', 'logseries', &
                                             'tiltedPS', 'genSibuya', 'negbin']
    character(len=18), parameter :: IN(13) = [character(len=18) :: 'logistic', 'neglogistic', &
      'HuslerReiss', 'mixed', 'asymlogistic', 'asymlogistic1', 'asymlogistic2', 'asymneglogistic', &
      'asymneglogistic1', 'asymneglogistic2', 'asymmixed', 'Khoudraji-HR', 'Khoudraji-mixed']
    character(len=24) :: g
    integer :: og, ie
    og = gen_og(fam); ie = gen_ie(fam); nm = ''
    if (.not. gen_valid(fam)) then
      nm = 'invalid'; return
    end if
    select case (og)
    case (0); g = ''
    case (1:7); g = GN(og)
    case default; g = trim(PN((og - 11) / 4 + 1)) // '.' // trim(PN(mod(og - 11, 4) + 1))
    end select
    if (ie == 0) then                                  ! Archimedean
      select case (og)
      case (1); nm = 'Gumbel'
      case (2); nm = 'Clayton'
      case (3); nm = 'Joe'
      case (4); nm = 'Frank'
      case (5); nm = 'BB9'
      case (6); nm = 'BB8'
      case (7); nm = 'BB10'
      case (12); nm = 'BB3'
      case (15); nm = 'BB1'
      case (16); nm = 'BB2'
      case (19); nm = 'BB6'
      case (20); nm = 'BB7'
      case default; nm = 'Arch(' // trim(g) // ')'
      end select
    else if (og == 0) then                             ! extreme value
      select case (ie)
      case (1); nm = 'Gumbel-EV'
      case (2); nm = 'Galambos'
      case (3); nm = 'HuslerReiss'
      case (4); nm = 'Tawn-mixed'
      case (5); nm = 'Tawn'
      case (6); nm = 'Tawn-type1'
      case (7); nm = 'Tawn-type2'
      case (8); nm = 'Joe-asymneglogistic'
      case (9); nm = 'asymneglogistic-type1'
      case (10); nm = 'asymneglogistic-type2'
      case (11); nm = 'Tawn-asymmixed'
      case default; nm = IN(ie)
      end select
    else if (ie == 2 .and. og == 2) then
      nm = 'BB4'
    else if (ie == 2 .and. og == 1) then
      nm = 'BB5'
    else
      nm = 'Archimax(' // trim(g) // '|' // trim(IN(ie)) // ')'
    end if
    if (gen_kh(fam)) nm = 'Khoudraji-' // trim(nm)
  end function gen_name

  pure real(dp) function pearson(n, x, y)
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n), y(n)
    real(dp) :: mx, my, sxy, sxx, syy
    if (n < 2) then
      pearson = 0.0_dp; return
    end if
    mx = sum(x) / n; my = sum(y) / n
    sxy = sum((x - mx) * (y - my)); sxx = sum((x - mx)**2); syy = sum((y - my)**2)
    pearson = sxy / max(sqrt(max(sxx * syy, 0.0_dp)), TINYV)
    pearson = min(max(pearson, -1.0_dp), 1.0_dp)
  end function pearson

  ! _estimate_tau_for_fit: 0 auto (normal-scores approx when n >= 800), 1 exact, 2 approx.
  real(dp) function tau_est(n, a, b, mode)
    integer, intent(in) :: n, mode
    real(dp), intent(in) :: a(n), b(n)
    real(dp), allocatable :: za(:), zb(:)
    integer :: i
    if (mode == 1 .or. (mode == 0 .and. n < 800)) then
      tau_est = tau_b(n, a, b); return
    end if
    allocate(za(n), zb(n))
    do i = 1, n
      za(i) = qnorm(clampu(a(i))); zb(i) = qnorm(clampu(b(i)))
    end do
    tau_est = (2.0_dp / PI) * asin(min(max(pearson(n, za, zb), -0.999999_dp), 0.999999_dp))
  end function tau_est

  ! ranks 1..n (ties broken by stable order), via merge sort
  subroutine rank_vec(n, x, r)
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n)
    real(dp), intent(out) :: r(n)
    integer, allocatable :: idx(:), buf(:)
    integer :: i
    allocate(idx(n), buf(n))
    do i = 1, n
      idx(i) = i
    end do
    call msort_val(idx, buf, 1, n)
    do i = 1, n
      r(idx(i)) = real(i, dp)
    end do
  contains
    recursive subroutine msort_val(ix, bf, lo, hi)
      integer, intent(inout) :: ix(:), bf(:)
      integer, intent(in) :: lo, hi
      integer :: mid, p, q, k
      if (hi <= lo) return
      mid = (lo + hi) / 2
      call msort_val(ix, bf, lo, mid)
      call msort_val(ix, bf, mid + 1, hi)
      p = lo; q = mid + 1; k = lo
      do while (p <= mid .and. q <= hi)
        if (x(ix(q)) < x(ix(p))) then
          bf(k) = ix(q); q = q + 1
        else
          bf(k) = ix(p); p = p + 1
        end if
        k = k + 1
      end do
      do while (p <= mid)
        bf(k) = ix(p); p = p + 1; k = k + 1
      end do
      do while (q <= hi)
        bf(k) = ix(q); q = q + 1; k = k + 1
      end do
      ix(lo:hi) = bf(lo:hi)
    end subroutine msort_val
  end subroutine rank_vec

  ! Tree criterion: 0 tau, 1 rho (Spearman), 2 joe, 3 hoeffd (binned, as torchvine),
  ! 4 Chatterjee xi (symmetric), 5 mutual information (KSG), 6 Ferreira tail coefficient.
  real(dp) function criterion(n, x, y, kind)
    integer, intent(in) :: n, kind
    real(dp), intent(in) :: x(n), y(n)
    real(dp), allocatable :: rx(:), ry(:)
    integer :: i
    real(dp) :: w, r
    if (n <= 10) then
      criterion = 0.0_dp; return
    end if
    select case (kind)
    case (1, 3)
      allocate(rx(n), ry(n))
      call rank_vec(n, x, rx); call rank_vec(n, y, ry)
      criterion = crit_ranks(n, rx, ry, kind, .false., x); return
    case (2)
      allocate(rx(n), ry(n))
      do i = 1, n
        rx(i) = qnorm(clampu(x(i))); ry(i) = qnorm(clampu(y(i)))
      end do
      r = pearson(n, rx, ry)
      w = -0.5_dp * log(1.0_dp - min(max(r * r, 0.0_dp), 1.0_dp - 1.0e-15_dp))
    case (4)
      w = xi_sym(n, x, y)
    case (5)
      w = mi_ksg(n, x, y)
    case (6)
      w = ferreira_tail(n, x, y)
    case default
      w = tau_b(n, x, y)
    end select
    if (w /= w) w = 0.0_dp
    criterion = abs(w)
  end function criterion

  ! the rank-based tree criteria from the ranks rx, ry of a pair (rank_vec): kind 1 Spearman's rho
  ! (weighted Pearson of the ranks when usew), kind 3 Hoeffding's D on the binned empirical copula.
  ! The vine ranks each pseudo-observation column once per tree and scores every candidate from them.
  real(dp) function crit_ranks(n, rx, ry, kind, usew, w) result(c)
    integer, intent(in) :: n, kind
    real(dp), intent(in) :: rx(n), ry(n), w(:)
    logical, intent(in) :: usew
    real(dp), allocatable :: cdf(:, :)
    integer, allocatable :: bx(:), by(:)
    integer :: i, j, k, nb
    real(dp) :: r, ux, uy
    c = 0.0_dp
    if (n <= 10) return
    if (kind == 1) then
      if (usew) then
        r = wpearson(n, rx, ry, w)
      else
        r = pearson(n, rx, ry)
      end if
    else
      if (n < 20) return
      nb = max(32, min(128, nint(sqrt(real(n, dp)))))
      allocate(cdf(nb, nb), bx(n), by(n))
      cdf = 0.0_dp
      do i = 1, n
        bx(i) = min(max(int((rx(i) - 0.5_dp) / n * nb), 0), nb - 1) + 1
        by(i) = min(max(int((ry(i) - 0.5_dp) / n * nb), 0), nb - 1) + 1
        cdf(bx(i), by(i)) = cdf(bx(i), by(i)) + 1.0_dp
      end do
      do j = 2, nb
        cdf(j, :) = cdf(j, :) + cdf(j - 1, :)
      end do
      do k = 2, nb
        cdf(:, k) = cdf(:, k) + cdf(:, k - 1)
      end do
      cdf = cdf / n
      r = 0.0_dp
      do i = 1, n
        ux = (rx(i) - 0.5_dp) / n; uy = (ry(i) - 0.5_dp) / n
        r = r + (cdf(bx(i), by(i)) - ux * uy)**2
      end do
      r = 30.0_dp * r / n
    end if
    if (r /= r) r = 0.0_dp
    c = abs(r)
  end function crit_ranks

  ! ---------------------------------------------------- per-edge family candidates
  ! ci: 1 selcrit (0 loglik,1 aic,2 bic,3 mbic), 2 method (0 mle,1 itau),
  !     3 preselect, 4 allow_rot, 5 aggressive_fast_mle, 6 tau_estimation
  ! Builds the candidate list (vinecopulib create_candidate_bicops + preselect)
  ! and each candidate's fit spec: method cm (-1 indep, 0 eval at x0, 1 MLE,
  ! 2 student with rho fixed), start cx0 and bounds clb/cub.
  subroutine edge_cands(n, a, b, nf, famset, ci, has_tau, tau_pre, nc, cf, cr, cm, cx0, clb, cub)
    integer, intent(in) :: n, nf, famset(nf), ci(6)
    real(dp), intent(in) :: a(n), b(n), tau_pre
    logical, intent(in) :: has_tau
    integer, intent(out) :: nc, cf(:), cr(:), cm(:)
    real(dp), intent(out) :: cx0(:, :), clb(:, :), cub(:, :)
    integer :: k, f, r, wr(2), i, m, j
    real(dp) :: tau0, tau, c1, c2, cdiff
    real(dp), allocatable :: z1(:), z2(:), q1(:), q2(:)
    logical :: keep, lt, ut, is9018, isbb, itau_ok
    if (has_tau) then
      tau0 = tau_pre
    else
      tau0 = tau_est(n, a, b, ci(6))
    end if
    if (tau0 > 0.0_dp) then
      wr = [0, 180]
    else
      wr = [90, 270]
    end if
    cdiff = 0.0_dp
    if (ci(3) /= 0) then
      allocate(z1(n), z2(n), q1(n), q2(n))
      do i = 1, n
        z1(i) = qnorm(clampu(a(i))); z2(i) = qnorm(clampu(b(i)))
      end do
      c1 = quad_cor(.true.); c2 = quad_cor(.false.)
      cdiff = c1 - c2
    end if
    nc = 0
    do k = 1, nf
      f = famset(k)
      itau_ok = (f >= 0 .and. f <= 6)
      if (ci(2) == 1 .and. .not. itau_ok) cycle
      do j = 1, 2
        if (.not. can_rotate(f)) then
          if (j == 2) exit
          r = 0
        else if (ci(4) /= 0) then
          r = wr(j)
        else
          if (j == 2 .or. tau0 <= 0.0_dp) exit
          r = 0
        end if
        if (ci(3) /= 0) then
          ! vinecopulib preselect_family
          isbb = (f >= 7 .and. f <= 10)
          if (.not. can_rotate(f)) then
            keep = .not. (abs(cdiff) > 0.3_dp .and. f == 5)
          else
            keep = .false.
            if (isbb .or. gen_is(f)) keep = (tau0 > 0.0_dp .and. (r == 0 .or. r == 180)) .or. &
                             (tau0 < 0.0_dp .and. (r == 90 .or. r == 270))
            is9018 = (r == 90 .or. r == 180)
            lt = (f == 3 .or. f == 7 .or. f == 9 .or. f == 11 .or. f == 17 .or. f == 18 .or. f == 22)
            ut = (f == 4 .or. f == 6 .or. isbb .or. f == 11 .or. f == 15 .or. f == 16 .or. f == 18 .or. f == 19 &
                  .or. f == 22 .or. f == 23) .or. (gen_is(f) .and. gen_ut(f))
            lt = lt .or. (gen_is(f) .and. gen_lt(f))
            if (cdiff > 0.05_dp) then
              keep = keep .or. (lt .and. is9018) .or. (ut .and. .not. is9018)
            else if (cdiff < -0.05_dp) then
              keep = keep .or. (lt .and. .not. is9018) .or. (ut .and. is9018)
            else
              keep = keep .or. (tau0 > 0.0_dp .and. (r == 0 .or. r == 180)) .or. &
                     (tau0 < 0.0_dp .and. (r == 90 .or. r == 270))
            end if
          end if
          if (.not. keep) cycle
        end if
        if (r == 90 .or. r == 270) then
          tau = -tau0
        else
          tau = tau0
        end if
        nc = nc + 1
        cf(nc) = f; cr(nc) = r
        call cand_spec(f, tau, ci, cm(nc), cx0(:, nc), clb(:, nc), cub(:, nc))
      end do
    end do
  contains
    real(dp) function quad_cor(first)
      logical, intent(in) :: first
      integer :: jj
      m = 0
      do jj = 1, n
        if (tau0 > 0.0_dp) then
          if (first) then
            keep = z1(jj) > 0.0_dp .and. z2(jj) > 0.0_dp
          else
            keep = z1(jj) < 0.0_dp .and. z2(jj) < 0.0_dp
          end if
        else
          if (first) then
            keep = z1(jj) < 0.0_dp .and. z2(jj) > 0.0_dp
          else
            keep = z1(jj) > 0.0_dp .and. z2(jj) < 0.0_dp
          end if
        end if
        if (keep) then
          m = m + 1; q1(m) = z1(jj); q2(m) = z2(jj)
        end if
      end do
      ! vinecopulib get_c1c2 uses the first count-1 rows
      if (m <= 1) then
        quad_cor = 0.0_dp
      else
        quad_cor = pearson(m - 1, q1(1:m - 1), q2(1:m - 1))
      end if
    end function quad_cor
  end subroutine edge_cands

  ! Fit spec of one candidate (ParBicop::fit logic): method cm (-1 indep, 0 eval
  ! at x0, 1 MLE, 2 student with rho fixed), start x0 and bounds.  tau is the
  ! Kendall tau of the ROTATED data.
  pure subroutine cand_spec(f, tau, ci, cm, x0, lb, ub)
    integer, intent(in) :: f, ci(6)
    real(dp), intent(in) :: tau
    integer, intent(out) :: cm
    real(dp), intent(out) :: x0(3), lb(3), ub(3)
    call par_bounds(f, lb, ub)
    x0 = 0.0_dp
    if (f == 0) then
      cm = -1
    else if (fam_is_np(f)) then
      cm = 3                       ! nonparametric (np_fit)
    else if (f >= 12) then
      x0 = 0.5_dp * (lb + ub)
      cm = 1
    else if (ci(2) == 1 .or. (ci(5) /= 0 .and. f >= 1 .and. f <= 6 .and. f /= 2)) then
      call tau2par(f, tau, x0)
      if (f == 2) then
        lb(1:2) = [x0(1), 2.01_dp]; ub(1:2) = [x0(1), 30.0_dp]; cm = 2
      else
        cm = 0
      end if
    else
      if (f == 11) then
        x0 = [0.5_dp, 0.5_dp, min(1.0_dp + 6.0_dp * abs(winsor(tau)), 60.0_dp)]
      else if (f >= 7) then
        x0 = lb + 0.1_dp
      else
        call tau2par(f, winsor(tau), x0)
      end if
      call adjust_bounds(f, tau, lb, ub)
      x0 = min(max(x0, lb), ub)
      cm = 1
    end if
  end subroutine cand_spec

  ! ======================================================================
  ! Nonparametric pair copulas (24 tll, 25 bernstein) stored in aux vectors,
  ! evaluation wrappers over every family, discrete margins, weighted Kendall
  ! tau and the extra tree criteria (Chatterjee xi, KSG mutual information,
  ! Ferreira tail coefficient).
  ! ======================================================================

  ! cell k (1..m-1) with g(k) <= x < g(k+1), clamped; g ascending
  pure integer function gcell(m, g, x)
    integer, intent(in) :: m
    real(dp), intent(in) :: g(m), x
    integer :: lo, hi, mid
    lo = 1; hi = m
    do while (hi - lo > 1)
      mid = (lo + hi) / 2
      if (g(mid) <= x) then
        lo = mid
      else
        hi = mid
      end if
    end do
    gcell = min(lo, m - 1)
  end function gcell

  ! integral from g(1) to x of the piecewise-linear function through (g(k), v(k))
  ! (vinecopulib InterpolationGrid::int_on_grid)
  pure real(dp) function int_grid(m, g, v, x)
    integer, intent(in) :: m
    real(dp), intent(in) :: g(m), v(m), x
    integer :: k
    real(dp) :: s, dx
    s = 0.0_dp
    do k = 1, m - 1
      if (x >= g(k + 1)) then
        s = s + 0.5_dp * (v(k) + v(k + 1)) * (g(k + 1) - g(k))
      else
        if (x > g(k)) then
          dx = x - g(k)
          s = s + (2.0_dp * v(k) + (v(k + 1) - v(k)) / (g(k + 1) - g(k)) * dx) * dx * 0.5_dp
        end if
        exit
      end if
    end do
    int_grid = s
  end function int_grid

  ! TLL: bilinear interpolation of the density grid; h-functions integrate the
  ! interpolated density along the free variable (clamped at 1e-4 and normalised),
  ! exactly as vinecopulib's InterpolationGrid::integrate_1d.
  pure subroutine tll_eval(ax, a, b, nh, lp, h1, h2)
    real(dp), intent(in) :: ax(:), a, b
    logical, intent(in) :: nh
    real(dp), intent(out) :: lp, h1, h2
    real(dp) :: g(TLLM), col(TLLM), x2x, xx1, y2y, yy1, dens
    integer :: i, j, k
    g = ax(4:3 + TLLM)
    i = gcell(TLLM, g, a); j = gcell(TLLM, g, b)
    x2x = g(i + 1) - a; xx1 = a - g(i); y2y = g(j + 1) - b; yy1 = b - g(j)
    dens = (tv(ax, i, j) * x2x * y2y + tv(ax, i + 1, j) * xx1 * y2y + tv(ax, i, j + 1) * x2x * yy1 &
            + tv(ax, i + 1, j + 1) * xx1 * yy1) / ((g(i + 1) - g(i)) * (g(j + 1) - g(j)))
    lp = log(max(dens, TINYV))
    h1 = b; h2 = a
    if (.not. nh) return
    do k = 1, TLLM
      col(k) = max((tv(ax, i, k) * x2x + tv(ax, i + 1, k) * xx1) / (g(i + 1) - g(i)), 1.0e-4_dp)
    end do
    h1 = min(max(int_grid(TLLM, g, col, b) / max(int_grid(TLLM, g, col, 1.0_dp), 1.0e-20_dp), 1.0e-10_dp), &
             1.0_dp - 1.0e-10_dp)
    do k = 1, TLLM
      col(k) = max((tv(ax, k, j) * y2y + tv(ax, k, j + 1) * yy1) / (g(j + 1) - g(j)), 1.0e-4_dp)
    end do
    h2 = min(max(int_grid(TLLM, g, col, a) / max(int_grid(TLLM, g, col, 1.0_dp), 1.0e-20_dp), 1.0e-10_dp), &
             1.0_dp - 1.0e-10_dp)
  end subroutine tll_eval

  pure real(dp) function tv(ax, p, q)
    real(dp), intent(in) :: ax(:)
    integer, intent(in) :: p, q
    tv = ax(3 + TLLM + p + (q - 1) * TLLM)
  end function tv

  ! exact inverse of the TLL h-function: the integrand is piecewise linear, so its
  ! integral is piecewise quadratic and each cell is solved in closed form.
  ! which = 1: v with h1(cnd, v) = w;  which = 2: v with h2(v, cnd) = w.
  pure real(dp) function tll_hinv(ax, cnd, w0, which)
    real(dp), intent(in) :: ax(:), cnd, w0
    integer, intent(in) :: which
    real(dp) :: g(TLLM), col(TLLM), x2x, xx1, tot, tgt, s, ci, rem, sl, dx, w
    integer :: i, k
    g = ax(4:3 + TLLM)
    w = min(max(w0, 1.0e-10_dp), 1.0_dp - 1.0e-10_dp)
    i = gcell(TLLM, g, cnd)
    x2x = g(i + 1) - cnd; xx1 = cnd - g(i)
    do k = 1, TLLM
      if (which == 1) then
        col(k) = max((tv(ax, i, k) * x2x + tv(ax, i + 1, k) * xx1) / (g(i + 1) - g(i)), 1.0e-4_dp)
      else
        col(k) = max((tv(ax, k, i) * x2x + tv(ax, k, i + 1) * xx1) / (g(i + 1) - g(i)), 1.0e-4_dp)
      end if
    end do
    tot = int_grid(TLLM, g, col, 1.0_dp)
    tgt = w * tot
    s = 0.0_dp
    tll_hinv = 1.0_dp
    do k = 1, TLLM - 1
      ci = 0.5_dp * (col(k) + col(k + 1)) * (g(k + 1) - g(k))
      if (s + ci >= tgt .or. k == TLLM - 1) then
        rem = max(tgt - s, 0.0_dp)
        sl = (col(k + 1) - col(k)) / (g(k + 1) - g(k))
        dx = 2.0_dp * rem / (col(k) + sqrt(max(col(k)**2 + 2.0_dp * sl * rem, 0.0_dp)))
        tll_hinv = g(k) + min(max(dx, 0.0_dp), g(k + 1) - g(k))
        exit
      end if
      s = s + ci
    end do
  end function tll_hinv

  ! binomial probabilities p(0:m) of Bin(m, x)
  pure subroutine binom_pmf(m, x, p)
    integer, intent(in) :: m
    real(dp), intent(in) :: x
    real(dp), intent(out) :: p(0:m)
    integer :: k
    real(dp) :: r
    if (x <= 0.5_dp) then
      p(0) = exp(m * log1p_(-x)); r = x / (1.0_dp - x)
      do k = 0, m - 1
        p(k + 1) = p(k) * r * real(m - k, dp) / real(k + 1, dp)
      end do
    else
      p(m) = exp(m * log(x)); r = (1.0_dp - x) / x
      do k = m, 1, -1
        p(k - 1) = p(k) * r * real(k, dp) / real(m - k + 1, dp)
      end do
    end if
  end subroutine binom_pmf

  ! Bernstein copula of order m with cell weights W (rows/cols sum to 1/m):
  ! c = m^2 sum W_ij b_{i,m-1}(a) b_{j,m-1}(b),
  ! h1 = m sum_i b_{i,m-1}(a) sum_j W_ij P(Bin(m, b) >= j)  (and symmetrically h2)
  pure subroutine bern_eval(ax, a, b, nh, lp, h1, h2)
    real(dp), intent(in) :: ax(:), a, b
    logical, intent(in) :: nh
    real(dp), intent(out) :: lp, h1, h2
    real(dp) :: pa(0:15), pb(0:15), qa(0:16), qb(0:16), ta(17), tb(17), s, dens
    integer :: m, i, j
    m = nint(ax(3))
    call binom_pmf(m - 1, a, pa(0:m - 1)); call binom_pmf(m - 1, b, pb(0:m - 1))
    dens = 0.0_dp
    do j = 1, m
      s = 0.0_dp
      do i = 1, m
        s = s + ax(3 + i + (j - 1) * m) * pa(i - 1)
      end do
      dens = dens + s * pb(j - 1)
    end do
    lp = log(max(dens * m * m, TINYV))
    h1 = b; h2 = a
    if (.not. nh) return
    call binom_pmf(m, a, qa(0:m)); call binom_pmf(m, b, qb(0:m))
    ta(m + 1) = 0.0_dp; tb(m + 1) = 0.0_dp
    do j = m, 1, -1
      ta(j) = ta(j + 1) + qa(j); tb(j) = tb(j + 1) + qb(j)
    end do
    h1 = 0.0_dp; h2 = 0.0_dp
    do j = 1, m
      do i = 1, m
        h1 = h1 + ax(3 + i + (j - 1) * m) * pa(i - 1) * tb(j)
        h2 = h2 + ax(3 + i + (j - 1) * m) * pb(j - 1) * ta(i)
      end do
    end do
    h1 = min(max(h1 * m, 0.0_dp), 1.0_dp); h2 = min(max(h2 * m, 0.0_dp), 1.0_dp)
  end subroutine bern_eval

  pure subroutine np_eval(ax, a, b, nh, lp, h1, h2)
    real(dp), intent(in) :: ax(:), a, b
    logical, intent(in) :: nh
    real(dp), intent(out) :: lp, h1, h2
    if (nint(ax(1)) == 24) then
      call tll_eval(ax, a, b, nh, lp, h1, h2)
    else
      call bern_eval(ax, a, b, nh, lp, h1, h2)
    end if
  end subroutine np_eval

  ! inverse h-function of a nonparametric pair copula (which as in tll_hinv)
  pure real(dp) function np_hinv(ax, cnd, w0, which)
    real(dp), intent(in) :: ax(:), cnd, w0
    integer, intent(in) :: which
    real(dp) :: lo, hi, v, vn, lp, h1, h2, f, w
    integer :: it
    if (nint(ax(1)) == 24) then
      np_hinv = tll_hinv(ax, cnd, w0, which); return
    end if
    ! bernstein: bracketed Newton with the density as the exact derivative
    w = clampu(w0); lo = 0.0_dp; hi = 1.0_dp; v = w
    do it = 1, 100
      if (which == 1) then
        call bern_eval(ax, cnd, clampu(v), .true., lp, h1, h2); f = h1 - w
      else
        call bern_eval(ax, clampu(v), cnd, .true., lp, h1, h2); f = h2 - w
      end if
      if (f < 0.0_dp) then
        lo = v
      else
        hi = v
      end if
      if (abs(f) <= 1.0e-15_dp .or. hi - lo <= 1.0e-15_dp) exit
      vn = v - f / exp(lp)
      if (.not. (vn > lo .and. vn < hi)) vn = 0.5_dp * (lo + hi)
      if (abs(vn - v) <= 1.0e-15_dp * max(v, 1.0e-300_dp)) then
        v = vn; exit
      end if
      v = vn
    end do
    np_hinv = v
  end function np_hinv

  ! transpose a nonparametric estimate (swap the roles of u1 and u2; edge flip)
  pure subroutine np_flip(ax)
    real(dp), intent(inout) :: ax(:)
    integer :: m, o, i, j
    real(dp) :: t
    m = nint(ax(3)); o = 3
    if (nint(ax(1)) == 24) o = 3 + TLLM
    do j = 1, m
      do i = j + 1, m
        t = ax(o + i + (j - 1) * m); ax(o + i + (j - 1) * m) = ax(o + j + (i - 1) * m); ax(o + j + (i - 1) * m) = t
      end do
    end do
  end subroutine np_flip

  ! any family, unrotated data (u1, u2); ax = aux vector (used by 24, 25 only)
  pure subroutine eval_x(fam, rot, pp, ax, u1, u2, nh, lp, h1, h2)
    integer, intent(in) :: fam, rot
    real(dp), intent(in) :: pp(NPP), ax(:), u1, u2
    logical, intent(in) :: nh
    real(dp), intent(out) :: lp, h1, h2
    if (fam == 0) then
      lp = 0.0_dp; h1 = u2; h2 = u1
    else if (fam_is_np(fam)) then
      call np_eval(ax, clampu(u1), clampu(u2), nh, lp, h1, h2)
    else
      call eval_rot(fam, rot, pp, u1, u2, nh, lp, h1, h2)
    end if
  end subroutine eval_x

  ! which = 1: v with h1(u1, v) = u2;  which = 2: v with h2(v, u2) = u1
  pure real(dp) function hinv_x(fam, rot, pp, ax, u1, u2, which)
    integer, intent(in) :: fam, rot, which
    real(dp), intent(in) :: pp(NPP), ax(:), u1, u2
    if (fam == 0) then
      hinv_x = merge(u2, u1, which == 1)
    else if (fam_is_np(fam)) then
      if (which == 1) then
        hinv_x = np_hinv(ax, clampu(u1), u2, 1)
      else
        hinv_x = np_hinv(ax, clampu(u2), u1, 2)
      end if
    else
      hinv_x = hinv_rot(fam, rot, pp, u1, u2, which)
    end if
  end function hinv_x

  ! Copula CDF C(a, b) of the unrotated family in closed form; ok = .false. where none
  ! exists (gaussian, student: bivariate normal / t integrals; tll).
  pure logical function has_cdf(fam)
    integer, intent(in) :: fam
    has_cdf = .not. (fam == 2)
  end function has_cdf

  pure real(dp) function cdf0(fam, pp, ax, a, b)
    integer, intent(in) :: fam
    real(dp), intent(in) :: pp(NPP), ax(:), a, b
    logical :: ok
    real(dp) :: th, de, t1, t2, ls, x, y, ua, ub, s, eta, za, zb, l1, m, la, lb, lab, t, aa, a1, a2, rp, d1, d2
    real(dp) :: qa(0:16), qb(0:16), ta(17), tb(17), rw(TLLM)
    integer :: i, j, mm
    ok = .true.; th = pp(1); de = pp(2)
    select case (fam)
    case (0)
      cdf0 = a * b
    case (3)
      if (th < 1.0e-10_dp) then
        cdf0 = a * b; return
      end if
      t1 = -th * log(a); t2 = -th * log(b)
      ls = max(t1, t2) + log(exp(t1 - max(t1, t2)) + exp(t2 - max(t1, t2)) - exp(-max(t1, t2)))
      cdf0 = exp(-ls / th)
    case (4)
      x = -log(a); y = -log(b); t1 = max(th, 1.0_dp) * log(x); t2 = max(th, 1.0_dp) * log(y)
      cdf0 = exp(-exp((max(t1, t2) + log1p_(exp(min(t1, t2) - max(t1, t2)))) / max(th, 1.0_dp)))
    case (5)
      if (abs(th) < 1.0e-10_dp) then
        cdf0 = a * b; return
      end if
      cdf0 = -log1p_(expm1_(-th * a) * expm1_(-th * b) / expm1_(-th)) / th
    case (6)
      ua = exp(max(th, 1.0_dp) * log1p_(-a)); ub = exp(max(th, 1.0_dp) * log1p_(-b))
      cdf0 = 1.0_dp - exp(log(ua + ub - ua * ub) / max(th, 1.0_dp))
    case (7:10)
      call arch(fam, th, de, a, t1, d1, d2); call arch(fam, th, de, b, t2, d1, d2)
      cdf0 = arch_psi(fam, th, de, t1 + t2)
    case (11, 15, 16, 19, 23)
      ! extreme value: C = exp(log(ab) A(t)), t = log b / log(ab)
      la = log(a); lb = log(b); lab = la + lb; t = lb / lab
      select case (fam)
      case (11)
        call tawn_pickands(pp, t, aa, a1, a2)
      case (15)
        aa = 1.0_dp - exp(-log(exp(-th * log(t)) + exp(-th * log1p_(-t))) / th)
      case (16)
        x = 1.0_dp / th + 0.5_dp * th * log(t / (1.0_dp - t)); y = 1.0_dp / th - 0.5_dp * th * log(t / (1.0_dp - t))
        aa = t * pnorm(x) + (1.0_dp - t) * pnorm(y)
      case (19)
        aa = (t**th + (1.0_dp - t)**th - exp(-log(exp(-th * de * log(t)) + exp(-th * de * log1p_(-t))) / de))**(1.0_dp / th)
      case default
        call tev_A(pp, t, aa)
      end select
      cdf0 = exp(lab * aa)
    case (12)
      eta = 1.0_dp + (th - 1.0_dp) * (a + b)
      cdf0 = 2.0_dp * th * a * b / (eta + sqrt(max(eta * eta - 4.0_dp * th * (th - 1.0_dp) * a * b, 0.0_dp)))
    case (13)
      cdf0 = a * b * (1.0_dp + th * (1.0_dp - a) * (1.0_dp - b))
    case (14)
      cdf0 = a * b / (1.0_dp - th * (1.0_dp - a) * (1.0_dp - b))
    case (17, 18)
      if (fam == 17) then
        za = de * expm1_(-th * log(a)); zb = de * expm1_(-th * log(b))
      else
        za = de * exp(th * log(-log(a))); zb = de * exp(th * log(-log(b)))
      end if
      if (za < 1.0_dp .and. zb < 1.0_dp) then             ! l1 = log(1 + expm1(za) + expm1(zb))
        l1 = log1p_(expm1_(za) + expm1_(zb))
      else
        l1 = max(za, zb) + log(1.0_dp + exp(min(za, zb) - max(za, zb)) - exp(-max(za, zb)))
      end if
      if (fam == 17) then
        cdf0 = exp(-log1p_(l1 / de) / th)
      else
        cdf0 = exp(-exp(log(l1 / de) / th))
      end if
    case (20)
      m = (1.0_dp / de - log(a))**th + (1.0_dp / de - log(b))**th - (1.0_dp / de)**th
      cdf0 = exp(1.0_dp / de - m**(1.0_dp / th))
    case (21)
      cdf0 = a * b * exp(-log(1.0_dp - de * (1.0_dp - a**th) * (1.0_dp - b**th)) / th)
    case (22)
      x = expm1_(-th * log(a)); y = expm1_(-th * log(b))
      s = exp(-log(exp(-de * log(x)) + exp(-de * log(y))) / de)
      cdf0 = exp(-log(1.0_dp + x + y - s) / th)
    case (1)
      cdf0 = bvn_lower(qnorm(a), qnorm(b), th)
    case (24)
      ! tll: exact integral of the bilinear density interpolant over [0,a] x [0,b]: for each
      ! grid row the integral over [0,b] is piecewise-linear in the first coordinate, so
      ! integrating those row integrals as a piecewise-linear function is exact
      ! (vinecopulib InterpolationGrid::integrate_2d, clipped to [1e-10, 1 - 1e-10])
      do i = 1, TLLM
        rw(i) = int_grid(TLLM, ax(4:3 + TLLM), ax(3 + TLLM + i:3 + TLLM + TLLM * TLLM:TLLM), b)
      end do
      cdf0 = min(max(int_grid(TLLM, ax(4:3 + TLLM), rw, a), 1.0e-10_dp), 1.0_dp - 1.0e-10_dp)
      return
    case (25)
      ! bernstein: C = sum W_ij P(Bin(m, a) >= i) P(Bin(m, b) >= j)
      mm = nint(ax(3))
      call binom_pmf(mm, a, qa(0:mm)); call binom_pmf(mm, b, qb(0:mm))
      ta(mm + 1) = 0.0_dp; tb(mm + 1) = 0.0_dp
      do j = mm, 1, -1
        ta(j) = ta(j + 1) + qa(j); tb(j) = tb(j + 1) + qb(j)
      end do
      cdf0 = 0.0_dp
      do j = 1, mm
        do i = 1, mm
          cdf0 = cdf0 + ax(3 + i + (j - 1) * mm) * ta(i) * tb(j)
        end do
      end do
    case default
      cdf0 = 0.0_dp; ok = .false.
      if (gen_is(fam)) then                       ! generic: C = Psi(l(Phi(a), Phi(b)))
        call gen_eval_rv(fam, pp(1:3), a, b, d1, d2, rp, t)
        cdf0 = exp(t); ok = .true.
      end if
    end select
    if (ok) cdf0 = min(max(cdf0, 0.0_dp), min(a, b))
    rp = 0.0_dp
  end function cdf0

  ! bivariate normal P(X < h, Y < k), correlation r: Genz (2004) BVND (Drezner-Wesolowsky
  ! with Gauss-Legendre 6/12/20 nodes by |r|, asymptotic expansion for |r| > 0.925),
  ! double-precision accurate (~1e-15)
  pure real(dp) function bvn_lower(h0, k0, r) result(p)
    real(dp), intent(in) :: h0, k0, r
    real(dp), parameter :: W6(3) = [0.1713244923791705_dp, 0.3607615730481384_dp, 0.4679139345726904_dp], &
      X6(3) = [0.9324695142031522_dp, 0.6612093864662647_dp, 0.2386191860831970_dp], &
      W12(6) = [0.04717533638651177_dp, 0.1069393259953183_dp, 0.1600783285433464_dp, 0.2031674267230659_dp, &
                0.2334925365383547_dp, 0.2491470458134029_dp], &
      X12(6) = [0.9815606342467191_dp, 0.9041172563704750_dp, 0.7699026741943050_dp, 0.5873179542866171_dp, &
                0.3678314989981802_dp, 0.1252334085114692_dp], &
      W20(10) = [0.01761400713915212_dp, 0.04060142980038694_dp, 0.06267204833410906_dp, &
                 0.08327674157670475_dp, 0.1019301198172404_dp, 0.1181945319615184_dp, 0.1316886384491766_dp, &
                 0.1420961093183821_dp, 0.1491729864726037_dp, 0.1527533871307259_dp], &
      X20(10) = [0.9931285991850949_dp, 0.9639719272779138_dp, 0.9122344282513259_dp, 0.8391169718222188_dp, &
                 0.7463319064601508_dp, 0.6360536807265150_dp, 0.5108670019508271_dp, 0.3737060887154196_dp, &
                 0.2277858511416451_dp, 0.07652652113349733_dp]
    real(dp) :: w(20), x(20), h, k, hk, bvn, hs, asr, sn, as, a, bs, c, d, b, sp, xs, rs, ep, tp, l
    integer :: ng, i
    tp = 2.0_dp * PI
    h = -h0; k = -k0            ! upper-orthant form: P(X > h, Y > k)
    if (abs(r) < 0.3_dp) then
      ng = 3; w(1:3) = W6; x(1:3) = X6
    else if (abs(r) < 0.75_dp) then
      ng = 6; w(1:6) = W12; x(1:6) = X12
    else
      ng = 10; w(1:10) = W20; x(1:10) = X20
    end if
    w(ng + 1:2 * ng) = w(1:ng); x(ng + 1:2 * ng) = 1.0_dp + x(1:ng); x(1:ng) = 1.0_dp - x(1:ng)
    hk = h * k; bvn = 0.0_dp
    if (abs(r) < 0.925_dp) then
      hs = 0.5_dp * (h * h + k * k); asr = 0.5_dp * asin(r)
      do i = 1, 2 * ng
        sn = sin(asr * x(i))
        bvn = bvn + w(i) * exp((sn * hk - hs) / (1.0_dp - sn * sn))
      end do
      bvn = bvn * asr / tp + pnorm(-h) * pnorm(-k)
    else
      if (r < 0.0_dp) then
        k = -k; hk = -hk
      end if
      if (abs(r) < 1.0_dp) then
        as = (1.0_dp - r) * (1.0_dp + r); a = sqrt(as); bs = (h - k)**2
        asr = -0.5_dp * (bs / as + hk); c = (4.0_dp - hk) / 8.0_dp; d = (12.0_dp - hk) / 80.0_dp
        if (asr > -100.0_dp) bvn = a * exp(asr) * (1.0_dp - c * (bs - as) * (1.0_dp - d * bs) / 3.0_dp + c * d * as * as)
        if (hk > -100.0_dp) then
          b = sqrt(bs); sp = sqrt(tp) * pnorm(-b / a)
          bvn = bvn - exp(-0.5_dp * hk) * sp * b * (1.0_dp - c * bs * (1.0_dp - d * bs) / 3.0_dp)
        end if
        a = 0.5_dp * a
        sp = 0.0_dp
        do i = 1, 2 * ng
          xs = (a * x(i))**2; asr = -0.5_dp * (bs / xs + hk)
          if (asr > -100.0_dp) then
            rs = sqrt(1.0_dp - xs)
            ep = exp(-0.5_dp * hk * xs / (1.0_dp + rs)**2) / rs
            sp = sp + w(i) * exp(asr) * ((1.0_dp + c * xs * (1.0_dp + 5.0_dp * d * xs)) - ep)
          end if
        end do
        bvn = (a * sp - bvn) / tp
      end if
      if (r > 0.0_dp) then
        bvn = bvn + pnorm(-max(h, k))
      else if (h >= k) then
        bvn = -bvn
      else
        if (h < 0.0_dp) then
          l = pnorm(k) - pnorm(h)
        else
          l = pnorm(-h) - pnorm(-k)
        end if
        bvn = l - bvn
      end if
    end if
    p = max(0.0_dp, min(1.0_dp, bvn))
  end function bvn_lower

  ! t-EV Pickands function A(t)
  pure subroutine tev_A(pp, t, aa)
    real(dp), intent(in) :: pp(NPP), t
    real(dp), intent(out) :: aa
    real(dp) :: r, nu, k, x1, x2
    r = min(max(pp(1), -0.999_dp), 0.999_dp); nu = min(max(pp(2), 1.0_dp), 60.0_dp)
    k = sqrt((nu + 1.0_dp) / ((1.0_dp - r) * (1.0_dp + r)))
    x1 = exp((log(t) - log1p_(-t)) / nu); x2 = exp((log1p_(-t) - log(t)) / nu)
    aa = t * pt_l(k * (x1 - r), nu + 1.0_dp, pp(8)) + (1.0_dp - t) * pt_l(k * (x2 - r), nu + 1.0_dp, pp(8))
  end subroutine tev_A

  ! CDF of the rotated copula at unrotated (u1, u2): 90: u2 - C(u2, 1-u1),
  ! 180: u1 + u2 - 1 + C(1-u1, 1-u2), 270: u1 - C(1-u2, u1)
  pure real(dp) function cdf_x(fam, rot, pp, ax, u1, u2)
    integer, intent(in) :: fam, rot
    real(dp), intent(in) :: pp(NPP), ax(:), u1, u2
    real(dp) :: a, b
    a = clampu(u1); b = clampu(u2)
    select case (rot)
    case (90)
      cdf_x = b - cdf0(fam, pp, ax, b, 1.0_dp - a)
    case (180)
      cdf_x = a + b - 1.0_dp + cdf0(fam, pp, ax, 1.0_dp - a, 1.0_dp - b)
    case (270)
      cdf_x = a - cdf0(fam, pp, ax, 1.0_dp - b, a)
    case default
      cdf_x = cdf0(fam, pp, ax, a, b)
    end select
  end function cdf_x

  ! copula CDF at n points (closed forms, see cdf0); student, which has none:
  ! C(u1, u2) = int_0^u1 h1(s, u2) ds by Gauss-Legendre (16 nodes on each of 8 panels)
  subroutine bicop_cdf_x(fam, rot, par, ax, n, u1, u2, out)
    integer, intent(in) :: fam, rot, n
    real(dp), intent(in) :: par(3), ax(:), u1(n), u2(n)
    real(dp), intent(out) :: out(n)
    real(dp), parameter :: gx(8) = [0.0950125098376374_dp, 0.2816035507792589_dp, 0.4580167776572274_dp, &
      0.6178762444026438_dp, 0.7554044083550030_dp, 0.8656312023878318_dp, 0.9445750230732326_dp, &
      0.9894009349916499_dp]
    real(dp), parameter :: gw(8) = [0.1894506104550685_dp, 0.1826034150449236_dp, 0.1691565193950025_dp, &
      0.1495959888165767_dp, 0.1246289712555339_dp, 0.0951585116824928_dp, 0.0622535239386479_dp, &
      0.0271524594117541_dp]
    real(dp) :: pp(NPP), s, c0, hw, x, l1, g1, g2, pl, ph
    integer :: i, ip, q, sg, np
    pp = prep_par(fam, par)
    !$omp parallel do private(i, s, c0, hw, x, l1, g1, g2, pl, ph, ip, q, sg, np) schedule(static)
    do i = 1, n
      if (has_cdf(fam)) then
        out(i) = cdf_x(fam, rot, pp, ax, u1(i), u2(i))
      else
        np = 8; s = 0.0_dp
        do ip = 1, np
          pl = clampu(u1(i)) * (ip - 1) / np; ph = clampu(u1(i)) * ip / np
          c0 = 0.5_dp * (pl + ph); hw = 0.5_dp * (ph - pl)
          do q = 1, 8
            do sg = -1, 1, 2
              x = c0 + sg * hw * gx(q)
              call eval_x(fam, rot, pp, ax, x, u2(i), .true., l1, g1, g2)
              s = s + gw(q) * hw * g1
            end do
          end do
        end do
        out(i) = min(max(s, 0.0_dp), 1.0_dp)
      end if
    end do
    !$omp end parallel do
  end subroutine bicop_cdf_x

  ! Pair copula with possibly discrete margins (vinecopulib var_types "d"): a, b are
  ! F(x) and am, bm = F(x-) of a discrete margin.  lp is the log density with respect
  ! to the product of the margins (C differences divided by the margin masses);
  ! h1 = F(b|a), h2 = F(a|b) and the minus versions h1m = F(b-|a), h2m = F(a-|b)
  ! feed the next tree.  For a discrete conditioning variable F(b|a) =
  ! (C(a,b) - C(am,b)) / (a - am) = mean of h1(s, b) over s in [am, a] (Gauss-Legendre).
  pure subroutine eval_disc(fam, rot, pp, ax, a, b, am, bm, da, db, lp, h1, h2, h1m, h2m, need_h)
    integer, intent(in) :: fam, rot
    real(dp), intent(in) :: pp(NPP), ax(:), a, b, am, bm
    logical, intent(in) :: da, db
    real(dp), intent(out) :: lp, h1, h2, h1m, h2m
    ! need_h = .false.: only lp is wanted (likelihood fits) -- the h-functions that only feed
    ! the next tree are skipped (lp itself is computed exactly as with need_h)
    logical, intent(in), optional :: need_h
    real(dp) :: t1, t2
    logical :: xa, xb, wh
    wh = .true.
    if (present(need_h)) wh = need_h
    xa = da .and. a - am > 1.0e-12_dp
    xb = db .and. b - bm > 1.0e-12_dp
    if (.not. (xa .or. xb)) then
      call eval_x(fam, rot, pp, ax, a, b, wh, lp, h1, h2)
      h1m = h1; h2m = h2
    else if (xa .and. .not. xb) then
      call eval_x(fam, rot, pp, ax, a, b, .true., lp, t1, h2)
      call eval_x(fam, rot, pp, ax, am, b, .true., lp, t1, h2m)
      lp = log(max((h2 - h2m) / (a - am), TINYV))
      if (.not. wh) then
        h1 = 0.0_dp; h1m = 0.0_dp; return
      end if
      h1 = avg_h(1, am, a, b); h1m = h1
    else if (xb .and. .not. xa) then
      call eval_x(fam, rot, pp, ax, a, b, .true., lp, h1, t2)
      call eval_x(fam, rot, pp, ax, a, bm, .true., lp, h1m, t2)
      lp = log(max((h1 - h1m) / (b - bm), TINYV))
      if (.not. wh) then
        h2 = 0.0_dp; h2m = 0.0_dp; return
      end if
      h2 = avg_h(2, bm, b, a); h2m = h2
    else
      h1 = avg_h(1, am, a, b); h1m = avg_h(1, am, a, bm)
      lp = log(max((h1 - h1m) / (b - bm), TINYV))
      if (.not. wh) then
        h2 = 0.0_dp; h2m = 0.0_dp; return
      end if
      h2 = avg_h(2, bm, b, a); h2m = avg_h(2, bm, b, am)
    end if
  contains
    ! mean over s in [lo, hi] of h1(s, y) (k = 1) or h2(y, s) (k = 2)
    pure real(dp) function avg_h(k, lo, hi, y)
      integer, intent(in) :: k
      real(dp), intent(in) :: lo, hi, y
      real(dp), parameter :: gx(8) = [0.0950125098376374_dp, 0.2816035507792589_dp, 0.4580167776572274_dp, &
        0.6178762444026438_dp, 0.7554044083550030_dp, 0.8656312023878318_dp, 0.9445750230732326_dp, &
        0.9894009349916499_dp]
      real(dp), parameter :: gw(8) = [0.1894506104550685_dp, 0.1826034150449236_dp, 0.1691565193950025_dp, &
        0.1495959888165767_dp, 0.1246289712555339_dp, 0.0951585116824928_dp, 0.0622535239386479_dp, &
        0.0271524594117541_dp]
      integer :: np, ip, q, sg
      real(dp) :: pl, ph, c0, hw, s, x, l1, g1, g2, ch, cl
      ! closed form (C(hi, y) - C(lo, y)) / (hi - lo) where the copula CDF has one; tiny
      ! cells (< 1e-3) keep the quadrature, which there is more accurate than a difference
      if (hi - lo >= 1.0e-3_dp .and. has_cdf(fam)) then
        if (k == 1) then
          ch = cdf_x(fam, rot, pp, ax, hi, y); cl = cdf_x(fam, rot, pp, ax, lo, y)
        else
          ch = cdf_x(fam, rot, pp, ax, y, hi); cl = cdf_x(fam, rot, pp, ax, y, lo)
        end if
        avg_h = min(max((ch - cl) / (hi - lo), 0.0_dp), 1.0_dp)
        return
      end if
      np = 1
      if (hi - lo > 0.02_dp) np = 4
      if (hi - lo > 0.2_dp) np = 8
      s = 0.0_dp
      do ip = 1, np
        pl = lo + (hi - lo) * (ip - 1) / np; ph = lo + (hi - lo) * ip / np
        c0 = 0.5_dp * (pl + ph); hw = 0.5_dp * (ph - pl)
        do q = 1, 8
          do sg = -1, 1, 2
            x = c0 + sg * hw * gx(q)
            if (k == 1) then
              call eval_x(fam, rot, pp, ax, x, y, .true., l1, g1, g2); s = s + gw(q) * hw * g1
            else
              call eval_x(fam, rot, pp, ax, y, x, .true., l1, g1, g2); s = s + gw(q) * hw * g2
            end if
          end do
        end do
      end do
      avg_h = min(max(s / (hi - lo), 0.0_dp), 1.0_dp)
    end function avg_h
  end subroutine eval_disc

  ! ----------------------------------------------------------- TLL estimation
  ! vinecopulib TllBicop::fit (nonparametric_method "constant", mult 1, 30 x 30
  ! normal-score grid), ported via torchvine's tll_fit: Gaussian-kernel density of
  ! the normal scores with bandwidth B = n^(-1/3) * cor-matrix * (|cor|/|mcor|)^(mcor/2)
  ! (mcor = ACE maximal correlation), c = f / (phi phi), margins normalised 3 times.
  ! npars = effective degrees of freedom: influence K(0) det(B^-1/2) / (n f) summed over
  ! the observations (interpolated from the grid).
  subroutine tll_fit(n, a, b, usew, w, ax, ll)
    integer, intent(in) :: n
    real(dp), intent(in) :: a(n), b(n), w(n)
    logical, intent(in) :: usew
    real(dp), intent(out) :: ax(:), ll
    integer, parameter :: m = TLLM
    real(dp), allocatable :: z1(:), z2(:), r1(:), r2(:), zd1(:), zd2(:), ww(:), p0(:), p1(:), e1(:)
    real(dp) :: g(m), zg(m), cor, mcor, base, scl, mult0, l11, l21, l22, det, ze1, ze2, s, dx2, dy2
    real(dp) :: f0(m, m), infl(m, m), v(m, m), rint, k0, ii, lp, h1, h2
    integer :: i, j, q, it
    allocate(z1(n), z2(n), r1(n), r2(n), zd1(n), zd2(n), ww(n), p0(n), p1(n), e1(n))
    ww = 1.0_dp
    if (usew) ww = w
    call rank_vec(n, a, r1); call rank_vec(n, b, r2)
    do i = 1, n
      z1(i) = qnorm(r1(i) / (n + 1.0_dp)); z2(i) = qnorm(r2(i) / (n + 1.0_dp))
    end do
    cor = min(max(wpearson(n, z1, z2, ww), -0.95_dp), 0.95_dp)
    call ace2(n, z1, z2, usew, ww, p0, p1)
    mcor = min(max(wpearson(n, p0, p1, ww), -0.99_dp), 0.99_dp)
    base = max(abs(cor) / max(abs(mcor), 1.0e-6_dp), 1.0e-12_dp)
    scl = base**(0.5_dp * mcor)
    mult0 = real(n, dp)**(-1.0_dp / 3.0_dp)
    l11 = sqrt(mult0 * scl); l21 = cor * mult0 * scl / l11; l22 = sqrt(mult0 * scl - l21 * l21)
    det = 1.0_dp / (l11 * l22)
    do i = 1, n
      zd1(i) = z1(i) / l11; zd2(i) = -l21 / (l11 * l22) * z1(i) + z2(i) / l22
    end do
    do i = 1, m
      zg(i) = -3.25_dp + 6.5_dp * (i - 1) / (m - 1)
      g(i) = pnorm(zg(i))
    end do
    k0 = 1.0_dp / (2.0_dp * PI)
    do i = 1, m
      ze1 = zg(i) / l11
      do q = 1, n
        dx2 = (ze1 - zd1(q))**2
        e1(q) = merge(0.0_dp, exp(-0.5_dp * dx2) * ww(q), dx2 > 80.0_dp)
      end do
      do j = 1, m
        ze2 = -l21 / (l11 * l22) * zg(i) + zg(j) / l22
        s = 0.0_dp
        do q = 1, n
          if (e1(q) == 0.0_dp) cycle
          dy2 = (ze2 - zd2(q))**2
          if (dy2 < 80.0_dp) s = s + e1(q) * exp(-0.5_dp * dy2)
        end do
        f0(i, j) = max(k0 * det * s / n, TINYV)
        infl(i, j) = k0 * det / f0(i, j) / n
        v(i, j) = f0(i, j) / (exp(-0.5_dp * (zg(i)**2 + zg(j)**2)) * k0)
      end do
    end do
    g(1) = 0.0_dp; g(m) = 1.0_dp
    do it = 1, 3
      do i = 1, m
        rint = max(int_grid(m, g, v(i, :), 1.0_dp), 1.0e-20_dp)
        v(i, :) = v(i, :) / rint
      end do
      do j = 1, m
        rint = max(int_grid(m, g, v(:, j), 1.0_dp), 1.0e-20_dp)
        v(:, j) = v(:, j) / rint
      end do
    end do
    ax = 0.0_dp
    ax(1) = 24.0_dp; ax(3) = m
    ax(4:3 + m) = g
    ax(4 + m:3 + m + m * m) = reshape(v, [m * m])
    ! effective degrees of freedom: influence grid interpolated at the observations
    s = 0.0_dp; ll = 0.0_dp
    do q = 1, n
      i = gcell(m, g, a(q)); j = gcell(m, g, b(q))
      ii = (infl(i, j) * (g(i + 1) - a(q)) * (g(j + 1) - b(q)) + infl(i + 1, j) * (a(q) - g(i)) * (g(j + 1) - b(q)) &
            + infl(i, j + 1) * (g(i + 1) - a(q)) * (b(q) - g(j)) + infl(i + 1, j + 1) * (a(q) - g(i)) * (b(q) - g(j))) &
           / ((g(i + 1) - g(i)) * (g(j + 1) - g(j)))
      s = s + ww(q) * ii
      call tll_eval(ax, a(q), b(q), .false., lp, h1, h2)
      ll = ll + ww(q) * lp
    end do
    ax(2) = max(1.0_dp, s)
  end subroutine tll_fit

  ! ACE maximal-correlation transforms (vinecopulib tools_stats::ace, torchvine port)
  subroutine ace2(n, x, y, usew, w, p0, p1)
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n), y(n), w(n)
    logical, intent(in) :: usew
    real(dp), intent(out) :: p0(n), p1(n)
    integer, allocatable :: o0(:), o1(:), i0(:), i1(:), bf(:)
    real(dp), allocatable :: t(:)
    integer :: i, wl, io, ii
    real(dp) :: sd
    allocate(o0(n), o1(n), i0(n), i1(n), bf(n), t(n))
    wl = ceiling(real(n, dp) / 5.0_dp)
    do i = 1, n
      o0(i) = i; o1(i) = i
    end do
    call msort_by(x, o0, bf, 1, n); call msort_by(y, o1, bf, 1, n)
    do i = 1, n
      i0(o0(i)) = i; i1(o1(i)) = i
    end do
    sd = sqrt(real(n, dp) * (n - 1.0_dp) / 12.0_dp)
    do i = 1, n
      p0(i) = ((i0(i) - 1) - ((n - 1.0_dp) / 2.0_dp - 1.0_dp)) / sd
      p1(i) = ((i1(i) - 1) - ((n - 1.0_dp) / 2.0_dp - 1.0_dp)) / sd
    end do
    if (usew) then
      p0 = p0 * w; p1 = p1 * w
    end if
    do io = 1, 15
      do ii = 1, 3
        call cef(p0 * w, o1, i1, t)
        t = t - sum(t) / n
        p1 = t / max(sqrt(sum(t * t) / (n - 1.0_dp)), TINYV)
      end do
      call cef(p1 * w, o0, i0, t)
      t = t - sum(t) / n
      p0 = t / max(sqrt(sum(t * t) / (n - 1.0_dp)), TINYV)
    end do
  contains
    ! conditional expectation by a moving-average window over the sorted order
    subroutine cef(v, ord, inv, out)
      real(dp), intent(in) :: v(n)
      integer, intent(in) :: ord(n), inv(n)
      real(dp), intent(out) :: out(n)
      real(dp) :: cs(0:n + 2 * wl), ys(n)
      integer :: k, kk
      cs(0) = 0.0_dp
      do k = 1, n + 2 * wl
        kk = min(max(k - wl, 1), n)
        cs(k) = cs(k - 1) + v(ord(kk))
      end do
      do k = 1, n
        ys(k) = (cs(k + 2 * wl) - cs(k - 1)) / (2 * wl + 1)
      end do
      if (wl < n) then
        ys(1:wl) = ys(wl + 1)
        ys(n - wl + 1:n) = ys(n - wl)
      end if
      do k = 1, n
        out(k) = ys(inv(k))
      end do
    end subroutine cef
  end subroutine ace2

  ! empirical Bernstein copula (Sancetta & Satchell 2004): cell masses of the data on an
  ! m x m grid, rebalanced (Sinkhorn) so every row/column holds exactly 1/m (a proper
  ! copula); order m = 1..15 chosen by AIC with (m-1)^2 parameters.
  subroutine bern_fit(n, a, b, usew, w, ax, ll)
    integer, intent(in) :: n
    real(dp), intent(in) :: a(n), b(n), w(n)
    logical, intent(in) :: usew
    real(dp), intent(out) :: ax(:), ll
    real(dp) :: cw(15, 15), tmp(NAUX), aic, best, l, wq, lp, h1, h2
    integer :: m, i, j, q, it
    best = huge(1.0_dp); ax = 0.0_dp; ll = 0.0_dp
    ax(1) = 25.0_dp; ax(2) = 0.0_dp; ax(3) = 1.0_dp; ax(4) = 1.0_dp
    best = 0.0_dp                                 ! m = 1: independence, AIC 0
    do m = 2, 15
      if ((m - 1)**2 >= n / 4) exit
      cw(1:m, 1:m) = 1.0e-9_dp / (m * m)
      do q = 1, n
        wq = 1.0_dp
        if (usew) wq = w(q)
        i = min(int(a(q) * m), m - 1) + 1; j = min(int(b(q) * m), m - 1) + 1
        cw(i, j) = cw(i, j) + wq
      end do
      cw(1:m, 1:m) = cw(1:m, 1:m) / sum(cw(1:m, 1:m))
      do it = 1, 50
        do i = 1, m
          cw(i, 1:m) = cw(i, 1:m) / (m * sum(cw(i, 1:m)))
        end do
        do j = 1, m
          cw(1:m, j) = cw(1:m, j) / (m * sum(cw(1:m, j)))
        end do
      end do
      tmp = 0.0_dp
      tmp(1) = 25.0_dp; tmp(2) = (m - 1)**2; tmp(3) = m
      tmp(4:3 + m * m) = reshape(cw(1:m, 1:m), [m * m])
      l = 0.0_dp
      do q = 1, n
        call bern_eval(tmp, a(q), b(q), .false., lp, h1, h2)
        if (usew) then
          l = l + w(q) * lp
        else
          l = l + lp
        end if
      end do
      aic = -2.0_dp * l + 2.0_dp * (m - 1)**2
      if (aic < best) then
        best = aic; ax = tmp; ll = l
      end if
    end do
  end subroutine bern_fit

  subroutine np_fit(fam, n, a, b, usew, w, ax, ll)
    integer, intent(in) :: fam, n
    real(dp), intent(in) :: a(n), b(n), w(n)
    logical, intent(in) :: usew
    real(dp), intent(out) :: ax(:), ll
    if (fam == 24) then
      call tll_fit(n, a, b, usew, w, ax, ll)
    else
      call bern_fit(n, a, b, usew, w, ax, ll)
    end if
  end subroutine np_fit

  ! ------------------------------------------------- weighted dependence measures
  pure real(dp) function wpearson(n, x, y, w)
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n), y(n), w(n)
    real(dp) :: sw, mx, my, cxy, vx, vy
    sw = max(sum(w), TINYV)
    mx = sum(w * x) / sw; my = sum(w * y) / sw
    cxy = sum(w * (x - mx) * (y - my)); vx = sum(w * (x - mx)**2); vy = sum(w * (y - my)**2)
    wpearson = cxy / max(sqrt(max(vx * vy, 0.0_dp)), TINYV)
    wpearson = min(max(wpearson, -1.0_dp), 1.0_dp)
  end function wpearson

  ! weighted Kendall tau-b, O(n log n): concordant/discordant pair weights w_i w_j
  ! counted with a Fenwick tree over y-ranks; reduces to tau_b for w = 1
  real(dp) function tau_w(n, x, y, w)
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n), y(n), w(n)
    integer, allocatable :: ix(:), bf(:), ry(:)
    real(dp), allocatable :: bit(:)
    integer :: i, j, k, nr
    real(dp) :: cc, dd, t0, tx, ty, gw, gw2, sall, sw
    allocate(ix(n), bf(n), ry(n))
    do i = 1, n
      ix(i) = i
    end do
    call msort_by(y, ix, bf, 1, n)
    nr = 0; ty = 0.0_dp; i = 1
    do while (i <= n)
      j = i; gw = 0.0_dp; gw2 = 0.0_dp
      nr = nr + 1
      do while (j <= n)
        if (y(ix(j)) /= y(ix(i))) exit
        ry(ix(j)) = nr; gw = gw + w(ix(j)); gw2 = gw2 + w(ix(j))**2
        j = j + 1
      end do
      ty = ty + 0.5_dp * (gw * gw - gw2)
      i = j
    end do
    do i = 1, n
      ix(i) = i
    end do
    call msort_by(x, ix, bf, 1, n)
    allocate(bit(nr)); bit = 0.0_dp
    cc = 0.0_dp; dd = 0.0_dp; tx = 0.0_dp; sall = 0.0_dp; i = 1
    do while (i <= n)
      j = i
      do while (j <= n)
        if (x(ix(j)) /= x(ix(i))) exit
        j = j + 1
      end do
      gw = 0.0_dp; gw2 = 0.0_dp
      do k = i, j - 1
        cc = cc + w(ix(k)) * bsum(ry(ix(k)) - 1)
        dd = dd + w(ix(k)) * (sall - bsum(ry(ix(k))))
        gw = gw + w(ix(k)); gw2 = gw2 + w(ix(k))**2
      end do
      do k = i, j - 1
        call badd(ry(ix(k)), w(ix(k)))
      end do
      sall = sall + gw
      tx = tx + 0.5_dp * (gw * gw - gw2)
      i = j
    end do
    sw = sum(w); t0 = 0.5_dp * (sw * sw - sum(w * w))
    tau_w = (cc - dd) / sqrt(max((t0 - tx) * (t0 - ty), TINYV))
  contains
    real(dp) function bsum(r0)
      integer, intent(in) :: r0
      integer :: r
      bsum = 0.0_dp; r = r0
      do while (r > 0)
        bsum = bsum + bit(r); r = r - iand(r, -r)
      end do
    end function bsum
    subroutine badd(r0, v)
      integer, intent(in) :: r0
      real(dp), intent(in) :: v
      integer :: r
      r = r0
      do while (r <= nr)
        bit(r) = bit(r) + v; r = r + iand(r, -r)
      end do
    end subroutine badd
  end function tau_w

  ! symmetric revised Chatterjee xi (torchvine bidep.chatterjee_xi, M = 1)
  real(dp) function xi_sym(n, x, y)
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n), y(n)
    real(dp), allocatable :: rx(:), ry(:)
    integer, allocatable :: px(:), py(:)
    integer :: i
    real(dp) :: s1, s2
    allocate(rx(n), ry(n), px(n), py(n))
    call rank_vec(n, x, rx); call rank_vec(n, y, ry)
    do i = 1, n
      px(nint(rx(i))) = i; py(nint(ry(i))) = i
    end do
    ! x-ranks in increasing-y order, and y-ranks in increasing-x order
    s1 = rx(py(n)); s2 = ry(px(n))
    do i = 1, n - 1
      s1 = s1 + min(rx(py(i)), rx(py(i + 1)))
      s2 = s2 + min(ry(px(i)), ry(px(i + 1)))
    end do
    xi_sym = -2.0_dp + 24.0_dp * max(s1, s2) / ((1.0_dp + n) * (2.0_dp + 4.0_dp * n))
  end function xi_sym

  ! KSG (Kraskov alg. 1) mutual information, k = 3, max-norm on standardised data
  ! (sklearn mutual_info_regression as used by torchvine, without its 1e-10 jitter)
  real(dp) function mi_ksg(n, x0, y0)
    integer, intent(in) :: n
    real(dp), intent(in) :: x0(n), y0(n)
    integer, parameter :: kk = 3
    real(dp), allocatable :: x(:), y(:), xs(:), ys(:), ysort(:)
    integer, allocatable :: ox(:), oy(:), bf(:)
    real(dp) :: dk(kk), dist, dxx, eps, sx, sy, s
    integer :: q, l, nx, ny, i
    allocate(x(n), y(n), xs(n), ys(n), ysort(n), ox(n), oy(n), bf(n))
    sx = sqrt(max(sum((x0 - sum(x0) / n)**2) / n, TINYV)); sy = sqrt(max(sum((y0 - sum(y0) / n)**2) / n, TINYV))
    x = x0 / sx; y = y0 / sy
    do i = 1, n
      ox(i) = i; oy(i) = i
    end do
    call msort_by(x, ox, bf, 1, n); call msort_by(y, oy, bf, 1, n)
    xs = x(ox); ys = y(ox); ysort = y(oy)
    s = 0.0_dp
    do q = 1, n
      dk = huge(1.0_dp)
      l = q - 1
      do while (l >= 1)
        dxx = xs(q) - xs(l)
        if (dxx >= dk(kk)) exit
        dist = max(dxx, abs(ys(l) - ys(q))); call ins(dist); l = l - 1
      end do
      l = q + 1
      do while (l <= n)
        dxx = xs(l) - xs(q)
        if (dxx >= dk(kk)) exit
        dist = max(dxx, abs(ys(l) - ys(q))); call ins(dist); l = l + 1
      end do
      eps = dk(kk)
      nx = nless(xs, xs(q) + eps) - nleq(xs, xs(q) - eps) - 1
      ny = nless(ysort, ys(q) + eps) - nleq(ysort, ys(q) - eps) - 1
      s = s + digamma(nx + 1.0_dp) + digamma(ny + 1.0_dp)
    end do
    mi_ksg = max(0.0_dp, digamma(real(n, dp)) + digamma(real(kk, dp)) - s / n)
  contains
    subroutine ins(v)
      real(dp), intent(in) :: v
      integer :: p
      if (v >= dk(kk)) return
      p = kk
      do while (p > 1)
        if (dk(p - 1) <= v) exit
        dk(p) = dk(p - 1); p = p - 1
      end do
      dk(p) = v
    end subroutine ins
    integer function nless(a, v)        ! #{a < v}, a sorted
      real(dp), intent(in) :: a(:), v
      integer :: lo, hi, mid
      lo = 0; hi = size(a)
      do while (lo < hi)
        mid = (lo + hi + 1) / 2
        if (a(mid) < v) then
          lo = mid
        else
          hi = mid - 1
        end if
      end do
      nless = lo
    end function nless
    integer function nleq(a, v)         ! #{a <= v}
      real(dp), intent(in) :: a(:), v
      integer :: lo, hi, mid
      lo = 0; hi = size(a)
      do while (lo < hi)
        mid = (lo + hi + 1) / 2
        if (a(mid) <= v) then
          lo = mid
        else
          hi = mid - 1
        end if
      end do
      nleq = lo
    end function nleq
  end function mi_ksg

  ! Ferreira (2013) tail-dependence coefficient, max over the four corners
  pure real(dp) function ferreira_tail(n, u, v)
    integer, intent(in) :: n
    real(dp), intent(in) :: u(n), v(n)
    real(dp) :: m(4)
    m(1) = sum(max(u, v)) / n
    m(2) = sum(max(v, 1.0_dp - u)) / n
    m(3) = sum(max(1.0_dp - u, 1.0_dp - v)) / n
    m(4) = sum(max(1.0_dp - v, u)) / n
    m = min(max(m, 0.5_dp), 2.0_dp / 3.0_dp)
    ferreira_tail = 3.0_dp - 1.0_dp / (1.0_dp - minval(m))
  end function ferreira_tail

  ! tree criterion with optional weights (tau, rho, joe weighted; others unweighted)
  real(dp) function criterion_w(n, x, y, kind, usew, w)
    integer, intent(in) :: n, kind
    real(dp), intent(in) :: x(n), y(n), w(n)
    logical, intent(in) :: usew
    real(dp), allocatable :: rx(:), ry(:)
    real(dp) :: r
    integer :: i
    if (.not. usew .or. kind > 2 .or. n <= 10) then
      criterion_w = criterion(n, x, y, kind); return
    end if
    allocate(rx(n), ry(n))
    select case (kind)
    case (1)
      call rank_vec(n, x, rx); call rank_vec(n, y, ry)
      criterion_w = crit_ranks(n, rx, ry, 1, .true., w); return
    case (2)
      do i = 1, n
        rx(i) = qnorm(clampu(x(i))); ry(i) = qnorm(clampu(y(i)))
      end do
      r = wpearson(n, rx, ry, w)
      r = -0.5_dp * log(1.0_dp - min(max(r * r, 0.0_dp), 1.0_dp - 1.0e-15_dp))
    case default
      r = tau_w(n, x, y, w)
    end select
    if (r /= r) r = 0.0_dp
    criterion_w = abs(r)
  end function criterion_w

  ! ------------------------------------------------- C API: single pair copula
  ! Fit a given family/rotation (u torch (n, 2)); ci as in edge_cands.
  subroutine fv_bicop_fit(n, u, fam, rot, ci, par, ll, nthreads) bind(C, name="fv_bicop_fit")
    integer(c_int), value :: n, fam, rot, nthreads
    real(c_double), intent(in) :: u(2, n)
    integer(c_int), intent(in) :: ci(6)
    real(c_double), intent(out) :: par(3), ll
    real(dp) :: wd(1), ax(NAUX)
    call bicop_fit_x(n, u, fam, rot, ci, .false., wd, par, ll, ax, nthreads)
  end subroutine fv_bicop_fit

  ! fit one family/rotation with optional weights; ax receives a nonparametric estimate
  subroutine bicop_fit_x(n, u, fam, rot, ci, usew, w, par, ll, ax, nthreads)
    integer, intent(in) :: n, fam, rot, nthreads, ci(6)
    real(dp), intent(in) :: u(2, n), w(:)
    logical, intent(in) :: usew
    real(dp), intent(out) :: par(3), ll
    real(dp), intent(inout) :: ax(:)
    real(dp), allocatable :: a(:), b(:)
    real(dp) :: x0(3), lb(3), ub(3), tau
    integer :: i, cm
    allocate(a(n), b(n))
    do i = 1, n
      call rotate(rot, clampu(u(1, i)), clampu(u(2, i)), a(i), b(i))
    end do
    if (usew) then
      tau = tau_w(n, a, b, w)
    else
      tau = tau_est(n, a, b, ci(6))
    end if
    call cand_spec(fam, tau, ci, cm, x0, lb, ub)
    if (cm < 0) then
      par = 0.0_dp; ll = 0.0_dp
    else
      call fit_one(n, u, fam, rot, cm, x0, lb, ub, usew, w, par, ll, max(1, nthreads), ax=ax)
    end if
  end subroutine bicop_fit_x

  ! Select family/rotation/parameters for one pair (vinecopulib Bicop::select).
  subroutine fv_bicop_select(n, u, nf, famset, ci, psi0, fam, rot, par, ll, nthreads) &
      bind(C, name="fv_bicop_select")
    integer(c_int), value :: n, nf, nthreads
    real(c_double), intent(in) :: u(2, n)
    integer(c_int), intent(in) :: famset(nf), ci(6)
    real(c_double), value :: psi0
    integer(c_int), intent(out) :: fam, rot
    real(c_double), intent(out) :: par(3), ll
    real(dp) :: wd(1), ax(NAUX)
    call bicop_select_x(n, u, nf, famset, ci, psi0, .false., wd, fam, rot, par, ll, ax, nthreads)
  end subroutine fv_bicop_select

  subroutine bicop_select_x(n, u, nf, famset, ci, psi0, usew, w, fam, rot, par, ll, ax, nthreads)
    integer, intent(in) :: n, nf, nthreads, famset(nf), ci(6)
    real(dp), intent(in) :: u(2, n), psi0, w(:)
    logical, intent(in) :: usew
    integer, intent(out) :: fam, rot
    real(dp), intent(out) :: par(3), ll
    real(dp), intent(inout) :: ax(:)
    integer :: nc, k, nouter, na
    integer, target :: pool(2)
    integer, allocatable :: cf(:), cr(:), cm(:)
    real(dp), allocatable :: cx0(:, :), clb(:, :), cub(:, :), cpar(:, :), cll(:), cbd(:), a(:), b(:), cax(:, :)
    real(dp) :: best, sc, eb
    allocate(cf(2 * nf + 2), cr(2 * nf + 2), cm(2 * nf + 2), cx0(3, 2 * nf + 2), clb(3, 2 * nf + 2), &
             cub(3, 2 * nf + 2), cpar(3, 2 * nf + 2), cll(2 * nf + 2), cbd(2 * nf + 2), a(n), b(n))
    na = merge(NAUX, 1, any(fam_is_np(famset)))
    allocate(cax(na, 2 * nf + 2)); cax = 0.0_dp
    do k = 1, n
      a(k) = clampu(u(1, k)); b(k) = clampu(u(2, k))
    end do
    if (usew) then
      call edge_cands(n, a, b, nf, famset, ci, .true., tau_w(n, a, b, w), nc, cf, cr, cm, cx0, clb, cub)
    else
      call edge_cands(n, a, b, nf, famset, ci, .false., 0.0_dp, nc, cf, cr, cm, cx0, clb, cub)
    end if
    nouter = max(1, min(nthreads, nc))
    !$ call omp_set_max_active_levels(2)
    ! racing (as in vine_select_x): a screened candidate is polished only if its bound can still
    ! beat the best score achieved so far
    best = -huge(1.0_dp)
    pool = 0
    !$omp parallel do private(k, sc, eb) num_threads(nouter) schedule(dynamic, 1)
    do k = 1, nc
      !$omp atomic
      pool(1) = pool(1) + 1
      !$omp atomic
      pool(2) = pool(2) + 1
      if (cm(k) < 0) then
        cpar(:, k) = 0.0_dp; cll(k) = 0.0_dp; cbd(k) = 0.0_dp
      else
        call fit_one(n, u, cf(k), cr(k), cm(k), cx0(:, k), clb(:, k), cub(:, k), usew, w, &
                     cpar(:, k), cll(k), nthreads, ax=cax(:, k), bnd=cbd(k), pool=pool)
        if (cbd(k) > cll(k)) then
          !$omp atomic read
          eb = best
          if (csc(k, cbd(k)) >= eb) &
            call fit_one(n, u, cf(k), cr(k), cm(k), cx0(:, k), clb(:, k), cub(:, k), usew, w, &
                         cpar(:, k), cll(k), nthreads, ax=cax(:, k), bnd=cbd(k), &
                         xr=[cpar(:, k)], tgt=eb - csc(k, 0.0_dp), pool=pool)
        end if
      end if
      if (finite_(cll(k))) then
        sc = csc(k, cll(k))
        !$omp critical (fv_best)
        best = max(best, sc)
        !$omp end critical (fv_best)
      end if
      !$omp atomic
      pool(1) = pool(1) - 1
      !$omp atomic
      pool(2) = pool(2) - 1
    end do
    !$omp end parallel do
    best = -huge(1.0_dp); fam = 0; rot = 0; par = 0.0_dp; ll = 0.0_dp
    do k = 1, nc
      if (.not. finite_(cll(k))) cycle
      sc = csc(k, cll(k))
      if (sc > best) then
        best = sc; fam = cf(k); rot = cr(k); par = cpar(:, k); ll = cll(k)
        if (fam_is_np(cf(k))) ax = cax(:, k)
      end if
    end do
  contains
    real(dp) function csc(k, v)
      integer, intent(in) :: k
      real(dp), intent(in) :: v
      if (fam_is_np(cf(k))) then
        csc = sel_score(ci(1), v, cf(k), n, psi0, cax(2, k))
      else
        csc = sel_score(ci(1), v, cf(k), n, psi0)
      end if
    end function csc
  end subroutine bicop_select_x

  ! relative fit cost class, used to schedule the slowest jobs first
  pure integer function job_cost(f, m)
    integer, intent(in) :: f, m
    if (m < 0) then
      job_cost = 0
    else if (f == 2 .or. f == 23 .or. fam_is_np(f)) then
      job_cost = 5
    else if (m == 0) then
      job_cost = 1
    else if (f == 10 .or. f == 11 .or. f == 8) then
      job_cost = 4
    else if (f >= 7) then
      job_cost = 3     ! BB, tawn and the extended families
    else
      job_cost = 2
    end if
  end function job_cost

  pure real(dp) function sel_score(selcrit, ll, f, n, psi0, kef)
    integer, intent(in) :: selcrit, f, n
    real(dp), intent(in) :: ll, psi0
    real(dp), intent(in), optional :: kef        ! effective npars (nonparametric families)
    real(dp) :: k
    k = npar_of(f)
    if (present(kef)) k = kef
    select case (selcrit)
    case (0); sel_score = ll
    case (1); sel_score = ll - k
    case (2); sel_score = ll - 0.5_dp * log(real(n, dp)) * k
    case default
      if (f == 0) then
        sel_score = ll + log(1.0_dp - psi0)
      else
        sel_score = ll - 0.5_dp * log(real(n, dp)) * k + log(psi0)
      end if
    end select
  end function sel_score

  ! ------------------------------------------------------- C API: vine select
  ! u: torch (n, d) -> Fortran u(d, n).  famset: family codes.
  ! ci(1..10): selcrit, method, preselect, allow_rot, aggressive, tau_est,
  !            tree_crit (0 tau,1 rho,2 joe,3 hoeffd), tree_algo (0 prim,1 kruskal),
  !            trunc_lvl (-1 = d-1), has_struct (fixed structure given)
  ! cr(1..2): threshold, psi0
  ! order_io / struct_in: fixed structure (natural-order labels) when ci(10)=1.
  ! Outputs (torch (d, d) [tree, edge] -> Fortran (edge, tree)):
  !   order_io(d), struct_out, fam_out, rot_out, par_out(2, d, d), ll_out(d, d); info=0 ok.
  subroutine fv_vine_select(n, d, u, nf, famset, ci, cr, order_io, struct_in, struct_out, &
                             fam_out, rot_out, par_out, ll_out, nthreads, info) &
                             bind(C, name="fv_vine_select")
    integer(c_int), value :: n, d, nf, nthreads
    real(c_double), intent(in) :: u(d, n), cr(2)
    integer(c_int), intent(in) :: famset(nf), ci(10), struct_in(d, d)
    integer(c_int), intent(inout) :: order_io(d)
    integer(c_int), intent(out) :: struct_out(d, d), fam_out(d, d), rot_out(d, d), info
    real(c_double), intent(out) :: par_out(3, d, d), ll_out(d, d)
    real(dp) :: wd(1), umd(1, 1), axd(1, d, d)
    logical :: dd(d)
    integer :: tr
    dd = .false.
    call vine_select_x(n, d, u, nf, famset, ci, cr, order_io, struct_in, struct_out, fam_out, rot_out, &
                       par_out, ll_out, nthreads, info, .false., wd, dd, umd, axd, .false., tr)
  end subroutine fv_vine_select

  ! Dissmann selection with observation weights (usew, w(n)), discrete variables
  ! (disc(d), um(d, n) = F(x-) in original order), nonparametric families (aux_out
  ! (NAUX, d, d) when the family set contains 24/25, else size 1) and mBICv-chosen
  ! truncation (seltrunc; the level found is returned in trunc_out).
  subroutine vine_select_x(n, d, u, nf, famset, ci, cr, order_io, struct_in, struct_out, &
                           fam_out, rot_out, par_out, ll_out, nthreads, info, usew, w, disc, um, aux_out, &
                           seltrunc, trunc_out, crit_out)
    integer, intent(in) :: n, d, nf, nthreads
    real(dp), intent(in) :: u(d, n), cr(2), w(:), um(:, :)
    integer, intent(in) :: famset(nf), ci(10), struct_in(d, d)
    integer, intent(inout) :: order_io(d)
    integer, intent(out) :: struct_out(d, d), fam_out(d, d), rot_out(d, d), info, trunc_out
    real(dp), intent(out) :: par_out(3, d, d), ll_out(d, d), aux_out(:, :, :)
    real(dp), intent(out), optional :: crit_out(d, d)     ! tree criterion of every edge (as par_out)
    logical, intent(in) :: usew, disc(d), seltrunc
    logical :: anyd, utau
    integer :: nax
    real(dp), allocatable :: colsm(:, :), vh1m(:, :), vh2m(:, :), pam(:, :), pbm(:, :), eh1m(:, :), eh2m(:, :)
    real(dp), allocatable :: E_aux(:, :, :), jaux(:, :, :)
    logical, allocatable :: vd1(:), vd2(:), eda(:), edb(:)
    real(dp) :: q1, q2, crk, cri, llt, npt, ps
    integer :: qn

    integer :: trunc, t, nv, i, j, k, ne, ncand, col, tt, pos, kk, e, f2
    integer, allocatable :: E_u(:, :), E_v(:, :), E_cond(:, :, :), E_all(:, :, :), E_nall(:, :)
    integer, allocatable :: E_cing(:, :, :), E_ncing(:, :), E_fam(:, :), E_rot(:, :)
    real(dp), allocatable :: E_par(:, :, :), E_ll(:, :), E_crit(:, :)
    logical, allocatable :: alive(:, :)
    real(dp), allocatable :: cols(:, :), vh1(:, :), vh2(:, :), eh1(:, :), eh2(:, :), pa(:, :), pb(:, :)
    integer, allocatable :: vprev(:, :), vall(:, :), vnall(:), vcond(:, :), vncond(:)
    integer, allocatable :: cu(:), cv(:), keepmask(:), mat(:, :), order0(:), ning(:), chk(:), mins(:, :)
    real(dp), allocatable :: cw(:), ccrit(:), taum(:, :)
    integer :: nning, nchk
    logical :: unknown, found
    real(dp) :: lpd, h1v, h2v, thr, psi0e
    integer :: civ(6), mc, nj, q, nouter
    integer, target :: pool(2)        ! racing pool: threads in use, fits running (grab / release)
    integer, allocatable :: jnc(:), jf(:, :), jr(:, :), jm(:, :), je(:), jk(:)
    real(dp), allocatable :: ptr(:, :, :), vr1(:, :), vr2(:, :)
    real(dp), allocatable :: jx0(:, :, :), jlb(:, :, :), jub(:, :, :), jpar(:, :, :), jll(:, :), jbd(:, :), ebest(:)
    real(dp), allocatable :: gsc(:, :)               ! generic pre-screen scores (selection criterion)
    logical, allocatable :: skipj(:, :)              ! candidate left out by the pre-screen
    integer, parameter :: KEEPG = 4                  ! generic candidates per edge that are really fitted
    real(dp) :: best, sc, eb, ppe(NPP)
    real(dp) :: tp(6), tq
    logical :: prof
    character(len=8) :: envv
    integer :: envst

    info = 0
    tp = 0.0_dp; tq = 0.0_dp
    anyd = any(disc)
    utau = .not. usew            ! unweighted: the fast tau-matrix path is available
    nax = size(aux_out, 1)
    trunc_out = 0
    call get_environment_variable("FORTVINE_PROFILE", envv, status=envst)
    prof = (envst == 0)
    thr = cr(1)
    trunc = ci(9)
    if (trunc < 0 .or. trunc > d - 1) trunc = d - 1
    unknown = ci(10) == 0
    allocate(cols(n, d))
    cols = transpose(u)
    do j = 1, d
      do i = 1, n
        cols(i, j) = clampu(cols(i, j))
      end do
    end do
    if (anyd) then
      allocate(colsm(n, d))
      do j = 1, d
        do i = 1, n
          if (disc(j)) then
            colsm(i, j) = clampu(um(j, i))
          else
            colsm(i, j) = cols(i, j)
          end if
        end do
      end do
    end if
    if (unknown) then
      do j = 1, d
        order_io(j) = j
      end do
    end if
    allocate(mins(d, d)); mins = 0
    if (.not. unknown) then
      ! min_array of the fixed structure (RVineStructure.compute_min_array)
      mins = struct_in
      do j = 1, d - 1
        do i = 2, min(d - j, trunc)
          mins(j, i) = min(struct_in(j, i), mins(j, i - 1))
        end do
      end do
    end if

    allocate(E_u(d, d), E_v(d, d), E_cond(2, d, d), E_all(d, d, d), E_nall(d, d), E_cing(d, d, d), &
             E_ncing(d, d), E_fam(d, d), E_rot(d, d), E_par(3, d, d), E_ll(d, d), E_crit(d, d), alive(d, d))
    allocate(E_aux(nax, d, d)); E_aux = 0.0_dp
    E_fam = 0; E_rot = 0; E_par = 0.0_dp; E_ll = 0.0_dp; alive = .false.; E_nall = 0; E_ncing = 0

    ! base: vertices = variables (natural order), prev = [root, i]
    nv = d
    allocate(vh1(n, d), vh2(n, d), vprev(2, d), vall(d, d), vnall(d), vcond(2, d), vncond(d))
    allocate(vd1(d), vd2(d), vh1m(merge(n, 1, anyd), d), vh2m(merge(n, 1, anyd), d))
    do i = 1, d
      vd1(i) = disc(order_io(i)); vd2(i) = vd1(i)
      if (anyd) then
        vh1m(:, i) = colsm(:, order_io(i)); vh2m(:, i) = vh1m(:, i)
      end if
      vh1(:, i) = cols(:, order_io(i)); vh2(:, i) = vh1(:, i)
      vprev(:, i) = [d, i - 1]
      vall(1, i) = order_io(i) - 1; vnall(i) = 1
      vcond(1, i) = order_io(i) - 1; vncond(i) = 1
    end do

    if (unknown .and. ci(7) == 0 .and. utau) then
      allocate(taum(d, d))
      call fv_tau_matrix(n, d, u, taum, nthreads)
    end if

    do t = 1, trunc
      !$ tq = omp_get_wtime()
      ! ranks of every vertex's two pseudo-observation columns, once per tree (rank-based criteria)
      if (allocated(vr1)) deallocate(vr1, vr2)
      if (ci(7) == 1 .or. ci(7) == 3) then
        allocate(vr1(n, nv), vr2(n, nv))
        !$omp parallel do private(i) num_threads(nthreads) schedule(dynamic, 1)
        do i = 1, nv
          call rank_vec(n, vh1(:, i), vr1(:, i))
          if (t == 1) then
            vr2(:, i) = vr1(:, i)
          else
            call rank_vec(n, vh2(:, i), vr2(:, i))
          end if
        end do
        !$omp end parallel do
      end if
      ! ---------- candidate edges among current vertices
      ncand = 0
      allocate(cu(nv * nv), cv(nv * nv), cw(nv * nv), ccrit(nv * nv), keepmask(nv * nv))
      if (unknown) then
        do i = 0, nv - 1
          do j = 0, i - 1
            if (common_nb(i + 1, j + 1) < 0) cycle
            ncand = ncand + 1; cu(ncand) = i; cv(ncand) = j
          end do
        end do
        !$omp parallel do private(k) num_threads(nthreads) schedule(dynamic, 1)
        do k = 1, ncand
          if (t == 1 .and. ci(7) == 0 .and. utau) then
            ccrit(k) = abs(taum(vall(1, cu(k) + 1) + 1, vall(1, cv(k) + 1) + 1))
          else
            ccrit(k) = edge_crit(cu(k) + 1, cv(k) + 1)
          end if
        end do
        !$omp end parallel do
        do k = 1, ncand
          if (ccrit(k) >= thr) then
            cw(k) = 1.0_dp - ccrit(k)
          else
            cw(k) = 1.0_dp
          end if
        end do
        if (nv > 2) then
          if (ci(8) == 1) then
            call mst_kruskal(nv, ncand, cu, cv, cw, keepmask)
          else
            call mst_prim(nv, ncand, cu, cv, cw, keepmask)
          end if
        else
          keepmask(1:ncand) = 1
        end if
      else
        do i = 0, nv - 2
          ncand = ncand + 1; cu(ncand) = i; cv(ncand) = mins(i + 1, t) - 1
        end do
        do k = 1, ncand
          ccrit(k) = edge_crit(cu(k) + 1, cv(k) + 1)
        end do
        keepmask(1:ncand) = 1
      end if

      ! ---------- kept edges, in insertion order
      ne = 0
      do k = 1, ncand
        if (keepmask(k) == 0) cycle
        ne = ne + 1
        E_u(ne, t) = cu(k); E_v(ne, t) = cv(k); E_crit(ne, t) = ccrit(k)
        alive(ne, t) = .true.
      end do
      if (ne /= nv - 1) then
        info = 2; return
      end if
      !$ tp(1) = tp(1) + omp_get_wtime() - tq; tq = omp_get_wtime()
      allocate(pa(n, ne), pb(n, ne), eh1(n, ne), eh2(n, ne), eda(ne), edb(ne))
      allocate(pam(merge(n, 1, anyd), ne), pbm(merge(n, 1, anyd), ne), eh1m(merge(n, 1, anyd), ne), &
               eh2m(merge(n, 1, anyd), ne))
      do e = 1, ne
        call edge_info(e, t)
      end do

      ! ---------- pair copulas: candidates per edge, then all (edge, family) fits in parallel
      psi0e = cr(2)
      civ = ci(1:6)
      if (ci(1) == 4) then
        civ(1) = 3; psi0e = cr(2)**t
      end if
      mc = 2 * nf + 2
      allocate(jnc(ne), jf(mc, ne), jr(mc, ne), jm(mc, ne), jx0(3, mc, ne), jlb(3, mc, ne), &
               jub(3, mc, ne), jpar(3, mc, ne), jll(mc, ne), jbd(mc, ne), ebest(ne), jaux(nax, mc, ne))
      !$omp parallel do private(e, k, q) num_threads(nthreads) schedule(dynamic, 1)
      do e = 1, ne
        if (E_crit(e, t) < thr) then
          jnc(e) = 0
        else
          if (usew) then
            call edge_cands(n, pa(:, e), pb(:, e), nf, famset, civ, .true., tau_w(n, pa(:, e), pb(:, e), w), &
                            jnc(e), jf(:, e), jr(:, e), jm(:, e), jx0(:, :, e), jlb(:, :, e), jub(:, :, e))
          else
            call edge_cands(n, pa(:, e), pb(:, e), nf, famset, civ, t == 1 .and. unknown .and. ci(7) == 0, &
                            taum_at(e), jnc(e), jf(:, e), jr(:, e), jm(:, e), jx0(:, :, e), jlb(:, :, e), &
                            jub(:, :, e))
          end if
          ! (nonparametric families on discrete pairs: fitted to mid-ranks, see fit_one_raw)
        end if
      end do
      !$omp end parallel do
      !$ tp(2) = tp(2) + omp_get_wtime() - tq; tq = omp_get_wtime()
      ! transforms of every continuous edge, shared by all its one-parameter candidates (mle1)
      allocate(ptr(n, 10, ne))
      !$omp parallel do private(e) num_threads(nthreads) schedule(dynamic, 1)
      do e = 1, ne
        if (jnc(e) > 0 .and. .not. (eda(e) .or. edb(e))) call pc_transforms(n, pa(:, e), pb(:, e), ptr(:, :, e))
      end do
      !$omp end parallel do
      ! generic families: every candidate is pre-screened on a ~NSUB subsample (one Newton fit),
      ! and only the KEEPG best per edge (by the selection criterion of the scaled log-likelihood)
      ! enter the full-data racing below -- most never touch all n rows
      allocate(skipj(size(jf, 1), ne)); skipj = .false.
      if (n > 3 * NSUB .and. any(gen_is(jf))) then
        allocate(gsc(size(jf, 1), ne)); gsc = -huge(1.0_dp)
        !$omp parallel do private(e, k, sc) num_threads(nthreads) schedule(dynamic, 1)
        do e = 1, ne
          if (eda(e) .or. edb(e)) cycle
          do k = 1, jnc(e)
            if (.not. gen_is(jf(k, e)) .or. jm(k, e) <= 0) cycle
            call gen_prescreen(n, pa(:, e), pb(:, e), jf(k, e), jr(k, e), jlb(:, k, e), jub(:, k, e), usew, w, sc)
            if (sc > -huge(1.0_dp)) gsc(k, e) = sel_score(civ(1), sc, jf(k, e), n, psi0e)
          end do
        end do
        !$omp end parallel do
        do e = 1, ne
          if (eda(e) .or. edb(e)) cycle
          q = 0
          do k = 1, jnc(e)
            if (.not. gen_is(jf(k, e)) .or. jm(k, e) <= 0) cycle
            ! rank of candidate k among the generic candidates of edge e
            if (count(gsc(1:jnc(e), e) > gsc(k, e)) >= KEEPG .or. gsc(k, e) <= -huge(1.0_dp)) then
              skipj(k, e) = .true.; jll(k, e) = -huge(1.0_dp); jpar(:, k, e) = 0.0_dp
            end if
          end do
        end do
      end if
      nj = sum(jnc) - count(skipj)                ! pre-screened-out candidates are not jobs
      allocate(je(max(nj, 1)), jk(max(nj, 1)))
      q = 0
      ! the cheap closed-form / one-parameter fits first (they set each edge's racing bar within
      ! a millisecond), then longest-processing-time-first: the slowest fits start first so no
      ! single late job stretches the tree's critical path
      do kk = 2, 0, -1
        do e = 1, ne
          do k = 1, jnc(e)
            if (job_cost(jf(k, e), jm(k, e)) /= kk) cycle
            if (skipj(k, e)) cycle
            q = q + 1; je(q) = e; jk(q) = k
          end do
        end do
      end do
      do kk = 5, 3, -1
        do e = 1, ne
          do k = 1, jnc(e)
            if (job_cost(jf(k, e), jm(k, e)) /= kk) cycle
            if (skipj(k, e)) cycle
            q = q + 1; je(q) = e; jk(q) = k
          end do
        end do
      end do
      nouter = max(1, min(nthreads, nj))
      pool = 0
      !$ call omp_set_max_active_levels(2)
      ! racing (no barrier): a screened candidate is polished only if its bound can still beat the
      ! best score its edge has achieved so far, and the polish gives up once it no longer can
      ! (ebest only rises and never exceeds the final best, so the selected candidate is always
      ! polished to the end)
      ebest = -huge(1.0_dp)
      !$omp parallel do private(q, e, k, sc, eb) num_threads(nouter) schedule(dynamic, 1)
      do q = 1, nj
        e = je(q); k = jk(q)
        !$omp atomic
        pool(1) = pool(1) + 1
        !$omp atomic
        pool(2) = pool(2) + 1
        call fit_job(q, nthreads, jbd(k, e))
        if (jbd(k, e) > jll(k, e)) then
          !$omp atomic read
          eb = ebest(e)
          if (jsc(k, e, jbd(k, e)) >= eb) &
            call fit_job(q, nthreads, jbd(k, e), xr=[jpar(:, k, e)], tgt=eb - jsc(k, e, 0.0_dp))
        end if
        if (finite_(jll(k, e))) then
          sc = jsc(k, e, jll(k, e))
          !$omp critical (fv_ebest)
          ebest(e) = max(ebest(e), sc)
          !$omp end critical (fv_ebest)
        end if
        !$omp atomic
        pool(1) = pool(1) - 1
        !$omp atomic
        pool(2) = pool(2) - 1
      end do
      !$omp end parallel do
      if (fv_abort) then                         ! host abort: the result is discarded anyway
        trunc = t - 1
        exit
      end if
      !$ tp(3) = tp(3) + omp_get_wtime() - tq; tq = omp_get_wtime()
      do e = 1, ne
        E_fam(e, t) = 0; E_rot(e, t) = 0; E_par(:, e, t) = 0.0_dp; E_ll(e, t) = 0.0_dp
        best = -huge(1.0_dp)
        do k = 1, jnc(e)
          if (.not. finite_(jll(k, e))) cycle
          if (fam_is_np(jf(k, e))) then
            sc = sel_score(civ(1), jll(k, e), jf(k, e), n, psi0e, jaux(2, k, e))
          else
            sc = sel_score(civ(1), jll(k, e), jf(k, e), n, psi0e)
          end if
          if (sc > best) then
            best = sc
            E_fam(e, t) = jf(k, e); E_rot(e, t) = jr(k, e)
            E_par(:, e, t) = jpar(:, k, e); E_ll(e, t) = jll(k, e)
            if (fam_is_np(jf(k, e))) E_aux(:, e, t) = jaux(:, k, e)
          end if
        end do
      end do
      ! mBICv-based truncation (vinecopulib select_trunc_lvl): keep tree t only if its
      ! mBICv contribution beats an all-independence tree t
      if (seltrunc) then
        llt = 0.0_dp; npt = 0.0_dp; qn = 0
        do e = 1, ne
          llt = llt + E_ll(e, t)
          if (fam_is_np(E_fam(e, t))) then
            npt = npt + E_aux(2, e, t)
          else
            npt = npt + npar_of(E_fam(e, t))
          end if
          if (E_fam(e, t) /= 0) qn = qn + 1
        end do
        ps = min(max(cr(2)**t, 1.0e-12_dp), 1.0_dp - 1.0e-12_dp)
        crk = -2.0_dp * llt + log(real(n, dp)) * npt - 2.0_dp * (qn * log(ps) + (ne - qn) * log(1.0_dp - ps))
        cri = -2.0_dp * ne * log(1.0_dp - ps)
        if (crk >= cri) then
          trunc = t - 1
          exit
        end if
      end if
      !$omp parallel do private(e, i, lpd, h1v, h2v, ppe, q1, q2) num_threads(nthreads) schedule(dynamic, 1)
      do e = 1, ne
        ppe = prep_par(E_fam(e, t), E_par(:, e, t))
        do i = 1, n
          if (eda(e) .or. edb(e)) then
            call eval_disc(E_fam(e, t), E_rot(e, t), ppe, E_aux(:, e, t), pa(i, e), pb(i, e), pam(i, e), &
                           pbm(i, e), eda(e), edb(e), lpd, h1v, h2v, q1, q2)
            eh1(i, e) = h1v; eh2(i, e) = h2v; eh1m(i, e) = q1; eh2m(i, e) = q2
          else if (E_fam(e, t) == 0) then
            eh1(i, e) = pb(i, e); eh2(i, e) = pa(i, e)
          else
            call eval_x(E_fam(e, t), E_rot(e, t), ppe, E_aux(:, e, t), &
                        pa(i, e), pb(i, e), .true., lpd, h1v, h2v)
            eh1(i, e) = h1v; eh2(i, e) = h2v
          end if
          if (anyd .and. .not. (eda(e) .or. edb(e))) then
            eh1m(i, e) = eh1(i, e); eh2m(i, e) = eh2(i, e)
          end if
        end do
      end do
      !$omp end parallel do
      deallocate(jnc, jf, jr, jm, jx0, jlb, jub, jpar, jll, jbd, ebest, je, jk, jaux, ptr)
      deallocate(skipj)
      if (allocated(gsc)) deallocate(gsc)
      !$ tp(4) = tp(4) + omp_get_wtime() - tq
#ifndef FORTVINE_NO_IO
      if (prof) print '(a,i3,a,i5,a,4f9.1)', '  [fortvine] tree', t, ' jobs', nj, &
        '  ms crit/info+cands/fits/h:', 1e3 * tp(1:4)
#endif
      tp(1:4) = 0.0_dp

      ! ---------- edges become the next tree's vertices
      deallocate(vh1, vh2, vprev, vall, vnall, vcond, vncond, vd1, vd2, vh1m, vh2m)
      nv = ne
      allocate(vh1(n, nv), vh2(n, nv), vprev(2, nv), vall(d, nv), vnall(nv), vcond(2, nv), vncond(nv))
      allocate(vd1(nv), vd2(nv), vh1m(merge(n, 1, anyd), nv), vh2m(merge(n, 1, anyd), nv))
      do e = 1, ne
        ! h1 = F(b | a) is discrete iff b is; h2 = F(a | b) iff a is
        vd1(e) = edb(e); vd2(e) = eda(e)
        if (anyd) then
          vh1m(:, e) = eh1m(:, e); vh2m(:, e) = eh2m(:, e)
        end if
        vh1(:, e) = eh1(:, e); vh2(:, e) = eh2(:, e)
        vprev(:, e) = [E_u(e, t), E_v(e, t)]
        vnall(e) = E_nall(e, t); vall(1:vnall(e), e) = E_all(1:vnall(e), e, t)
        vcond(:, e) = E_cond(:, e, t); vncond(e) = 2
      end do
      deallocate(cu, cv, cw, ccrit, keepmask, pa, pb, eh1, eh2, eda, edb, pam, pbm, eh1m, eh2m)
    end do
    trunc_out = trunc

    ! ================= finalize: trees -> natural-order R-vine matrix
    trunc_out = trunc
    struct_out = 0; fam_out = 0; rot_out = 0; par_out = 0.0_dp; ll_out = 0.0_dp; aux_out = 0.0_dp
    if (present(crit_out)) crit_out = 0.0_dp
    if (trunc == 0) then
      do j = 1, d
        order_io(j) = j
      end do
      return
    end if
    allocate(mat(d, d), order0(d), ning(d), chk(d))
    mat = -1
    do col = 0, d - 2
      tt = max(min(trunc, d - 1 - col), 1)
      found = .false.
      do e = 1, d - tt
        if (.not. alive(e, tt)) cycle
        if (min(degree(E_u(e, tt), tt), degree(E_v(e, tt), tt)) > 1) cycle
        found = .true.; exit
      end do
      if (.not. found) then
        info = 3; return
      end if
      if (degree(E_v(e, tt), tt) == 1) then
        pos = 1
      else
        pos = 0
      end if
      if (pos == 1) call flip(e, tt)
      order0(col + 1) = E_cond(pos + 1, e, tt)
      call put(tt - 1, col, e, tt, E_cond(2 - pos, e, tt))
      nning = E_ncing(e, tt); ning(1:nning) = E_cing(1:nning, e, tt)
      alive(e, tt) = .false.
      do kk = 1, tt - 1
        nchk = nning + 1
        chk(1) = order0(col + 1); chk(2:nchk) = ning(1:nning)
        found = .false.
        do f2 = 1, d - (tt - kk)
          if (.not. alive(f2, tt - kk)) cycle
          if (same_set(E_all(1:E_nall(f2, tt - kk), f2, tt - kk), E_nall(f2, tt - kk), chk, nchk)) then
            found = .true.; exit
          end if
        end do
        if (.not. found) then
          info = 4; return
        end if
        if (order0(col + 1) == E_cond(2, f2, tt - kk)) then
          pos = 1
        else
          pos = 0
        end if
        if (pos == 1) call flip(f2, tt - kk)
        call put(tt - kk - 1, col, f2, tt - kk, E_cond(2 - pos, f2, tt - kk))
        nning = E_ncing(f2, tt - kk); ning(1:nning) = E_cing(1:nning, f2, tt - kk)
        alive(f2, tt - kk) = .false.
      end do
    end do
    order0(d) = mat(d - 1, 1)
    ! 1-based labels, then relabel to natural order: new = argsort(order)[old-1] + 1
    do j = 1, d
      order_io(j) = order0(j) + 1
    end do
    if (allocated(keepmask)) deallocate(keepmask)
    allocate(keepmask(d))   ! reuse as argsort(order)
    do j = 1, d
      keepmask(order_io(j)) = j
    end do
    do t = 1, min(d - 1, trunc)
      do j = 1, d - t
        struct_out(j, t) = keepmask(mat(j, t) + 1)
      end do
    end do
  contains
    ! one (edge, family) fit of the tree's job list on nin inner threads
    subroutine fit_job(q, nin, bnd, xr, tgt)
      integer, intent(in) :: q, nin
      real(dp), intent(out), optional :: bnd        ! racing screen (fit_one); = jll for a final fit
      real(dp), intent(in), optional :: xr(3), tgt  ! polish from the screen point, give up below tgt
      real(dp), allocatable :: ju(:, :), jum(:, :)
      logical :: edisc, okm
      real(dp) :: tjob
      if (jm(jk(q), je(q)) == 1 .and. any(jf(jk(q), je(q)) == [1, 3, 4, 5, 6]) .and. .not. (eda(je(q)) .or. edb(je(q)))) then
        call mle1(jf(jk(q), je(q)), jr(jk(q), je(q)), n, pa(:, je(q)), pb(:, je(q)), ptr(:, :, je(q)), usew, w, &
                  jx0(:, jk(q), je(q)), jlb(:, jk(q), je(q)), jub(:, jk(q), je(q)), jpar(:, jk(q), je(q)), jll(jk(q), je(q)), okm)
        if (okm) then
          if (present(bnd)) bnd = jll(jk(q), je(q))
          return
        end if
      end if
      if (present(bnd)) bnd = 0.0_dp
      if (jm(jk(q), je(q)) < 0) then
        jpar(:, jk(q), je(q)) = 0.0_dp; jll(jk(q), je(q)) = 0.0_dp
      else
        allocate(ju(2, n))
        ju(1, :) = pa(:, je(q)); ju(2, :) = pb(:, je(q))
        edisc = eda(je(q)) .or. edb(je(q))
        if (edisc) then
          allocate(jum(2, n))
          jum(1, :) = pam(:, je(q)); jum(2, :) = pbm(:, je(q))
          call fit_one(n, ju, jf(jk(q), je(q)), jr(jk(q), je(q)), jm(jk(q), je(q)), &
                       jx0(:, jk(q), je(q)), jlb(:, jk(q), je(q)), jub(:, jk(q), je(q)), &
                       usew, w, jpar(:, jk(q), je(q)), jll(jk(q), je(q)), nin, &
                       ax=jaux(:, jk(q), je(q)), um=jum, dsc=[eda(je(q)), edb(je(q))], pool=pool)
          if (present(bnd)) bnd = jll(jk(q), je(q))
          deallocate(jum)
        else
          !$ tjob = omp_get_wtime()
          call fit_one(n, ju, jf(jk(q), je(q)), jr(jk(q), je(q)), jm(jk(q), je(q)), &
                       jx0(:, jk(q), je(q)), jlb(:, jk(q), je(q)), jub(:, jk(q), je(q)), &
                       usew, w, jpar(:, jk(q), je(q)), jll(jk(q), je(q)), nin, &
                       ax=jaux(:, jk(q), je(q)), bnd=bnd, xr=xr, tgt=tgt, pool=pool)
#ifndef FORTVINE_NO_IO
          !$ if (prof .and. omp_get_wtime() - tjob > 0.2_dp) print '(a,i3,a,i3,a,f8.3,a,3f9.4)', '  [fortvine] slow job fam', &
          !$   jf(jk(q), je(q)), ' rot', jr(jk(q), je(q)), ' s', omp_get_wtime() - tjob, ' par', jpar(:, jk(q), je(q))
#endif
        end if
        deallocate(ju)
      end if
    end subroutine fit_job

    real(dp) function jsc(k, e, v)
      integer, intent(in) :: k, e
      real(dp), intent(in) :: v
      if (fam_is_np(jf(k, e))) then
        jsc = sel_score(civ(1), v, jf(k, e), n, psi0e, jaux(2, k, e))
      else
        jsc = sel_score(civ(1), v, jf(k, e), n, psi0e)
      end if
    end function jsc

    real(dp) function taum_at(e0)
      integer, intent(in) :: e0
      taum_at = 0.0_dp
      if (t == 1 .and. unknown .and. ci(7) == 0 .and. utau) taum_at = taum(E_cond(1, e0, t) + 1, E_cond(2, e0, t) + 1)
    end function taum_at

    integer function common_nb(i0, i1)
      integer, intent(in) :: i0, i1
      integer :: q, s
      common_nb = -1
      do q = 1, 2
        do s = 1, 2
          if (vprev(q, i0) == vprev(s, i1)) then
            common_nb = vprev(q, i0); return
          end if
        end do
      end do
    end function common_nb

    ! pair data for vertices (i0, i1), 1-based (+ minus data and discreteness flags)
    subroutine pc_data(i0, i1, xa, xb, xam, xbm, fa, fb)
      integer, intent(in) :: i0, i1
      real(dp), intent(out) :: xa(n), xb(n)
      real(dp), intent(out), optional :: xam(:), xbm(:)
      logical, intent(out), optional :: fa, fb
      integer :: cn
      logical :: s0, s1
      cn = common_nb(i0, i1)
      s0 = vprev(1, i0) == cn; s1 = vprev(1, i1) == cn
      if (s0) then
        xa = vh1(:, i0)
      else
        xa = vh2(:, i0)
      end if
      if (s1) then
        xb = vh1(:, i1)
      else
        xb = vh2(:, i1)
      end if
      if (.not. present(fa)) return
      fa = merge(vd1(i0), vd2(i0), s0); fb = merge(vd1(i1), vd2(i1), s1)
      if (anyd) then
        if (s0) then
          xam = vh1m(:, i0)
        else
          xam = vh2m(:, i0)
        end if
        if (s1) then
          xbm = vh1m(:, i1)
        else
          xbm = vh2m(:, i1)
        end if
      end if
    end subroutine pc_data

    real(dp) function edge_crit(i0, i1)
      integer, intent(in) :: i0, i1
      real(dp), allocatable :: xa(:), xb(:)
      integer :: cn
      allocate(xa(n), xb(n))
      if (allocated(vr1)) then
        cn = common_nb(i0, i1)
        if (vprev(1, i0) == cn) then
          xa = vr1(:, i0)
        else
          xa = vr2(:, i0)
        end if
        if (vprev(1, i1) == cn) then
          xb = vr1(:, i1)
        else
          xb = vr2(:, i1)
        end if
        edge_crit = crit_ranks(n, xa, xb, ci(7), usew .and. ci(7) == 1, w)
        return
      end if
      call pc_data(i0, i1, xa, xb)
      edge_crit = criterion_w(n, xa, xb, ci(7), usew, w)
    end function edge_crit

    subroutine edge_info(e0, t0)
      integer, intent(in) :: e0, t0
      integer :: i0, i1, q, m, c(d)
      logical :: inb
      i0 = E_u(e0, t0) + 1; i1 = E_v(e0, t0) + 1
      call pc_data(i0, i1, pa(:, e0), pb(:, e0), pam(:, e0), pbm(:, e0), eda(e0), edb(e0))
      ! conditioned = sym-diff (ordered), conditioning = intersection (v0 order)
      m = 0
      do q = 1, vnall(i0)
        if (.not. any(vall(1:vnall(i1), i1) == vall(q, i0))) then
          m = m + 1; c(m) = vall(q, i0)
        end if
      end do
      do q = 1, vnall(i1)
        if (.not. any(vall(1:vnall(i0), i0) == vall(q, i1))) then
          m = m + 1; c(m) = vall(q, i1)
        end if
      end do
      E_cond(1:2, e0, t0) = c(1:2)
      E_ncing(e0, t0) = 0
      do q = 1, vnall(i0)
        inb = any(vall(1:vnall(i1), i1) == vall(q, i0))
        if (inb) then
          E_ncing(e0, t0) = E_ncing(e0, t0) + 1
          E_cing(E_ncing(e0, t0), e0, t0) = vall(q, i0)
        end if
      end do
      E_nall(e0, t0) = 2 + E_ncing(e0, t0)
      E_all(1:2, e0, t0) = c(1:2)
      E_all(3:E_nall(e0, t0), e0, t0) = E_cing(1:E_ncing(e0, t0), e0, t0)
    end subroutine edge_info

    integer function degree(vx, t0)
      integer, intent(in) :: vx, t0
      integer :: q
      degree = 0
      do q = 1, d - t0
        if (alive(q, t0) .and. (E_u(q, t0) == vx .or. E_v(q, t0) == vx)) degree = degree + 1
      end do
    end function degree

    subroutine flip(e0, t0)
      integer, intent(in) :: e0, t0
      real(dp) :: tmp
      if (E_fam(e0, t0) == 11) then   ! TawnBicop::flip swaps psi1, psi2
        tmp = E_par(1, e0, t0); E_par(1, e0, t0) = E_par(2, e0, t0); E_par(2, e0, t0) = tmp
      end if
      if (fam_is_np(E_fam(e0, t0))) call np_flip(E_aux(:, e0, t0))   ! transpose the estimate
      if (gen_is(E_fam(e0, t0))) call gen_flip(E_fam(e0, t0), E_par(:, e0, t0))   ! asymmetric generic
      if (E_rot(e0, t0) == 90) then
        E_rot(e0, t0) = 270
      else if (E_rot(e0, t0) == 270) then
        E_rot(e0, t0) = 90
      end if
    end subroutine flip

    ! mat[row][col] = label; also copy the pair copula into the output slot
    subroutine put(row, cl, e0, t0, label)
      integer, intent(in) :: row, cl, e0, t0, label
      mat(cl + 1, row + 1) = label
      fam_out(cl + 1, row + 1) = E_fam(e0, t0)
      rot_out(cl + 1, row + 1) = E_rot(e0, t0)
      par_out(:, cl + 1, row + 1) = E_par(:, e0, t0)
      ll_out(cl + 1, row + 1) = E_ll(e0, t0)
      aux_out(:, cl + 1, row + 1) = E_aux(:, e0, t0)
      if (present(crit_out)) crit_out(cl + 1, row + 1) = E_crit(e0, t0)
    end subroutine put

    logical function same_set(s1, m1, s2, m2)
      integer, intent(in) :: m1, m2, s1(m1), s2(m2)
      integer :: q
      same_set = (m1 == m2)
      if (.not. same_set) return
      do q = 1, m1
        if (.not. any(s2 == s1(q))) then
          same_set = .false.; return
        end if
      end do
    end function same_set
  end subroutine vine_select_x

  ! Prim as in torchvine (_mst_prim): deterministic tie-breaking by vertex index.
  subroutine mst_prim(nv, ne, eu, ev, ew, keep)
    integer, intent(in) :: nv, ne, eu(ne), ev(ne)
    real(dp), intent(in) :: ew(ne)
    integer, intent(out) :: keep(ne)
    real(dp), allocatable :: wm(:, :), bw(:)
    integer, allocatable :: bp(:)
    logical, allocatable :: intree(:)
    integer :: it, v, i, q, a, b
    real(dp) :: wv
    allocate(wm(nv, nv), bw(nv), bp(nv), intree(nv))
    wm = huge(1.0_dp)
    do q = 1, ne
      wm(eu(q) + 1, ev(q) + 1) = ew(q); wm(ev(q) + 1, eu(q) + 1) = ew(q)
    end do
    bw = huge(1.0_dp); bp = -1; intree = .false.
    bw(1) = 0.0_dp
    do it = 1, nv
      v = -1; wv = huge(1.0_dp)
      do i = 1, nv
        if ((.not. intree(i)) .and. bw(i) < wv) then
          v = i; wv = bw(i)
        end if
      end do
      if (v < 0) exit
      intree(v) = .true.
      do i = 1, nv
        if (intree(i) .or. wm(v, i) == huge(1.0_dp)) cycle
        if (wm(v, i) < bw(i) - 1.0e-15_dp .or. (abs(wm(v, i) - bw(i)) <= 1.0e-15_dp .and. &
            (bp(i) < 0 .or. v - 1 < bp(i)))) then
          bw(i) = wm(v, i); bp(i) = v - 1
        end if
      end do
    end do
    keep = 0
    do q = 1, ne
      a = min(eu(q), ev(q)); b = max(eu(q), ev(q))
      if (bp(b + 1) == a .or. bp(a + 1) == b) keep(q) = 1
    end do
  end subroutine mst_prim

  subroutine mst_kruskal(nv, ne, eu, ev, ew, keep)
    integer, intent(in) :: nv, ne, eu(ne), ev(ne)
    real(dp), intent(in) :: ew(ne)
    integer, intent(out) :: keep(ne)
    integer, allocatable :: idx(:), par(:)
    integer :: i, j, k, ra, rb, got
    allocate(idx(ne), par(nv))
    do i = 1, ne
      idx(i) = i
    end do
    ! insertion sort by (weight, min, max) -- ne is O(d^2)
    do i = 2, ne
      k = idx(i); j = i - 1
      do while (j >= 1)
        if (.not. klt(k, idx(j))) exit
        idx(j + 1) = idx(j); j = j - 1
      end do
      idx(j + 1) = k
    end do
    do i = 1, nv
      par(i) = i
    end do
    keep = 0; got = 0
    do i = 1, ne
      ra = root(eu(idx(i)) + 1); rb = root(ev(idx(i)) + 1)
      if (ra /= rb) then
        par(ra) = rb; keep(idx(i)) = 1; got = got + 1
        if (got == nv - 1) exit
      end if
    end do
  contains
    logical function klt(p, q)
      integer, intent(in) :: p, q
      if (ew(p) /= ew(q)) then
        klt = ew(p) < ew(q)
      else if (min(eu(p), ev(p)) /= min(eu(q), ev(q))) then
        klt = min(eu(p), ev(p)) < min(eu(q), ev(q))
      else
        klt = max(eu(p), ev(p)) < max(eu(q), ev(q))
      end if
    end function klt
    integer function root(x0)
      integer, intent(in) :: x0
      integer :: x
      x = x0
      do while (par(x) /= x)
        par(x) = par(par(x)); x = par(x)
      end do
      root = x
    end function root
  end subroutine mst_kruskal

end module fortvine
