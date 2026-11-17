! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Shahab Aldin Shojaeezadeh. Part of fortvine (see LICENSE). If you use it, please cite:
! Shojaeezadeh, S. A. (2026). fortvine: a complete compendium to work on vine copula modeling.
! fortvine_api.f90 -- typed Fortran API for fortvine.
!
!   use fortvine_api
!   type(fv_vinecop) :: vc
!   type(fv_controls) :: ctl
!   call vc%select(u, ctl)            ! u(n, d) pseudo-observations
!   call vc%pdf(u, dens)              ! dens(n)
!
! Data layout is Fortran-natural: u(n, d), columns = variables.
module fortvine_api
  use iso_c_binding, only: c_int, c_double
  use fortvine, only: NAUX, fam_is_np, gen_name, gen_codes, gen_is, gen_valid, gen_dup, gen_npar, gen_flip, &
    gen_exch, gen_lt, gen_ut, vine_pdf_x, vine_pdf_disc, vine_pdf_disc_x, vine_rosenblatt_x, vine_inverse_rosenblatt_x, &
    vine_select_x, bicop_eval_x, bicop_fit_x, bicop_select_x, fv_kendall_tau_c, fv_tau_matrix, fv_pseudo_obs, &
    criterion_w, par_bounds, vine_edge_uv, qnorm, fv_fitcache_enable, fv_fitcache_clear, fv_fitcache_disable, fv_fitcache_stats, &
    fv_fitcache_capacity, &
    vine_rosenblatt_disc, bicop_cdf_x
  use iso_fortran_env, only: int64
  !$ use omp_lib
  implicit none
  private
  integer, parameter :: dp = kind(1.0d0)

  ! families
  integer, parameter, public :: FV_INDEP = 0, FV_GAUSSIAN = 1, FV_STUDENT = 2, FV_CLAYTON = 3, &
    FV_GUMBEL = 4, FV_FRANK = 5, FV_JOE = 6, FV_BB1 = 7, FV_BB6 = 8, FV_BB7 = 9, FV_BB8 = 10, FV_TAWN = 11, &
    FV_PLACKETT = 12, FV_FGM = 13, FV_AMH = 14, FV_GALAMBOS = 15, FV_HUSLER_REISS = 16, FV_BB2 = 17, &
    FV_BB3 = 18, FV_BB5 = 19, FV_BB9 = 20, FV_BB10 = 21, FV_BB4 = 22, FV_TEV = 23, FV_TLL = 24, &
    FV_BERNSTEIN = 25
  integer, parameter, public :: FV_NCLASSIC = 26   ! hand-coded families 0..25
  integer, parameter, public :: FV_NFAM = 192      ! capacity of fv_controls%family_set (generic families: 10000+)
  character(len=12), parameter, public :: FV_FAMILY_NAME(0:25) = [character(len=12) :: 'indep', 'gaussian', &
    'student', 'clayton', 'gumbel', 'frank', 'joe', 'bb1', 'bb6', 'bb7', 'bb8', 'tawn', 'plackett', 'fgm', &
    'amh', 'galambos', 'husler_reiss', 'bb2', 'bb3', 'bb5', 'bb9', 'bb10', 'bb4', 'tev', 'tll', 'bernstein']
  ! vinecopulib's parametric set, every parametric family, the nonparametric ones
  integer, parameter, public :: FV_PARAMETRIC(12) = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
  integer, parameter, public :: FV_ALL(24) = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, &
    18, 19, 20, 21, 22, 23]
  integer, parameter, public :: FV_NONPAR(2) = [24, 25]
  integer, parameter, public :: FV_ONE_PAR(5) = [1, 3, 4, 5, 6]
  ! selection criteria, methods, tree criteria / algorithms, tau estimation
  integer, parameter, public :: FV_LOGLIK = 0, FV_AIC = 1, FV_BIC = 2, FV_MBIC = 3, FV_MBICV = 4
  integer, parameter, public :: FV_MLE = 0, FV_ITAU = 1
  integer, parameter, public :: FV_TREE_TAU = 0, FV_TREE_RHO = 1, FV_TREE_JOE = 2, FV_TREE_HOEFFD = 3, &
    FV_TREE_XI = 4, FV_TREE_MI = 5, FV_TREE_FERREIRA = 6
  integer, parameter, public :: FV_PRIM = 0, FV_KRUSKAL = 1
  integer, parameter, public :: FV_TAU_AUTO = 0, FV_TAU_EXACT = 1, FV_TAU_APPROX = 2

  type, public :: fv_controls
    integer :: family_set(FV_NFAM) = [FV_ALL, FV_NONPAR, spread(0, 1, FV_NFAM - FV_NCLASSIC)]   ! first nfam used
    integer :: nfam = size(FV_ALL) - 1   ! default: every parametric family except the slow-to-fit tev;
                                         ! nfam = size(FV_ALL) adds it, nfam = FV_NCLASSIC adds tll and bernstein; fv_generic_families() the generic ones
    integer :: selcrit = FV_BIC
    integer :: method = FV_MLE
    integer :: tree_crit = FV_TREE_TAU
    integer :: tree_algo = FV_PRIM
    integer :: tau_mode = FV_TAU_EXACT
    integer :: trunc_lvl = -1           ! -1 = d - 1
    integer :: nthreads = 0             ! 0 = FORTVINE_THREADS or half the logical CPUs
    logical :: preselect = .true.
    logical :: allow_rotations = .true.
    logical :: fast_mle = .false.       ! .true.: one-parameter families by tau inversion only
    real(dp) :: threshold = 0.0_dp
    real(dp) :: psi0 = 0.9_dp
    logical :: select_trunc = .false.   ! truncation level chosen by mBICv (vinecopulib select_trunc_lvl)
    logical :: select_threshold = .false.  ! threshold chosen by mBICv (vinecopulib select_threshold)
  end type fv_controls

  type, public :: fv_bicop
    integer :: family = FV_INDEP
    integer :: rotation = 0
    real(dp) :: par(3) = 0.0_dp
    real(dp) :: loglik = 0.0_dp
    real(dp), allocatable :: aux(:)           ! nonparametric estimate (tll, bernstein)
  contains
    procedure :: pdf => bc_pdf
    procedure :: logpdf => bc_logpdf
    procedure :: hfunc1 => bc_hfunc1
    procedure :: hfunc2 => bc_hfunc2
    procedure :: hinv1 => bc_hinv1
    procedure :: hinv2 => bc_hinv2
    procedure :: cdf => bc_cdf
    procedure :: fit => bc_fit
    procedure :: select => bc_select
    procedure :: npars => bc_npars
    procedure :: npars_eff => bc_npars_eff
  end type fv_bicop

  type, public :: fv_vinecop
    integer :: d = 0
    integer :: trunc = 0
    integer, allocatable :: order(:)          ! structure order (1-based labels)
    integer, allocatable :: st(:, :)          ! natural-order struct array (edge, tree)
    integer, allocatable :: fam(:, :), rot(:, :)
    real(dp), allocatable :: par(:, :, :)     ! (3, edge, tree)
    real(dp), allocatable :: edge_ll(:, :)
    real(dp), allocatable :: aux(:, :, :)     ! (NAUX, edge, tree) when nonparametric families are used
    logical, allocatable :: disc(:)           ! discrete variables (original labels), set by select
    real(dp) :: threshold = 0.0_dp            ! threshold used by the last select (chosen if select_threshold)
    ! non-simplified (relaxed) edges, set by relax: parameters piecewise linear in
    ! w = sum_j wa(j, r) qnorm(u_j) (the edge's conditioning variables) between anchors wk(:, r)
    integer :: nrel = 0
    integer, allocatable :: erel(:, :)        ! (2, nrel): edge, tree
    real(dp), allocatable :: wa(:, :)         ! (d, nrel) loadings, original labels
    real(dp), allocatable :: wk(:, :)         ! (NK, nrel) anchor values of w (increasing)
    real(dp), allocatable :: pk(:, :, :)      ! (3, NK, nrel) pair-copula parameters at the anchors
  contains
    procedure :: init => vc_init
    procedure :: set_dvine => vc_set_dvine
    procedure :: set_cvine => vc_set_cvine
    procedure :: select => vc_select
    procedure :: pdf => vc_pdf
    procedure :: pdf_grid => vc_pdf_grid
    procedure :: loglik => vc_loglik
    procedure :: rosenblatt => vc_rosenblatt
    procedure :: inverse_rosenblatt => vc_inverse_rosenblatt
    procedure :: simulate => vc_simulate
    procedure :: npars => vc_npars
    procedure :: npars_eff => vc_npars_eff
    procedure :: matrix => vc_matrix
    procedure :: print => vc_print
    procedure :: relax => vc_relax
    procedure :: unrelax => vc_unrelax
    procedure :: rel_params => vc_rel_params
  end type fv_vinecop

  public :: fv_threads, fv_kendall_tau, fv_tau_mat, fv_pseudo_obs_nd, fv_dependence, NAUX, fv_par_bounds, &
    fam_is_np, gen_is, gen_valid, gen_dup, gen_npar, gen_flip, gen_exch, gen_lt, gen_ut
  public :: fv_fam_name, fv_generic_families
  ! pair-fit memo cache shared by every select/fit (exact reuse across tuning trials)
  public :: fv_fitcache_enable, fv_fitcache_clear, fv_fitcache_disable, fv_fitcache_stats, fv_fitcache_capacity

contains

  ! name of any family code: the hand-coded families 0..25, then the generic ones (gen_name)
  pure function fv_fam_name(fam) result(nm)
    integer, intent(in) :: fam
    character(len=48) :: nm
    if (fam >= 0 .and. fam < FV_NCLASSIC) then
      nm = FV_FAMILY_NAME(fam)
    else
      nm = gen_name(fam)
    end if
  end function fv_fam_name

  ! every valid generic family code; new_only (default .true.) leaves out the codes that
  ! duplicate a hand-coded family or another generic code
  pure function fv_generic_families(new_only) result(cs)
    logical, intent(in), optional :: new_only
    integer, allocatable :: cs(:)
    logical :: n1
    n1 = .true.
    if (present(new_only)) n1 = new_only
    cs = gen_codes(n1)
  end function fv_generic_families

  ! counter-based uniform in (0, 1): lowbias32 integer hash of (seed, i, j), 32-bit
  ! arithmetic kept inside int64 (no overflow); independent of thread count and order
  pure real(dp) function hash_unif(seed, i, j) result(r)
    integer, intent(in) :: seed, i, j
    integer(int64), parameter :: M32 = 4294967295_int64
    integer(int64) :: x
    x = iand(int(seed, int64) * 1664525_int64 + int(i, int64) * 22695477_int64 + int(j, int64) * 69069_int64 &
             + 1013904223_int64, M32)
    x = mix(x); x = mix(ieor(x, int(i, int64))); x = mix(ieor(x, int(j, int64)))
    r = (real(x, dp) + 0.5_dp) / 4294967296.0_dp
  contains
    pure integer(int64) function mix(y0) result(y)
      integer(int64), intent(in) :: y0
      y = y0
      y = ieor(y, shiftr(y, 16)); y = mul32(y, 2146121005_int64)
      y = ieor(y, shiftr(y, 15)); y = mul32(y, 2221713035_int64)
      y = ieor(y, shiftr(y, 16))
    end function mix
    pure integer(int64) function mul32(a, b)        ! (a b) mod 2^32 for a, b < 2^32
      integer(int64), intent(in) :: a, b
      mul32 = iand(a * iand(b, 65535_int64) + iand(a * shiftr(b, 16), 65535_int64) * 65536_int64, M32)
    end function mul32
  end function hash_unif

  ! parameter box of a pair-copula family (the MLE search box); indep / tll / bernstein: lb = ub = 0
  pure subroutine fv_par_bounds(fam, lb, ub)
    integer, intent(in) :: fam
    real(dp), intent(out) :: lb(3), ub(3)
    call par_bounds(fam, lb, ub)
  end subroutine fv_par_bounds

  ! ------------------------------------------------------------ utilities
  integer function fv_threads(req)
    integer, intent(in), optional :: req
    character(len=16) :: env
    integer :: st, v
    fv_threads = 1
    if (present(req)) then
      if (req > 0) then
        fv_threads = req; return
      end if
    end if
    call get_environment_variable("FORTVINE_THREADS", env, status=st)
    if (st == 0) then
      read(env, *, iostat=st) v
      if (st == 0 .and. v > 0) then
        fv_threads = v; return
      end if
    end if
    !$ fv_threads = max(1, omp_get_num_procs() / 2)
  end function fv_threads

  pure function ci6(c) result(ci)
    type(fv_controls), intent(in) :: c
    integer(c_int) :: ci(6)
    ci = [c%selcrit, c%method, merge(1, 0, c%preselect), merge(1, 0, c%allow_rotations), &
          merge(1, 0, c%fast_mle), c%tau_mode]
  end function ci6

  real(dp) function fv_kendall_tau(x, y)
    real(dp), intent(in) :: x(:), y(:)
    fv_kendall_tau = fv_kendall_tau_c(size(x), x, y)
  end function fv_kendall_tau

  ! |dependence| between x and y as used for the tree weights: kind = FV_TREE_* (tau, rho,
  ! joe, hoeffd, xi, mi, ferreira); optional weights (tau, rho, joe)
  real(dp) function fv_dependence(x, y, kind, weights)
    real(dp), intent(in) :: x(:), y(:)
    integer, intent(in) :: kind
    real(dp), intent(in), optional :: weights(:)
    real(dp) :: wd(1)
    if (present(weights)) then
      fv_dependence = criterion_w(size(x), x, y, kind, .true., weights)
    else
      fv_dependence = criterion_w(size(x), x, y, kind, .false., wd)
    end if
  end function fv_dependence

  ! Kendall tau matrix of u(n, d)
  subroutine fv_tau_mat(u, tm, nthreads)
    real(dp), intent(in) :: u(:, :)
    real(dp), intent(out) :: tm(:, :)
    integer, intent(in), optional :: nthreads
    real(dp), allocatable :: ut(:, :)
    ut = transpose(u)
    call fv_tau_matrix(size(u, 1), size(u, 2), ut, tm, fv_threads(nthreads))
  end subroutine fv_tau_mat

  ! pseudo-observations (average ranks / (n+1)) of x(n, d)
  subroutine fv_pseudo_obs_nd(x, u, nthreads)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(out) :: u(:, :)
    integer, intent(in), optional :: nthreads
    real(dp), allocatable :: xt(:, :), ut(:, :)
    xt = transpose(x)
    allocate(ut(size(x, 2), size(x, 1)))
    call fv_pseudo_obs(size(x, 1), size(x, 2), xt, ut, fv_threads(nthreads))
    u = transpose(ut)
  end subroutine fv_pseudo_obs_nd

  ! ------------------------------------------------------------ fv_bicop
  subroutine bc_eval(self, u1, u2, what, out)
    class(fv_bicop), intent(in) :: self
    real(dp), intent(in) :: u1(:), u2(:)
    integer, intent(in) :: what
    real(dp), intent(out) :: out(:)
    real(dp), allocatable :: u(:, :)
    real(dp) :: axd(1)
    allocate(u(2, size(u1)))
    u(1, :) = u1; u(2, :) = u2
    if (fam_is_np(self%family)) then
      call bicop_eval_x(self%family, self%rotation, self%par, self%aux, size(u1), u, what, out, fv_threads())
    else
      call bicop_eval_x(self%family, self%rotation, self%par, axd, size(u1), u, what, out, fv_threads())
    end if
  end subroutine bc_eval

  subroutine bc_pdf(self, u1, u2, out)
    class(fv_bicop), intent(in) :: self
    real(dp), intent(in) :: u1(:), u2(:)
    real(dp), intent(out) :: out(:)
    call bc_eval(self, u1, u2, 1, out)
  end subroutine bc_pdf
  subroutine bc_logpdf(self, u1, u2, out)
    class(fv_bicop), intent(in) :: self
    real(dp), intent(in) :: u1(:), u2(:)
    real(dp), intent(out) :: out(:)
    call bc_eval(self, u1, u2, 0, out)
  end subroutine bc_logpdf
  ! copula distribution function C(u1, u2) (closed form for every family except gaussian /
  ! student: Gauss-Legendre over the exact h-function in normal scores; tll: exact integral of the interpolant)
  subroutine bc_cdf(self, u1, u2, out)
    class(fv_bicop), intent(in) :: self
    real(dp), intent(in) :: u1(:), u2(:)
    real(dp), intent(out) :: out(:)
    real(dp) :: axd(1)
    if (allocated(self%aux)) then
      call bicop_cdf_x(self%family, self%rotation, self%par, self%aux, size(u1), u1, u2, out)
    else
      call bicop_cdf_x(self%family, self%rotation, self%par, axd, size(u1), u1, u2, out)
    end if
  end subroutine bc_cdf

  subroutine bc_hfunc1(self, u1, u2, out)
    class(fv_bicop), intent(in) :: self
    real(dp), intent(in) :: u1(:), u2(:)
    real(dp), intent(out) :: out(:)
    call bc_eval(self, u1, u2, 2, out)
  end subroutine bc_hfunc1
  subroutine bc_hfunc2(self, u1, u2, out)
    class(fv_bicop), intent(in) :: self
    real(dp), intent(in) :: u1(:), u2(:)
    real(dp), intent(out) :: out(:)
    call bc_eval(self, u1, u2, 3, out)
  end subroutine bc_hfunc2
  subroutine bc_hinv1(self, u1, u2, out)
    class(fv_bicop), intent(in) :: self
    real(dp), intent(in) :: u1(:), u2(:)
    real(dp), intent(out) :: out(:)
    call bc_eval(self, u1, u2, 4, out)
  end subroutine bc_hinv1
  subroutine bc_hinv2(self, u1, u2, out)
    class(fv_bicop), intent(in) :: self
    real(dp), intent(in) :: u1(:), u2(:)
    real(dp), intent(out) :: out(:)
    call bc_eval(self, u1, u2, 5, out)
  end subroutine bc_hinv2

  integer function bc_npars(self)
    class(fv_bicop), intent(in) :: self
    select case (self%family)
    case (FV_INDEP); bc_npars = 0
    case (FV_STUDENT, FV_BB1, FV_BB6, FV_BB7, FV_BB8, FV_BB2, FV_BB3, FV_BB5, FV_BB9, FV_BB10, FV_BB4, FV_TEV)
      bc_npars = 2
    case (FV_TAWN); bc_npars = 3
    case (FV_TLL, FV_BERNSTEIN); bc_npars = nint(self%npars_eff())
    case default
      bc_npars = 1
      if (gen_is(self%family)) bc_npars = gen_npar(self%family)
    end select
  end function bc_npars

  ! number of parameters; effective degrees of freedom for tll / bernstein
  real(dp) function bc_npars_eff(self)
    class(fv_bicop), intent(in) :: self
    if (fam_is_np(self%family)) then
      bc_npars_eff = 0.0_dp
      if (allocated(self%aux)) bc_npars_eff = self%aux(2)
    else
      bc_npars_eff = self%npars()
    end if
  end function bc_npars_eff

  ! fit the current family/rotation to (u1, u2)
  subroutine bc_fit(self, u1, u2, ctl, weights)
    class(fv_bicop), intent(inout) :: self
    real(dp), intent(in) :: u1(:), u2(:)
    type(fv_controls), intent(in), optional :: ctl
    real(dp), intent(in), optional :: weights(:)
    type(fv_controls) :: c
    real(dp), allocatable :: u(:, :)
    real(dp) :: ll, wd(1)
    if (present(ctl)) c = ctl
    allocate(u(2, size(u1)))
    u(1, :) = u1; u(2, :) = u2
    if (allocated(self%aux)) deallocate(self%aux)
    allocate(self%aux(merge(NAUX, 1, fam_is_np(self%family)))); self%aux = 0.0_dp
    if (present(weights)) then
      call bicop_fit_x(size(u1), u, self%family, self%rotation, ci6(c), .true., weights, self%par, ll, self%aux, &
                       fv_threads(c%nthreads))
    else
      call bicop_fit_x(size(u1), u, self%family, self%rotation, ci6(c), .false., wd, self%par, ll, self%aux, &
                       fv_threads(c%nthreads))
    end if
    self%loglik = ll
  end subroutine bc_fit

  ! select family, rotation and parameters (vinecopulib Bicop::select)
  subroutine bc_select(self, u1, u2, ctl, weights)
    class(fv_bicop), intent(inout) :: self
    real(dp), intent(in) :: u1(:), u2(:)
    type(fv_controls), intent(in), optional :: ctl
    real(dp), intent(in), optional :: weights(:)
    type(fv_controls) :: c
    real(dp), allocatable :: u(:, :)
    real(dp) :: ll, wd(1)
    integer :: f, r
    if (present(ctl)) c = ctl
    allocate(u(2, size(u1)))
    u(1, :) = u1; u(2, :) = u2
    if (allocated(self%aux)) deallocate(self%aux)
    allocate(self%aux(merge(NAUX, 1, any(fam_is_np(c%family_set(1:c%nfam)))))); self%aux = 0.0_dp
    if (present(weights)) then
      call bicop_select_x(size(u1), u, c%nfam, c%family_set(1:c%nfam), ci6(c), c%psi0, .true., weights, f, r, &
                          self%par, ll, self%aux, fv_threads(c%nthreads))
    else
      call bicop_select_x(size(u1), u, c%nfam, c%family_set(1:c%nfam), ci6(c), c%psi0, .false., wd, f, r, &
                          self%par, ll, self%aux, fv_threads(c%nthreads))
    end if
    self%family = f; self%rotation = r; self%loglik = ll
  end subroutine bc_select

  ! ------------------------------------------------------------ fv_vinecop
  subroutine vc_init(self, d)
    class(fv_vinecop), intent(inout) :: self
    integer, intent(in) :: d
    integer :: j
    self%d = d; self%trunc = 0
    call self%unrelax()
    if (allocated(self%order)) deallocate(self%order, self%st, self%fam, self%rot, self%par, self%edge_ll)
    if (allocated(self%aux)) deallocate(self%aux)
    if (allocated(self%disc)) deallocate(self%disc)
    allocate(self%order(d), self%st(d, d), self%fam(d, d), self%rot(d, d), self%par(3, d, d), self%edge_ll(d, d))
    allocate(self%aux(1, d, d), self%disc(d))
    self%order = [(j, j = 1, d)]
    self%st = 0; self%fam = 0; self%rot = 0; self%par = 0.0_dp; self%edge_ll = 0.0_dp; self%aux = 0.0_dp
    self%disc = .false.
  end subroutine vc_init

  ! fixed D-vine / C-vine structure for `order` (families chosen by select)
  subroutine vc_set_dvine(self, order, trunc)
    class(fv_vinecop), intent(inout) :: self
    integer, intent(in) :: order(:)
    integer, intent(in), optional :: trunc
    integer :: t, e, d
    d = size(order)
    call self%init(d)
    self%order = order
    self%trunc = d - 1
    if (present(trunc)) self%trunc = max(0, min(d - 1, trunc))
    do t = 1, self%trunc
      do e = 1, d - t
        self%st(e, t) = t + e
      end do
    end do
  end subroutine vc_set_dvine

  subroutine vc_set_cvine(self, order, trunc)
    class(fv_vinecop), intent(inout) :: self
    integer, intent(in) :: order(:)
    integer, intent(in), optional :: trunc
    integer :: t, d
    d = size(order)
    call self%init(d)
    self%order = order
    self%trunc = d - 1
    if (present(trunc)) self%trunc = max(0, min(d - 1, trunc))
    do t = 1, self%trunc
      self%st(1:d - t, t) = d - t + 1
    end do
  end subroutine vc_set_cvine

  ! Select structure (unless fixed) + families + parameters. u(n, d).
  ! fixed defaults to .true. when a structure is already set (set_dvine/set_cvine).
  ! weights(n): observation weights.  discrete(d) + u_minus(n, d): discrete variables
  ! (vinecopulib var_types "d"): u = F(x), u_minus = F(x-) for the discrete columns.
  ! select (vinecopulib Vinecop::select).  ctl%select_threshold: the threshold is chosen by
  ! mBICv as in vinecopulib's sparse_select_all_trees -- thresholds 1, then each time the
  ! ceil(5 %)-th largest criterion among the currently thresholded edges, until mBICv stops
  ! improving or the threshold drops below 0.01; the mBICv of a pass is accumulated over
  ! the trees up to the last check point (taken whenever more than 10 % of the pairs have
  ! been added since the previous one), as in vinecopulib.  Passes share the pair fits
  ! through the fit cache (enabled for the search if the caller has not enabled it).
  subroutine vc_select(self, u, ctl, fixed, info, weights, discrete, u_minus)
    class(fv_vinecop), intent(inout) :: self
    real(dp), intent(in) :: u(:, :)
    type(fv_controls), intent(in), optional :: ctl
    logical, intent(in), optional :: fixed
    integer, intent(out), optional :: info
    real(dp), intent(in), optional :: weights(:), u_minus(:, :)
    logical, intent(in), optional :: discrete(:)
    type(fv_controls) :: c
    type(fv_vinecop) :: cur, best
    real(dp), allocatable :: crit(:, :), tc(:)
    real(dp) :: thr, mb, mbt, mopt, chg, tot, neff, ps, ll, np, mt(size(u, 2))
    integer :: inf, d, t, e, m, nn, k
    logical :: owncache, prof, sel_tr
    character(len=8) :: ev
    integer :: es
    call get_environment_variable("FORTVINE_PROFILE", ev, status=es); prof = es == 0
    call self%unrelax()
    if (present(ctl)) c = ctl
    if (.not. c%select_threshold) then
      call vc_select_once(self, u, c, fixed, info, weights, discrete, u_minus)
      self%threshold = c%threshold
      return
    end if
    d = size(u, 2)
    if (self%d /= d) call self%init(d)
    owncache = fv_fitcache_capacity() == 0
    if (owncache) call fv_fitcache_enable(max(4096, 8 * d * d * c%nfam))
    neff = size(u, 1)
    if (present(weights)) neff = sum(weights)**2 / max(sum(weights**2), tiny(1.0_dp))
    c%select_threshold = .false.
    sel_tr = c%select_trunc
    c%select_trunc = .false.        ! passes fit every tree; truncation is decided below
    thr = 1.0_dp; mopt = huge(1.0_dp)
    allocate(crit(d, d))
    do
      select type (self)
      type is (fv_vinecop)
        cur = self
      end select
      c%threshold = thr
      call vc_select_once(cur, u, c, fixed, inf, weights, discrete, u_minus, crit)
      if (inf /= 0) exit
      tot = 0.5_dp * d * (d - 1); chg = 0.0_dp; mbt = 0.0_dp; mb = 0.0_dp
      do t = 1, d - 1
        ll = 0.0_dp; np = 0.0_dp; nn = 0
        if (t <= cur%trunc) then
          do e = 1, d - t
            ll = ll + cur%edge_ll(e, t)
            if (cur%fam(e, t) /= 0) nn = nn + 1
            if (fam_is_np(cur%fam(e, t))) then
              np = np + cur%aux(2, e, t)
            else if (cur%fam(e, t) /= 0) then
              np = np + npar_fam(cur%fam(e, t))
            end if
          end do
        end if
        ps = min(max(c%psi0**t, 1.0e-300_dp), 1.0_dp - epsilon(1.0_dp))
        ! mBICv of tree t (Nagler et al. 2019): prior psi0^t for a non-independence pair
        mt(t) = -2.0_dp * ll + log(neff) * np - 2.0_dp * (nn * log(ps) + (d - t - nn) * log(1.0_dp - ps))
        mbt = mbt + mt(t)
        chg = chg + (d - t)
        if (chg / tot > 0.1_dp) then
          chg = 0.0_dp
          if (sel_tr .and. t > 1 .and. mbt >= mb) then
            ! vinecopulib's truncation rollback: drop tree t, then earlier trees that do
            ! not lower the mBICv; the pass keeps the previous check point's mBICv
            mbt = mbt - mt(t); k = t - 1
            do while (k > 1)
              if (mt(k) <= 0.0_dp) exit
              mbt = mbt - mt(k); k = k - 1
            end do
            cur%trunc = min(cur%trunc, k)
            exit
          end if
          mb = mbt
        end if
      end do
#ifndef FORTVINE_NO_IO
      if (prof) print '(a,f12.8,a,f14.4,a,f14.4,a,i3)', '  [fortvine] threshold', thr, '  mbicv', mb, &
        '  loglik', cur%loglik(u, u_minus), '  trunc', cur%trunc
#endif
      if (mb >= mopt) exit
      best = cur; best%threshold = thr; mopt = mb
      if (thr < 0.01_dp) exit
      m = 0
      if (allocated(tc)) deallocate(tc)
      allocate(tc(d * d))
      do t = 1, cur%trunc
        do e = 1, d - t
          if (crit(e, t) < thr) then
            m = m + 1; tc(m) = crit(e, t)
          end if
        end do
      end do
      if (m == 0) exit
      call sort_desc(tc(1:m))
      k = max(1, ceiling(0.05_dp * m))
      thr = tc(k)
    end do
    if (owncache) call fv_fitcache_disable()
    if (present(info)) info = 0
    if (best%d == d) then
      select type (self)
      type is (fv_vinecop)
        self = best
      end select
    else
      call vc_select_once(self, u, c, fixed, info, weights, discrete, u_minus)
    end if
  contains
    integer function npar_fam(f)
      integer, intent(in) :: f
      type(fv_bicop) :: b
      b%family = f
      npar_fam = b%npars()
    end function npar_fam
    pure subroutine sort_desc(x)
      real(dp), intent(inout) :: x(:)
      integer :: i, j
      real(dp) :: v
      do i = 2, size(x)
        v = x(i); j = i - 1
        do while (j >= 1)
          if (x(j) >= v) exit
          x(j + 1) = x(j); j = j - 1
        end do
        x(j + 1) = v
      end do
    end subroutine sort_desc
  end subroutine vc_select

  subroutine vc_select_once(self, u, c, fixed, info, weights, discrete, u_minus, crit)
    class(fv_vinecop), intent(inout) :: self
    real(dp), intent(in) :: u(:, :)
    logical, intent(in), optional :: fixed
    integer, intent(out), optional :: info
    real(dp), intent(in), optional :: weights(:), u_minus(:, :)
    logical, intent(in), optional :: discrete(:)
    type(fv_controls), intent(in) :: c
    real(dp), intent(out), optional :: crit(:, :)
    real(dp), allocatable :: ut(:, :), umt(:, :), aux(:, :, :)
    real(dp) :: wd(1)
    integer :: ci(10), inf, st_out(size(u, 2), size(u, 2)), tout
    logical :: fx, dsc(size(u, 2))
    integer :: n, d, T
    n = size(u, 1); d = size(u, 2)
    if (self%d /= d) call self%init(d)
    fx = self%trunc > 0
    if (present(fixed)) fx = fixed
    T = d - 1
    if (c%trunc_lvl >= 0) T = min(d - 1, c%trunc_lvl)
    if (fx) T = min(T, self%trunc)
    ci(1:6) = ci6(c)
    ci(7:10) = [c%tree_crit, c%tree_algo, T, merge(1, 0, fx)]
    ut = transpose(u)
    dsc = .false.
    if (present(discrete)) dsc = discrete
    if (any(dsc)) then
      umt = transpose(u_minus)
    else
      allocate(umt(1, 1)); umt = 0.0_dp
    end if
    allocate(aux(merge(NAUX, 1, any(fam_is_np(c%family_set(1:c%nfam)))), d, d))
    if (present(weights)) then
      call vine_select_x(n, d, ut, c%nfam, c%family_set(1:c%nfam), ci, [c%threshold, c%psi0], self%order, &
                         self%st, st_out, self%fam, self%rot, self%par, self%edge_ll, fv_threads(c%nthreads), inf, &
                         .true., weights, dsc, umt, aux, c%select_trunc, tout, crit)
    else
      call vine_select_x(n, d, ut, c%nfam, c%family_set(1:c%nfam), ci, [c%threshold, c%psi0], self%order, &
                         self%st, st_out, self%fam, self%rot, self%par, self%edge_ll, fv_threads(c%nthreads), inf, &
                         .false., wd, dsc, umt, aux, c%select_trunc, tout, crit)
    end if
    if (present(info)) info = inf
    if (inf /= 0) return
    self%st = st_out
    self%trunc = tout
    call move_alloc(aux, self%aux)
    self%disc = dsc
  end subroutine vc_select_once

  ! density on a prediction grid: ug(d, bins, n) where only variable `tgt`
  ! (1-based) varies across bins; out(bins, n).  Edges not touching `tgt` are
  ! evaluated once per sample.
  ! ugm(d, bins, n): F(x-) grid for discrete variables (needed when the vine has any)
  subroutine vc_pdf_grid(self, ug, tgt, out, nthreads, ugm)
    class(fv_vinecop), intent(in) :: self
    real(dp), intent(in) :: ug(:, :, :)
    integer, intent(in) :: tgt
    real(dp), intent(out) :: out(:, :)
    integer, intent(in), optional :: nthreads
    real(dp), intent(in), optional :: ugm(:, :, :)
    if (any(self%disc)) then
      call vine_pdf_disc_x(size(ug, 3), size(ug, 2), self%d, tgt - 1, ug, ugm, self%disc, self%order, self%trunc, &
                           self%st, self%fam, self%rot, self%par, size(self%aux, 1), self%aux, out, fv_threads(nthreads))
    else
      ! (relaxed edges: erel is an absent argument while unallocated)
      call vine_pdf_x(size(ug, 3), size(ug, 2), self%d, tgt - 1, ug, self%order, self%trunc, self%st, self%fam, &
                      self%rot, self%par, size(self%aux, 1), self%aux, out, fv_threads(nthreads), self%erel, &
                      self%rel_params(ug(:, 1, :)))
    end if
  end subroutine vc_pdf_grid

  ! u_minus(n, d) = F(x-) for discrete variables (needed when the vine has any)
  subroutine vc_pdf(self, u, out, nthreads, u_minus)
    class(fv_vinecop), intent(in) :: self
    real(dp), intent(in) :: u(:, :)
    real(dp), intent(out) :: out(:)
    integer, intent(in), optional :: nthreads
    real(dp), intent(in), optional :: u_minus(:, :)
    real(dp), allocatable :: ut(:, :, :)
    if (any(self%disc)) then
      call vine_pdf_disc(size(u, 1), self%d, transpose(u), transpose(u_minus), self%disc, self%order, self%trunc, &
                         self%st, self%fam, self%rot, self%par, size(self%aux, 1), self%aux, out, fv_threads(nthreads))
      return
    end if
    allocate(ut(self%d, 1, size(u, 1)))
    ut(:, 1, :) = transpose(u)
    call vine_pdf_x(size(u, 1), 1, self%d, -1, ut, self%order, self%trunc, self%st, self%fam, self%rot, &
                    self%par, size(self%aux, 1), self%aux, out, fv_threads(nthreads), self%erel, &
                    self%rel_params(ut(:, 1, :)))
  end subroutine vc_pdf

  real(dp) function vc_loglik(self, u, u_minus)
    class(fv_vinecop), intent(in) :: self
    real(dp), intent(in) :: u(:, :)
    real(dp), intent(in), optional :: u_minus(:, :)
    real(dp), allocatable :: p(:)
    allocate(p(size(u, 1)))
    call self%pdf(u, p, u_minus=u_minus)
    vc_loglik = sum(log(p))
  end function vc_loglik

  ! Rosenblatt transform out(n, d).  With discrete variables (vinecopulib convention)
  ! u_minus(n, d) = F(x-) is required; out is Brockwell's randomised transform
  ! W F(x_k|...) + (1 - W) F(x_k-|...) with independent uniforms W (reproducible from
  ! seed, default 1), which is iid uniform under the model; randomize = .false. returns
  ! the conditional CDF F(x_k|...) itself, and out_minus(n, d) receives F(x_k-|...).
  subroutine vc_rosenblatt(self, u, out, nthreads, u_minus, out_minus, randomize, seed)
    class(fv_vinecop), intent(in) :: self
    real(dp), intent(in) :: u(:, :)
    real(dp), intent(out) :: out(:, :)
    integer, intent(in), optional :: nthreads, seed
    real(dp), intent(in), optional :: u_minus(:, :)
    real(dp), intent(out), optional :: out_minus(:, :)
    logical, intent(in), optional :: randomize
    real(dp), allocatable :: ut(:, :), ot(:, :), om(:, :)
    integer :: i, j, sd
    logical :: rnd
    if (any(self%disc)) then
#ifndef FORTVINE_NO_IO
      if (.not. present(u_minus)) error stop 'fv_vinecop%rosenblatt: u_minus is required for discrete variables'
#else
      if (.not. present(u_minus)) then      ! ponytail: no stop under R (CRAN); -1 is an impossible value
        out = -1.0_dp
        return
      end if
#endif
      allocate(ot(self%d, size(u, 1)), om(self%d, size(u, 1)))
      call vine_rosenblatt_disc(size(u, 1), self%d, transpose(u), transpose(u_minus), self%disc, self%order, &
                                self%trunc, self%st, self%fam, self%rot, self%par, size(self%aux, 1), self%aux, &
                                ot, om, fv_threads(nthreads))
      if (present(out_minus)) out_minus = transpose(om)
      rnd = .true.
      if (present(randomize)) rnd = randomize
      if (rnd) then
        sd = 1
        if (present(seed)) sd = seed
        do i = 1, size(u, 1)
          do j = 1, self%d
            if (self%disc(j)) ot(j, i) = om(j, i) + hash_unif(sd, i, j) * (ot(j, i) - om(j, i))
          end do
        end do
      end if
      out = transpose(ot)
      return
    end if
    ut = transpose(u)
    allocate(ot(self%d, size(u, 1)))
    call vine_rosenblatt_x(size(u, 1), self%d, ut, self%order, self%trunc, self%st, self%fam, self%rot, &
                           self%par, size(self%aux, 1), self%aux, ot, fv_threads(nthreads))
    out = transpose(ot)
  end subroutine vc_rosenblatt

  subroutine vc_inverse_rosenblatt(self, u, out, nthreads)
    class(fv_vinecop), intent(in) :: self
    real(dp), intent(in) :: u(:, :)
    real(dp), intent(out) :: out(:, :)
    integer, intent(in), optional :: nthreads
    real(dp), allocatable :: ut(:, :), ot(:, :)
    ut = transpose(u)
    allocate(ot(self%d, size(u, 1)))
    call vine_inverse_rosenblatt_x(size(u, 1), self%d, ut, self%order, self%trunc, self%st, self%fam, &
                                   self%rot, self%par, size(self%aux, 1), self%aux, ot, fv_threads(nthreads))
    out = transpose(ot)
  end subroutine vc_inverse_rosenblatt

  ! n draws from the vine, out(n, d); seed > 0 makes it reproducible
  subroutine vc_simulate(self, n, out, seed, nthreads)
    class(fv_vinecop), intent(in) :: self
    integer, intent(in) :: n
    real(dp), intent(out) :: out(:, :)
    integer, intent(in), optional :: seed, nthreads
    real(dp), allocatable :: w(:, :)
    integer :: ns, k
    integer, allocatable :: sd(:)
    if (present(seed)) then
      call random_seed(size=ns)
      allocate(sd(ns))
      sd = [(seed + 7919 * k, k = 1, ns)]
      call random_seed(put=sd)
    end if
    allocate(w(n, self%d))
    call random_number(w)
    call self%inverse_rosenblatt(w, out, nthreads)
  end subroutine vc_simulate

  integer function vc_npars(self)
    class(fv_vinecop), intent(in) :: self
    vc_npars = nint(self%npars_eff())
  end function vc_npars

  ! number of parameters, with effective degrees of freedom for tll / bernstein edges
  real(dp) function vc_npars_eff(self)
    class(fv_vinecop), intent(in) :: self
    type(fv_bicop) :: b
    integer :: t, e
    vc_npars_eff = 0.0_dp
    do t = 1, self%trunc
      do e = 1, self%d - t
        if (fam_is_np(self%fam(e, t))) then
          vc_npars_eff = vc_npars_eff + self%aux(2, e, t)
        else
          b%family = self%fam(e, t)
          vc_npars_eff = vc_npars_eff + b%npars()
        end if
      end do
    end do
    do e = 1, self%nrel                      ! relaxed edges: one parameter set per anchor
      b%family = self%fam(self%erel(1, e), self%erel(2, e))
      vc_npars_eff = vc_npars_eff + (size(self%wk, 1) - 1) * b%npars()
    end do
  end function vc_npars_eff

  ! vinecopulib R-vine matrix in ORIGINAL variable labels
  function vc_matrix(self) result(m)
    class(fv_vinecop), intent(in) :: self
    integer :: m(self%d, self%d), t, e, j
    m = 0
    do j = 1, self%d
      m(self%d - j + 1, j) = self%order(j)
    end do
    do t = 1, self%trunc
      do e = 1, self%d - t
        m(t, e) = self%order(self%st(e, t))
      end do
    end do
  end function vc_matrix

  subroutine vc_print(self, unit)
    class(fv_vinecop), intent(in) :: self
    integer, intent(in), optional :: unit
    integer :: u, t, e
    u = 6
    if (present(unit)) u = unit
#ifndef FORTVINE_NO_IO
    write(u, '(a,i0,a,i0)') 'fortvine vine copula: d = ', self%d, ', trunc = ', self%trunc
    write(u, '(a,*(1x,i0))') 'order:', self%order
    do t = 1, self%trunc
      do e = 1, self%d - t
        write(u, '(2x,a,i0,a,i0,a,a12,a,i3,a,3(1x,es11.4))') 'tree ', t, ' edge ', e, ': ', &
          trim(fv_fam_name(self%fam(e, t))), ' rot ', self%rot(e, t), ' par', self%par(:, e, t)
      end do
    end do
    do e = 1, self%nrel
      write(u, '(2x,a,i0,a,i0,a,*(1x,es10.3))') 'non-simplified tree ', self%erel(2, e), ' edge ', self%erel(1, e), &
        ': first parameter at the w anchors', self%pk(1, :, e)
    end do
#endif
  end subroutine vc_print


  subroutine vc_unrelax(self)
    class(fv_vinecop), intent(inout) :: self
    self%nrel = 0
    if (allocated(self%erel)) deallocate(self%erel, self%wa, self%wk, self%pk)
  end subroutine vc_unrelax

  ! parameters of the relaxed edges at the samples uc(d, n) (original order), prel(3, nrel, n)
  function vc_rel_params(self, uc) result(prel)
    class(fv_vinecop), intent(in) :: self
    real(dp), intent(in) :: uc(:, :)
    real(dp) :: prel(3, self%nrel, size(uc, 2))
    integer :: i, r, j, k, nk
    real(dp) :: w, a
    if (self%nrel == 0) return             ! simplified vine: prel is empty (wk unallocated)
    nk = size(self%wk, 1)
    do i = 1, size(uc, 2)
      do r = 1, self%nrel
        w = 0.0_dp
        do j = 1, self%d
          if (self%wa(j, r) /= 0.0_dp) w = w + self%wa(j, r) * qnorm(min(max(uc(j, i), 1e-10_dp), 1.0_dp - 1e-10_dp))
        end do
        if (w <= self%wk(1, r)) then
          prel(:, r, i) = self%pk(:, 1, r)
        else if (w >= self%wk(nk, r)) then
          prel(:, r, i) = self%pk(:, nk, r)
        else
          k = 1
          do while (self%wk(k + 1, r) < w)
            k = k + 1
          end do
          a = (w - self%wk(k, r)) / max(self%wk(k + 1, r) - self%wk(k, r), 1e-12_dp)
          prel(:, r, i) = (1.0_dp - a) * self%pk(:, k, r) + a * self%pk(:, k + 1, r)
        end if
      end do
    end do
  end function vc_rel_params

  ! Non-simplified pair copulas (local likelihood, Acar, Craiu & Yao 2011).  Every parametric edge of
  ! tree >= 2 whose conditioned pair contains variable tgt (0: any) and whose conditioning set D does
  ! not may get parameters that vary with w = the first principal component of the normal scores of
  ! D: kernel-weighted fits of its family at NK anchors in the rank of w, linear in w in between.
  ! mode 1 relaxes only the edges whose constant-conditional-correlation test rejects at 5 % (after
  ! Kurz & Spanhel 2022: Fisher z of the normal-score correlation over the quartiles of w, null by
  ! circular shifts of w, which keeps serial dependence); mode 2 relaxes them all.  Later trees are
  ! refitted (same families) on the inputs the relaxed edges produce.  u(n, d) in original order.
  subroutine vc_relax(self, u, tgt, mode, ctl, weights, nrel_out)
    class(fv_vinecop), intent(inout) :: self
    real(dp), intent(in) :: u(:, :)
    integer, intent(in) :: tgt, mode
    type(fv_controls), intent(in) :: ctl
    real(dp), intent(in), optional :: weights(:)
    integer, intent(out), optional :: nrel_out
    integer, parameter :: NK = 5, NBS = 99
    real(dp), parameter :: AK(NK) = [0.1_dp, 0.3_dp, 0.5_dp, 0.7_dp, 0.9_dp], HK = 0.15_dp
    integer :: n, d, t, e, k, nr, nth, q, i
    real(dp) :: lo(1)
    integer, allocatable :: er(:, :), DD(:), idx(:), rk(:)
    real(dp), allocatable :: ut(:, :), v1(:, :), v2(:, :), wa(:, :), wk(:, :), pk(:, :, :), a(:), w(:), kw(:), &
                             ow(:), prel(:, :, :)
    type(fv_bicop) :: b
    logical :: tg, changed
    call self%unrelax()
    if (present(nrel_out)) nrel_out = 0
    if (mode <= 0 .or. self%trunc < 2) return
    if (allocated(self%disc)) then
      if (any(self%disc)) return
    end if
    n = size(u, 1); d = self%d; nth = fv_threads(ctl%nthreads)
    ut = transpose(u)
    allocate(v1(d, n), v2(d, n), er(2, d * d), wa(d, d * d), wk(NK, d * d), pk(3, NK, d * d), kw(n), ow(n), w(n), &
             rk(n))
    ow = 1.0_dp
    if (present(weights)) ow = weights
    wa = 0.0_dp; nr = 0
    do t = 2, self%trunc
      changed = nr > 0
      prel = rel_at(nr)
      call vine_edge_uv(n, d, ut, self%order, self%trunc, self%st, self%fam, self%rot, self%par, size(self%aux, 1), &
                        self%aux, t, v1, v2, nth, er(:, 1:nr), prel)
      do e = 1, d - t
        DD = [(self%order(self%st(e, q)), q = 1, t - 1)]
        tg = (tgt == 0 .or. self%order(e) == tgt .or. self%order(self%st(e, t)) == tgt) .and. .not. any(DD == tgt) &
             .and. .not. fam_is_np(self%fam(e, t))
        if (tg) then
          call pc1(DD)
          if (mode /= 2) tg = ccc_p(v1(e, :), v2(e, :)) < 0.05_dp
        end if
        if (tg) then
          nr = nr + 1
          if (self%fam(e, t) == FV_INDEP) then
            self%fam(e, t) = FV_GAUSSIAN; self%rot(e, t) = 0
          end if
          er(:, nr) = [e, t]; wa(:, nr) = 0.0_dp; wa(DD, nr) = a
          b%family = self%fam(e, t); b%rotation = self%rot(e, t)
          do k = 1, NK
            wk(k, nr) = w(idx(max(1, min(n, nint(AK(k) * n)))))
            kw = ow * exp(-0.5_dp * ((real(rk, dp) / n - AK(k)) / HK)**2)
            call b%fit(v1(e, :), v2(e, :), ctl, kw)
            pk(:, k, nr) = b%par
          end do
          self%par(:, e, t) = pk(:, (NK + 1) / 2, nr)
          prel = rel_at(nr)                          ! its log-likelihood at the per-sample parameters
          self%edge_ll(e, t) = 0.0_dp
          do i = 1, n
            call bicop_eval_x(self%fam(e, t), self%rot(e, t), prel(:, nr, i), self%aux(:, e, t), 1, &
                              reshape([v1(e, i), v2(e, i)], [2, 1]), 0, lo, 1)
            self%edge_ll(e, t) = self%edge_ll(e, t) + ow(i) * lo(1)
          end do
        else if (changed .and. self%fam(e, t) /= FV_INDEP) then
          b%family = self%fam(e, t); b%rotation = self%rot(e, t)
          call b%fit(v1(e, :), v2(e, :), ctl, ow)
          self%par(:, e, t) = b%par; self%edge_ll(e, t) = b%loglik
          if (fam_is_np(self%fam(e, t))) self%aux(:, e, t) = b%aux(1:size(self%aux, 1))
        end if
      end do
    end do
    if (nr > 0) then
      self%nrel = nr
      self%erel = er(:, 1:nr); self%wa = wa(:, 1:nr); self%wk = wk(:, 1:nr); self%pk = pk(:, :, 1:nr)
    end if
    if (present(nrel_out)) nrel_out = nr
  contains
    ! parameters of the edges relaxed so far, at the training samples
    function rel_at(m) result(p)
      integer, intent(in) :: m
      real(dp), allocatable :: p(:, :, :)
      type(fv_vinecop) :: tmp
      allocate(p(3, m, n))
      if (m == 0) return
      tmp%d = d; tmp%nrel = m
      tmp%erel = er(:, 1:m); tmp%wa = wa(:, 1:m); tmp%wk = wk(:, 1:m); tmp%pk = pk(:, :, 1:m)
      p = tmp%rel_params(ut)
    end function rel_at

    ! w = first principal component of the normal scores of the columns cs (loadings a, on the raw
    ! normal scores: their mean is ~0 and sd ~1), its order idx and ranks rk
    subroutine pc1(cs)
      integer, intent(in) :: cs(:)
      real(dp), allocatable :: z(:, :), cm(:, :), an(:)
      integer :: j, l, it, nd
      nd = size(cs)
      allocate(z(n, nd))
      do j = 1, nd
        do l = 1, n
          z(l, j) = qnorm(min(max(u(l, cs(j)), 1e-10_dp), 1.0_dp - 1e-10_dp))
        end do
      end do
      cm = matmul(transpose(z), z) / n
      a = [(1.0_dp / sqrt(real(nd, dp)), j = 1, nd)]
      do it = 1, 200                                 ! power iteration
        an = matmul(cm, a)
        an = an / max(sqrt(sum(an**2)), 1e-300_dp)
        if (maxval(abs(an - a)) < 1e-10_dp) then
          a = an; exit
        end if
        a = an
      end do
      w = matmul(z, a)
      idx = argsort(w)
      do l = 1, n
        rk(idx(l)) = l
      end do
    end subroutine pc1

    ! p-value of the CCC test of the pair (x1, x2) against the quartiles of w (circular-shift null)
    real(dp) function ccc_p(x1, x2) result(pv)
      real(dp), intent(in) :: x1(:), x2(:)
      real(dp), allocatable :: z1(:), z2(:)
      integer, allocatable :: g(:)
      real(dp) :: t0
      integer :: j, cnt
      allocate(z1(n), z2(n), g(n))
      do j = 1, n
        z1(j) = qnorm(min(max(x1(j), 1e-10_dp), 1.0_dp - 1e-10_dp))
        z2(j) = qnorm(min(max(x2(j), 1e-10_dp), 1.0_dp - 1e-10_dp))
        g(j) = min(4, 1 + (4 * (rk(j) - 1)) / n)
      end do
      t0 = cstat(g, z1, z2); cnt = 0
      do j = 1, NBS
        if (cstat(cshift(g, (j * n) / (NBS + 1)), z1, z2) >= t0) cnt = cnt + 1
      end do
      pv = real(1 + cnt, dp) / real(NBS + 1, dp)
    end function ccc_p

    ! Fisher-z heterogeneity of the normal-score correlation over the groups gg (1..4)
    real(dp) function cstat(gg, z1, z2)
        integer, intent(in) :: gg(:)
        real(dp), intent(in) :: z1(:), z2(:)
        real(dp) :: sn(4), s1(4), s2(4), s11(4), s22(4), s12(4), r, zf(4), zb, ww
        integer :: l, h
        sn = 0; s1 = 0; s2 = 0; s11 = 0; s22 = 0; s12 = 0
        do l = 1, n
          h = gg(l)
          sn(h) = sn(h) + 1; s1(h) = s1(h) + z1(l); s2(h) = s2(h) + z2(l)
          s11(h) = s11(h) + z1(l)**2; s22(h) = s22(h) + z2(l)**2; s12(h) = s12(h) + z1(l) * z2(l)
        end do
        zb = 0; ww = 0
        do h = 1, 4
          r = (s12(h) - s1(h) * s2(h) / sn(h)) / &
              sqrt(max((s11(h) - s1(h)**2 / sn(h)) * (s22(h) - s2(h)**2 / sn(h)), 1e-300_dp))
          zf(h) = atanh(max(-0.999999_dp, min(0.999999_dp, r)))
          zb = zb + (sn(h) - 3) * zf(h); ww = ww + (sn(h) - 3)
        end do
        zb = zb / ww
        cstat = sum((sn - 3) * (zf - zb)**2)
    end function cstat
  end subroutine vc_relax

  ! ascending order of x (bottom-up merge sort, stable)
  function argsort(x) result(ix)
    real(dp), intent(in) :: x(:)
    integer, allocatable :: ix(:), tmp(:)
    integer :: n, w, lo, mid, hi, i, j, k
    n = size(x)
    ix = [(i, i = 1, n)]
    allocate(tmp(n))
    w = 1
    do while (w < n)
      do lo = 1, n, 2 * w
        mid = min(lo + w - 1, n); hi = min(lo + 2 * w - 1, n)
        i = lo; j = mid + 1; k = lo
        do while (i <= mid .and. j <= hi)
          if (x(ix(j)) < x(ix(i))) then
            tmp(k) = ix(j); j = j + 1
          else
            tmp(k) = ix(i); i = i + 1
          end if
          k = k + 1
        end do
        do while (i <= mid)
          tmp(k) = ix(i); i = i + 1; k = k + 1
        end do
        do while (j <= hi)
          tmp(k) = ix(j); j = j + 1; k = k + 1
        end do
      end do
      ix = tmp
      w = 2 * w
    end do
  end function argsort

end module fortvine_api
