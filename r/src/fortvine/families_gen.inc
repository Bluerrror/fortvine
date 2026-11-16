! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Shahab Aldin Shojaeezadeh. Part of fortvine (see LICENSE). If you use it, please cite:
! Shojaeezadeh, S. A. (2026). fortvine: a complete compendium to work on vine copula modeling.
! families_gen.inc -- the generic pair-copula engine, written once and compiled from
! fortvine.f90 for real(dp) (REALMODE) and for every dual type VT (exact parameter scores).
!
! One formula covers Archimedean, extreme-value, Archimax and mixture-of-max-id copulas
! (Joe 2014, sections 3.5 and 3.18; Caperaa, Fougeres & Genest 2000):
!
!     C(u, v) = Psi( l(Phi(u), Phi(v)) ),     Phi = Psi^-1,
!
!   Psi  outer generator: none (Psi = exp(-s): C is the inner copula itself), a one- or
!        two-parameter Laplace-transform family, or a composition psi(-log zeta(s))
!        ("psi stopped zeta", Joe 2014 appendix A.1.2);
!   l    log-scale inner max-id copula: l = s + t (independence: C is Archimedean) or
!        l(s, t) = (s + t) A(t / (s + t)) with a Pickands function A (C is Archimax; with
!        Psi = exp(-s) it is extreme-value).
!
! With S = Phi(u), T = Phi(v), L = l(S, T):
!     h1 = dC/du = Psi'(L) l_s / Psi'(S),   h2 = Psi'(L) l_t / Psi'(T),
!     c  = [Psi''(L) l_s l_t + (-Psi'(L)) (-l_st)] / (Psi'(S) Psi'(T)).
! Psi'' >= 0, Psi' < 0, l_s, l_t >= 0 and l_st <= 0, so both density terms are >= 0 and
! everything is evaluated in log space: each generator returns log Psi, log(-Psi'),
! log Psi''; each inner returns l, log l_s, log l_t and -l_st (linear, may be 0).  Every
! argument s travels with log s, because Phi(u) of a gamma-inner composition (BB2, BB3, ...)
! is exp(delta x) and overflows long before log Phi(u) does.
! A Khoudraji asymmetrisation C_K(u, v) = u^(1-a) v^(1-b) C(u^a, v^b) wraps any member.
! Family codes, parameter layout and bounds: gen_* in fortvine.f90.
  pure subroutine GEN_NAME(fam, P, u, v, lp, lh1, lh2, lc)
    integer, intent(in) :: fam
    VT, intent(in) :: P(3)
    real(dp), intent(in) :: u, v          ! rotated and clamped to (0, 1)
    VT, intent(out) :: lp, lh1, lh2, lc  ! log c, log h1, log h2, log C
    VT :: Q(2), R(3), ka, kb, x, y, lx, ly, lp0, lh10, lh20, lc0, br
    integer :: og, ie, no, ni
    real(dp) :: lu, lv
    og = gen_og(fam); ie = gen_ie(fam); no = gen_npo(og); ni = gen_npi(ie)
    Q = VCONST(0.0_dp); R = VCONST(0.0_dp)
    if (no > 0) Q(1:no) = P(1:no)
    if (ni > 0) R(1:ni) = P(no + 1:no + ni)
    lu = log(u); lv = log(v)
    if (.not. gen_kh(fam)) then
      call amx(VCONST(u), VCONST(lu), VCONST(log1p_(-u)), VCONST(v), VCONST(lv), VCONST(log1p_(-v)), lp, lh1, lh2, lc)
      return
    end if
    ! Khoudraji: c_K = u^-a v^-b C0 [(1-a)(1-b) + (1-a) b y d2 + (1-b) a x d1 + a b x y c0 / C0]
    ! with x = u^a, y = v^b and d1 = h1_0 / C0, d2 = h2_0 / C0 (all four terms >= 0)
    ka = P(no + ni + 1); kb = P(no + ni + 2)
    lx = ka * lu; ly = kb * lv; x = fexp(lx); y = fexp(ly)
    call amx(x, lx, l1me(lx), y, ly, l1me(ly), lp0, lh10, lh20, lc0)
    br = (1.0_dp - ka) * (1.0_dp - kb) + (1.0_dp - ka) * kb * fexp(ly + lh20 - lc0) &
         + (1.0_dp - kb) * ka * fexp(lx + lh10 - lc0) + ka * kb * fexp(lx + ly + lp0 - lc0)
    lp = -(ka * lu) - kb * lv + lc0 + flog(br)
    lh1 = -(ka * lu) + (1.0_dp - kb) * lv + lc0 + flog((1.0_dp - ka) + ka * fexp(lx + lh10 - lc0))
    lh2 = (1.0_dp - ka) * lu - kb * lv + lc0 + flog((1.0_dp - kb) + kb * fexp(ly + lh20 - lc0))
    lc = (1.0_dp - ka) * lu + (1.0_dp - kb) * lv + lc0
  contains

    ! log(1 - exp(x)), x < 0, accurate at both ends (Maechler 2012)
    pure function l1me(x) result(r)
      VT, intent(in) :: x
      VT :: r
      if (fval(x) > -0.6931471805599453_dp) then
        r = flog(-fexpm1(x))
      else
        r = flog1p(-fexp(x))
      end if
    end function l1me

    ! log Phi(z); below z = -30 the asymptotic (Mills ratio) series, where Phi underflows near -38
    pure function lpn(z) result(r)
      VT, intent(in) :: z
      VT :: r, q
      if (fval(z) > -30.0_dp) then
        r = flog(fpnorm(z))
      else
        q = 1.0_dp / (z * z)
        r = -(0.5_dp * z * z) - 0.9189385332046728_dp - flog(-z) &
            + flog1p(q * (-1.0_dp + q * (3.0_dp + q * (-15.0_dp + q * (105.0_dp - 945.0_dp * q)))))
      end if
    end function lpn

    ! log(1 - (1-e)^(1/d)) with e = e^-s, w = 1-e (lw = log w): for tiny e it is
    ! log(e/d) + log1p((1-1/d) e/2) + O(e^2), where log(1 - w^(1/d)) would underflow to log 0
    pure function sib_l1(d, s, e, lw) result(r)
      VT, intent(in) :: d, s, e, lw
      VT :: r
      if (fval(e) < 1.0e-8_dp) then
        r = -s - flog(d) + flog1p(0.5_dp * (1.0_dp - 1.0_dp / d) * e)
      else
        r = l1me(lw / d)
      end if
    end function sib_l1

    ! logarithmic series pieces at s: z = (1-e^-d) e^-s, log z, log(1-z) and log L, L = -log(1-z);
    ! 1 - z = (1 - e^-s) + e^-(d+s) has no cancellation near z = 1, and log L = log z + log1p(z/2)
    ! for tiny z
    pure subroutine lsr(d, s, lz, z, l1z, lll)
      VT, intent(in) :: d, s
      VT, intent(out) :: lz, z, l1z, lll
      lz = l1me(-d) - s; z = fexp(lz)
      if (fval(z) > 0.5_dp) then
        l1z = flse(l1me(-s), -d - s)
      else
        l1z = flog1p(-z)
      end if
      if (fval(z) < 1.0e-8_dp) then
        lll = lz + flog1p(0.5_dp * z)
      else
        lll = flog(-l1z)
      end if
    end subroutine lsr

    ! log(1 + s) from s and log s (s may have overflowed)
    pure function lop(s, ls) result(r)
      VT, intent(in) :: s, ls
      VT :: r
      if (fval(ls) > 0.0_dp) then
        r = ls + flog1p(fexp(-ls))
      else
        r = flog1p(s)
      end if
    end function lop

    ! s = expm1(t) and log s for t > 0, without overflow
    pure subroutine em1(t, s, ls)
      VT, intent(in) :: t
      VT, intent(out) :: s, ls
      s = fexpm1(t)
      if (fval(t) > 30.0_dp) then
        ls = t + l1me(-t)
      else
        ls = flog(s)
      end if
    end subroutine em1

    ! the master formula at (x, y) with log x, log(1-x), ... supplied
    pure subroutine amx(x, lx, l1x, y, ly, l1y, lp, lh1, lh2, lc)
      VT, intent(in) :: x, lx, l1x, y, ly, l1y
      VT, intent(out) :: lp, lh1, lh2, lc
      VT :: S, T, lS, lT, L, lL, lls, llt, lmst, d0, d2, sS, sT, pL0, pL1, pL2
      call gphi(x, lx, l1x, S, lS); call gphi(y, ly, l1y, T, lT)
      call gpsi(S, lS, d0, sS, d2, .false.); call gpsi(T, lT, d0, sT, d2, .false.)   ! log(-Psi') only
      if (ie == 0) then
        lL = flse(lS, lT); L = fexp(lL)
        lls = VCONST(0.0_dp); llt = VCONST(0.0_dp)
        call gpsi(L, lL, pL0, pL1, pL2, .true.)
        lp = pL2 + lls + llt - sS - sT
      else
        call inner(lS, lT, L, lL, lls, llt, lmst)  ! lmst = log(-l_st)
        call gpsi(L, lL, pL0, pL1, pL2, .true.)
        lp = flse(pL2 + lls + llt, pL1 + lmst) - sS - sT   ! both density terms in log space
      end if
      lh1 = pL1 + lls - sS
      lh2 = pL1 + llt - sT
      lc = pL0
    end subroutine amx

    ! ------------------------------------------------------------------ outer generator
    ! log Psi, log(-Psi'), log Psi'' at s > 0 (ls = log s)
    pure subroutine gpsi(s, ls, l0, l1, l2, full)
      VT, intent(in) :: s, ls
      VT, intent(out) :: l0, l1, l2
      logical, intent(in) :: full           ! .false.: l1 only matters (l2 is left at the outer part)
      VT :: g, lgg, lg1, mg2, a0, a1, a2
      integer :: po, pin
      select case (og)
      case (0)                                     ! Psi = exp(-s)
        l0 = -s; l1 = -s; l2 = -s
      case (1:7)
        call prim(og, Q(1), Q(2), s, ls, l0, l1, l2)
      case default                                 ! composition psi(g(s)), g = -log zeta
        po = (og - 11) / 4 + 1; pin = mod(og - 11, 4) + 1
        call bern(pin, Q(2), s, ls, g, lgg, lg1, mg2, full)
        call prim(po, Q(1), Q(2), g, lgg, a0, a1, a2)
        l0 = a0; l1 = a1 + lg1; l2 = a2
        if (full) l2 = a2 + 2.0_dp * lg1 + flog1p(fexp(a1 - a2 - 2.0_dp * lg1) * mg2)
      end select
    end subroutine gpsi

    ! s = Phi(x) = Psi^-1(x) and log s, x in (0, 1)
    pure subroutine gphi(x, lx, l1x, s, ls)
      VT, intent(in) :: x, lx, l1x
      VT, intent(out) :: s, ls
      VT :: xs, lxs
      integer :: po, pin
      select case (og)
      case (0)
        s = -lx; ls = flog(s)
      case (1:7)
        call primphi(og, Q(1), Q(2), x, lx, l1x, s, ls)
      case default
        po = (og - 11) / 4 + 1; pin = mod(og - 11, 4) + 1
        call primphi(po, Q(1), Q(2), x, lx, l1x, xs, lxs)
        call berninv(pin, Q(2), xs, lxs, s, ls)
      end select
    end subroutine gphi

    ! Laplace-transform families (Joe 2014 appendix A): 1 positive stable, 2 gamma, 3 Sibuya,
    ! 4 logarithmic series (one parameter th); 5 exponentially tilted positive stable (th, de),
    ! 6 generalised Sibuya (th, de), 7 shifted negative binomial (th, pi = de)
    pure subroutine prim(k, th, de, s, ls, l0, l1, l2)
      integer, intent(in) :: k
      VT, intent(in) :: th, de, s, ls
      VT, intent(out) :: l0, l1, l2
      VT :: lq, q, l1s, lw, lk, lz, z, l1z, r, lr, e, ld
      select case (k)
      case (1)                     ! exp(-s^(1/th))
        lq = ls / th; q = fexp(lq)
        l0 = -q; l1 = lq - flog(th) - ls - q
        l2 = -q + lq - flog(th) - 2.0_dp * ls + flog(q / th + 1.0_dp - 1.0_dp / th)
      case (2)                     ! (1+s)^(-1/th)
        l1s = lop(s, ls)
        l0 = -(l1s / th); l1 = -flog(th) - (1.0_dp / th + 1.0_dp) * l1s
        l2 = -flog(th) + flog(1.0_dp / th + 1.0_dp) - (1.0_dp / th + 2.0_dp) * l1s
      case (3)                     ! 1 - (1 - e^-s)^(1/th)
        e = fexp(-s); lw = l1me(-s)
        l0 = sib_l1(th, s, e, lw); l1 = -flog(th) + (1.0_dp / th - 1.0_dp) * lw - s
        l2 = -flog(th) + (1.0_dp / th - 2.0_dp) * lw - s + flog1p(-(e / th))
      case (4)                     ! -log(1 - (1 - e^-th) e^-s) / th
        call lsr(th, s, lz, z, l1z, lk)
        l0 = lk - flog(th); l1 = lz - l1z - flog(th); l2 = lz - 2.0_dp * l1z - flog(th)
      case (5)                     ! exp(1/de - (de^-th + s)^(1/th))
        r = fexp(-(th * flog(de))) + s; lr = flog(r); lq = lr / th; q = fexp(lq)
        l0 = 1.0_dp / de - q; l1 = l0 + lq - flog(th) - lr
        l2 = l0 + lq - flog(th) - 2.0_dp * lr + flog(q / th + 1.0_dp - 1.0_dp / th)
      case (6)                     ! [1 - (1 - eta e^-s)^(1/th)] / de,  eta = 1 - (1-de)^th
        lz = l1me(th * flog1p(-de)) - s; z = fexp(lz); l1z = flog1p(-z)
        l0 = l1me(l1z / th) - flog(de)
        l1 = -flog(de) - flog(th) + (1.0_dp / th - 1.0_dp) * l1z + lz
        l2 = -flog(de * th) + (1.0_dp / th - 2.0_dp) * l1z + lz + flog1p(-(z / th))
      case default                 ! [(1-pi) / (e^s - pi)]^(1/th)
        e = fexp(-s); ld = flog1p(-(de * e))
        l0 = (flog1p(-de) - s - ld) / th
        l1 = l0 - flog(th) - ld
        l2 = l0 - 2.0_dp * ld + flog(1.0_dp / (th * th) + de * e / th)
      end select
    end subroutine prim

    pure subroutine primphi(k, th, de, x, lx, l1x, s, ls)
      integer, intent(in) :: k
      VT, intent(in) :: th, de, x, lx, l1x
      VT, intent(out) :: s, ls
      VT :: r
      select case (k)
      case (1); ls = th * flog(-lx); s = fexp(ls)
      case (2); call em1(-(th * lx), s, ls)
      case (3); s = -l1me(th * l1x); ls = flog(s)
      case (4)                     ! -log[(1 - e^-th x)/(1 - e^-th)], = -log1p(-r) near x = 1
        r = fexp(-th - l1me(-th)) * fexpm1(th * fexp(l1x))
        if (fval(r) < 0.5_dp) then
          s = -flog1p(-r)
        else
          s = l1me(-th) - l1me(-(th * x))
        end if
        ls = flog(s)
      case (5); s = fexp(-(th * flog(de))) * fexpm1(th * flog1p(-(de * lx))); ls = flog(s)
      case (6); s = l1me(th * flog1p(-de)) - l1me(th * flog1p(-(de * x))); ls = flog(s)
      case default; s = flog1p((1.0_dp - de) * fexpm1(-(th * lx))); ls = flog(s)
      end select
    end subroutine primphi

    ! inner Bernstein function g = -log zeta of an infinitely divisible LT zeta(.; d):
    ! g, log g, log g', -g'' (>= 0, linear: it can be exactly 0)
    pure subroutine bern(k, d, s, ls, g, lgg, lg1, mg2, full)
      integer, intent(in) :: k
      VT, intent(in) :: d, s, ls
      VT, intent(out) :: g, lgg, lg1, mg2
      logical, intent(in) :: full           ! .false.: -g'' not needed
      VT :: l1s, e, lw, lzt, bb, dl, ek, dk, lk, lz, z, l1z, ll, lmz, zk, tau
      integer :: j
      select case (k)
      case (1)                     ! positive stable: s^(1/d)
        lgg = ls / d; g = fexp(lgg)
        lg1 = -flog(d) + (1.0_dp / d - 1.0_dp) * ls
        mg2 = (1.0_dp / d) * (1.0_dp - 1.0_dp / d) * fexp((1.0_dp / d - 2.0_dp) * ls)
      case (2)                     ! gamma: log(1+s)/d
        l1s = lop(s, ls); g = l1s / d; lgg = flog(l1s) - flog(d)
        lg1 = -flog(d) - l1s; mg2 = fexp(-flog(d) - 2.0_dp * l1s)
      case (3)                     ! Sibuya: -log(1 - (1-e^-s)^(1/d))
        e = fexp(-s); lw = l1me(-s); lzt = sib_l1(d, s, e, lw); g = -lzt; lgg = flog(g)
        lg1 = -flog(d) + (1.0_dp / d - 1.0_dp) * lw - s - lzt
        ! -g'' = d^-1 w^(1/d-2) e B / (1-w^(1/d))^2 with B = 1 - e/d - w^(1/d) >= 0;
        ! for small e, B = (1-e/d)(1 - exp(D)) with D summed so the O(e) terms cancel exactly
        mg2 = VCONST(0.0_dp)
        if (.not. full) return
        if (fval(e) < 1.0e-3_dp) then         ! direct form below: relative error < 1e-10 for e >= 1e-3
          dl = VCONST(0.0_dp); ek = e; dk = 1.0_dp / d
          do j = 2, 40
            ek = ek * e; dk = dk / d
            dl = dl - (ek / real(j, dp)) * (1.0_dp / d - dk)
            if (fval(ek) < 1.0e-17_dp) exit
          end do
          bb = fexp(flog1p(-(e / d))) * (-fexpm1(dl))
        else
          bb = 1.0_dp - e / d - fexp(lw / d)
        end if
        if (fval(e) < 1.0e-8_dp) then          ! -g'' = (1 - 1/d) e / 2 + O(e^2): exp(s) * 0 otherwise
          mg2 = 0.5_dp * (1.0_dp - 1.0_dp / d) * e
        else
          mg2 = fexp(-flog(d) + (1.0_dp / d - 2.0_dp) * lw - s - 2.0_dp * lzt) * bb
        end if
      case default                 ! logarithmic series: -log(-log(1 - k e^-s)/d), k = 1-e^-d
        call lsr(d, s, lz, z, l1z, lk); ll = -l1z         ! lk = log L here
        ! near s = 0, L = d - log1p(tau) with tau = expm1(d) (1 - e^-s): g = -log1p(-log1p(tau)/d)
        ! without the cancellation of log d - log L (g ~ g'(0) s there)
        tau = flog1p(fexpm1(d) * (-fexpm1(-s)))
        if (fval(tau) < 0.5_dp * fval(d)) then
          lk = flog(d - tau); g = -flog1p(-(tau / d))
        else
          g = flog(d) - lk
        end if
        lgg = flog(g)
        lg1 = lz - l1z - lk
        mg2 = VCONST(0.0_dp)
        if (.not. full) return
        if (fval(z) < 1.0e-3_dp) then            ! L - z = sum_{j>=2} z^j / j
          lmz = VCONST(0.0_dp); zk = z
          do j = 2, 60
            zk = zk * z; lmz = lmz + zk / real(j, dp)
            if (fval(zk) < 1.0e-17_dp * fval(lmz)) exit
          end do
        else
          lmz = ll - z
        end if
        if (fval(z) < 1.0e-8_dp) then          ! -g'' = z / 2 + O(z^2)
          mg2 = 0.5_dp * z
        else
          mg2 = fexp(lz - 2.0_dp * l1z - 2.0_dp * lk) * lmz
        end if
      end select
    end subroutine bern

    pure subroutine berninv(k, d, xg, lxg, s, ls)
      integer, intent(in) :: k
      VT, intent(in) :: d, xg, lxg
      VT, intent(out) :: s, ls
      VT :: e, a
      select case (k)
      case (1); ls = d * lxg; s = fexp(ls)
      case (2); call em1(d * xg, s, ls)
      case (3)                     ! -log(1 - (1 - e^-x)^d) = x - log d + (d-1) e^-x / 2 + ...
        e = fexp(-xg)
        if (fval(e) < 1.0e-8_dp) then
          s = xg - flog(d) + 0.5_dp * (d - 1.0_dp) * e
        else
          s = -l1me(d * l1me(-xg))
        end if
        ls = flog(s)
      case default                 ! log(1-e^-d) - log(1 - exp(-d e^-x))
        a = flog(d) - xg
        if (fval(a) < -20.0_dp) then
          s = l1me(-d) - (a - 0.5_dp * fexp(a))
        else
          s = l1me(-d) - l1me(-fexp(a))
        end if
        ls = flog(s)
      end select
    end subroutine berninv

    ! ------------------------------------------------------------------ inner (Pickands)
    ! l = (S+T) A(w), w = T/(S+T), from log S and log T: log w and log(1-w) are exact even
    ! where w itself underflows (S and T of very different size: gamma-inner compositions)
    pure subroutine inner(lS, lT, L, lL, lls, llt, lmst)
      VT, intent(in) :: lS, lT
      VT, intent(out) :: L, lL, lls, llt, lmst
      VT :: lst, lw, l1w, w, aa, la2
      lst = flse(lS, lT); lw = lT - lst; l1w = lS - lst; w = fexp(lw)
      call pick(ie, R, w, lw, l1w, aa, lls, llt, la2)
      lL = lst + flog(aa); L = fexp(lL)
      lmst = lw + l1w + la2 - lst                     ! log(-l_st) = log(w (1-w) A2 / (S+T))
    end subroutine inner

    ! Pickands function at w (w = share of the second margin, as fortvine's Tawn / Galambos):
    ! A, log l_s = log(A - w A1), log l_t = log(A + (1-w) A1) and log A2 (A1, A2 = first and
    ! second derivative).  The partials are written per family without the cancellation of
    ! A - w A1 (it is ~0 at strong dependence)
    pure recursive subroutine pick(k, Rp, w, lw, l1w, aa, lsd, ltd, la2)
      integer, intent(in) :: k
      VT, intent(in) :: Rp(3), w, lw, l1w
      VT, intent(out) :: aa, lsd, ltd, la2
      VT :: z1, z2, ph, one, lv, vv, lzz, l1zz, zz, b0, bls, blt, lb2
      one = VCONST(1.0_dp)
      select case (k)
      case (1); call asl(one, one, Rp(1), w, lw, l1w, aa, lsd, ltd, la2)          ! logistic (Gumbel)
      case (2); call anl(one, one, Rp(1), w, lw, l1w, aa, lsd, ltd, la2)          ! negative logistic
      case (3)                                                                    ! Husler-Reiss
        z1 = 1.0_dp / Rp(1) + 0.5_dp * Rp(1) * (lw - l1w)
        z2 = 1.0_dp / Rp(1) - 0.5_dp * Rp(1) * (lw - l1w)
        aa = w * fpnorm(z1) + (1.0_dp - w) * fpnorm(z2)
        lsd = lpn(z2); ltd = lpn(z1)                                              ! l_s = Phi(z2), l_t = Phi(z1)
        la2 = flog(Rp(1)) - 0.5_dp * z1 * z1 - 0.9189385332046728_dp - 0.6931471805599453_dp - lw - 2.0_dp * l1w
      case (4)                                                                    ! Tawn symmetric mixed
        aa = 1.0_dp - Rp(1) * fexp(lw + l1w)
        lsd = flog1p(-(Rp(1) * fexp(2.0_dp * lw))); ltd = flog1p(-(Rp(1) * fexp(2.0_dp * l1w)))
        la2 = flog(2.0_dp * Rp(1))
      case (5); call asl(Rp(1), Rp(2), Rp(3), w, lw, l1w, aa, lsd, ltd, la2)      ! asymmetric logistic
      case (6); call asl(Rp(2), one, Rp(1), w, lw, l1w, aa, lsd, ltd, la2)        ! Tawn type 1 (th, psi1)
      case (7); call asl(one, Rp(2), Rp(1), w, lw, l1w, aa, lsd, ltd, la2)        ! Tawn type 2 (th, psi2)
      case (8); call anl(Rp(1), Rp(2), Rp(3), w, lw, l1w, aa, lsd, ltd, la2)      ! asym. negative logistic
      case (9); call anl(Rp(2), one, Rp(1), w, lw, l1w, aa, lsd, ltd, la2)        ! type 1 (th, psi1)
      case (10); call anl(one, Rp(2), Rp(1), w, lw, l1w, aa, lsd, ltd, la2)       ! type 2 (th, psi2)
      case (11)                                                                   ! asymmetric mixed
        ph = -(Rp(1) / 3.0_dp) + Rp(2) * ((1.0_dp - Rp(1)) / 2.0_dp + Rp(1) / 3.0_dp)
        aa = 1.0_dp - (Rp(1) + ph) * w + Rp(1) * w * w + ph * w * w * w
        lsd = flog(1.0_dp - Rp(1) * w * w - 2.0_dp * ph * w * w * w)
        ltd = flog(1.0_dp - Rp(1) - ph + 2.0_dp * Rp(1) * w + (3.0_dp * ph - Rp(1)) * w * w - 2.0_dp * ph * w * w * w)
        la2 = flog(2.0_dp * Rp(1) + 6.0_dp * ph * w)
      case default                                                                ! 12/13: Khoudraji of HR / mixed
        ! l_K(x, y) = (1-a) x + (1-b) y + l0(a x, b y): A_K(w) = (1-a)(1-w) + (1-b) w + v A0(z),
        ! v = a(1-w) + b w, z = b w / v; l_K,s = (1-a) + a l0_s, l_K,t = (1-b) + b l0_t;
        ! A_K2 = a^2 b^2 A0_2(z) / v^3
        lv = flse(flog(Rp(2)) + l1w, flog(Rp(3)) + lw); vv = fexp(lv)
        lzz = flog(Rp(3)) + lw - lv; l1zz = flog(Rp(2)) + l1w - lv; zz = fexp(lzz)
        call pick(merge(3, 4, k == 12), Rp, zz, lzz, l1zz, b0, bls, blt, lb2)
        aa = (1.0_dp - Rp(2)) * (1.0_dp - w) + (1.0_dp - Rp(3)) * w + vv * b0
        lsd = flog((1.0_dp - Rp(2)) + Rp(2) * fexp(bls)); ltd = flog((1.0_dp - Rp(3)) + Rp(3) * fexp(blt))
        la2 = 2.0_dp * (flog(Rp(2)) + flog(Rp(3))) + lb2 - 3.0_dp * lv
      end select
    end subroutine pick

    ! asymmetric logistic: A = (1-p1)(1-w) + (1-p2)w + M, M = [(p1(1-w))^th + (p2 w)^th]^(1/th);
    ! with sp + sq = 1 the two shares of the bracket: l_s = (1-p1) + M sp/(1-w),
    ! l_t = (1-p2) + M sq/w, A2 = (th-1) M sp sq (p1 p2 / (p q))^2
    pure subroutine asl(p1, p2, th, w, lw, l1w, aa, lsd, ltd, la2)
      VT, intent(in) :: p1, p2, th, w, lw, l1w
      VT, intent(out) :: aa, lsd, ltd, la2
      VT :: lp_, lq_, lpp, lqq, lm, lmm
      lp_ = flog(p1) + l1w; lq_ = flog(p2) + lw
      lpp = th * lp_; lqq = th * lq_; lm = flse(lpp, lqq); lmm = lm / th
      aa = (1.0_dp - p1) * (1.0_dp - w) + (1.0_dp - p2) * w + fexp(lmm)
      lsd = flog((1.0_dp - p1) + fexp(lmm + lpp - lm - l1w))
      ltd = flog((1.0_dp - p2) + fexp(lmm + lqq - lm - lw))
      if (fval(th) > 1.0_dp) then
        la2 = flog(th - 1.0_dp) + lmm + (lpp - lm) + (lqq - lm) + 2.0_dp * (flog(p1) + flog(p2) - lp_ - lq_)
      else
        la2 = VCONST(-745.0_dp)                     ! th = 1: independence part, A2 = 0
      end if
    end subroutine asl

    ! asymmetric negative logistic: A = 1 - M, M = [(p1(1-w))^-th + (p2 w)^-th]^(-1/th); with
    ! R = (p/q)^th: l_s = 1 - p1 (1+R)^(-1-1/th), l_t = 1 - p2 (1+1/R)^(-1-1/th),
    ! A2 = (1+th) M rp rq (p1 p2 / (p q))^2
    pure subroutine anl(p1, p2, th, w, lw, l1w, aa, lsd, ltd, la2)
      VT, intent(in) :: p1, p2, th, w, lw, l1w
      VT, intent(out) :: aa, lsd, ltd, la2
      VT :: lp_, lq_, lpp, lqq, lm, lmm, lr
      lp_ = flog(p1) + l1w; lq_ = flog(p2) + lw
      lpp = -(th * lp_); lqq = -(th * lq_); lm = flse(lpp, lqq); lmm = -(lm / th)
      aa = 1.0_dp - fexp(lmm)
      lr = th * (lp_ - lq_)
      lsd = flog((1.0_dp - p1) - p1 * fexpm1(-(1.0_dp + 1.0_dp / th) * flse(VCONST(0.0_dp), lr)))
      ltd = flog((1.0_dp - p2) - p2 * fexpm1(-(1.0_dp + 1.0_dp / th) * flse(VCONST(0.0_dp), -lr)))
      la2 = flog(1.0_dp + th) + lmm + (lpp - lm) + (lqq - lm) + 2.0_dp * (flog(p1) + flog(p2) - lp_ - lq_)
    end subroutine anl
  end subroutine GEN_NAME
