! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Shahab Aldin Shojaeezadeh. Part of fortvine (see LICENSE). If you use it, please cite:
! Shojaeezadeh, S. A. (2026). fortvine: a complete compendium to work on vine copula modeling.
! families_ext.inc -- extended pair-copula families, written ONCE in generic
! form and compiled twice from fortvine.f90 (cpp):
!   VT = real(dp)  + REALMODE : log-density and both h-functions (evaluation)
!   VT = type(hd)             : log-density with exact parameter scores (MLE)
! Inputs (a, b) are rotated and clamped to (0,1); P(1:3) are the parameters.
!
! 12 plackett  (th > 0)            13 fgm  (th in [-1,1])     14 amh (th in [-1,1))
! 15 galambos  (th > 0, EV)        16 husler_reiss (lam > 0, EV)
! 17 bb2  (th > 0, de > 0)         18 bb3  (th >= 1, de > 0)  19 bb5 (th >= 1, de > 0, EV)
! 20 bb9  (th >= 1, de > 0)        21 bb10 (th > 0, pi in [0,1])
! 22 bb4  (th > 0, de > 0)
!
! Archimedean families use log c = log psi2(s) + log|phi1(a)| + log|phi1(b)| and
! h1 = psi1(s) phi1(a) with s = phi(a) + phi(b) (C itself is never formed);
! extreme-value families use the Pickands function A(t) with t = log b / log(ab).
  pure subroutine EXT_NAME(fam, P, a, b, la0, lb0, nh, lp, h1, h2)
    integer, intent(in) :: fam
    logical, intent(in) :: nh          ! compute h1, h2 (REALMODE only)
    VT, intent(in) :: P(3)
    real(dp), intent(in) :: a, b, la0, lb0   ! la0, lb0 = log a, log b (parameter-free, precomputed)
    VT, intent(out) :: lp
    real(dp), intent(out) :: h1, h2
    VT :: th, de, s, d, q, r, l1, za, zb, lpa, lpb, lps1, lps2, aa, a1, a2, g, g1, g2, bb, b1, b2
    VT :: m, rp, rq, lpt, lqt, lm, mt, z1, z2, ua, ub, brk, sq
    real(dp) :: la, lb, lab, t, xa, xb
    h1 = b; h2 = a
    th = P(1); de = P(2)
    select case (fam)

    case (12)   ! Plackett
      s = 1.0_dp + (th - 1.0_dp) * (a + b)
      d = s * s - 4.0_dp * th * (th - 1.0_dp) * a * b
      lp = flog(th) + flog(1.0_dp + (th - 1.0_dp) * (a + b - 2.0_dp * a * b)) - 1.5_dp * flog(d)
#ifdef REALMODE
      if (nh) then
        sq = sqrt(d)
        h1 = 0.5_dp * (1.0_dp - (s - 2.0_dp * th * b) / sq)
        h2 = 0.5_dp * (1.0_dp - (s - 2.0_dp * th * a) / sq)
      end if
#endif

    case (13)   ! Farlie-Gumbel-Morgenstern
      lp = flog1p(th * (1.0_dp - 2.0_dp * a) * (1.0_dp - 2.0_dp * b))
#ifdef REALMODE
      if (nh) then
        h1 = b * (1.0_dp + th * (1.0_dp - b) * (1.0_dp - 2.0_dp * a))
        h2 = a * (1.0_dp + th * (1.0_dp - a) * (1.0_dp - 2.0_dp * b))
      end if
#endif

    case (14)   ! Ali-Mikhail-Haq
      q = 1.0_dp - th * (1.0_dp - a) * (1.0_dp - b)
      lp = flog(1.0_dp + th * ((1.0_dp + a) * (1.0_dp + b) - 3.0_dp) + th * th * (1.0_dp - a) * (1.0_dp - b)) &
           - 3.0_dp * flog(q)
#ifdef REALMODE
      if (nh) then
        h1 = b * (1.0_dp - th * (1.0_dp - b)) / (q * q)
        h2 = a * (1.0_dp - th * (1.0_dp - a)) / (q * q)
      end if
#endif

    case (15, 16, 19)   ! extreme-value families
      la = la0; lb = lb0; lab = la + lb
      t = lb / lab
      if (fam == 15) then
        call galambos_A(th, t, aa, a1, a2)
      else if (fam == 16) then
        ! Husler-Reiss: z1,2 = 1/lam +- (lam/2) log(t/(1-t)); A' = Phi(z1) - Phi(z2) because
        ! t phi(z1) = (1-t) phi(z2); A'' = lam phi(z1) / (2 t (1-t)^2)
        z1 = 1.0_dp / th + 0.5_dp * th * log(t / (1.0_dp - t))
        z2 = 1.0_dp / th - 0.5_dp * th * log(t / (1.0_dp - t))
        aa = t * fpnorm(z1) + (1.0_dp - t) * fpnorm(z2)
        a1 = fpnorm(z1) - fpnorm(z2)
        a2 = th * fdnorm(z1) / (2.0_dp * t * (1.0_dp - t) * (1.0_dp - t))
      else
        ! BB5: A = (t^th + (1-t)^th - g^th)^(1/th) with g = (t^-ga + (1-t)^-ga)^(-1/ga),
        ! ga = th*de.  Every g-term is G = g^th = exp(-log m / de) times a bounded factor
        ! (K1 = g'/g, K2 = g''/g from the Galambos ratios), so nothing hits 0^negative.
        lpt = -((th * de) * log(t)); lqt = -((th * de) * log(1.0_dp - t))
        lm = flse(lpt, lqt)
        rp = fexp(lpt - lm); rq = fexp(lqt - lm)
        g = fexp(-(lm / de))
        g1 = rp / t - rq / (1.0_dp - t)
        g2 = -(1.0_dp + th * de) * rp * rq / (t * (1.0_dp - t))**2
        bb = t**th + (1.0_dp - t)**th - g
        b1 = th * (t**(th - 1.0_dp) - (1.0_dp - t)**(th - 1.0_dp) - g * g1)
        b2 = th * (th - 1.0_dp) * (t**(th - 2.0_dp) + (1.0_dp - t)**(th - 2.0_dp) - g * g1 * g1) - th * g * g2
        aa = bb**(1.0_dp / th)
        a1 = aa * b1 / (th * bb)
        a2 = aa * ((1.0_dp / th - 1.0_dp) * b1 * b1 / (th * bb * bb) + b2 / (th * bb))
      end if
      lp = lab * aa - lab + flog(aa * aa + (1.0_dp - 2.0_dp * t) * a1 * aa - (1.0_dp - t) * t * (a1 * a1 + a2 / lab))
#ifdef REALMODE
      if (nh) then
        h1 = exp(lab * aa) * (aa - t * a1) / a
        h2 = exp(lab * aa) * (aa + (1.0_dp - t) * a1) / b
      end if
#endif

    case (17)   ! BB2: phi = exp(de (t^-th - 1)) - 1, psi = (1 + log1p(s)/de)^(-1/th)
      xa = la0; xb = lb0
      za = de * fexpm1(-(th * xa)); zb = de * fexpm1(-(th * xb))
      lpa = flog(de * th) - (th + 1.0_dp) * xa + za
      lpb = flog(de * th) - (th + 1.0_dp) * xb + zb
      l1 = log1p_es(za, zb)
      r = 1.0_dp + l1 / de
      lps2 = -flog(th * de) + (-1.0_dp / th - 2.0_dp) * flog(r) - 2.0_dp * l1 + flog((1.0_dp + 1.0_dp / th) / de + r)
      lps1 = -flog(th * de) + (-1.0_dp / th - 1.0_dp) * flog(r) - l1
      lp = lps2 + lpa + lpb
#ifdef REALMODE
      if (nh) then
        h1 = exp(lps1 + lpa); h2 = exp(lps1 + lpb)
      end if
#endif

    case (18)   ! BB3: phi = exp(de (-log t)^th) - 1, psi = exp(-(log1p(s)/de)^(1/th))
      xa = -la0; xb = -lb0
      za = de * fexp(th * log(xa)); zb = de * fexp(th * log(xb))
      lpa = za + flog(de * th) + (th - 1.0_dp) * log(xa) + xa
      lpb = zb + flog(de * th) + (th - 1.0_dp) * log(xb) + xb
      l1 = log1p_es(za, zb)
      q = fexp(flog(l1 / de) / th)
      lps2 = -q + flog(q) - 2.0_dp * flog(th) - 2.0_dp * flog(l1) - 2.0_dp * l1 &
             + flog(q - 1.0_dp + th * (1.0_dp + l1))
      lps1 = -q + flog(q) - flog(th) - flog(l1) - l1
      lp = lps2 + lpa + lpb
#ifdef REALMODE
      if (nh) then
        h1 = exp(lps1 + lpa); h2 = exp(lps1 + lpb)
      end if
#endif

    case (20)   ! BB9 (Crowder): phi = (1/de - log t)^th - de^-th, psi = exp(-(s + de^-th)^(1/th) + 1/de)
      xa = -la0; xb = -lb0
      lpa = flog(th) + (th - 1.0_dp) * flog(1.0_dp / de + xa) + xa
      lpb = flog(th) + (th - 1.0_dp) * flog(1.0_dp / de + xb) + xb
      m = (1.0_dp / de + xa)**th + (1.0_dp / de + xb)**th - (1.0_dp / de)**th
      q = m**(1.0_dp / th)
      lps2 = -q + 1.0_dp / de + flog(q) + flog(q - 1.0_dp + th) - 2.0_dp * flog(th) - 2.0_dp * flog(m)
      lps1 = -q + 1.0_dp / de + flog(q) - flog(th) - flog(m)
      lp = lps2 + lpa + lpb
#ifdef REALMODE
      if (nh) then
        h1 = exp(lps1 + lpa); h2 = exp(lps1 + lpb)
      end if
#endif

    case (21)   ! BB10: C = uv [1 - pi (1-u^th)(1-v^th)]^(-1/th);  de = pi
      ua = fexp(th * la0); ub = fexp(th * lb0)
      d = 1.0_dp - de * (1.0_dp - ua) * (1.0_dp - ub)
      brk = (1.0_dp - de * (1.0_dp - ub) + de * th * ub) * d &
            - (1.0_dp + th) * de * (1.0_dp - de * (1.0_dp - ub)) * ub * (1.0_dp - ua)
      lp = (-1.0_dp / th - 2.0_dp) * flog(d) + flog(brk)
#ifdef REALMODE
      if (nh) then
        h1 = b * (1.0_dp - de * (1.0_dp - ub)) * exp((-1.0_dp / th - 1.0_dp) * log(d))
        h2 = a * (1.0_dp - de * (1.0_dp - ua)) * exp((-1.0_dp / th - 1.0_dp) * log(d))
      end if
#endif

    case (22)   ! BB4 (Joe 1997): x = a^-th - 1, y = b^-th - 1, S = (x^-de + y^-de)^(-1/de),
      ! W = 1 + x + y - S, C = W^(-1/th); p = (S/x)^(de+1), q = (S/y)^(de+1) (both < 1):
      ! h1 = W^(-1/th-1) (1-p) a^(-th-1),
      ! c = (ab)^(-th-1) W^(-1/th-2) [(1+th)(1-p)(1-q) + th (1+de) W p q / S]
      la = la0; lb = lb0
      za = fexpm1(-(th * la)); zb = fexpm1(-(th * lb))
      lpa = flog(za); lpb = flog(zb)
      l1 = -(flse(-(de * lpa), -(de * lpb)) / de)           ! log S
      rp = fexp((de + 1.0_dp) * (l1 - lpa)); rq = fexp((de + 1.0_dp) * (l1 - lpb))
      ua = -fexpm1((de + 1.0_dp) * (l1 - lpa)); ub = -fexpm1((de + 1.0_dp) * (l1 - lpb))   ! 1-p, 1-q
      g = 1.0_dp + za + zb - fexp(l1)
      lp = (-1.0_dp / th - 2.0_dp) * flog(g) - (th + 1.0_dp) * (la + lb) &
           + flog((1.0_dp + th) * ua * ub + th * (1.0_dp + de) * g * rp * rq * fexp(-l1))
#ifdef REALMODE
      if (nh) then
        h1 = exp((-1.0_dp / th - 1.0_dp) * log(g) + log(ua) - (th + 1.0_dp) * la)
        h2 = exp((-1.0_dp / th - 1.0_dp) * log(g) + log(ub) - (th + 1.0_dp) * lb)
      end if
#endif

    case default
      lp = 0.0_dp * th
    end select
#ifdef REALMODE
      if (nh) then
      h1 = min(max(h1, 0.0_dp), 1.0_dp)
      h2 = min(max(h2, 0.0_dp), 1.0_dp)
      end if
#endif
  contains
    ! log(1 + s) with s = expm1(za) + expm1(zb), za, zb >= 0, without overflow
    ! or cancellation
    pure function log1p_es(za, zb) result(l)
      VT, intent(in) :: za, zb
      VT :: l
      if (fval(za) < 1.0_dp .and. fval(zb) < 1.0_dp) then
        l = flog1p(fexpm1(za) + fexpm1(zb))
      else if (fval(za) >= fval(zb)) then
        l = za + flog(1.0_dp + fexp(zb - za) - fexp(-za))
      else
        l = zb + flog(1.0_dp + fexp(za - zb) - fexp(-zb))
      end if
    end function log1p_es

    ! Galambos Pickands function in log space: A = 1 - m^(-1/th), m = t^-th + (1-t)^-th,
    ! A' = m^(-1/th) (rq/(1-t) - rp/t), A'' = (1+th) m^(-1/th) rp rq / (t(1-t))^2,
    ! rp = t^-th/m, rq = (1-t)^-th/m  (rp + rq = 1, so no cancellation anywhere)
    pure subroutine galambos_A(th, t, aa, a1, a2)
      VT, intent(in) :: th
      real(dp), intent(in) :: t
      VT, intent(out) :: aa, a1, a2
      VT :: lpt, lqt, lm, mt, rp, rq
      lpt = -(th * log(t)); lqt = -(th * log(1.0_dp - t))
      lm = flse(lpt, lqt)
      mt = fexp(-(lm / th))
      rp = fexp(lpt - lm); rq = fexp(lqt - lm)
      aa = 1.0_dp - mt
      a1 = mt * (rq / (1.0_dp - t) - rp / t)
      a2 = (1.0_dp + th) * mt * rp * rq / (t * (1.0_dp - t))**2
    end subroutine galambos_A
  end subroutine EXT_NAME
