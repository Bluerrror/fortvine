! families_core.inc -- log-density of the core pair-copula families 1, 3-11 (12-22 via
! families_ext.inc) written once in generic form and compiled from fortvine.f90 for each
! dual type VT: type(hd) (3 score slots), hd1 (1), hd2 (2).
  pure function LPDF_NAME(fam, P, a, b, x, y) result(lp)
    integer, intent(in) :: fam
    VT, intent(in) :: P(3)
    real(dp), intent(in) :: a, b, x, y
    VT :: lp, r, s, t1, t2, ls, st, eu, ev, ed, den, ua, ub, pa, pb, da, db, dum, cc, dc, ddc
    VT :: aa, a1, a2, tm, lx2, ly2, hi2
    VT :: gq1, gq2, gq3
    VT :: ta, tb, lta, ltb, e1, e2, w, ra, rb, lra, lrb, z, lz, lq, eta, l1da, l1db, ga, gb, l1ga, l1gb, emz
    real(dp) :: lx, ly, la, lb, m, tt, lab
    select case (fam)
    case (1)
      r = P(1)
      s = 1.0_dp - r * r
      lp = -0.5_dp * flog(s) - (r * r * (x * x + y * y) - 2.0_dp * x * y * r) / (2.0_dp * s)
    case (3)
      la = log(a); lb = log(b)
      t1 = -(P(1) * la); t2 = -(P(1) * lb)
      if (fval(t1) >= fval(t2)) then
        ls = t1 + flog(1.0_dp + fexp(t2 - t1) - fexp(-t1))
      else
        ls = t2 + flog(1.0_dp + fexp(t1 - t2) - fexp(-t2))
      end if
      lp = flog1p(P(1)) - (1.0_dp + P(1)) * (la + lb) - (2.0_dp + 1.0_dp / P(1)) * ls
    case (4)
      lx = log(x); ly = log(y)
      t1 = P(1) * lx; t2 = P(1) * ly
      if (fval(t1) >= fval(t2)) then
        ls = t1 + flog1p(fexp(t2 - t1))
      else
        ls = t2 + flog1p(fexp(t1 - t2))
      end if
      st = fexp(ls / P(1))
      lp = -st + (P(1) - 1.0_dp) * (lx + ly) + (x + y) + (1.0_dp / P(1) - 2.0_dp) * ls &
           + flog(P(1) - 1.0_dp + st)
    case (5)
      eu = fexpm1(-(P(1) * a)); ev = fexpm1(-(P(1) * b)); ed = fexpm1(-P(1))
      den = ed + eu * ev
      lp = flog(-(P(1) * ed) * (eu + 1.0_dp) * (ev + 1.0_dp)) - 2.0_dp * flog(fabs(den))
    case (6)
      la = log(1.0_dp - a); lb = log(1.0_dp - b)
      ua = fexp(P(1) * la); ub = fexp(P(1) * lb)
      s = ua + ub - ua * ub
      lp = (1.0_dp / P(1) - 2.0_dp) * flog(s) + (P(1) - 1.0_dp) * (la + lb) + flog(P(1) - 1.0_dp + s)
    case (7)   ! BB1   (x = log a, y = log b)
      ta = fexpm1(-(P(1) * x)); tb = fexpm1(-(P(1) * y))
      lta = flog(ta); ltb = flog(tb)
      e1 = P(2) * lta; e2 = P(2) * ltb
      ls = flse(e1, e2)
      w = fexp(ls / P(2))
      lp = flog(P(1) * (P(2) - 1.0_dp) + (P(1) * P(2) + 1.0_dp) * w) &
           + (-1.0_dp / P(1) - 2.0_dp) * flog1p(w) + (1.0_dp / P(2) - 2.0_dp) * ls &
           - (P(1) + 1.0_dp) * (x + y) + (P(2) - 1.0_dp) * (lta + ltb)
    case (8)   ! BB6   (x = log(1-a), y = log(1-b))
      ua = fexp(P(1) * x); ub = fexp(P(1) * y)
      ra = -flog1p(-ua); rb = -flog1p(-ub)
      lra = flog(ra); lrb = flog(rb)
      ls = flse(P(2) * lra, P(2) * lrb)
      w = fexp(ls / P(2))
      z = -fexpm1(-w); lz = flog(z)
      lp = flog(P(1)) + (1.0_dp / P(1) - 2.0_dp) * lz - w + flog(w) - 2.0_dp * ls &
           + flog(z * (w + P(2) - 1.0_dp) + (1.0_dp - 1.0_dp / P(1)) * fexp(-w) * w) &
           + (P(2) - 1.0_dp) * (lra + lrb) + (P(1) - 1.0_dp) * (x + y) + ra + rb
    case (9)   ! BB7   (x = log(1-a), y = log(1-b))
      ua = fexp(P(1) * x); ub = fexp(P(1) * y)
      ra = flog1p(-ua); rb = flog1p(-ub)
      ls = flog1p(fexpm1(-(P(2) * ra)) + fexpm1(-(P(2) * rb)))
      lq = -(ls / P(2))
      z = -fexpm1(lq); lz = flog(z)
      lp = flog(P(1)) + (1.0_dp / P(1) - 2.0_dp) * lz + lq - 2.0_dp * ls &
           + flog((1.0_dp - 1.0_dp / P(1)) * fexp(lq) + z * (1.0_dp + P(2))) &
           - (P(2) + 1.0_dp) * (ra + rb) + (P(1) - 1.0_dp) * (x + y)
    case (10)  ! BB8   (x = a, y = b)
      eta = -fexpm1(P(1) * flog1p(-P(2)))
      l1da = flog1p(-(P(2) * a)); l1db = flog1p(-(P(2) * b))
      ga = fexp(P(1) * l1da); gb = fexp(P(1) * l1db)
      l1ga = flog1p(-ga); l1gb = flog1p(-gb)
      emz = fexp(l1ga + l1gb) / eta
      lz = flog(1.0_dp - emz)
      lp = flog(P(2) * P(1)) + (1.0_dp / P(1) - 2.0_dp) * lz + flog(emz) + flog(1.0_dp - emz / P(1)) &
           + (P(1) - 1.0_dp) * (l1da + l1db) - l1ga - l1gb
    case (11)
      lab = x + y
      tt = y / lab
      tm = (P(2) * tt)**P(3) + (P(1) * (1.0_dp - tt))**P(3)
      aa = (1.0_dp - P(1)) * (1.0_dp - tt) + (1.0_dp - P(2)) * tt + tm**(1.0_dp / P(3))
      a1 = P(1) - P(2) + tm**(1.0_dp / P(3) - 1.0_dp) * (P(2) * (P(2) * tt)**(P(3) - 1.0_dp) &
           - P(1) * (P(1) * (1.0_dp - tt))**(P(3) - 1.0_dp))
      lx2 = P(3) * flog(P(2) * tt); ly2 = P(3) * flog(P(1) * (1.0_dp - tt))
      if (fval(lx2) >= fval(ly2)) then
        hi2 = lx2 + flog1p(fexp(ly2 - lx2))
      else
        hi2 = ly2 + flog1p(fexp(lx2 - ly2))
      end if
      a2 = fexp(flog(P(3) - 1.0_dp) + (1.0_dp / P(3) - 2.0_dp) * hi2 + lx2 + ly2 &
                - 2.0_dp * log(tt) - 2.0_dp * log1p_(-tt))
      lp = lab * aa - lab + flog(aa * aa + (1.0_dp - 2.0_dp * tt) * a1 * aa &
           - (1.0_dp - tt) * tt * (a1 * a1 + a2 / lab))
    case (12:22)
      call EXT_NAME(fam, P, a, b, x, y, .false., lp, lx, ly)   ! x, y = log a, log b (prep_xy)
    case default
      if (gen_is(fam)) then                                  ! generic families (families_gen.inc)
        call GEN_NAME(fam, P, a, b, lp, gq1, gq2, gq3)
      else
        lp = VCONST(0.0_dp)
      end if
    end select
  end function LPDF_NAME
