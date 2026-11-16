! fortvine_ad.f90 -- dual numbers: value and exact gradient with respect to up to 3
! copula parameters (forward-mode differentiation).
!
! Every operation applies the analytic chain rule, so per-observation scores
! d log c / d theta are exact to rounding -- no finite differences.  The MLE in
! fortvine.f90 uses them for Fisher scoring (BHHH): curvature = sum of score outer
! products, which avoids carrying a Hessian through every operation.
module fortvine_ad
  implicit none
  private
  integer, parameter :: dp = kind(1.0d0)

  type, public :: hd
    real(dp) :: v = 0.0_dp
    real(dp) :: g(3) = 0.0_dp
  end type hd

  public :: operator(+), operator(-), operator(*), operator(/), operator(**)
  public :: hexp, hlog, hlog1p, hexpm1, hsqrt, hconst, hvar, hunary, habs
  ! generic math over real(dp) and dual numbers: family formulas are written once
  ! (src/families_ext.inc) and compiled for both value evaluation and exact scores
  public :: fval, flog, fexp, flog1p, fexpm1, fsqrt, fabs, fpnorm, fdnorm, flse

  interface fval
    module procedure hval, rval
  end interface
  interface flog
    module procedure hlog, rlog
  end interface
  interface fexp
    module procedure hexp, rexp
  end interface
  interface flog1p
    module procedure hlog1p, rlog1p
  end interface
  interface fexpm1
    module procedure hexpm1, rexpm1
  end interface
  interface fsqrt
    module procedure hsqrt, rsqrt
  end interface
  interface fabs
    module procedure habs, rabs
  end interface
  interface fpnorm
    module procedure hpnorm, rpnorm
  end interface
  interface fdnorm
    module procedure hdnorm, rdnorm
  end interface
  interface flse
    module procedure hlse, rlse
  end interface

  interface operator(+)
    module procedure add_hh, add_hr, add_rh
  end interface
  interface operator(-)
    module procedure sub_hh, sub_hr, sub_rh, neg_h
  end interface
  interface operator(*)
    module procedure mul_hh, mul_hr, mul_rh
  end interface
  interface operator(/)
    module procedure div_hh, div_hr, div_rh
  end interface
  interface operator(**)
    module procedure pow_hr, pow_hh, pow_rh
  end interface

contains

  elemental function hconst(x) result(r)
    real(dp), intent(in) :: x
    type(hd) :: r
    r%v = x
  end function hconst

  pure function hvar(x, i) result(r)
    real(dp), intent(in) :: x
    integer, intent(in) :: i
    type(hd) :: r
    r%v = x
    r%g(i) = 1.0_dp
  end function hvar

  ! f(x) given f, f', f'' at x%v
  pure function hunary(x, f0, f1, f2) result(r)
    type(hd), intent(in) :: x
    real(dp), intent(in) :: f0, f1, f2
    type(hd) :: r
    r%v = f0
    r%g = f1 * x%g
  end function hunary

  pure function add_hh(a, b) result(r)
    type(hd), intent(in) :: a, b
    type(hd) :: r
    r%v = a%v + b%v; r%g = a%g + b%g
  end function add_hh
  pure function add_hr(a, b) result(r)
    type(hd), intent(in) :: a
    real(dp), intent(in) :: b
    type(hd) :: r
    r = a; r%v = a%v + b
  end function add_hr
  pure function add_rh(a, b) result(r)
    real(dp), intent(in) :: a
    type(hd), intent(in) :: b
    type(hd) :: r
    r = b; r%v = a + b%v
  end function add_rh

  pure function sub_hh(a, b) result(r)
    type(hd), intent(in) :: a, b
    type(hd) :: r
    r%v = a%v - b%v; r%g = a%g - b%g
  end function sub_hh
  pure function sub_hr(a, b) result(r)
    type(hd), intent(in) :: a
    real(dp), intent(in) :: b
    type(hd) :: r
    r = a; r%v = a%v - b
  end function sub_hr
  pure function sub_rh(a, b) result(r)
    real(dp), intent(in) :: a
    type(hd), intent(in) :: b
    type(hd) :: r
    r%v = a - b%v; r%g = -b%g
  end function sub_rh
  pure function neg_h(a) result(r)
    type(hd), intent(in) :: a
    type(hd) :: r
    r%v = -a%v; r%g = -a%g
  end function neg_h

  pure function mul_hh(a, b) result(r)
    type(hd), intent(in) :: a, b
    type(hd) :: r
    r%v = a%v * b%v
    r%g = a%v * b%g + b%v * a%g
  end function mul_hh
  pure function mul_hr(a, b) result(r)
    type(hd), intent(in) :: a
    real(dp), intent(in) :: b
    type(hd) :: r
    r%v = a%v * b; r%g = a%g * b
  end function mul_hr
  pure function mul_rh(a, b) result(r)
    real(dp), intent(in) :: a
    type(hd), intent(in) :: b
    type(hd) :: r
    r%v = a * b%v; r%g = a * b%g
  end function mul_rh

  pure function hinv_(x) result(r)
    type(hd), intent(in) :: x
    type(hd) :: r
    real(dp) :: iv
    iv = 1.0_dp / x%v
    r = hunary(x, iv, -iv * iv, 2.0_dp * iv * iv * iv)
  end function hinv_
  pure function div_hh(a, b) result(r)
    type(hd), intent(in) :: a, b
    type(hd) :: r
    r = mul_hh(a, hinv_(b))
  end function div_hh
  pure function div_hr(a, b) result(r)
    type(hd), intent(in) :: a
    real(dp), intent(in) :: b
    type(hd) :: r
    r = mul_hr(a, 1.0_dp / b)
  end function div_hr
  pure function div_rh(a, b) result(r)
    real(dp), intent(in) :: a
    type(hd), intent(in) :: b
    type(hd) :: r
    r = mul_rh(a, hinv_(b))
  end function div_rh

  pure function hexp(x) result(r)
    type(hd), intent(in) :: x
    type(hd) :: r
    real(dp) :: e
    e = exp(x%v)
    r = hunary(x, e, e, e)
  end function hexp
  pure function hlog(x) result(r)
    type(hd), intent(in) :: x
    type(hd) :: r
    r = hunary(x, log(x%v), 1.0_dp / x%v, -1.0_dp / (x%v * x%v))
  end function hlog
  pure function hlog1p(x) result(r)
    type(hd), intent(in) :: x
    type(hd) :: r
    real(dp) :: u, l
    u = 1.0_dp + x%v
    if (u == 1.0_dp) then
      l = x%v
    else
      l = log(u) * x%v / (u - 1.0_dp)
    end if
    r = hunary(x, l, 1.0_dp / u, -1.0_dp / (u * u))
  end function hlog1p
  pure function hexpm1(x) result(r)
    type(hd), intent(in) :: x
    type(hd) :: r
    real(dp) :: e, m
    e = exp(x%v)
    if (abs(x%v) > 0.5_dp) then
      m = e - 1.0_dp
    else if (e == 1.0_dp) then
      m = x%v
    else
      m = (e - 1.0_dp) * x%v / log(e)
    end if
    r = hunary(x, m, e, e)
  end function hexpm1
  pure function hsqrt(x) result(r)
    type(hd), intent(in) :: x
    type(hd) :: r
    real(dp) :: s
    s = sqrt(x%v)
    r = hunary(x, s, 0.5_dp / s, -0.25_dp / (s * x%v))
  end function hsqrt
  pure function habs(x) result(r)
    type(hd), intent(in) :: x
    type(hd) :: r
    if (x%v < 0.0_dp) then
      r = neg_h(x)
    else
      r = x
    end if
  end function habs

  ! x**p, real p (x > 0)
  pure function pow_hr(x, p) result(r)
    type(hd), intent(in) :: x
    real(dp), intent(in) :: p
    type(hd) :: r
    real(dp) :: f0
    f0 = x%v**p
    r = hunary(x, f0, p * f0 / x%v, p * (p - 1.0_dp) * f0 / (x%v * x%v))
  end function pow_hr
  ! x**y = exp(y log x)
  pure function pow_hh(x, y) result(r)
    type(hd), intent(in) :: x, y
    type(hd) :: r
    r = hexp(mul_hh(y, hlog(x)))
  end function pow_hh
  ! c**y = exp(y log c), c > 0
  pure function pow_rh(c, y) result(r)
    real(dp), intent(in) :: c
    type(hd), intent(in) :: y
    type(hd) :: r
    r = hexp(mul_hr(y, log(c)))
  end function pow_rh

  ! ---- real counterparts and extra dual functions for the generic interfaces
  elemental real(dp) function rval(x)
    real(dp), intent(in) :: x
    rval = x
  end function rval
  elemental real(dp) function hval(x)
    type(hd), intent(in) :: x
    hval = x%v
  end function hval
  elemental real(dp) function rlog(x)
    real(dp), intent(in) :: x
    rlog = log(x)
  end function rlog
  elemental real(dp) function rexp(x)
    real(dp), intent(in) :: x
    rexp = exp(x)
  end function rexp
  elemental real(dp) function rlog1p(x)
    real(dp), intent(in) :: x
    real(dp) :: u
    u = 1.0_dp + x
    if (u == 1.0_dp) then
      rlog1p = x
    else
      rlog1p = log(u) * x / (u - 1.0_dp)
    end if
  end function rlog1p
  elemental real(dp) function rexpm1(x)
    real(dp), intent(in) :: x
    real(dp) :: u
    u = exp(x)
    if (abs(x) > 0.5_dp) then
      rexpm1 = u - 1.0_dp
    else if (u == 1.0_dp) then
      rexpm1 = x
    else
      rexpm1 = (u - 1.0_dp) * x / log(u)
    end if
  end function rexpm1
  elemental real(dp) function rsqrt(x)
    real(dp), intent(in) :: x
    rsqrt = sqrt(x)
  end function rsqrt
  elemental real(dp) function rabs(x)
    real(dp), intent(in) :: x
    rabs = abs(x)
  end function rabs
  elemental real(dp) function rpnorm(x)
    real(dp), intent(in) :: x
    rpnorm = 0.5_dp * erfc(-x * 0.70710678118654752440_dp)
  end function rpnorm
  elemental real(dp) function rdnorm(x)
    real(dp), intent(in) :: x
    rdnorm = 0.39894228040143267794_dp * exp(-0.5_dp * x * x)
  end function rdnorm
  pure function hpnorm(x) result(r)
    type(hd), intent(in) :: x
    type(hd) :: r
    real(dp) :: f
    f = rdnorm(x%v)
    r = hunary(x, rpnorm(x%v), f, -x%v * f)
  end function hpnorm
  pure function hdnorm(x) result(r)
    type(hd), intent(in) :: x
    type(hd) :: r
    real(dp) :: f
    f = rdnorm(x%v)
    r = hunary(x, f, -x%v * f, (x%v * x%v - 1.0_dp) * f)
  end function hdnorm
  ! log(exp(p) + exp(q))
  elemental real(dp) function rlse(p, q)
    real(dp), intent(in) :: p, q
    rlse = max(p, q) + rlog1p(exp(min(p, q) - max(p, q)))
  end function rlse
  pure function hlse(p, q) result(r)
    type(hd), intent(in) :: p, q
    type(hd) :: r
    if (p%v >= q%v) then
      r = p + hlog1p(hexp(q - p))
    else
      r = q + hlog1p(hexp(p - q))
    end if
  end function hlse

end module fortvine_ad

! hd1 / hd2: width-specialised dual numbers for the scores of one- and two-parameter
! families (fortvine_adk.inc, compiled once per width)
#define AD_MOD fortvine_ad1
#define AD_T hd1
#define AD_N 1
#define AD_VAR h1var
#define AD_CONST h1const
#include "fortvine_adk.inc"
#undef AD_MOD
#undef AD_T
#undef AD_N
#undef AD_VAR
#undef AD_CONST
#define AD_MOD fortvine_ad2
#define AD_T hd2
#define AD_N 2
#define AD_VAR h2var
#define AD_CONST h2const
#include "fortvine_adk.inc"
