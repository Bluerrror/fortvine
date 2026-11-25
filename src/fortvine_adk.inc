! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Shahab Aldin Shojaeezadeh. Part of fortvine (see LICENSE). If you use it, please cite:
! Shojaeezadeh, S. A. (2026). fortvine: a complete compendium to work on vine copula modeling.
! fortvine_adk.inc -- dual numbers of a fixed width AD_N (1, 2 or 3 copula parameters): value
! and exact gradient, compiled once per width by fortvine_ad.f90 (types hd1, hd2, hd).
module AD_MOD
  use fv_util, only: pnorm, dnorm
  implicit none
  private
  integer, parameter :: dp = kind(1.0d0)

  type, public :: AD_T
    real(dp) :: v = 0.0_dp
    real(dp) :: g(AD_N) = 0.0_dp
  end type AD_T

  public :: operator(+), operator(-), operator(*), operator(/), operator(**)
  public :: AD_VAR, AD_CONST
  public :: fval, flog, fexp, flog1p, fexpm1, fsqrt, fabs, fpnorm, fdnorm, flse

  interface fval
    module procedure kval
  end interface
  interface flog
    module procedure klog
  end interface
  interface fexp
    module procedure kexp
  end interface
  interface flog1p
    module procedure klog1p
  end interface
  interface fexpm1
    module procedure kexpm1
  end interface
  interface fsqrt
    module procedure ksqrt
  end interface
  interface fabs
    module procedure kabs
  end interface
  interface fpnorm
    module procedure kpnorm
  end interface
  interface fdnorm
    module procedure kdnorm
  end interface
  interface flse
    module procedure klse
  end interface
  interface operator(+)
    module procedure add_kk, add_kr, add_rk
  end interface
  interface operator(-)
    module procedure sub_kk, sub_kr, sub_rk, neg_k
  end interface
  interface operator(*)
    module procedure mul_kk, mul_kr, mul_rk
  end interface
  interface operator(/)
    module procedure div_kk, div_kr, div_rk
  end interface
  interface operator(**)
    module procedure pow_kr, pow_kk, pow_rk
  end interface

contains

  elemental function AD_CONST(x) result(r)
    real(dp), intent(in) :: x
    type(AD_T) :: r
    r%v = x
  end function AD_CONST

  pure function AD_VAR(x, i) result(r)
    real(dp), intent(in) :: x
    integer, intent(in) :: i
    type(AD_T) :: r
    r%v = x
    if (i <= AD_N) r%g(i) = 1.0_dp
  end function AD_VAR

  ! f(x) given f, f', f'' at x%v
  pure function kun(x, f0, f1, f2) result(r)
    type(AD_T), intent(in) :: x
    real(dp), intent(in) :: f0, f1, f2
    type(AD_T) :: r
    r%v = f0
    r%g = f1 * x%g
  end function kun

  pure function add_kk(a, b) result(r)
    type(AD_T), intent(in) :: a, b
    type(AD_T) :: r
    r%v = a%v + b%v; r%g = a%g + b%g
  end function add_kk
  pure function add_kr(a, b) result(r)
    type(AD_T), intent(in) :: a
    real(dp), intent(in) :: b
    type(AD_T) :: r
    r = a; r%v = a%v + b
  end function add_kr
  pure function add_rk(a, b) result(r)
    real(dp), intent(in) :: a
    type(AD_T), intent(in) :: b
    type(AD_T) :: r
    r = b; r%v = a + b%v
  end function add_rk
  pure function sub_kk(a, b) result(r)
    type(AD_T), intent(in) :: a, b
    type(AD_T) :: r
    r%v = a%v - b%v; r%g = a%g - b%g
  end function sub_kk
  pure function sub_kr(a, b) result(r)
    type(AD_T), intent(in) :: a
    real(dp), intent(in) :: b
    type(AD_T) :: r
    r = a; r%v = a%v - b
  end function sub_kr
  pure function sub_rk(a, b) result(r)
    real(dp), intent(in) :: a
    type(AD_T), intent(in) :: b
    type(AD_T) :: r
    r%v = a - b%v; r%g = -b%g
  end function sub_rk
  pure function neg_k(a) result(r)
    type(AD_T), intent(in) :: a
    type(AD_T) :: r
    r%v = -a%v; r%g = -a%g
  end function neg_k
  pure function mul_kk(a, b) result(r)
    type(AD_T), intent(in) :: a, b
    type(AD_T) :: r
    r%v = a%v * b%v
    r%g = a%v * b%g + b%v * a%g
  end function mul_kk
  pure function mul_kr(a, b) result(r)
    type(AD_T), intent(in) :: a
    real(dp), intent(in) :: b
    type(AD_T) :: r
    r%v = a%v * b; r%g = a%g * b
  end function mul_kr
  pure function mul_rk(a, b) result(r)
    real(dp), intent(in) :: a
    type(AD_T), intent(in) :: b
    type(AD_T) :: r
    r%v = a * b%v; r%g = a * b%g
  end function mul_rk
  pure function kinv(x) result(r)
    type(AD_T), intent(in) :: x
    type(AD_T) :: r
    real(dp) :: iv
    iv = 1.0_dp / x%v
    r = kun(x, iv, -iv * iv, 2.0_dp * iv * iv * iv)
  end function kinv
  pure function div_kk(a, b) result(r)
    type(AD_T), intent(in) :: a, b
    type(AD_T) :: r
    r = mul_kk(a, kinv(b))
  end function div_kk
  pure function div_kr(a, b) result(r)
    type(AD_T), intent(in) :: a
    real(dp), intent(in) :: b
    type(AD_T) :: r
    r = mul_kr(a, 1.0_dp / b)
  end function div_kr
  pure function div_rk(a, b) result(r)
    real(dp), intent(in) :: a
    type(AD_T), intent(in) :: b
    type(AD_T) :: r
    r = mul_rk(a, kinv(b))
  end function div_rk

  pure function kexp(x) result(r)
    type(AD_T), intent(in) :: x
    type(AD_T) :: r
    real(dp) :: e
    e = exp(x%v)
    r = kun(x, e, e, e)
  end function kexp
  pure function klog(x) result(r)
    type(AD_T), intent(in) :: x
    type(AD_T) :: r
    r = kun(x, log(x%v), 1.0_dp / x%v, -1.0_dp / (x%v * x%v))
  end function klog
  pure function klog1p(x) result(r)
    type(AD_T), intent(in) :: x
    type(AD_T) :: r
    real(dp) :: u, l
    u = 1.0_dp + x%v
    if (u == 1.0_dp) then
      l = x%v
    else
      l = log(u) * x%v / (u - 1.0_dp)
    end if
    r = kun(x, l, 1.0_dp / u, -1.0_dp / (u * u))
  end function klog1p
  pure function kexpm1(x) result(r)
    type(AD_T), intent(in) :: x
    type(AD_T) :: r
    real(dp) :: e, m
    e = exp(x%v)
    if (abs(x%v) > 0.5_dp) then
      m = e - 1.0_dp
    else if (e == 1.0_dp) then
      m = x%v
    else
      m = (e - 1.0_dp) * x%v / log(e)
    end if
    r = kun(x, m, e, e)
  end function kexpm1
  pure function ksqrt(x) result(r)
    type(AD_T), intent(in) :: x
    type(AD_T) :: r
    real(dp) :: s
    s = sqrt(x%v)
    r = kun(x, s, 0.5_dp / s, -0.25_dp / (s * x%v))
  end function ksqrt
  pure function kabs(x) result(r)
    type(AD_T), intent(in) :: x
    type(AD_T) :: r
    if (x%v < 0.0_dp) then
      r = neg_k(x)
    else
      r = x
    end if
  end function kabs
  ! x**p, real p (x > 0)
  pure function pow_kr(x, p) result(r)
    type(AD_T), intent(in) :: x
    real(dp), intent(in) :: p
    type(AD_T) :: r
    real(dp) :: f0
    f0 = x%v**p
    r = kun(x, f0, p * f0 / x%v, p * (p - 1.0_dp) * f0 / (x%v * x%v))
  end function pow_kr
  pure function pow_kk(x, y) result(r)
    type(AD_T), intent(in) :: x, y
    type(AD_T) :: r
    r = kexp(mul_kk(y, klog(x)))
  end function pow_kk
  pure function pow_rk(c, y) result(r)
    real(dp), intent(in) :: c
    type(AD_T), intent(in) :: y
    type(AD_T) :: r
    r = kexp(mul_kr(y, log(c)))
  end function pow_rk
  elemental real(dp) function kval(x)
    type(AD_T), intent(in) :: x
    kval = x%v
  end function kval
  pure function kpnorm(x) result(r)
    type(AD_T), intent(in) :: x
    type(AD_T) :: r
    real(dp) :: f
    f = dnorm(x%v)
    r = kun(x, pnorm(x%v), f, -x%v * f)
  end function kpnorm
  pure function kdnorm(x) result(r)
    type(AD_T), intent(in) :: x
    type(AD_T) :: r
    real(dp) :: f
    f = dnorm(x%v)
    r = kun(x, f, -x%v * f, (x%v * x%v - 1.0_dp) * f)
  end function kdnorm
  ! log(exp(p) + exp(q))
  pure function klse(p, q) result(r)
    type(AD_T), intent(in) :: p, q
    type(AD_T) :: r
    if (p%v >= q%v) then
      r = p + klog1p(kexp(q - p))
    else
      r = q + klog1p(kexp(p - q))
    end if
  end function klse
end module AD_MOD
