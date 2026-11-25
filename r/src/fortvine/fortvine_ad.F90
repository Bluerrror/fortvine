! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Shahab Aldin Shojaeezadeh. Part of fortvine (see LICENSE). If you use it, please cite:
! Shojaeezadeh, S. A. (2026). fortvine: a complete compendium to work on vine copula modeling.
! fortvine_ad.f90 -- dual numbers: value and exact gradient with respect to up to 3
! copula parameters (forward-mode differentiation).  The dual types are compiled from one source,
! fortvine_adk.inc, once per width: hd (3 slots, any family), hd1 and hd2 (one- and two-parameter
! families: their scores evaluate 1.1-2x faster without the unused slots).  fortvine_ad adds the
! real(dp) specifics of the generic functions, so a family source compiles for real values as well.
!
! Every operation applies the analytic chain rule, so per-observation scores
! d log c / d theta are exact to rounding -- no finite differences.  The MLE in
! fortvine.f90 uses them for Fisher scoring (BHHH): curvature = sum of score outer
! products, which avoids carrying a Hessian through every operation.
#define AD_MOD fortvine_adh
#define AD_T hd
#define AD_N 3
#define AD_VAR hvar
#define AD_CONST hconst
#include "fortvine_adk.inc"
#undef AD_MOD
#undef AD_T
#undef AD_N
#undef AD_VAR
#undef AD_CONST
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
#undef AD_MOD
#undef AD_T
#undef AD_N
#undef AD_VAR
#undef AD_CONST

module fortvine_ad
  use fv_util, only: log1p_, expm1_, pnorm, dnorm
  use fortvine_adh
  implicit none
  private
  integer, parameter :: dp = kind(1.0d0)

  public :: hd, hvar, hconst
  public :: operator(+), operator(-), operator(*), operator(/), operator(**)
  public :: fval, flog, fexp, flog1p, fexpm1, fsqrt, fabs, fpnorm, fdnorm, flse

  ! the real(dp) specifics, added to fortvine_adh's generic functions
  interface fval
    module procedure rval
  end interface
  interface flog
    module procedure rlog
  end interface
  interface fexp
    module procedure rexp
  end interface
  interface flog1p
    module procedure log1p_
  end interface
  interface fexpm1
    module procedure expm1_
  end interface
  interface fsqrt
    module procedure rsqrt
  end interface
  interface fabs
    module procedure rabs
  end interface
  interface fpnorm
    module procedure pnorm
  end interface
  interface fdnorm
    module procedure dnorm
  end interface
  interface flse
    module procedure rlse
  end interface

contains

  elemental real(dp) function rval(x)
    real(dp), intent(in) :: x
    rval = x
  end function rval
  elemental real(dp) function rlog(x)
    real(dp), intent(in) :: x
    rlog = log(x)
  end function rlog
  elemental real(dp) function rexp(x)
    real(dp), intent(in) :: x
    rexp = exp(x)
  end function rexp
  elemental real(dp) function rsqrt(x)
    real(dp), intent(in) :: x
    rsqrt = sqrt(x)
  end function rsqrt
  elemental real(dp) function rabs(x)
    real(dp), intent(in) :: x
    rabs = abs(x)
  end function rabs
  elemental real(dp) function rlse(p, q)
    real(dp), intent(in) :: p, q
    rlse = max(p, q) + log1p_(exp(min(p, q) - max(p, q)))
  end function rlse

end module fortvine_ad
