! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Shahab Aldin Shojaeezadeh. Part of fortvine (see LICENSE). If you use it, please cite:
! Shojaeezadeh, S. A. (2026). fortvine: a complete compendium to work on vine copula modeling.
! fv_util.f90 -- the one copy of the small numerics every fortvine module shares: log1p / expm1, the standard
! normal (pdf, cdf, quantile), digamma / trigamma, stable sorting, binary search, median, NSE / KGE and
! upper-casing.  No dependencies; every other module uses it.
module fv_util
  use iso_c_binding, only: c_double
  implicit none
  private
  integer, parameter :: dp = c_double
  real(dp), parameter :: SQRT2 = 1.4142135623730951_dp
  public :: log1p_, expm1_, pnorm, dnorm, qnorm, digamma, trigamma, sort_ab, argsort, sort, median, lower_bound, &
            upper_bound, nse, kge, upper

  interface sort
    module procedure sort_r, sort_i
  end interface sort

contains

  ! log(1 + x) without losing x when 1 + x rounds (Goldberg's correction)
  elemental real(dp) function log1p_(x)
    !DIR$ ATTRIBUTES FORCEINLINE :: log1p_
    real(dp), intent(in) :: x
    real(dp) :: u
    u = 1.0_dp + x
    if (u == 1.0_dp) then
      log1p_ = x
    else
      log1p_ = log(u) * x / (u - 1.0_dp)
    end if
  end function log1p_

  ! exp(x) - 1, accurate near 0 (Kahan's (u - 1) x / log u)
  elemental real(dp) function expm1_(x)
    !DIR$ ATTRIBUTES FORCEINLINE :: expm1_
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

  elemental real(dp) function pnorm(x)
    !DIR$ ATTRIBUTES FORCEINLINE :: pnorm
    real(dp), intent(in) :: x
    pnorm = 0.5_dp * erfc(-x / SQRT2)
  end function pnorm

  elemental real(dp) function dnorm(x)
    !DIR$ ATTRIBUTES FORCEINLINE :: dnorm
    real(dp), intent(in) :: x
    dnorm = exp(-0.5_dp * x * x) * 0.3989422804014327_dp
  end function dnorm

  ! standard normal quantile: Wichura (1988) AS241 PPND16, ~1e-16 relative; -huge / +huge at p <= 0 / p >= 1
  elemental real(dp) function qnorm(p)
    real(dp), intent(in) :: p
    real(dp) :: q, r, v
    if (p <= 0.0_dp) then
      qnorm = -huge(1.0_dp); return
    else if (p >= 1.0_dp) then
      qnorm = huge(1.0_dp); return
    end if
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

  ! digamma psi(x), x > 0: recurrence up to 12, then the asymptotic series (~1e-14)
  elemental real(dp) function digamma(x0)
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

  ! trigamma psi'(x), x > 0 (recurrence + asymptotic series)
  elemental real(dp) function trigamma(x0)
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

  ! stable ascending sort by (a, b) lexicographically, permuting a, b and ix together; cnt grows by the number
  ! of inversions (pairs i < j with (a, b)(j) < (a, b)(i)).  Insertion sort of runs of RUN, then bottom-up merges on
  ! contiguous keys.  A stable sort's order and inversion count are unique, so any stable sort gives the same.
  pure subroutine sort_ab(n, a, b, ix, cnt)
    integer, intent(in) :: n
    real(dp), intent(inout) :: a(n), b(n)
    integer, intent(inout) :: ix(n)
    integer(8), intent(inout) :: cnt
    integer, parameter :: RUN = 32
    real(dp), allocatable :: a2(:), b2(:)
    integer, allocatable :: i2(:)
    real(dp) :: ta, tb
    integer :: lo, mid, hi, i, j, p, q, k, ti, w
    do lo = 1, n, RUN
      hi = min(lo + RUN - 1, n)
      do i = lo + 1, hi
        ta = a(i); tb = b(i); ti = ix(i); j = i - 1
        do while (j >= lo)
          if (.not. (ta < a(j) .or. (ta == a(j) .and. tb < b(j)))) exit
          a(j + 1) = a(j); b(j + 1) = b(j); ix(j + 1) = ix(j); j = j - 1
        end do
        cnt = cnt + (i - 1 - j)
        a(j + 1) = ta; b(j + 1) = tb; ix(j + 1) = ti
      end do
    end do
    if (n <= RUN) return
    allocate(a2(n), b2(n), i2(n))
    w = RUN
    do while (w < n)
      do lo = 1, n, 2 * w
        mid = min(lo + w - 1, n); hi = min(lo + 2 * w - 1, n)
        p = lo; q = mid + 1; k = lo
        do while (p <= mid .and. q <= hi)
          if (a(q) < a(p) .or. (a(q) == a(p) .and. b(q) < b(p))) then
            a2(k) = a(q); b2(k) = b(q); i2(k) = ix(q); q = q + 1
            cnt = cnt + (mid - p + 1)
          else
            a2(k) = a(p); b2(k) = b(p); i2(k) = ix(p); p = p + 1
          end if
          k = k + 1
        end do
        a2(k:k + mid - p) = a(p:mid); b2(k:k + mid - p) = b(p:mid); i2(k:k + mid - p) = ix(p:mid)
        k = k + mid - p + 1
        a2(k:hi) = a(q:hi); b2(k:hi) = b(q:hi); i2(k:hi) = ix(q:hi)
      end do
      a = a2; b = b2; ix = i2
      w = 2 * w
    end do
  end subroutine sort_ab

  ! the stable ascending order of x (ties keep their order): x(argsort(x)) is sorted
  pure function argsort(x) result(ix)
    real(dp), intent(in) :: x(:)
    integer :: ix(size(x))
    real(dp) :: a(size(x)), z(size(x))
    integer(8) :: cnt
    integer :: i
    a = x; z = 0.0_dp; cnt = 0
    ix = [(i, i = 1, size(x))]
    call sort_ab(size(x), a, z, ix, cnt)
  end function argsort

  ! sort x in place, ascending (descending = .true.: descending)
  pure subroutine sort_r(x, descending)
    real(dp), intent(inout) :: x(:)
    logical, intent(in), optional :: descending
    real(dp) :: z(size(x))
    integer :: ix(size(x)), i
    integer(8) :: cnt
    z = 0.0_dp; cnt = 0
    ix = [(i, i = 1, size(x))]
    call sort_ab(size(x), x, z, ix, cnt)
    if (present(descending)) then
      if (descending) x = x(size(x):1:-1)
    end if
  end subroutine sort_r

  pure subroutine sort_i(x, descending)
    integer, intent(inout) :: x(:)
    logical, intent(in), optional :: descending
    real(dp) :: a(size(x))
    a = real(x, dp)
    x = x(argsort(a))
    if (present(descending)) then
      if (descending) x = x(size(x):1:-1)
    end if
  end subroutine sort_i

  ! median: the middle value, or the mean of the two middle values for an even count
  pure real(dp) function median(x)
    real(dp), intent(in) :: x(:)
    real(dp) :: s(size(x))
    integer :: n
    n = size(x)
    if (n == 0) then
      median = 0.0_dp; return
    end if
    s = x; call sort_r(s)
    if (mod(n, 2) == 1) then
      median = s((n + 1) / 2)
    else
      median = 0.5_dp * (s(n / 2) + s(n / 2 + 1))
    end if
  end function median

  ! first k with xs(k) >= x in sorted xs (searchsorted 'left'), 1..n+1
  pure integer function lower_bound(xs, x)
    real(dp), intent(in) :: xs(:), x
    integer :: lo, hi, mid
    lo = 1; hi = size(xs) + 1
    do while (lo < hi)
      mid = (lo + hi) / 2
      if (xs(mid) < x) then
        lo = mid + 1
      else
        hi = mid
      end if
    end do
    lower_bound = lo
  end function lower_bound

  ! first k with xs(k) > x in sorted xs (searchsorted 'right'), 1..n+1
  pure integer function upper_bound(xs, x)
    real(dp), intent(in) :: xs(:), x
    integer :: lo, hi, mid
    lo = 1; hi = size(xs) + 1
    do while (lo < hi)
      mid = (lo + hi) / 2
      if (xs(mid) <= x) then
        lo = mid + 1
      else
        hi = mid
      end if
    end do
    upper_bound = lo
  end function upper_bound

  ! Nash-Sutcliffe efficiency of sim against obs
  pure real(dp) function nse(obs, sim)
    real(dp), intent(in) :: obs(:), sim(:)
    nse = 1.0_dp - sum((sim - obs)**2) / max(sum((obs - sum(obs) / size(obs))**2), 1e-300_dp)
  end function nse

  ! Kling-Gupta efficiency (Gupta et al. 2009) of sim against obs
  pure real(dp) function kge(obs, sim)
    real(dp), intent(in) :: obs(:), sim(:)
    real(dp) :: mo, ms, so, ss, r
    integer :: n
    n = size(obs)
    mo = sum(obs) / n; ms = sum(sim) / n
    so = sqrt(sum((obs - mo)**2) / (n - 1)); ss = sqrt(sum((sim - ms)**2) / (n - 1))
    r = sum((obs - mo) * (sim - ms)) / sqrt(max(sum((obs - mo)**2) * sum((sim - ms)**2), 1e-300_dp))
    kge = 1.0_dp - sqrt((r - 1)**2 + (ss / max(so, 1e-300_dp) - 1)**2 + (ms / max(mo, 1e-300_dp) - 1)**2)
  end function kge

  ! s with a-z upper-cased (same length)
  pure function upper(s) result(r)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: r
    integer :: i, c
    r = s
    do i = 1, len(s)
      c = iachar(s(i:i))
      if (c >= 97 .and. c <= 122) r(i:i) = achar(c - 32)
    end do
  end function upper

end module fv_util
