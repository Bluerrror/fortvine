! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Shahab Aldin Shojaeezadeh. Part of fortvine (see LICENSE). If you use it, please cite:
! Shojaeezadeh, S. A. (2026). fortvine: a complete compendium to work on vine copula modeling.
! newton1.inc -- included in the CONTAINS part of mle1 and mle1d, so newton1 calls the host's
! internal `pass` directly. (Passing an internal procedure as an argument makes gfortran build
! a trampoline, and a shared library with a trampoline needs an executable stack, which
! modern glibc refuses to load.)
    ! Newton's method on [lb, ub] from the closed-form start x0; each pass gives the log-likelihood f,
    ! score g and curvature h.  The iteration stops when the quadratic model predicts a remaining gain
    ! g^2 / 2|h| below 1e-6 log-likelihood units (the exact MLE to that accuracy), so a start inside
    ! the Le Cam one-step neighbourhood costs two passes.  Non-concave points step uphill; any step that
    ! lowers the likelihood is halved; the box is respected.
    subroutine newton1(x0, lb, ub, th, f, ok)
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
