! Per-family kernel speed, single thread: ns per evaluation on 1e6 points.
program famspeed
  use fortvine_api
  use fortvine, only: fv_bicop_eval
  implicit none
  integer, parameter :: dp = kind(1.0d0), n = 1000000
  real(dp), allocatable :: u1(:), u2(:), out(:), uu(:, :)
  type(fv_bicop) :: b
  integer :: k, what
  integer(8) :: t0, t1, rate
  real(dp) :: ns(6)
  character(len=8) :: names(6) = [character(len=8) :: 'logpdf', 'pdf', 'hfunc1', 'hfunc2', 'hinv1', 'hinv2']
  allocate(u1(n), u2(n), out(n))
  call random_number(u1); call random_number(u2)
  u1 = 0.001_dp + 0.998_dp * u1; u2 = 0.001_dp + 0.998_dp * u2
  allocate(uu(2, n)); uu(1, :) = u1; uu(2, :) = u2
  print '(a13,6a10)', 'family', names
  do k = 1, 23
    call setfam(k, b)
    do what = 0, 5
      call system_clock(t0, rate)
      call fv_bicop_eval(b%family, 0, b%par, n, uu, what, out, 1)
      call system_clock(t1)
      ns(what + 1) = 1.0e9_dp * real(t1 - t0, dp) / rate / n
    end do
    print '(a13,6f10.1)', trim(FV_FAMILY_NAME(b%family)), ns
  end do
contains
  include 'setfam.inc'

end program famspeed
