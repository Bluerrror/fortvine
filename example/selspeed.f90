! Selection cost by candidate family: every family (both rotations of the tau sign) fitted to
! data of every truth family, n = 5000, one thread.  Prints ms per candidate summed over truths.
program selspeed
  use fortvine_api
  implicit none
  integer, parameter :: dp = kind(1.0d0), n = 5000, nt = 23
  type(fv_bicop) :: truth, est
  type(fv_controls) :: ctl
  real(dp) :: w1(n), w2(n), u2(n), tf(0:25), tsel(nt), llsum(0:25)
  integer :: k, f, j, r(2)
  integer(8) :: t0, t1, rate
  ctl%nthreads = 1
  tf = 0; llsum = 0
  call random_seed(put=[(12345 + 7 * k, k = 1, 64)])
  do k = 1, nt
    call setfam(k, truth)
    call random_number(w1); call random_number(w2)
    call truth%hinv1(w1, w2, u2)
    do f = 1, 25
      r = [0, 180]
      do j = 1, 2
        if (j == 2 .and. any(f == [1, 2, 5, 12, 13, 14, 24, 25])) exit
        est%family = f; est%rotation = r(j)
        call system_clock(t0, rate)
        call est%fit(w1, u2, ctl)
        call system_clock(t1)
        tf(f) = tf(f) + 1e3_dp * real(t1 - t0, dp) / rate
      end do
    end do
    ctl%nfam = FV_NCLASSIC
    call system_clock(t0, rate)
    call est%select(w1, u2, ctl)
    call system_clock(t1)
    tsel(k) = 1e3_dp * real(t1 - t0, dp) / rate
    print '(a13,a,f9.1,a,a13,i4,3f10.5,f14.5)', trim(FV_FAMILY_NAME(truth%family)), '  select(26) ms', tsel(k), '  ->', &
          trim(FV_FAMILY_NAME(est%family)), est%rotation, est%par, est%loglik
  end do
  print '(a)', '--- ms per candidate family, summed over 23 truths (both rotations)'
  do f = 1, 25
    print '(a13,f10.1)', trim(FV_FAMILY_NAME(f)), tf(f)
  end do
  print '(a,f10.1,a,f10.1)', 'total fit ms', sum(tf), '   total select ms', sum(tsel)
contains
  include 'setfam.inc'
end program selspeed
