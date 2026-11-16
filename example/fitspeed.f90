! Time a single MLE fit per family (n = 3000, one thread) and count likelihood evaluations.
program fitspeed
  use fortvine_api
  implicit none
  integer, parameter :: dp = kind(1.0d0), n = 3000
  type(fv_bicop) :: truth, est
  type(fv_controls) :: ctl
  real(dp) :: w1(n), w2(n), u2(n)
  integer :: k
  integer(8) :: t0, t1, rate
  ctl%nthreads = 1
  call random_seed(put=[(12345 + 7 * k, k = 1, 64)])
  do k = 1, 23
    call setfam(k, truth)
    call random_number(w1); call random_number(w2)
    call truth%hinv1(w1, w2, u2)
    est%family = truth%family; est%rotation = 0
    call system_clock(t0, rate)
    call est%fit(w1, u2, ctl)
    call system_clock(t1)
    print '(a13,f9.2,a,3f9.4,a,f12.4)', trim(FV_FAMILY_NAME(truth%family)), 1e3 * real(t1 - t0, dp) / rate, ' ms  par', est%par, '  ll', est%loglik
  end do
contains
  include 'setfam.inc'
end program fitspeed
