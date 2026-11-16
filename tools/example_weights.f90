! example_weights -- data for the figures of docs/examples/weights-discrete.md
! (tools/example_figures.py weights).  Re-runs both parts of example/ex05_weights_discrete.f90 with
! the same seeds and writes, into the current directory:
!   w_sample.txt   u1 u2 weight  (two-regime Clayton sample, exponential time-decay weights)
!   w_profile.txt  theta  l(theta)/n  l_w(theta)/sum(w)   Clayton log-likelihood profiles
!   d_data.txt     z1 x2 z3 u2 u2minus   (D-vine sample with variable 2 a Poisson(2) count;
!                  z = normal scores of the continuous variables)
!   d_ros.txt      randomised Rosenblatt residuals r1 r2 r3 (seed 1), then the non-randomised
!                  F(x|.) of the three columns and F(x-|.) of the three columns
! stdout: theta estimates (equal, decay weights, each regime alone).
program example_weights
  use fortvine_api
  implicit none
  integer, parameter :: dp = kind(1.0d0), n = 3000, d = 3, NT = 221
  real(dp), parameter :: LAM = 2.0_dp
  type(fv_bicop) :: cop, fit
  type(fv_vinecop) :: truth, vd
  real(dp), allocatable :: u1(:), u2(:), w(:), wt(:), lp(:)
  real(dp), allocatable :: z(:, :), u(:, :), um(:, :), r(:, :), rp(:, :), rm(:, :)
  real(dp) :: th(NT), pr(NT, 2)
  integer :: i, k, un

  ! ---------------------------------------------------------------- (1) weights, = ex05
  allocate(u1(n), u2(n), w(n), wt(n), lp(n))
  call seed_rng(11)
  call random_number(u1); call random_number(w)
  cop%family = FV_CLAYTON
  cop%par(1) = 1.0_dp; call cop%hinv1(u1(1:n/2), w(1:n/2), u2(1:n/2))
  cop%par(1) = 4.0_dp; call cop%hinv1(u1(n/2+1:), w(n/2+1:), u2(n/2+1:))
  wt = [(0.998_dp**(n - i), i = 1, n)]
  fit%family = FV_CLAYTON
  call fit%fit(u1, u2);                     print '(a,f9.5)', 'theta_equal ', fit%par(1)
  call fit%fit(u1, u2, weights=wt);         print '(a,f9.5)', 'theta_decay ', fit%par(1)
  call fit%fit(u1(1:n/2), u2(1:n/2));       print '(a,f9.5)', 'theta_old   ', fit%par(1)
  call fit%fit(u1(n/2+1:), u2(n/2+1:));     print '(a,f9.5)', 'theta_new   ', fit%par(1)
  open(newunit=un, file='w_sample.txt', status='replace')
  do i = 1, n
    write(un, '(3es24.16)') u1(i), u2(i), wt(i)
  end do
  close(un)
  th = [(0.5_dp + 6.0_dp * (i - 1) / (NT - 1), i = 1, NT)]
  do k = 1, NT
    cop%par(1) = th(k)
    call cop%logpdf(u1, u2, lp)
    pr(k, :) = [sum(lp) / n, sum(wt * lp) / sum(wt)]
  end do
  open(newunit=un, file='w_profile.txt', status='replace')
  do k = 1, NT
    write(un, '(3es20.11)') th(k), pr(k, :)
  end do
  close(un)

  ! ---------------------------------------------------------------- (2) a discrete variable, = ex05
  allocate(z(n, d), u(n, d), um(n, d), r(n, d), rp(n, d), rm(n, d))
  call truth%set_dvine([1, 2, 3])
  truth%fam(1:2, 1) = [FV_GAUSSIAN, FV_CLAYTON]
  truth%par(1, 1, 1) = 0.7_dp; truth%par(1, 2, 1) = 2.0_dp
  truth%fam(1, 2) = FV_FRANK; truth%par(1, 1, 2) = 2.0_dp
  call truth%simulate(n, z, seed=5)
  u = z; um = z
  do i = 1, n
    k = qpois(z(i, 2))
    u(i, 2) = ppois(k); um(i, 2) = ppois(k - 1)
  end do
  call vd%set_dvine([1, 2, 3])
  call vd%select(u, discrete=[.false., .true., .false.], u_minus=um)
  call vd%rosenblatt(u, r, u_minus=um, seed=1)
  call vd%rosenblatt(u, rp, u_minus=um, randomize=.false., out_minus=rm)
  open(newunit=un, file='d_data.txt', status='replace')
  do i = 1, n
    write(un, '(es24.16,i4,3es24.16)') qn(u(i, 1)), qpois(z(i, 2)), qn(u(i, 3)), u(i, 2), um(i, 2)
  end do
  close(un)
  open(newunit=un, file='d_ros.txt', status='replace')
  do i = 1, n
    write(un, '(9es24.16)') r(i, :), rp(i, :), rm(i, :)
  end do
  close(un)
  print '(a,3f8.4)', 'rosenblatt means', sum(r, 1) / n

contains
  subroutine seed_rng(seed)
    integer, intent(in) :: seed
    integer :: ns, j
    integer, allocatable :: s(:)
    call random_seed(size=ns)
    s = [(seed + 7919 * j, j = 1, ns)]
    call random_seed(put=s)
  end subroutine seed_rng

  pure real(dp) function ppois(k)
    integer, intent(in) :: k
    real(dp) :: p
    integer :: j
    ppois = 0.0_dp
    if (k < 0) return
    p = exp(-LAM); ppois = p
    do j = 1, k
      p = p * LAM / j; ppois = ppois + p
    end do
  end function ppois

  pure integer function qpois(p)
    real(dp), intent(in) :: p
    qpois = 0
    do while (ppois(qpois) < p)
      qpois = qpois + 1
    end do
  end function qpois

  ! standard normal quantile by bisection on erfc (only for the plot coordinates)
  pure real(dp) function qn(p)
    real(dp), intent(in) :: p
    real(dp) :: lo, hi
    integer :: j
    lo = -9.0_dp; hi = 9.0_dp
    do j = 1, 80
      qn = 0.5_dp * (lo + hi)
      if (0.5_dp * erfc(-qn / sqrt(2.0_dp)) < p) then
        lo = qn
      else
        hi = qn
      end if
    end do
  end function qn
end program example_weights
