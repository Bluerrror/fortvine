! example_tails -- data for the figure of docs/examples/tails.md (tools/example_figures.py tails).
! The five copulas of example/ex04_tails.f90 plus a Gaussian copula with the same Kendall tau as the
! Gumbel one (rho = sin(pi tau / 2), tau = 1/2).  Writes, into the current directory:
!   tails_curve.txt  q, then lambda_L(q) = C(q,q)/q and lambda_U(q) = (2q - 1 + C(1-q,1-q))/q for the
!                    six copulas (12 columns), q log-spaced from 1e-6 to 0.5
!   tails_emp.txt    q, then the empirical lambda_L(q), lambda_U(q) of the 200 000-point samples of
!                    ex04 (same seeds) for the five copulas (10 columns)
program example_tails
  use fortvine_api
  implicit none
  integer, parameter :: dp = kind(1.0d0), n = 200000, nc = 6, NQ = 121, NE = 9
  real(dp), parameter :: QE(NE) = [0.001_dp, 0.002_dp, 0.005_dp, 0.01_dp, 0.02_dp, 0.05_dp, 0.1_dp, 0.2_dp, 0.4_dp]
  type(fv_bicop) :: cop(nc)
  real(dp), allocatable :: u1(:), w(:), u2(:)
  real(dp) :: q(NQ), cl(NQ), cu(NQ), lam(NQ, 2 * nc), emp(NE, 10)
  integer :: k, i, un

  allocate(u1(n), w(n), u2(n))
  cop(1)%family = FV_GUMBEL;  cop(1)%par = [2.0_dp, 0.0_dp, 0.0_dp]
  cop(2)%family = FV_CLAYTON; cop(2)%par = [2.0_dp, 0.0_dp, 0.0_dp]
  cop(3)%family = FV_CLAYTON; cop(3)%par = [2.0_dp, 0.0_dp, 0.0_dp]; cop(3)%rotation = 180
  cop(4)%family = FV_TAWN;    cop(4)%par = [0.4_dp, 0.9_dp, 3.0_dp]
  cop(5)%family = 20300;      cop(5)%par = [3.0_dp, 0.5_dp, 0.9_dp]
  cop(6)%family = FV_GAUSSIAN; cop(6)%par = [sin(acos(-1.0_dp) / 4.0_dp), 0.0_dp, 0.0_dp]

  q = [(10.0_dp**(-6.0_dp + (6.0_dp - log10(2.0_dp)) * (i - 1) / (NQ - 1)), i = 1, NQ)]
  do k = 1, nc
    call cop(k)%cdf(q, q, cl)
    call cop(k)%cdf(1.0_dp - q, 1.0_dp - q, cu)
    lam(:, 2 * k - 1) = cl / q
    lam(:, 2 * k) = (2.0_dp * q - 1.0_dp + cu) / q
  end do
  open(newunit=un, file='tails_curve.txt', status='replace')
  do i = 1, NQ
    write(un, '(*(es18.9))') q(i), lam(i, :)
  end do
  close(un)

  do k = 1, 5
    call simulate(cop(k), 1000 + k)                  ! = ex04
    do i = 1, NE
      emp(i, 2 * k - 1) = count(u1 < QE(i) .and. u2 < QE(i)) / (n * QE(i))
      emp(i, 2 * k) = count(u1 > 1 - QE(i) .and. u2 > 1 - QE(i)) / (n * QE(i))
    end do
  end do
  open(newunit=un, file='tails_emp.txt', status='replace')
  do i = 1, NE
    write(un, '(*(es18.9))') QE(i), emp(i, :)
  end do
  close(un)
  print '(a,2f8.4)', 'gaussian lambda_U at 1e-2, 1e-4:', lam(minloc(abs(q - 1e-2_dp), 1), 12), &
    lam(minloc(abs(q - 1e-4_dp), 1), 12)

contains
  subroutine simulate(c, seed)
    type(fv_bicop), intent(in) :: c
    integer, intent(in) :: seed
    integer :: ns, j
    integer, allocatable :: s(:)
    call random_seed(size=ns)
    s = [(seed + 7919 * j, j = 1, ns)]
    call random_seed(put=s)
    call random_number(u1); call random_number(w)
    call c%hinv1(u1, w, u2)
  end subroutine simulate
end program example_tails
