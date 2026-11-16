! ex04_tails -- tail dependence of five pair copulas, three ways: the limit from the family
! formula, the exact coefficient at a finite level from the copula CDF, and the empirical
! coefficient of a simulated sample.  Then: which family does BIC select from each sample?
!
!   lower tail  lambda_L(q) = C(q, q) / q                      -> lambda_L as q -> 0
!   upper tail  lambda_U(q) = (1 - 2 q + C(q, q)) / (1 - q)    -> lambda_U as q -> 1
program ex04_tails
  use fortvine_api
  implicit none
  integer, parameter :: dp = kind(1.0d0), n = 200000, m = 2000, nc = 5
  type(fv_bicop) :: cop(nc), chosen
  type(fv_controls) :: ctl
  real(dp), allocatable :: u1(:), w(:), u2(:)   ! large samples: heap, not stack (-fopenmp)
  real(dp) :: lim(2, nc), th, a, b
  integer, allocatable :: gen(:)
  integer :: k

  allocate(u1(n), w(n), u2(n))
  ! the five copulas (par = family parameters, unused slots 0)
  cop(1)%family = FV_GUMBEL;  cop(1)%par = [2.0_dp, 0.0_dp, 0.0_dp]
  cop(2)%family = FV_CLAYTON; cop(2)%par = [2.0_dp, 0.0_dp, 0.0_dp]
  cop(3)%family = FV_CLAYTON; cop(3)%par = [2.0_dp, 0.0_dp, 0.0_dp]; cop(3)%rotation = 180
  cop(4)%family = FV_TAWN;    cop(4)%par = [0.4_dp, 0.9_dp, 3.0_dp]     ! psi1, psi2, theta
  cop(5)%family = 20300;      cop(5)%par = [3.0_dp, 0.5_dp, 0.9_dp]     ! Khoudraji-Joe: theta, a, b

  ! limits (lambda_L, lambda_U) from the family formulas
  th = 2.0_dp
  lim(:, 1) = [0.0_dp, 2.0_dp - 2.0_dp**(1.0_dp / th)]                 ! Gumbel
  lim(:, 2) = [2.0_dp**(-1.0_dp / th), 0.0_dp]                         ! Clayton
  lim(:, 3) = [0.0_dp, 2.0_dp**(-1.0_dp / th)]                         ! survival Clayton
  a = 0.4_dp; b = 0.9_dp; th = 3.0_dp
  lim(:, 4) = [0.0_dp, a + b - (a**th + b**th)**(1.0_dp / th)]         ! Tawn
  a = 0.5_dp; b = 0.9_dp; th = 3.0_dp
  lim(:, 5) = [0.0_dp, a + b - (a**th + b**th)**(1.0_dp / th)]         ! Khoudraji-Joe

  print '(a)', '                        |      lower tail lambda_L(q)         |      upper tail lambda_U(q)'
  print '(a)', 'copula              rot |      limit  q=1e-2  q=1e-4  sample |      limit  q=1e-2  q=1e-4  sample'
  do k = 1, nc
    call simulate(cop(k), 1000 + k)
    print '(a20,i3,a,4f8.3,a,4f8.3)', fv_fam_name(cop(k)%family), cop(k)%rotation, ' |   ', &
      lim(1, k), lam_lower(cop(k), 1.0e-2_dp), lam_lower(cop(k), 1.0e-4_dp), &
      count(u1 < 1.0e-2_dp .and. u2 < 1.0e-2_dp) / (n * 1.0e-2_dp), ' |   ', &
      lim(2, k), lam_upper(cop(k), 1.0e-2_dp), lam_upper(cop(k), 1.0e-4_dp), &
      count(u1 > 1 - 1.0e-2_dp .and. u2 > 1 - 1.0e-2_dp) / (n * 1.0e-2_dp)
  end do

  ! BIC selection from the first m rows of each sample: default family set (hand-coded
  ! parametric families with rotations), then the same plus every generic family
  gen = fv_generic_families()
  print '(/,a)', 'BIC selection from 2000 observations:'
  do k = 1, nc
    call simulate(cop(k), 1000 + k)
    call chosen%select(u1(1:m), u2(1:m))
    write(*, '(2x,a14,a,a10,a,i3)', advance='no') fv_fam_name(cop(k)%family), ' default set: ', &
      fv_fam_name(chosen%family), ' rot', chosen%rotation
    ctl%nfam = 24 + size(gen); ctl%family_set = 0
    ctl%family_set(1:24) = FV_ALL; ctl%family_set(25:ctl%nfam) = gen
    call chosen%select(u1(1:m), u2(1:m), ctl)
    print '(a,a14,a,i3)', '   + generic families: ', fv_fam_name(chosen%family), ' rot', chosen%rotation
  end do

contains

  ! n draws of (u1, u2) from c by the inverse Rosenblatt transform, reproducible from seed
  subroutine simulate(c, seed)
    type(fv_bicop), intent(in) :: c
    integer, intent(in) :: seed
    integer :: ns, j
    integer, allocatable :: s(:)
    call random_seed(size=ns)
    s = [(seed + 7919 * j, j = 1, ns)]
    call random_seed(put=s)
    call random_number(u1); call random_number(w)
    call c%hinv1(u1, w, u2)               ! u2 = h1^-1(w | u1)
  end subroutine simulate

  real(dp) function lam_lower(c, q)
    type(fv_bicop), intent(in) :: c
    real(dp), intent(in) :: q
    real(dp) :: cq(1)
    call c%cdf([q], [q], cq)
    lam_lower = cq(1) / q
  end function lam_lower

  real(dp) function lam_upper(c, q)       ! q is the distance from 1
    type(fv_bicop), intent(in) :: c
    real(dp), intent(in) :: q
    real(dp) :: cq(1)
    call c%cdf([1 - q], [1 - q], cq)
    lam_upper = (2 * q - 1 + cq(1)) / q
  end function lam_upper
end program ex04_tails
