! ex05_weights_discrete -- (1) observation weights: a time-decay weighted fit that follows a
! change in dependence, and integer weights = duplicated rows; (2) a vine with a discrete
! (Poisson count) variable, passed as u = F(x) and u_minus = F(x-) with discrete = .true.
program ex05_weights_discrete
  use fortvine_api
  implicit none
  integer, parameter :: dp = kind(1.0d0), n = 3000, d = 3
  real(dp), parameter :: LAM = 2.0_dp                ! Poisson mean of the count variable
  type(fv_bicop) :: cop, fit
  type(fv_vinecop) :: truth, vd, vn
  real(dp), allocatable :: u1(:), u2(:), w(:), wt(:), x1(:), x2(:)
  real(dp), allocatable :: z(:, :), u(:, :), um(:, :), umid(:, :), r(:, :), dens(:)
  integer :: i, k, cnt

  ! ---------------------------------------------------------------- (1) weights
  allocate(u1(n), u2(n), w(n), wt(n), x1(2 * n), x2(2 * n))
  call seed_rng(11)
  call random_number(u1); call random_number(w)
  cop%family = FV_CLAYTON
  cop%par(1) = 1.0_dp; call cop%hinv1(u1(1:n/2), w(1:n/2), u2(1:n/2))           ! old regime: tau 1/3
  cop%par(1) = 4.0_dp; call cop%hinv1(u1(n/2+1:), w(n/2+1:), u2(n/2+1:))       ! recent regime: tau 2/3

  fit%family = FV_CLAYTON
  call fit%fit(u1, u2)
  print '(a,f7.3)', 'Clayton fit, equal weights:          theta =', fit%par(1)
  wt = [(0.998_dp**(n - i), i = 1, n)]                ! exponential time decay, half-life ~350 rows
  call fit%fit(u1, u2, weights=wt)
  print '(a,f7.3,a,f7.1,a)', 'Clayton fit, time-decay weights:     theta =', fit%par(1), &
    '   (effective n =', sum(wt)**2 / sum(wt**2), ')'

  ! integer weights are the same as repeating rows
  call random_number(wt)
  wt = merge(2.0_dp, 1.0_dp, wt > 0.5_dp)
  k = 0
  do i = 1, n
    do cnt = 1, nint(wt(i))
      k = k + 1; x1(k) = u1(i); x2(k) = u2(i)
    end do
  end do
  call fit%fit(u1, u2, weights=wt)
  print '(a,f12.8,a,f12.4)', 'weights 1/2:     theta =', fit%par(1), '  loglik =', fit%loglik
  call fit%fit(x1(1:k), x2(1:k))
  print '(a,f12.8,a,f12.4)', 'duplicated rows: theta =', fit%par(1), '  loglik =', fit%loglik

  ! ---------------------------------------------------------------- (2) a discrete variable
  ! D-vine 1-2-3: gaussian(0.7) and clayton(2) in tree 1, frank(2) in tree 2
  allocate(z(n, d), u(n, d), um(n, d), umid(n, d), r(n, d), dens(n))
  call truth%set_dvine([1, 2, 3])
  truth%fam(1:2, 1) = [FV_GAUSSIAN, FV_CLAYTON]
  truth%par(1, 1, 1) = 0.7_dp; truth%par(1, 2, 1) = 2.0_dp
  truth%fam(1, 2) = FV_FRANK; truth%par(1, 1, 2) = 2.0_dp
  call truth%simulate(n, z, seed=5)

  ! variable 2 becomes a Poisson(2) count x = F^-1(z); the copula data are F(x) and F(x-)
  u = z; um = z                                      ! continuous columns: u_minus is not used
  do i = 1, n
    k = qpois(z(i, 2))
    u(i, 2) = ppois(k); um(i, 2) = ppois(k - 1)
  end do
  print '(/,a,i0,a,i0,a,f6.3)', 'count variable: ', count(um(:, 2) == 0.0_dp), ' zeros, maximum ', &
    maxval([(qpois(z(i, 2)), i = 1, n)]), ', mean ', sum([(qpois(z(i, 2)), i = 1, n)]) / real(n, dp)

  ! (a) the discrete-aware fit on the fixed D-vine structure
  call vd%set_dvine([1, 2, 3])
  call vd%select(u, discrete=[.false., .true., .false.], u_minus=um)
  ! (b) a naive fit that treats the count as continuous (mid-point of its CDF jump)
  umid = u; umid(:, 2) = 0.5_dp * (u(:, 2) + um(:, 2))
  call vn%set_dvine([1, 2, 3])
  call vn%select(umid)
  ! select may store the same D-vine in the reverse order, so edges are found by their pair
  print '(a)', 'tree-1 pair      truth             discrete-aware    naive'
  do k = 1, 2
    print '(a,i0,a,i0,a,3(a10,f7.3,1x))', '  (', k, ',', k + 1, ')    ', &
      fv_fam_name(truth%fam(edge(truth, k, k + 1), 1)), truth%par(1, edge(truth, k, k + 1), 1), &
      fv_fam_name(vd%fam(edge(vd, k, k + 1), 1)), vd%par(1, edge(vd, k, k + 1), 1), &
      fv_fam_name(vn%fam(edge(vn, k, k + 1), 1)), vn%par(1, edge(vn, k, k + 1), 1)
  end do

  ! log-likelihood of the mixed data (density in x1, x3 times probability mass in x2)
  truth%disc = [.false., .true., .false.]
  print '(a,f10.2,a,f10.2)', 'loglik of the mixed data: truth', truth%loglik(u, um), &
    '   discrete-aware fit', vd%loglik(u, um)
  call vd%pdf(u(1:3, :), dens(1:3), u_minus=um(1:3, :))
  print '(a,3es11.3)', 'first three mixed densities:', dens(1:3)

  ! randomised Rosenblatt transform (Brockwell 2007): iid uniform when the model is right
  call vd%rosenblatt(u, r, u_minus=um, seed=1)
  print '(a,3f8.4)', 'Rosenblatt column means (1/2):       ', sum(r, 1) / n
  print '(a,3f8.4)', 'Rosenblatt column variances (1/12):  ', &
    [(sum((r(:, k) - sum(r(:, k)) / n)**2) / n, k = 1, d)]
  print '(a,3f8.4)', 'Kendall tau (1,2) (1,3) (2,3), ~0:   ', fv_kendall_tau(r(:, 1), r(:, 2)), &
    fv_kendall_tau(r(:, 1), r(:, 3)), fv_kendall_tau(r(:, 2), r(:, 3))

contains

  ! index of the tree-1 edge joining variables a and b: edge e joins m(d-e+1, e) and m(1, e)
  integer function edge(v, a, b)
    type(fv_vinecop), intent(in) :: v
    integer, intent(in) :: a, b
    integer :: m(v%d, v%d), e
    m = v%matrix()
    do e = 1, v%d - 1
      edge = e
      if (minval([m(v%d - e + 1, e), m(1, e)]) == min(a, b) .and. maxval([m(v%d - e + 1, e), m(1, e)]) == max(a, b)) &
        return
    end do
    edge = 0
  end function edge

  subroutine seed_rng(seed)
    integer, intent(in) :: seed
    integer :: ns, j
    integer, allocatable :: s(:)
    call random_seed(size=ns)
    s = [(seed + 7919 * j, j = 1, ns)]
    call random_seed(put=s)
  end subroutine seed_rng

  ! Poisson(LAM) distribution function F(k) and quantile min{k : F(k) >= p}
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
end program ex05_weights_discrete
