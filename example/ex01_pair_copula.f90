! ex01_pair_copula -- one bivariate copula: evaluate, simulate, fit a known family, select a family.
program ex01_pair_copula
  use fortvine_api
  implicit none
  integer :: ns
  integer, allocatable :: seed(:)
  integer, parameter :: dp = kind(1.0d0), n = 2000
  type(fv_bicop) :: truth, fitted, chosen
  real(dp) :: u1(n), w(n), u2(n), dens(3)

  ! a Clayton copula with theta = 2 (Kendall's tau = theta / (theta + 2) = 0.5)
  truth%family = FV_CLAYTON; truth%par = [2.0_dp, 0.0_dp, 0.0_dp]

  ! density at three points
  call truth%pdf([0.1_dp, 0.5_dp, 0.9_dp], [0.1_dp, 0.5_dp, 0.9_dp], dens)
  print '(a,3f10.4)', 'Clayton(2) density on the diagonal:', dens

  ! simulate by inverse Rosenblatt: u2 = h^-1(w | u1) with independent uniforms u1, w
  call random_seed(size=ns); allocate(seed(ns)); seed = 20261116; call random_seed(put=seed)   ! reproducible sample
  call random_number(u1); call random_number(w)
  call truth%hinv1(u1, w, u2)
  print '(a,f7.3)', 'empirical Kendall tau:', fv_kendall_tau(u1, u2)

  ! maximum likelihood fit of a given family
  fitted%family = FV_CLAYTON
  call fitted%fit(u1, u2)
  print '(a,f8.4,a,f10.2)', 'fitted theta:', fitted%par(1), '   log-likelihood:', fitted%loglik

  ! family selection by BIC over the default family set (with rotations)
  call chosen%select(u1, u2)
  print '(a,a,a,i0,a,3f8.3)', 'selected: ', trim(fv_fam_name(chosen%family)), '  rotation ', chosen%rotation, &
    '  par', chosen%par
end program ex01_pair_copula
