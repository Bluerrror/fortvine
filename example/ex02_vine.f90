! ex02_vine -- a 5-dimensional vine copula: simulate from a known D-vine, select a new vine from
! the sample (structure, families, parameters), then density, Rosenblatt transform and simulation.
program ex02_vine
  use fortvine_api
  implicit none
  integer, parameter :: dp = kind(1.0d0), n = 3000, d = 5
  type(fv_vinecop) :: truth, vc
  type(fv_controls) :: ctl
  real(dp) :: u(n, d), w(n, d), dens(n), sim(10, d)
  integer :: e, t

  ! a D-vine 1-2-3-4-5 with Gumbel pairs in tree 1 and Frank pairs in tree 2, independence above
  call truth%set_dvine([1, 2, 3, 4, 5])
  do e = 1, d - 1
    truth%fam(e, 1) = FV_GUMBEL; truth%par(1, e, 1) = 2.0_dp
  end do
  do e = 1, d - 2
    truth%fam(e, 2) = FV_FRANK; truth%par(1, e, 2) = 3.0_dp
  end do
  call truth%simulate(n, u, seed=42)

  ! select a regular vine: Dissmann's algorithm, all parametric families, BIC
  ctl%nthreads = 4
  call vc%select(u, ctl)
  call vc%print()
  print '(a,f12.2,a,i0)', 'log-likelihood:', vc%loglik(u), '   parameters:', vc%npars()

  ! density, Rosenblatt transform (independent uniforms if the model is right) and new samples
  call vc%pdf(u, dens)
  call vc%rosenblatt(u, w)
  print '(a,f8.4)', 'Kendall tau of the first two Rosenblatt columns (should be ~0):', &
    fv_kendall_tau(w(:, 1), w(:, 2))
  call vc%simulate(10, sim, seed=7)
  print '(a,5f8.4)', 'first simulated row:', sim(1, :)
  do t = 1, 2
    print '(a,i0,a,*(1x,a))', 'tree ', t, ' families:', (trim(fv_fam_name(vc%fam(e, t))), e = 1, d - t)
  end do
end program ex02_vine
