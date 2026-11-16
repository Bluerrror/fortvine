! ex03_generic_families -- the generic engine (codes 10000+): list the families, use an Archimax
! copula and an asymmetric Khoudraji copula, and select among classic + generic families.
program ex03_generic_families
  use fortvine_api
  implicit none
  integer :: ns
  integer, allocatable :: seed(:)
  integer, parameter :: dp = kind(1.0d0), n = 3000
  type(fv_bicop) :: truth, fitted, chosen
  type(fv_controls) :: ctl
  integer, allocatable :: codes(:)
  real(dp) :: u1(n), w(n), u2(n), lo(3), hi(3)
  integer :: k

  codes = fv_generic_families()                      ! the 85 new families (duplicates excluded)
  print '(a,i0,a)', 'generic families: ', size(codes), '; the first five:'
  do k = 1, 5
    call fv_par_bounds(codes(k), lo, hi)
    print '(i7,2x,a40,a,i0,a)', codes(k), fv_fam_name(codes(k)), '  (', gen_npar(codes(k)), ' parameters)'
  end do

  ! Archimax(gamma | Husler-Reiss) = code 10203: gamma (Clayton) generator, Husler-Reiss Pickands
  truth%family = 10203; truth%par = [1.0_dp, 1.5_dp, 0.0_dp]
  call random_seed(size=ns); allocate(seed(ns)); seed = 20261116; call random_seed(put=seed)   ! reproducible sample
  call random_number(u1); call random_number(w)
  call truth%hinv1(u1, w, u2)
  fitted%family = 10203
  call fitted%fit(u1, u2)
  print '(a,a,a,2f8.3)', 'fitted ', trim(fv_fam_name(fitted%family)), ':', fitted%par(1:2)

  ! Khoudraji-Clayton = code 20200: asymmetric, C(u,v) = u^(1-a) v^(1-b) Clayton(u^a, v^b)
  truth%family = 20200; truth%par = [2.0_dp, 0.6_dp, 0.95_dp]
  call truth%hinv1(u1, w, u2)
  fitted%family = 20200
  call fitted%fit(u1, u2)
  print '(a,3f8.3)', 'fitted Khoudraji-Clayton (theta, a, b):', fitted%par

  ! selection among every parametric family + every generic family
  ctl%nfam = 24 + size(codes); ctl%family_set = 0
  ctl%family_set(1:24) = FV_ALL; ctl%family_set(25:ctl%nfam) = codes
  call chosen%select(u1, u2, ctl)
  print '(a,a,a,i0)', 'selected among all families: ', trim(fv_fam_name(chosen%family)), '  rotation ', chosen%rotation
end program ex03_generic_families
