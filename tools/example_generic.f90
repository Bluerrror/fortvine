! example_generic -- data for the figures of docs/examples/generic.md (tools/example_figures.py generic).
! Re-runs the two samples of example/ex03_generic_families.f90 (same seed), refits them and writes,
! into the current directory:
!   gen_grid.txt      4 blocks of 121 x 121 densities on the normal-score grid: Archimax 10203 true,
!                     fitted; Khoudraji-Clayton 20200 true, fitted
!   gen_pickands.txt  w  A(w) for 8 curves (columns): the inner Pickands function of 10203 recovered
!                     from its CDF, A(w) = phi(C(psi(1-w), psi(w))) with the gamma generator, at the
!                     true and the fitted parameters; the Husler-Reiss family 16 (lambda = 1.5) from
!                     A(w) = -log C(exp(-(1-w)), exp(-w)); and the five generic extreme-value codes
!                     listed by ex03 (10004, 10006, 10007, 10008, 10009)
! stdout: fitted parameters and the parameters of the extreme-value curves.
program example_generic
  use fortvine_api
  implicit none
  integer, parameter :: dp = kind(1.0d0), n = 3000, NG = 121, NW = 201
  integer :: ns, i, j, k, un
  integer, allocatable :: seed(:)
  type(fv_bicop) :: truth, fitted, ev(5), hr
  real(dp) :: u1(n), w(n), u2(n), z(NG), gu(NG*NG), gv(NG*NG), gd(NG*NG)
  real(dp) :: ww(NW), a(NW, 8), pu(NW), pv(NW), c(NW), th_arch(2)

  call random_seed(size=ns); allocate(seed(ns)); seed = 20261116; call random_seed(put=seed)   ! = ex03
  call random_number(u1); call random_number(w)

  do i = 1, NG
    z(i) = -3.0_dp + 6.0_dp * (i - 1) / (NG - 1)
  end do
  do j = 1, NG
    do i = 1, NG
      gu(i + (j - 1) * NG) = pn(z(i)); gv(i + (j - 1) * NG) = pn(z(j))
    end do
  end do
  open(newunit=un, file='gen_grid.txt', status='replace')

  ! Archimax(gamma | Husler-Reiss), code 10203
  truth%family = 10203; truth%par = [1.0_dp, 1.5_dp, 0.0_dp]
  call truth%hinv1(u1, w, u2)
  fitted%family = 10203
  call fitted%fit(u1, u2)
  print '(a,3f10.5)', 'archimax_hat', fitted%par
  th_arch = [truth%par(1), fitted%par(1)]
  call grid(truth); call grid(fitted)
  ww = [((i - 1.0_dp) / (NW - 1), i = 1, NW)]
  ww(1) = 1.0e-6_dp; ww(NW) = 1.0_dp - 1.0e-6_dp
  do k = 1, 2
    pu = (2.0_dp - ww)**(-1.0_dp / th_arch(k)); pv = (1.0_dp + ww)**(-1.0_dp / th_arch(k))
    if (k == 1) call truth%cdf(pu, pv, c)
    if (k == 2) call fitted%cdf(pu, pv, c)
    a(:, k) = c**(-th_arch(k)) - 1.0_dp
  end do

  ! Khoudraji-Clayton, code 20200 (same uniforms, as in ex03)
  truth%family = 20200; truth%par = [2.0_dp, 0.6_dp, 0.95_dp]
  call truth%hinv1(u1, w, u2)
  fitted%family = 20200
  call fitted%fit(u1, u2)
  print '(a,3f10.5)', 'khoudraji_hat', fitted%par
  call grid(truth); call grid(fitted)
  close(un)

  ! Pickands functions of extreme-value copulas: A(w) = -log C(exp(-(1-w)), exp(-w))
  pu = exp(-(1.0_dp - ww)); pv = exp(-ww)
  hr%family = FV_HUSLER_REISS; hr%par = [1.5_dp, 0.0_dp, 0.0_dp]
  call hr%cdf(pu, pv, c); a(:, 3) = -log(c)
  ev(1)%family = 10004; ev(1)%par = [0.8_dp, 0.0_dp, 0.0_dp]     ! Tawn-mixed: theta
  ev(2)%family = 10006; ev(2)%par = [3.0_dp, 0.5_dp, 0.0_dp]     ! Tawn-type1: theta, psi
  ev(3)%family = 10007; ev(3)%par = [3.0_dp, 0.5_dp, 0.0_dp]     ! Tawn-type2: theta, psi
  ev(4)%family = 10008; ev(4)%par = [0.6_dp, 0.9_dp, 2.0_dp]     ! asym. neg. logistic: psi1, psi2, theta
  ev(5)%family = 10009; ev(5)%par = [2.0_dp, 0.5_dp, 0.0_dp]     ! its type 1: theta, psi
  do k = 1, 5
    call ev(k)%cdf(pu, pv, c); a(:, 3 + k) = -log(c)
    print '(a,i6,3f8.3)', 'ev', ev(k)%family, ev(k)%par
  end do
  open(newunit=un, file='gen_pickands.txt', status='replace')
  do i = 1, NW
    write(un, '(9es18.9)') ww(i), a(i, :)
  end do
  close(un)

contains
  subroutine grid(cop)
    type(fv_bicop), intent(in) :: cop
    integer :: ii, jj
    call cop%pdf(gu, gv, gd)
    do jj = 1, NG
      write(un, '(*(1x,es12.4e3))') (gd(ii + (jj - 1) * NG), ii = 1, NG)
    end do
  end subroutine grid
  real(dp) function pn(x)
    real(dp), intent(in) :: x
    pn = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function pn
end program example_generic
