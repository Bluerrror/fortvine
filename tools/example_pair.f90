! example_pair -- data for the figures of docs/examples/pair-copula.md (tools/example_figures.py pair).
! Re-runs the sample of example/ex01_pair_copula.f90 (same seed), refits Clayton and writes, into the
! current directory:
!   pair_sample.txt  u1 u2 of the 2000 draws
!   pair_grid.txt    true and fitted density on the 121 x 121 normal-score grid (two blocks)
!   pair_h.txt       u  v  h1(v|u)  c(u,v)   for u = 0.1, 0.5, 0.9 on a v grid
!   pair_hinv.txt    u  w  h1^-1(w|u)
! stdout: fitted theta and the round-trip error max |h1^-1(h1(v|u)|u) - v|.
program example_pair
  use fortvine_api
  implicit none
  integer, parameter :: dp = kind(1.0d0), n = 2000, NG = 121, NV = 401, NW = 19
  integer :: ns, i, j, k, un
  integer, allocatable :: seed(:)
  type(fv_bicop) :: truth, fitted
  real(dp) :: u1(n), w(n), u2(n), z(NG), gu(NG*NG), gv(NG*NG), gd(NG*NG)
  real(dp) :: v(NV), uu(NV), h(NV), back(NV), cd(NV), ww(NW), uw(NW), hi(NW), err
  real(dp), parameter :: UC(3) = [0.1_dp, 0.5_dp, 0.9_dp]

  truth%family = FV_CLAYTON; truth%par = [2.0_dp, 0.0_dp, 0.0_dp]
  call random_seed(size=ns); allocate(seed(ns)); seed = 20261116; call random_seed(put=seed)   ! = ex01
  call random_number(u1); call random_number(w)
  call truth%hinv1(u1, w, u2)
  fitted%family = FV_CLAYTON
  call fitted%fit(u1, u2)
  print '(a,f10.5,a,f10.3)', 'theta_hat', fitted%par(1), '  loglik', fitted%loglik

  open(newunit=un, file='pair_sample.txt', status='replace')
  do i = 1, n
    write(un, '(2es24.16)') u1(i), u2(i)
  end do
  close(un)

  do i = 1, NG
    z(i) = -3.0_dp + 6.0_dp * (i - 1) / (NG - 1)
  end do
  do j = 1, NG
    do i = 1, NG
      gu(i + (j - 1) * NG) = pn(z(i)); gv(i + (j - 1) * NG) = pn(z(j))
    end do
  end do
  open(newunit=un, file='pair_grid.txt', status='replace')
  call truth%pdf(gu, gv, gd)
  do j = 1, NG
    write(un, '(*(1x,es12.4e3))') (gd(i + (j - 1) * NG), i = 1, NG)
  end do
  call fitted%pdf(gu, gv, gd)
  do j = 1, NG
    write(un, '(*(1x,es12.4e3))') (gd(i + (j - 1) * NG), i = 1, NG)
  end do
  close(un)

  v = [((i - 0.5_dp) / NV, i = 1, NV)]
  ww = [((i * 0.05_dp), i = 1, NW)]
  err = 0.0_dp
  open(newunit=un, file='pair_h.txt', status='replace')
  open(newunit=j, file='pair_hinv.txt', status='replace')
  do k = 1, 3
    uu = UC(k)
    call truth%hfunc1(uu, v, h)
    call truth%hinv1(uu, h, back)
    call truth%pdf(uu, v, cd)
    err = max(err, maxval(abs(back - v)))
    do i = 1, NV
      write(un, '(4es24.16)') UC(k), v(i), h(i), cd(i)
    end do
    uw = UC(k)
    call truth%hinv1(uw, ww, hi)
    do i = 1, NW
      write(j, '(3es24.16)') UC(k), ww(i), hi(i)
    end do
  end do
  close(un); close(j)
  print '(a,es10.2)', 'roundtrip', err
contains
  real(dp) function pn(x)
    real(dp), intent(in) :: x
    pn = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function pn
end program example_pair
