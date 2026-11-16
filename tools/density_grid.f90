! density_grid -- density grids, Kendall tau (1 - 4 int h1 h2) and tail proxies for the documentation figures
! stdin lines: fam rot p1 p2 p3 k lo hi target outfile   (k = 0: no tuning)
program density_grid
  use fortvine_api
  implicit none
  integer, parameter :: dp = kind(1.0d0), NG = 121
  type(fv_bicop) :: cop
  integer :: fam, rot, k, ios, it, i, j
  real(dp) :: p(3), lo, hi, tgt, tl, th, tm, tau, z(NG), uu(NG*NG), vv(NG*NG), dd(NG*NG)
  real(dp) :: lamL, lamU, t1(1), t2(1), c1(1)
  character(len=256) :: outf
  do
    read(*, *, iostat=ios) fam, rot, p, k, lo, hi, tgt, outf
    if (ios /= 0) exit
    cop%family = fam; cop%rotation = rot; cop%par = p
    if (k > 0) then
      cop%par(k) = lo; tl = ktau(cop, 160)
      cop%par(k) = hi; th = ktau(cop, 160)
      do it = 1, 50
        cop%par(k) = 0.5_dp * (lo + hi); tm = ktau(cop, 160)
        if ((tm - tgt) * (th - tl) < 0.0_dp) then
          lo = cop%par(k); tl = tm
        else
          hi = cop%par(k); th = tm
        end if
        if (hi - lo < 1.0e-6_dp * max(1.0_dp, abs(hi))) exit
      end do
    end if
    tau = ktau(cop, 500)
    t1 = 1.0e-5_dp; call cop%cdf(t1, t1, c1); lamL = c1(1) / t1(1)
    t2 = 1.0_dp - 1.0e-5_dp; call cop%cdf(t2, t2, c1); lamU = (1.0_dp - 2.0_dp * t2(1) + c1(1)) / (1.0_dp - t2(1))
    write(*, '(i6,1x,i4,3(1x,es14.6),3(1x,f9.5))') fam, rot, cop%par, tau, lamL, lamU
    if (trim(outf) /= '-') then
      do i = 1, NG
        z(i) = -3.0_dp + 6.0_dp * (i - 1) / (NG - 1)
      end do
      do j = 1, NG
        do i = 1, NG
          uu(i + (j - 1) * NG) = pn(z(i)); vv(i + (j - 1) * NG) = pn(z(j))
        end do
      end do
      call cop%pdf(uu, vv, dd)
      open(10, file=trim(outf), status='replace')
      write(10, '(a,i6,a,i4,a,3es14.6,a,f9.5)') '# fam ', fam, ' rot ', rot, ' par ', cop%par, ' tau ', tau
      do j = 1, NG
        write(10, '(*(1x,es12.4e3))') (dd(i + (j - 1) * NG), i = 1, NG)
      end do
      close(10)
    end if
  end do
contains
  real(dp) function pn(x)
    real(dp), intent(in) :: x
    pn = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function pn
  real(dp) function ktau(c, n)
    type(fv_bicop), intent(in) :: c
    integer, intent(in) :: n
    real(dp) :: a(n*n), b(n*n), h1(n*n), h2(n*n)
    integer :: ii, jj
    do jj = 1, n
      do ii = 1, n
        a(ii + (jj - 1) * n) = (ii - 0.5_dp) / n; b(ii + (jj - 1) * n) = (jj - 0.5_dp) / n
      end do
    end do
    call c%hfunc1(a, b, h1); call c%hfunc2(a, b, h2)
    ktau = 1.0_dp - 4.0_dp * sum(h1 * h2) / real(n, dp)**2
  end function ktau
end program density_grid
