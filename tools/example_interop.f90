! example_interop -- data for the figure of docs/examples/interop.md (tools/example_figures.py interop).
! Evaluates the density, h-function and inverse h-function of twelve pair copulas at 500 fixed points
! through the Fortran API (fv_bicop) and writes, into the current directory:
!   interop_specs.txt    fam rot p1 p2 p3 of the twelve copulas
!   interop_points.txt   u v of the 500 points (an R2 low-discrepancy sequence, 17 digits)
!   interop_fortran.txt  per point: pdf, h1, hinv1 of copula 1, then of copula 2, ... (36 columns)
! The same specs and points are then evaluated from C (tools/example_interop.c) and from Python
! through ctypes, so the three sets of numbers can be compared.
program example_interop
  use fortvine_api
  implicit none
  integer, parameter :: dp = kind(1.0d0), m = 500, nc = 12
  integer :: fam(nc), rot(nc), k, i, un
  real(dp) :: par(3, nc), u(m), v(m), out(m, 3 * nc)
  type(fv_bicop) :: c
  fam = [1, 2, 3, 4, 5, 6, 7, 3, 11, 16, 10203, 20300]
  rot = [0, 0, 0, 0, 0, 0, 0, 90, 0, 0, 0, 0]
  par = reshape([0.5_dp, 0.0_dp, 0.0_dp,  0.5_dp, 4.0_dp, 0.0_dp,  2.0_dp, 0.0_dp, 0.0_dp, &
                 2.0_dp, 0.0_dp, 0.0_dp,  5.0_dp, 0.0_dp, 0.0_dp,  2.0_dp, 0.0_dp, 0.0_dp, &
                 1.0_dp, 1.5_dp, 0.0_dp,  2.0_dp, 0.0_dp, 0.0_dp,  0.4_dp, 0.9_dp, 3.0_dp, &
                 1.5_dp, 0.0_dp, 0.0_dp,  1.0_dp, 1.5_dp, 0.0_dp,  3.0_dp, 0.5_dp, 0.9_dp], [3, nc])
  do i = 1, m       ! R2 sequence (Roberts 2018): fractional parts of i * (1/g, 1/g^2), g = 1.3247...
    u(i) = modulo(0.5_dp + i * 0.7548776662466927_dp, 1.0_dp)
    v(i) = modulo(0.5_dp + i * 0.5698402909980532_dp, 1.0_dp)
  end do
  do k = 1, nc
    c%family = fam(k); c%rotation = rot(k); c%par = par(:, k)
    call c%pdf(u, v, out(:, 3 * k - 2))
    call c%hfunc1(u, v, out(:, 3 * k - 1))
    call c%hinv1(u, v, out(:, 3 * k))
  end do
  open(newunit=un, file='interop_specs.txt', status='replace')
  do k = 1, nc
    write(un, '(i6,i5,3es24.16)') fam(k), rot(k), par(:, k)
  end do
  close(un)
  open(newunit=un, file='interop_points.txt', status='replace')
  do i = 1, m
    write(un, '(2es25.17)') u(i), v(i)
  end do
  close(un)
  open(newunit=un, file='interop_fortran.txt', status='replace')
  do i = 1, m
    write(un, '(*(es25.17))') out(i, :)
  end do
  close(un)
end program example_interop
