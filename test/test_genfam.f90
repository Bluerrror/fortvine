! Tests of the generic pair-copula engine (src/families_gen.inc):
!  1. codes that duplicate a hand-coded family reproduce its log-density and h-functions
!  2. every valid generic code is self-consistent on a grid: dh1/dv = c, dh2/du = c, dC/du = h1
!     (finite differences, two interior parameter points each; non-finite values at the far
!     corners of strongly dependent members are tolerated, wrong finite values are not).
!     Draws with Blomqvist beta = 4 C(1/2, 1/2) - 1 > 0.9 are skipped: near comonotonicity the
!     density lives within rounding of the diagonal and double precision runs out (the
!     hand-coded BB families too); that range is a known limit, not a tested one
!  3. MLE recovers simulated parameters (log-likelihood at the fit >= at the truth)
program test_genfam
  use fortvine_api
  implicit none
  integer, parameter :: dp = kind(1.0d0), NG = 9
  real(dp) :: ug(NG), u1(NG * NG), u2(NG * NG)
  integer :: i, j, nfail = 0, nskip = 0
  integer, allocatable :: cs(:)
  ug = [0.003_dp, 0.05_dp, 0.17_dp, 0.33_dp, 0.5_dp, 0.66_dp, 0.83_dp, 0.95_dp, 0.997_dp]
  do i = 1, NG
    do j = 1, NG
      u1((i - 1) * NG + j) = ug(i); u2((i - 1) * NG + j) = ug(j)
    end do
  end do
  call same(10100, [2.3_dp, 0.0_dp, 0.0_dp], 4, [2.3_dp, 0.0_dp, 0.0_dp])          ! Gumbel
  call same(10200, [1.7_dp, 0.0_dp, 0.0_dp], 3, [1.7_dp, 0.0_dp, 0.0_dp])          ! Clayton
  call same(10300, [2.6_dp, 0.0_dp, 0.0_dp], 6, [2.6_dp, 0.0_dp, 0.0_dp])          ! Joe
  call same(10400, [6.0_dp, 0.0_dp, 0.0_dp], 5, [6.0_dp, 0.0_dp, 0.0_dp])          ! Frank
  call same(10500, [2.0_dp, 1.5_dp, 0.0_dp], 20, [2.0_dp, 1.5_dp, 0.0_dp])         ! BB9
  call same(10600, [3.0_dp, 0.7_dp, 0.0_dp], 10, [3.0_dp, 0.7_dp, 0.0_dp])         ! BB8
  call same(10700, [1.5_dp, 0.6_dp, 0.0_dp], 21, [1.5_dp, 0.6_dp, 0.0_dp])         ! BB10
  call same(11500, [0.8_dp, 1.6_dp, 0.0_dp], 7, [0.8_dp, 1.6_dp, 0.0_dp])          ! BB1
  call same(11600, [1.2_dp, 0.9_dp, 0.0_dp], 17, [1.2_dp, 0.9_dp, 0.0_dp])         ! BB2
  call same(11200, [1.4_dp, 0.8_dp, 0.0_dp], 18, [1.4_dp, 0.8_dp, 0.0_dp])         ! BB3
  call same(11900, [1.5_dp, 1.4_dp, 0.0_dp], 8, [1.5_dp, 1.4_dp, 0.0_dp])          ! BB6
  call same(12000, [1.6_dp, 1.1_dp, 0.0_dp], 9, [1.6_dp, 1.1_dp, 0.0_dp])          ! BB7
  call same(10002, [1.3_dp, 0.0_dp, 0.0_dp], 15, [1.3_dp, 0.0_dp, 0.0_dp])         ! Galambos
  call same(10005, [0.6_dp, 0.9_dp, 2.5_dp], 11, [0.6_dp, 0.9_dp, 2.5_dp])         ! Tawn
  call same(10202, [1.1_dp, 1.3_dp, 0.0_dp], 22, [1.1_dp, 1.3_dp, 0.0_dp])         ! BB4
  call same(10102, [1.4_dp, 1.2_dp, 0.0_dp], 19, [1.4_dp, 1.2_dp, 0.0_dp])         ! BB5
  call same(20001, [2.5_dp, 0.6_dp, 0.9_dp], 11, [0.6_dp, 0.9_dp, 2.5_dp])         ! Khoudraji-Gumbel = Tawn
  call same(20003, [1.8_dp, 0.6_dp, 0.9_dp], 10012, [1.8_dp, 0.6_dp, 0.9_dp])      ! wrapper = Pickands form
  cs = fv_generic_families(.false.)
  call random_seed()
  do i = 1, size(cs)
    do j = 0, 270, 90
      call consist(cs(i), j)
    end do
  end do
  print '(a,i0,a)', '      ', nskip, ' near-(counter)monotone parameter draws skipped (|beta| > 0.9)'
  call check(nfail == 0, 'every generic family x 4 rotations: dh/du = c and dC/du = h on a 9 x 9 grid')
  call mle(10012, [2.0_dp, 0.5_dp, 0.9_dp])
  call mle(10203, [1.0_dp, 1.5_dp, 0.0_dp])
  call mle(11700, [1.0_dp, 2.0_dp, 0.0_dp])
  call mle(20400, [8.0_dp, 0.6_dp, 0.9_dp])
  call flipchk()
  call check(size(fv_generic_families()) > 80 .and. all(gen_valid(fv_generic_families())), 'generic family list')
  if (nfail == 0) then
    print '(a)', 'ALL GENERIC FAMILY TESTS PASSED'
  else
    print '(i0,a)', nfail, ' GENERIC FAMILY TEST(S) FAILED'
    error stop 1
  end if
contains
  subroutine check(ok, msg)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: msg
    if (ok) then
      print '(2a)', 'ok    ', msg
    else
      print '(2a)', 'FAIL  ', msg
      nfail = nfail + 1
    end if
  end subroutine check

  subroutine same(fg, pg, fh, ph)
    integer, intent(in) :: fg, fh
    real(dp), intent(in) :: pg(3), ph(3)
    type(fv_bicop) :: a, b
    real(dp), dimension(NG * NG) :: la, lb, h1a, h1b, h2a, h2b
    a%family = fg; a%par = pg; b%family = fh; b%par = ph
    call a%logpdf(u1, u2, la); call b%logpdf(u1, u2, lb)
    call a%hfunc1(u1, u2, h1a); call b%hfunc1(u1, u2, h1b)
    call a%hfunc2(u1, u2, h2a); call b%hfunc2(u1, u2, h2b)
    call check(maxval(abs(la - lb)) < 1e-6_dp .and. max(maxval(abs(h1a - h1b)), maxval(abs(h2a - h2b))) < 1e-9_dp, &
               'generic ' // trim(fv_fam_name(fg)) // ' reproduces the hand-coded family')
  end subroutine same

  subroutine consist(f, rot)
    integer, intent(in) :: f, rot
    type(fv_bicop) :: b
    real(dp) :: lo(3), hi(3), r(3), e1, e2, e3
    real(dp), dimension(NG * NG) :: c, h1, hp, hm, cp, cm, lc, du, dv, hh(1)
    real(dp), parameter :: EPS = 1e-6_dp
    integer :: q, np
    logical :: fin(NG * NG)
    np = gen_npar(f)
    call fv_par_bounds(f, lo, hi)
    du = EPS * u1 * (1 - u1); dv = EPS * u2 * (1 - u2)
    do q = 1, 2
      call random_number(r)
      b%family = f; b%rotation = rot; b%par = 0.0_dp
      b%par(1:np) = lo(1:np) + (hi(1:np) - lo(1:np)) * (0.1_dp + 0.5_dp * r(1:np))**merge(2.0_dp, 1.0_dp, hi(1:np) > 2)
      call b%cdf([0.5_dp], [0.5_dp], hh)
      if (abs(4.0_dp * hh(1) - 1.0_dp) > 0.9_dp) then   ! |Blomqvist beta|: rotations 90 / 270 are negative
        nskip = nskip + 1; cycle
      end if
      call b%logpdf(u1, u2, lc); c = exp(lc); fin = lc == lc
      call b%hfunc1(u1, u2, h1)
      call b%hfunc1(u1, u2 + dv, hp); call b%hfunc1(u1, u2 - dv, hm)
      e1 = maxval(abs((hp - hm) / (2 * dv) - c) / (1 + c), mask=fin)
      call b%hfunc2(u1 + du, u2, hp); call b%hfunc2(u1 - du, u2, hm)
      e2 = maxval(abs((hp - hm) / (2 * du) - c) / (1 + c), mask=fin)
      call b%cdf(u1 + du, u2, cp); call b%cdf(u1 - du, u2, cm)
      e3 = maxval(abs((cp - cm) / (2 * du) - h1), mask=fin)
      if (max(e1, e2, e3) > 1e-3_dp .or. count(fin) < NG * NG / 2) then
        nfail = nfail + 1
        print '(a,i6,a,i4,1x,a,3f9.4,3es10.2)', 'FAIL  consistency ', f, ' rot', rot, trim(fv_fam_name(f)), b%par, e1, e2, e3
      end if
    end do
  end subroutine consist

  ! the vine edge flip of every asymmetric code is its transpose: C_flip(u, v) = C(v, u)
  subroutine flipchk()
    type(fv_bicop) :: a, b
    real(dp) :: lo(3), hi(3), ca(NG * NG), cb(NG * NG), e
    integer :: q, f, np, r
    e = 0.0_dp
    do q = 1, size(cs)
      f = cs(q)
      if (gen_exch(f)) cycle
      np = gen_npar(f); call fv_par_bounds(f, lo, hi)
      do r = 0, 270, 90                     ! flip of a rotated edge: transpose + rotation 90 <-> 270
        a%family = f; a%rotation = r; a%par = 0.0_dp; a%par(1:np) = lo(1:np) + 0.37_dp * (hi(1:np) - lo(1:np))**0.7_dp
        b = a; call gen_flip(b%family, b%par)
        if (r == 90) b%rotation = 270
        if (r == 270) b%rotation = 90
        call a%logpdf(u2, u1, ca); call b%logpdf(u1, u2, cb)
        e = max(e, maxval(abs(ca - cb), mask=ca == ca .and. abs(ca) < 600))
      end do
    end do
    call check(e < 1e-10_dp, 'edge flip of every asymmetric generic family x 4 rotations = its transpose')
  end subroutine flipchk

  subroutine mle(f, p)
    integer, intent(in) :: f
    real(dp), intent(in) :: p(3)
    integer, parameter :: n = 3000
    type(fv_bicop) :: t, e
    real(dp) :: u(n), w(n), v(n), lp(n)
    t%family = f; t%par = p
    call random_number(u); call random_number(w)
    call t%hinv1(u, w, v); call t%logpdf(u, v, lp)
    e%family = f; call e%fit(u, v)
    call check(e%loglik >= sum(lp) - 1e-6_dp * n, 'MLE of ' // trim(fv_fam_name(f)) // ' reaches the truth''s likelihood')
  end subroutine mle
end program test_genfam
