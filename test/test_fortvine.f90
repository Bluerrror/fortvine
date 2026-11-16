! Property tests for fortvine (no reference library needed):
!  1. hinv inverts hfunc for every family x rotation
!  2. h1 equals the integral of the density (h1 = dC/du1 => dh1/du2 = c)
!  3. MLE recovers parameters from simulated data
!  4. vine: inverse_rosenblatt(rosenblatt(u)) = u, pdf_grid = pdf, select beats truth-free baseline
!  5. Kendall tau equals the O(n^2) definition (with ties)
!  6. nonparametric families (tll, bernstein): hinv, h = integral of c, fit quality
!  7. weights: integer weights = duplicated rows; unit weights = unweighted
!  8. discrete margins: pair densities equal brute-force integrals of c; select
!  9. mBICv truncation; xi / mutual information / Ferreira tree criteria
! 10. copula CDFs (bivariate normal, tll, bernstein) = integrals of the density; discrete
!     Rosenblatt (exact 2-d case, randomised PIT uniform, F >= F-); trivariate discrete
!     density = exact rectangle probability; tll / bernstein on discrete pairs;
!     select_threshold; fit cache; t-EV MLE optimality
! 11. non-simplified edges (relax): the CCC test rejects a varying conditional correlation and
!     not a constant one; the relaxed vine predicts held-out data better
program test_fortvine
  use fortvine_api
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer :: nfail = 0

  call test_hinv_and_density()
  call test_fit_recovery()
  call test_vine()
  call test_tau()
  call test_nonpar()
  call test_weights()
  call test_discrete()
  call test_trunc_and_criteria()
  call test_core_extensions()
  call test_nonsimplified()
  if (nfail == 0) then
    print '(a)', 'ALL TESTS PASSED'
  else
    print '(a,i0,a)', 'FAILED: ', nfail, ' check(s)'
    error stop 1
  end if

contains

  subroutine check(ok, msg)
    logical, intent(in) :: ok
    character(*), intent(in) :: msg
    if (.not. ok) then
      nfail = nfail + 1
      print '(a,a)', 'FAIL: ', msg
    end if
  end subroutine check

  subroutine family_case(k, b)
    integer, intent(in) :: k
    type(fv_bicop), intent(out) :: b
    select case (k)
    case (1); b%family = FV_GAUSSIAN; b%par = [0.6_dp, 0.0_dp, 0.0_dp]
    case (2); b%family = FV_STUDENT; b%par = [0.5_dp, 4.0_dp, 0.0_dp]
    case (3); b%family = FV_CLAYTON; b%par = [2.5_dp, 0.0_dp, 0.0_dp]
    case (4); b%family = FV_GUMBEL; b%par = [2.2_dp, 0.0_dp, 0.0_dp]
    case (5); b%family = FV_FRANK; b%par = [6.0_dp, 0.0_dp, 0.0_dp]
    case (6); b%family = FV_JOE; b%par = [2.8_dp, 0.0_dp, 0.0_dp]
    case (7); b%family = FV_BB1; b%par = [0.8_dp, 1.6_dp, 0.0_dp]
    case (8); b%family = FV_BB6; b%par = [1.5_dp, 1.8_dp, 0.0_dp]
    case (9); b%family = FV_BB7; b%par = [1.6_dp, 1.2_dp, 0.0_dp]
    case (10); b%family = FV_BB8; b%par = [3.0_dp, 0.7_dp, 0.0_dp]
    case (11); b%family = FV_TAWN; b%par = [0.4_dp, 0.9_dp, 3.0_dp]
    case (12); b%family = FV_PLACKETT; b%par = [6.0_dp, 0.0_dp, 0.0_dp]
    case (13); b%family = FV_FGM; b%par = [0.7_dp, 0.0_dp, 0.0_dp]
    case (14); b%family = FV_AMH; b%par = [0.8_dp, 0.0_dp, 0.0_dp]
    case (15); b%family = FV_GALAMBOS; b%par = [1.4_dp, 0.0_dp, 0.0_dp]
    case (16); b%family = FV_HUSLER_REISS; b%par = [1.8_dp, 0.0_dp, 0.0_dp]
    case (17); b%family = FV_BB2; b%par = [0.8_dp, 1.2_dp, 0.0_dp]
    case (18); b%family = FV_BB3; b%par = [1.6_dp, 0.8_dp, 0.0_dp]
    case (19); b%family = FV_BB5; b%par = [1.8_dp, 0.9_dp, 0.0_dp]
    case (20); b%family = FV_BB9; b%par = [2.0_dp, 0.7_dp, 0.0_dp]
    case (21); b%family = FV_BB10; b%par = [2.0_dp, 0.8_dp, 0.0_dp]
    case (22); b%family = FV_BB4; b%par = [1.2_dp, 1.5_dp, 0.0_dp]
    case (23); b%family = FV_TEV; b%par = [0.5_dp, 4.0_dp, 0.0_dp]
    end select
  end subroutine family_case

  subroutine test_hinv_and_density()
    integer, parameter :: n = 400, m = 4000
    type(fv_bicop) :: b
    real(dp) :: u1(n), u2(n), h(n), back(n), hres(n), g(m), c(m), uu(m), hh(1), err, interr
    integer :: k, r, i, rots(4)
    character(len=200) :: msg
    rots = [0, 90, 180, 270]
    call random_seed()
    call random_number(u1); call random_number(u2)
    u1 = 0.001_dp + 0.998_dp * u1; u2 = 0.001_dp + 0.998_dp * u2
    do k = 1, 23
      do r = 1, 4
        call family_case(k, b)
        if ((b%family == FV_GAUSSIAN .or. b%family == FV_STUDENT .or. b%family == FV_FRANK .or. &
             b%family == FV_PLACKETT .or. b%family == FV_FGM .or. b%family == FV_AMH) .and. r > 1) cycle
        b%rotation = rots(r)
        call b%hfunc1(u1, u2, h); call b%hinv1(u1, h, back)
        ! residual |h(hinv(h)) - h|: the forward error |back - u2| is ill-conditioned
        ! wherever the density is ~0 (strong tail dependence), for any method
        call b%hfunc1(u1, back, hres)
        err = maxval(abs(hres - h))
        call b%hfunc2(u1, u2, h); call b%hinv2(h, u2, back)
        call b%hfunc2(back, u2, hres)
        err = max(err, maxval(abs(hres - h)))
        write(msg, '(a,a,a,i0,a,es9.2)') 'hinv o hfunc ', trim(FV_FAMILY_NAME(b%family)), ' rot ', rots(r), ' resid ', err
        call check(err < 2.0e-10_dp, msg)     ! inputs are clamped to [1e-10, 1-1e-10]
        ! h1(0.37, v) = integral_0^v c(0.37, t) dt  (midpoint rule)
        do i = 1, m
          g(i) = (i - 0.5_dp) / m * 0.8_dp
        end do
        uu = 0.37_dp
        call b%pdf(uu, g, c)
        interr = 0.0_dp
        call b%hfunc1([0.37_dp], [0.8_dp], hh)
        interr = abs(sum(c) * 0.8_dp / m - (hh(1) - hfunc_at0(b)))
        write(msg, '(a,a,a,i0)') 'h1 = integral of c ', trim(FV_FAMILY_NAME(b%family)), ' rot ', rots(r)
        call check(interr < 2.0e-4_dp, msg)
        block
          real(dp) :: hp(1), hm(1), cc(1), fd
          call b%hfunc2([0.37_dp + 1.0e-5_dp], [0.61_dp], hp)
          call b%hfunc2([0.37_dp - 1.0e-5_dp], [0.61_dp], hm)
          call b%pdf([0.37_dp], [0.61_dp], cc)
          fd = (hp(1) - hm(1)) / 2.0e-5_dp
          write(msg, '(a,a,a,i0)') 'dh2/du = c ', trim(FV_FAMILY_NAME(b%family)), ' rot ', rots(r)
          call check(abs(fd - cc(1)) < 1.0e-5_dp * max(1.0_dp, cc(1)), msg)
        end block
        block
          real(dp) :: h0(1), h1v(1)
          call b%hfunc1([0.37_dp], [1.0_dp - 1.0e-12_dp], h1v)
          call b%hfunc1([0.37_dp], [1.0e-12_dp], h0)
          write(msg, '(a,a,a,i0)') 'h1 limits 0 and 1 ', trim(FV_FAMILY_NAME(b%family)), ' rot ', rots(r)
          call check(h1v(1) > 1.0_dp - 1.0e-6_dp .and. h0(1) < 1.0e-6_dp, msg)
        end block
      end do
    end do
  end subroutine test_hinv_and_density

  real(dp) function hfunc_at0(b)
    type(fv_bicop), intent(in) :: b
    real(dp) :: h(1)
    call b%hfunc1([0.37_dp], [1.0e-12_dp], h)
    hfunc_at0 = h(1)
  end function hfunc_at0

  subroutine test_fit_recovery()
    integer, parameter :: n = 20000
    type(fv_bicop) :: truth, est
    real(dp) :: w1(n), w2(n), u2(n), lpv(n)
    integer :: k, j
    character(len=200) :: msg
    do k = 1, 23
      call family_case(k, truth)
      call random_number(w1); call random_number(w2)
      call truth%hinv1(w1, w2, u2)
      est%family = truth%family; est%rotation = 0
      call est%fit(w1, u2)
      call truth%logpdf(w1, u2, lpv)
      ! MLE property: the optimum is at least as likely as the true parameters
      write(msg, '(a,a,a,2f14.4)') 'MLE optimum ', trim(FV_FAMILY_NAME(truth%family)), ' ll est/truth', est%loglik, sum(lpv)
      if (truth%family == FV_TEV) then
        ! t-EV is fitted derivative-free (Nelder-Mead): optimum to 0.05 log-lik units
        call check(est%loglik >= sum(lpv) - 0.05_dp, msg)
      else
        call check(est%loglik >= sum(lpv) - 1.0e-6_dp * abs(sum(lpv)), msg)
      end if
      ! BB5, BB9, t-EV: parameters are weakly identified (near-flat ridge), so only the
      ! optimum check above applies to them
      if (truth%family == FV_BB5 .or. truth%family == FV_BB9 .or. truth%family == FV_TEV) cycle
      do j = 1, truth%npars()
        write(msg, '(a,a,a,i0,a,2f9.4)') 'MLE recovery ', trim(FV_FAMILY_NAME(truth%family)), ' par ', j, ' est/true', &
          est%par(j), truth%par(j)
        call check(abs(est%par(j) - truth%par(j)) < 0.3_dp * max(1.0_dp, abs(truth%par(j))), msg)
      end do
    end do
  end subroutine test_fit_recovery

  subroutine test_vine()
    integer, parameter :: n = 3000, d = 6, bins = 50
    type(fv_vinecop) :: truth, vc
    type(fv_controls) :: ctl
    real(dp) :: u(n, d), w(n, d), back(n, d), p(n), pg(bins, 20), flat(20 * bins, d), pf(20 * bins)
    real(dp) :: ug(d, bins, 20)
    integer :: t, e, i, k
    call truth%set_dvine([3, 1, 5, 2, 6, 4])
    do t = 1, d - 1
      do e = 1, d - t
        truth%fam(e, t) = merge(FV_CLAYTON, FV_GUMBEL, mod(e + t, 2) == 0)
        truth%par(1, e, t) = merge(2.0_dp, 1.6_dp, mod(e + t, 2) == 0) / t
        if (truth%fam(e, t) == FV_GUMBEL) truth%par(1, e, t) = max(1.0_dp + 0.6_dp / t, 1.0_dp)
      end do
    end do
    truth%fam(1, 1) = FV_TAWN; truth%rot(1, 1) = 90; truth%par(:, 1, 1) = [0.5_dp, 0.9_dp, 3.0_dp]
    truth%fam(2, 1) = FV_STUDENT; truth%par(:, 2, 1) = [0.5_dp, 5.0_dp, 0.0_dp]
    call truth%simulate(n, u, seed=42)
    call truth%rosenblatt(u, w)
    call truth%inverse_rosenblatt(w, back)
    call check(maxval(abs(back - u)) < 1.0e-7_dp, 'inverse_rosenblatt(rosenblatt(u)) = u')
    ! families + parameters on the true structure must reach the truth's likelihood
    call vc%set_dvine([3, 1, 5, 2, 6, 4])
    call vc%select(u, ctl)
    print '(a,2f12.2,2i5)', '  fixed-structure loglik selected / truth, npars:', vc%loglik(u), truth%loglik(u), &
      vc%npars(), truth%npars()
    if (vc%loglik(u) < truth%loglik(u) - 2.0_dp * truth%npars()) call vc%print()
    call check(vc%loglik(u) >= truth%loglik(u) - 2.0_dp * truth%npars(), &
               'family/parameter selection on the true structure reaches the true loglik')
    call vc%init(d)
    call vc%select(u, ctl)                      ! full Dissmann selection (structure too)
    ! pdf_grid vs pdf on the flattened grid
    do i = 1, 20
      do k = 1, bins
        ug(:, k, i) = u(i, :)
        ug(4, k, i) = (k - 0.5_dp) / bins
        flat((i - 1) * bins + k, :) = ug(:, k, i)
      end do
    end do
    call vc%pdf_grid(ug, 4, pg)
    call vc%pdf(flat, pf)
    call check(maxval(abs(reshape(pg, [20 * bins]) - pf) / pf) < 1.0e-12_dp, 'pdf_grid = pdf')
    call vc%pdf(u, p)
    call check(all(p > 0.0_dp), 'pdf positive')
  end subroutine test_vine

  subroutine test_tau()
    integer, parameter :: n = 700
    real(dp) :: x(n), y(n), num, tx, ty, s
    integer :: i, j
    call random_number(x); call random_number(y)
    x = real(nint(10 * x), dp); y = real(nint(10 * (0.5_dp * x / 10 + 0.5_dp * y)), dp)  ! many ties
    num = 0; tx = 0; ty = 0
    do i = 1, n
      do j = i + 1, n
        s = sign(1.0_dp, x(i) - x(j)) * merge(0.0_dp, 1.0_dp, x(i) == x(j)) &
            * sign(1.0_dp, y(i) - y(j)) * merge(0.0_dp, 1.0_dp, y(i) == y(j))
        num = num + s
        if (x(i) /= x(j)) tx = tx + 1
        if (y(i) /= y(j)) ty = ty + 1
      end do
    end do
    call check(abs(fv_kendall_tau(x, y) - num / sqrt(tx * ty)) < 1.0e-12_dp, 'Kendall tau-b with ties')
  end subroutine test_tau

  ! ---------------------------------------------------------------- tll / bernstein
  subroutine test_nonpar()
    integer, parameter :: n = 2000, m = 4000
    type(fv_bicop) :: truth, b
    type(fv_controls) :: ctl
    real(dp) :: w1(n), w2(n), u2(n), h(n), back(n), hres(n), lpv(n), g(m), c(m), uu(m), hh(2), err, interr
    integer :: k, i
    character(len=200) :: msg
    truth%family = FV_GUMBEL; truth%par = [1.8_dp, 0.0_dp, 0.0_dp]
    call random_number(w1); call random_number(w2)
    call truth%hinv1(w1, w2, u2)
    call truth%logpdf(w1, u2, lpv)
    do k = FV_TLL, FV_BERNSTEIN
      b%family = k; b%rotation = 0
      call b%fit(w1, u2)
      write(msg, '(a,a,a,2f10.2,a,f8.2)') 'nonparametric fit ', trim(FV_FAMILY_NAME(k)), ' ll / truth', b%loglik, &
        sum(lpv), ' df', b%npars_eff()
      ! a smoothed estimate pays for its bias, but must capture most of the dependence
      call check(b%loglik > 0.75_dp * sum(lpv) .and. b%npars_eff() > 1.0_dp, msg)
      call b%hfunc1(w1, u2, h); call b%hinv1(w1, h, back); call b%hfunc1(w1, back, hres)
      err = maxval(abs(hres - h))
      call b%hfunc2(w1, u2, h); call b%hinv2(h, u2, back); call b%hfunc2(back, u2, hres)
      err = max(err, maxval(abs(hres - h)))
      write(msg, '(a,a,a,es9.2)') 'hinv o hfunc ', trim(FV_FAMILY_NAME(k)), ' resid ', err
      call check(err < 1.0e-8_dp, msg)
      do i = 1, m
        g(i) = (i - 0.5_dp) / m * 0.8_dp
      end do
      uu = 0.37_dp
      call b%pdf(uu, g, c)
      call b%hfunc1([0.37_dp, 0.37_dp], [0.8_dp, 1.0e-12_dp], hh)
      interr = abs(sum(c) * 0.8_dp / m - (hh(1) - hh(2)))
      ! tll: h integrates the (clamped, renormalised) interpolant, as vinecopulib does
      write(msg, '(a,a,a,es9.2)') 'h1 = integral of c ', trim(FV_FAMILY_NAME(k)), ' err ', interr
      call check(interr < merge(2.0e-2_dp, 1.0e-6_dp, k == FV_TLL), msg)
    end do
    ! selection among parametric + nonparametric candidates
    ctl%nfam = 4; ctl%family_set(1:4) = [FV_INDEP, FV_GUMBEL, FV_TLL, FV_BERNSTEIN]
    call b%select(w1, u2, ctl)
    call check(b%family == FV_GUMBEL, 'select prefers the true parametric family over tll/bernstein (BIC)')
  end subroutine test_nonpar

  ! ---------------------------------------------------------------- weights
  subroutine test_weights()
    integer, parameter :: n = 1500, d = 4
    type(fv_bicop) :: truth, b1, b2
    type(fv_vinecop) :: v1, v2
    type(fv_controls) :: ctl
    real(dp) :: w1(n), w2(n), u2(n), wt(n), x1(2 * n), x2(2 * n), u(n, d), one(n)
    integer :: i, k, fams(2)
    character(len=200) :: msg
    truth%family = FV_CLAYTON; truth%par = [2.0_dp, 0.0_dp, 0.0_dp]
    call random_number(w1); call random_number(w2)
    call truth%hinv1(w1, w2, u2)
    call random_number(wt)
    wt = merge(2.0_dp, 1.0_dp, wt > 0.6_dp)
    k = 0
    do i = 1, n
      k = k + 1; x1(k) = w1(i); x2(k) = u2(i)
      if (wt(i) > 1.5_dp) then
        k = k + 1; x1(k) = w1(i); x2(k) = u2(i)
      end if
    end do
    call check(abs(fv_dependence(w1, u2, FV_TREE_TAU, wt) - fv_kendall_tau(x1(1:k), x2(1:k))) < 1.0e-12_dp, &
               'weighted tau = tau of duplicated rows')
    fams = [FV_CLAYTON, FV_BB1]
    do i = 1, 2
      b1%family = fams(i); b2%family = fams(i)
      call b1%fit(w1, u2, weights=wt)
      call b2%fit(x1(1:k), x2(1:k))
      write(msg, '(a,a,4f12.5)') 'weighted MLE = MLE on duplicated rows ', trim(FV_FAMILY_NAME(fams(i))), &
        b1%par(1:2), b2%par(1:2)
      call check(maxval(abs(b1%par - b2%par)) < 1.0e-3_dp .and. abs(b1%loglik - b2%loglik) < 1.0e-6_dp * abs(b2%loglik), msg)
    end do
    ! unit weights reproduce the unweighted vine
    call v1%set_dvine([1, 2, 3, 4])
    v1%fam(1:3, 1) = FV_GUMBEL; v1%par(1, 1:3, 1) = 1.7_dp
    v1%fam(1:2, 2) = FV_FRANK; v1%par(1, 1:2, 2) = 3.0_dp
    call v1%simulate(n, u, seed=7)
    one = 1.0_dp
    ctl%nfam = 12; ctl%family_set(1:12) = FV_PARAMETRIC
    call v1%init(d); call v1%select(u, ctl)
    call v2%init(d); call v2%select(u, ctl, weights=one)
    call check(all(v1%fam == v2%fam) .and. all(v1%st == v2%st) .and. abs(v1%loglik(u) - v2%loglik(u)) < 1.0e-6_dp, &
               'unit weights reproduce the unweighted selection')
  end subroutine test_weights

  ! ---------------------------------------------------------------- discrete margins
  subroutine test_discrete()
    integer, parameter :: m = 300, n = 2000
    type(fv_vinecop) :: v, sel
    type(fv_bicop) :: b
    type(fv_controls) :: ctl
    real(dp) :: uu(1, 2), um(1, 2), p(1), s, a, am, bb, bm, x(m), cc(m), y(m), ref
    real(dp) :: u(n, 3), umn(n, 3), z(n, 3)
    integer :: k, i, j, cse, r, rots(4), kk
    character(len=200) :: msg
    rots = [0, 90, 180, 270]
    ! every family x rotation: closed-form copula CDF differences (and the quadrature
    ! fallback for gaussian / student) against brute-force integrals of the density
    do kk = 1, 25
      k = kk
      if (kk == 24) k = -1          ! bernstein, fitted below
      if (kk == 25) exit
      do r = 1, 4
      if (k > 0) then
        call family_case(k, b)
      else
        call bern_case(b)
      end if
      if (r > 1 .and. (b%family <= FV_STUDENT .or. b%family == FV_FRANK .or. (b%family >= FV_PLACKETT .and. &
          b%family <= FV_AMH) .or. b%family >= FV_TLL)) cycle
      b%rotation = rots(r)
      call v%set_dvine([1, 2])
      v%fam(1, 1) = b%family; v%rot(1, 1) = b%rotation; v%par(:, 1, 1) = b%par
      if (b%family >= FV_TLL) then
        deallocate(v%aux); allocate(v%aux(NAUX, 2, 2)); v%aux = 0.0_dp; v%aux(:, 1, 1) = b%aux
      end if
      a = 0.62_dp; am = 0.31_dp; bb = 0.47_dp; bm = 0.29_dp
      do cse = 1, 3
        v%disc = [cse /= 2, cse /= 1]          ! 1: (d, c)  2: (c, d)  3: (d, d)
        uu(1, :) = [a, bb]; um(1, :) = [am, bm]
        call v%pdf(uu, p, u_minus=um)
        ! brute force: midpoint integral of the density over the discrete cell(s)
        do i = 1, m
          x(i) = (i - 0.5_dp) / m
        end do
        select case (cse)
        case (1)
          y = bb
          call b%pdf(am + (a - am) * x, y, cc); ref = sum(cc) / m
        case (2)
          y = a
          call b%pdf(y, bm + (bb - bm) * x, cc); ref = sum(cc) / m
        case (3)
          ref = 0.0_dp
          do j = 1, m
            y = bm + (bb - bm) * (j - 0.5_dp) / m
            call b%pdf(am + (a - am) * x, y, cc); ref = ref + sum(cc)
          end do
          ref = ref / (real(m, dp) * m)
        end select
        write(msg, '(a,a,a,i0,a,i0,2f12.7)') 'discrete pair density ', trim(FV_FAMILY_NAME(b%family)), ' rot ', &
          rots(r), ' case ', cse, p(1), ref
        call check(abs(p(1) - ref) < 2.0e-5_dp * max(1.0_dp, ref), msg)
      end do
      end do
    end do
    ! selection with a discrete variable: the selected vine is at least as likely as the truth
    call v%set_dvine([1, 2, 3])
    v%fam(1:2, 1) = [FV_GAUSSIAN, FV_CLAYTON]; v%par(1, 1, 1) = 0.7_dp; v%par(1, 2, 1) = 2.0_dp
    v%fam(1, 2) = FV_FRANK; v%par(1, 1, 2) = 2.0_dp
    call v%simulate(n, z, seed=3)
    u = z; umn = z
    do i = 1, n                              ! variable 2 -> 5 equiprobable categories
      u(i, 2) = (floor(z(i, 2) * 5) + 1) / 5.0_dp
      umn(i, 2) = floor(z(i, 2) * 5) / 5.0_dp
    end do
    ctl%nfam = 12; ctl%family_set(1:12) = FV_PARAMETRIC
    call sel%init(3)
    call sel%select(u, ctl, discrete=[.false., .true., .false.], u_minus=umn)
    v%disc = [.false., .true., .false.]
    s = v%loglik(u, umn)
    write(msg, '(a,2f12.3)') 'discrete select loglik >= truth - npars ', sel%loglik(u, umn), s
    call check(sel%loglik(u, umn) >= s - real(v%npars(), dp), msg)
  end subroutine test_discrete

  subroutine bern_case(b)
    type(fv_bicop), intent(out) :: b
    type(fv_bicop) :: t
    real(dp) :: w1(3000), w2(3000), u2(3000)
    t%family = FV_CLAYTON; t%par = [2.0_dp, 0.0_dp, 0.0_dp]
    call random_number(w1); call random_number(w2)
    call t%hinv1(w1, w2, u2)
    b%family = FV_BERNSTEIN
    call b%fit(w1, u2)
  end subroutine bern_case

  ! ---------------------------------------------------------------- truncation, criteria
  subroutine test_trunc_and_criteria()
    integer, parameter :: n = 3000, d = 6
    type(fv_vinecop) :: v, sel
    type(fv_controls) :: ctl
    real(dp) :: u(n, d), z(n, 2), r, mi
    integer :: k, info
    character(len=200) :: msg
    ! only tree 1 carries dependence -> mBICv truncates after tree 1
    call v%set_dvine([1, 2, 3, 4, 5, 6])
    v%fam(1:5, 1) = FV_GAUSSIAN; v%par(1, 1:5, 1) = 0.6_dp
    call v%simulate(n, u, seed=11)
    ctl%nfam = 12; ctl%family_set(1:12) = FV_PARAMETRIC
    ctl%selcrit = FV_MBICV; ctl%select_trunc = .true.
    call sel%init(d); call sel%select(u, ctl)
    write(msg, '(a,i0)') 'mBICv truncation level (truth 1): ', sel%trunc
    call check(sel%trunc == 1, msg)
    ! dependence in every tree -> no truncation
    do k = 2, 5
      v%fam(1:d - k, k) = FV_CLAYTON; v%par(1, 1:d - k, k) = 1.5_dp
    end do
    call v%simulate(n, u, seed=12)
    call sel%init(d); call sel%select(u, ctl)
    write(msg, '(a,i0)') 'mBICv truncation level (truth 5): ', sel%trunc
    call check(sel%trunc == d - 1, msg)
    ! new tree criteria give valid vines at least as good (within 5%) as the tau-based one
    ! (this truth is hard for every criterion: strong higher trees make non-neighbours
    ! more dependent than the tree-1 pairs)
    ctl%selcrit = FV_BIC; ctl%select_trunc = .false.
    call sel%init(d); call sel%select(u, ctl)
    mi = sel%loglik(u)
    do k = FV_TREE_XI, FV_TREE_FERREIRA
      ctl%tree_crit = k
      call sel%init(d); call sel%select(u, ctl, info=info)
      write(msg, '(a,i0,a,f12.2,a,f12.2)') 'tree criterion ', k, ' loglik ', sel%loglik(u), ' tau-based ', mi
      call check(info == 0 .and. sel%loglik(u) > 0.95_dp * mi, msg)
    end do
    ! KSG mutual information of a Gaussian pair vs -log(1 - rho^2)/2
    r = 0.7_dp
    call random_number(z)
    z = min(max(z, 1.0e-12_dp), 1.0_dp - 1.0e-12_dp)
    z(:, 2) = r * qn(z(:, 1)) + sqrt(1 - r * r) * qn(z(:, 2)); z(:, 1) = qn(z(:, 1))
    mi = fv_dependence(z(:, 1), z(:, 2), FV_TREE_MI)
    write(msg, '(a,2f9.4)') 'KSG mutual information vs exact ', mi, -0.5_dp * log(1 - r * r)
    call check(abs(mi + 0.5_dp * log(1 - r * r)) < 0.04_dp, msg)
    call check(fv_dependence(z(:, 1), z(:, 2), FV_TREE_XI) > 0.2_dp, 'Chatterjee xi detects dependence')
  end subroutine test_trunc_and_criteria

  ! ---------------------------------------------------------------- core extensions
  subroutine test_core_extensions()
    use iso_fortran_env, only: int64
    integer, parameter :: n = 2000, m = 1000
    type(fv_vinecop) :: v, sel, s2
    type(fv_bicop) :: b, g
    type(fv_controls) :: ctl
    real(dp) :: z(n, 3), u(n, 3), um(n, 3), w(n, 3), wm(n, 3), w2(n, 3), h(n), hm(n), x(m), cc(m), y(m)
    real(dp) :: ref, cv(1), a, bb, p(1), r, pc, s, e1, mx, ks, ll0, lp(n), cr, c4(4)
    real(dp) :: uu(1, 3), uum(1, 3), d1, d3
    integer :: i, j, k, fams(3)
    integer(int64) :: hits, miss, miss0
    character(len=200) :: msg
    ! -- copula CDFs against midpoint double integrals of the density over [0,a] x [0,b]
    do i = 1, m
      x(i) = (i - 0.5_dp) / m
    end do
    fams = [FV_GAUSSIAN, FV_TLL, FV_BERNSTEIN]
    a = 0.62_dp; bb = 0.47_dp
    do k = 1, 3
      if (fams(k) == FV_GAUSSIAN) then
        b = fv_bicop(); b%family = FV_GAUSSIAN; b%par = [0.6_dp, 0.0_dp, 0.0_dp]
      else if (fams(k) == FV_BERNSTEIN) then
        call bern_case(b)
      else
        call tll_case(b)
      end if
      ref = 0.0_dp
      do j = 1, m
        y = bb * (j - 0.5_dp) / m
        call b%pdf(a * x, y, cc); ref = ref + sum(cc)
      end do
      ref = ref * a * bb / (real(m, dp) * m)
      call b%cdf([a], [bb], cv)
      write(msg, '(a,a,2f14.9)') 'copula cdf = integral of the density: ', trim(FV_FAMILY_NAME(b%family)), cv(1), ref
      call check(abs(cv(1) - ref) < 2.0e-6_dp, msg)
    end do
    b = fv_bicop(); b%family = FV_GAUSSIAN; b%par = [0.6_dp, 0.0_dp, 0.0_dp]
    call b%cdf([0.5_dp], [0.5_dp], cv)
    e1 = cv(1) - (0.25_dp + asin(0.6_dp) / (4.0_dp * acos(0.0_dp)))
    write(msg, '(a,es10.2)') 'bivariate normal C(1/2, 1/2) = 1/4 + asin(r)/(2 pi), err ', e1
    call check(abs(e1) < 1.0e-14_dp, msg)
    ! -- tll on a discrete pair: the (d, d) density is the exact integral of the interpolant
    !    over the cell (differences of the copula CDF)
    call tll_case(b)
    call v%set_dvine([1, 2])
    v%fam(1, 1) = FV_TLL
    deallocate(v%aux); allocate(v%aux(NAUX, 2, 2)); v%aux = 0.0_dp; v%aux(:, 1, 1) = b%aux
    v%disc = [.true., .true.]
    uu(1, 1:2) = [0.62_dp, 0.47_dp]; uum(1, 1:2) = [0.31_dp, 0.29_dp]
    call v%pdf(uu(:, 1:2), p, u_minus=uum(:, 1:2))
    ref = 0.0_dp
    do j = 1, m
      y = 0.29_dp + 0.18_dp * (j - 0.5_dp) / m
      call b%pdf(0.31_dp + 0.31_dp * x, y, cc); ref = ref + sum(cc)
    end do
    ref = ref / (real(m, dp) * m)
    write(msg, '(a,2f12.7)') 'tll discrete (d, d) pair density = cell integral ', p(1), ref
    call check(abs(p(1) - ref) < 1.0e-6_dp, msg)
    ! -- discrete Rosenblatt: data from a vine with a 5-category variable 2
    call v%set_dvine([1, 2, 3])
    v%fam(1:2, 1) = [FV_GAUSSIAN, FV_CLAYTON]; v%par(1, 1, 1) = 0.7_dp; v%par(1, 2, 1) = 2.0_dp
    v%fam(1, 2) = FV_FRANK; v%par(1, 1, 2) = 2.0_dp
    call v%simulate(n, z, seed=21)
    u = z; um = z
    do i = 1, n
      u(i, 2) = (floor(z(i, 2) * 5) + 1) / 5.0_dp; um(i, 2) = floor(z(i, 2) * 5) / 5.0_dp
    end do
    v%disc = [.false., .true., .false.]
    call v%rosenblatt(u, w, u_minus=um, randomize=.false., out_minus=wm)
    call check(all(w >= wm - 1.0e-15_dp) .and. all(abs(w(:, [1, 3]) - wm(:, [1, 3])) < 1.0e-15_dp), &
               'discrete rosenblatt: F >= F- (equal for continuous variables)')
    call v%rosenblatt(u, w2, u_minus=um, seed=5)
    call check(all(w2 >= wm - 1.0e-15_dp .and. w2 <= w + 1.0e-15_dp), 'randomised PIT lies between F- and F')
    mx = 0.0_dp
    do k = 1, 3
      h = w2(:, k); call sortv(h)
      ks = 0.0_dp
      do i = 1, n
        ks = max(ks, abs(h(i) - (i - 0.5_dp) / n))
      end do
      mx = max(mx, ks)
    end do
    s = max(abs(corr(w2(:, 1), w2(:, 2))), abs(corr(w2(:, 1), w2(:, 3))), abs(corr(w2(:, 2), w2(:, 3))))
    write(msg, '(a,2f9.4)') 'randomised discrete PIT is iid uniform: max KS, max |corr| ', mx, s
    call check(mx < 1.63_dp / sqrt(real(n, dp)) .and. s < 4.0_dp / sqrt(real(n, dp)), msg)
    ! exact two-variable case: the conditioned column is h1 at F(x2) and at F(x2-)
    call s2%set_dvine([1, 2]); s2%fam(1, 1) = FV_GAUSSIAN; s2%par(1, 1, 1) = 0.7_dp; s2%disc = [.false., .true.]
    call s2%rosenblatt(u(:, 1:2), w(:, 1:2), u_minus=um(:, 1:2), randomize=.false., out_minus=wm(:, 1:2))
    b = fv_bicop(); b%family = FV_GAUSSIAN; b%par = [0.7_dp, 0.0_dp, 0.0_dp]
    call b%hfunc1(u(:, 1), u(:, 2), h); call b%hfunc1(u(:, 1), max(um(:, 2), 1.0e-10_dp), hm)
    h = min(max(h, 1.0e-10_dp), 1.0_dp - 1.0e-10_dp); hm = min(max(hm, 1.0e-10_dp), 1.0_dp - 1.0e-10_dp)
    e1 = max(maxval(abs(w(:, 1) - u(:, 1))), maxval(abs(w(:, 2) - h)), maxval(abs(wm(:, 2) - hm)))
    if (maxval(abs(w(:, 1) - u(:, 1))) > 1.0e-12_dp) then
      ! structure ordered the other way: variable 2 unconditioned, variable 1 given 2
      e1 = maxval(abs(w(:, 2) - u(:, 2)))
    end if
    write(msg, '(a,es10.2)') 'discrete rosenblatt (c, d) = h-functions at F and F-: err ', e1
    call check(e1 < 1.0e-9_dp, msg)
    ! -- trivariate gaussian (d, c, d): vine density = exact rectangle probability / masses
    r = 0.5_dp; pc = 0.3_dp
    call v%set_dvine([1, 2, 3])
    v%fam = 0; v%fam(1:2, 1) = FV_GAUSSIAN; v%par = 0.0_dp; v%par(1, 1:2, 1) = r
    v%fam(1, 2) = FV_GAUSSIAN; v%par(1, 1, 2) = pc
    v%disc = [.true., .false., .true.]
    g = fv_bicop(); g%family = FV_GAUSSIAN; g%par = [pc, 0.0_dp, 0.0_dp]
    mx = 0.0_dp
    do i = 1, 40
      uu(1, :) = [0.15_dp + 0.02_dp * i, 0.05_dp + 0.022_dp * i, 0.9_dp - 0.015_dp * i]
      uum(1, :) = [uu(1, 1) - 0.12_dp, uu(1, 2), uu(1, 3) - 0.2_dp]
      call v%pdf(uu, p, u_minus=uum)
      s = qn(uu(1, 2)); e1 = sqrt(1.0_dp - r * r)
      d1 = 0.5_dp * erfc(-((qn(uu(1, 1)) - r * s) / e1) / sqrt(2.0_dp))
      a = 0.5_dp * erfc(-((qn(uum(1, 1)) - r * s) / e1) / sqrt(2.0_dp))
      d3 = 0.5_dp * erfc(-((qn(uu(1, 3)) - r * s) / e1) / sqrt(2.0_dp))
      bb = 0.5_dp * erfc(-((qn(uum(1, 3)) - r * s) / e1) / sqrt(2.0_dp))
      call g%cdf([d1, a, d1, a], [d3, d3, bb, bb], c4)
      ref = (c4(1) - c4(2) - c4(3) + c4(4)) / (0.12_dp * 0.2_dp)
      mx = max(mx, abs(p(1) / ref - 1.0_dp))
    end do
    write(msg, '(a,es10.2)') 'trivariate gaussian (d, c, d) density = exact rectangle probability, rel err ', mx
    call check(mx < 1.0e-8_dp, msg)
    ! -- tll / bernstein selected on discrete pairs
    call v%set_dvine([1, 2, 3])
    v%fam(1:2, 1) = [FV_GAUSSIAN, FV_CLAYTON]; v%par(1, 1, 1) = 0.7_dp; v%par(1, 2, 1) = 2.0_dp
    v%fam(1, 2) = FV_FRANK; v%par(1, 1, 2) = 2.0_dp
    v%disc = [.false., .true., .false.]
    ll0 = v%loglik(u, um)
    do k = FV_TLL, FV_BERNSTEIN
      ctl = fv_controls(); ctl%nfam = 1; ctl%family_set(1) = k
      call sel%init(3); call sel%select(u, ctl, discrete=[.false., .true., .false.], u_minus=um)
      write(msg, '(a,a,2f12.3)') 'nonparametric on discrete pairs (loglik, truth): ', trim(FV_FAMILY_NAME(k)), &
        sel%loglik(u, um), ll0
      call check(sel%fam(1, 1) == k .and. sel%loglik(u, um) > 0.7_dp * ll0, msg)
    end do
    ! -- select_threshold: the chosen threshold refitted as a fixed threshold gives the same vine
    ctl = fv_controls(); ctl%nfam = 12; ctl%family_set(1:12) = FV_PARAMETRIC; ctl%selcrit = FV_MBICV
    ctl%select_threshold = .true.
    call v%set_dvine([1, 2, 3])
    v%fam(1:2, 1) = [FV_GAUSSIAN, FV_CLAYTON]; v%par(1, 1, 1) = 0.7_dp; v%par(1, 2, 1) = 0.1_dp
    v%fam(1, 2) = FV_FRANK; v%par(1, 1, 2) = 0.3_dp; v%disc = .false.
    call v%simulate(n, z, seed=22)
    call sel%init(3); call sel%select(z, ctl)
    cr = sel%threshold
    ctl%select_threshold = .false.; ctl%threshold = cr
    call s2%init(3); call s2%select(z, ctl)
    write(msg, '(a,f9.5,2f12.4)') 'select_threshold: threshold, loglik, refit at that threshold ', cr, &
      sel%loglik(z), s2%loglik(z)
    call check(cr > 0.0_dp .and. cr <= 1.0_dp .and. abs(sel%loglik(z) - s2%loglik(z)) < 1.0e-9_dp, msg)
    ! -- fit cache: a repeated selection is served from the cache (no new misses) and is identical
    ! (one thread: with several, racing decisions depend on which fits finish first)
    ctl = fv_controls(); ctl%nthreads = 1
    call fv_fitcache_enable(4096)
    call sel%init(3); call sel%select(z, ctl)
    call fv_fitcache_stats(hits, miss0)
    call s2%init(3); call s2%select(z, ctl)
    call fv_fitcache_stats(hits, miss)
    call fv_fitcache_disable()
    write(msg, '(a,3i6)') 'fit cache: second selection all hits (hits, misses 1st, misses both) ', hits, miss0, miss
    call check(hits > 0 .and. miss == miss0 .and. all(sel%fam == s2%fam) .and. all(sel%par == s2%par), msg)
    ! -- t-EV: the joint (rho, nu) fit is a local maximum of the likelihood
    call v%set_dvine([1, 2]); v%fam(1, 1) = FV_TEV; v%par(:, 1, 1) = [0.6_dp, 5.0_dp, 0.0_dp]; v%disc = .false.
    call v%simulate(n, z(:, 1:2), seed=23)
    b = fv_bicop(); b%family = FV_TEV; call b%fit(z(:, 1), z(:, 2))
    g = b; mx = -huge(1.0_dp)
    do i = -1, 1
      do j = -1, 1
        if (i == 0 .and. j == 0) cycle
        g%par(1) = b%par(1) + 2.0e-3_dp * i; g%par(2) = b%par(2) * (1.0_dp + 1.0e-2_dp * j)
        call g%logpdf(z(:, 1), z(:, 2), lp); mx = max(mx, sum(lp))
      end do
    end do
    call b%logpdf(z(:, 1), z(:, 2), lp)
    write(msg, '(a,2f10.4,2f12.5)') 't-EV fit (rho, nu), loglik vs best neighbour ', b%par(1:2), sum(lp), mx
    call check(sum(lp) >= mx .and. abs(sum(lp) - b%loglik) < 1.0e-8_dp * abs(sum(lp)), msg)
  end subroutine test_core_extensions

  subroutine tll_case(b)
    type(fv_bicop), intent(out) :: b
    type(fv_bicop) :: t
    real(dp) :: w1(3000), w2(3000), u2(3000)
    t%family = FV_CLAYTON; t%par = [2.0_dp, 0.0_dp, 0.0_dp]
    call random_number(w1); call random_number(w2)
    call t%hinv1(w1, w2, u2)
    b%family = FV_TLL
    call b%fit(w1, u2)
  end subroutine tll_case

  subroutine sortv(a)
    real(dp), intent(inout) :: a(:)
    integer :: i, j
    real(dp) :: t
    do i = 2, size(a)
      t = a(i); j = i - 1
      do while (j >= 1)
        if (a(j) <= t) exit
        a(j + 1) = a(j); j = j - 1
      end do
      a(j + 1) = t
    end do
  end subroutine sortv

  ! y = rho(x1) x2 + sqrt(1 - rho^2) e (normal scores), rho = 0.8 tanh(z1): the copula of (y, x2)
  ! given x1 varies with x1.  D-vine y - x1 - x2: tree-2 edge (y, x2 | x1) is the target edge.
  subroutine test_nonsimplified()
    integer, parameter :: n = 4000
    real(dp) :: u(n, 3), ut(n, 3), llr, lls
    type(fv_vinecop) :: v, vs
    type(fv_controls) :: c
    integer :: nr, k, nrej
    print '(a)', '== non-simplified edges'
    c%nfam = 2; c%family_set(1:2) = [FV_INDEP, FV_GAUSSIAN]
    do k = 1, 2
      call gen_ns(u, k == 1)
      call v%set_dvine([3, 1, 2])
      call v%select(u, c)
      vs = v
      call v%relax(u, 3, 1, c, nrel_out=nr)
      if (k == 1) then
        call check(nr == 1, 'CCC test rejects a varying conditional correlation (1 edge relaxed)')
        call gen_ns(ut, .true.)
        llr = v%loglik(ut); lls = vs%loglik(ut)
        print '(a,2f12.2)', '  held-out loglik relaxed / simplified', llr, lls
        call check(llr > lls + 50.0_dp, 'relaxed vine predicts held-out data better')
        call check(v%npars() > vs%npars(), 'relaxed edges count their anchor parameters')
      end if
    end do
    nrej = 0                                  ! constant conditional correlation: ~5 % rejections
    do k = 1, 20
      call gen_ns(u, .false.)
      call v%set_dvine([3, 1, 2])
      call v%select(u, c)
      call v%relax(u, 3, 1, c, nrel_out=nr)
      nrej = nrej + nr
    end do
    print '(a,i0,a)', '  constant correlation: relaxed in ', nrej, ' of 20 samples'
    call check(nrej <= 4, 'CCC test holds its level on a simplified vine')
  end subroutine test_nonsimplified


  subroutine gen_ns(uu, varying)
      real(dp), intent(out) :: uu(:, :)
      logical, intent(in) :: varying
      real(dp) :: z(size(uu, 1), 3), e(size(uu, 1)), rh(size(uu, 1)), w1(size(uu, 1)), w2(size(uu, 1))
      integer :: j
      do j = 1, 3
        call random_number(w1); call random_number(w2)
        z(:, j) = sqrt(-2.0_dp * log(max(w1, 1e-300_dp))) * cos(6.283185307179586_dp * w2)
      end do
      e = z(:, 3)
      if (varying) then
        rh = 0.8_dp * tanh(z(:, 1))
      else
        rh = 0.5_dp
      end if
      z(:, 3) = rh * z(:, 2) + sqrt(1.0_dp - rh**2) * e
      uu = 0.5_dp * erfc(-z / sqrt(2.0_dp))
  end subroutine gen_ns

  real(dp) function corr(a, b)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: ma, mb
    ma = sum(a) / size(a); mb = sum(b) / size(b)
    corr = sum((a - ma) * (b - mb)) / sqrt(sum((a - ma)**2) * sum((b - mb)**2))
  end function corr

  elemental real(dp) function qn(p)          ! probit by bisection on erfc (test-only reference)
    real(dp), intent(in) :: p
    real(dp) :: lo, hi, mid
    integer :: it
    lo = -9.0_dp; hi = 9.0_dp
    do it = 1, 60
      mid = 0.5_dp * (lo + hi)
      if (0.5_dp * erfc(-mid / sqrt(2.0_dp)) < p) then
        lo = mid
      else
        hi = mid
      end if
    end do
    qn = 0.5_dp * (lo + hi)
  end function qn

end program test_fortvine
