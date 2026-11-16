! fv_bench -- fortvine side of the cross-library benchmark (bench/run_bench.py drives it).
!   fv_bench eval    fam p1 p2 p3 n reps          ns per evaluation: logpdf, hfunc1, hinv1
!   fv_bench fit     data.csv fam reps            ms per MLE fit of a given family, log-likelihood
!   fv_bench select  data.csv reps                ms per family selection (common family set), result
!   fv_bench vine    data.csv nthreads reps       s per vine selection, log-likelihood, parameters,
!                                                 then s for pdf of all rows and for simulating n rows
! Data: CSV without header, n rows x d columns of pseudo-observations in (0, 1).
! Every line of output is "key=value" pairs, one result per line.
program fv_bench
  use fortvine_api
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  ! the families every compared library has: indep, gaussian, student, clayton, gumbel, frank,
  ! joe, bb1, bb6, bb7, bb8, tawn
  integer, parameter :: COMMON(12) = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
  character(len=512) :: mode, a
  real(dp), allocatable :: u(:, :)
  call get_command_argument(1, mode)
  select case (trim(mode))
  case ('eval'); call do_eval()
  case ('fit'); call do_fit()
  case ('select'); call do_select()
  case ('vine'); call do_vine()
  case default
    print '(a)', 'usage: fv_bench eval|fit|select|vine ...'
  end select
contains
  real(dp) function now()
    integer(8) :: c, r
    call system_clock(c, r); now = real(c, dp) / real(r, dp)
  end function now

  integer function iarg(k)
    integer, intent(in) :: k
    call get_command_argument(k, a); read(a, *) iarg
  end function iarg

  real(dp) function rarg(k)
    integer, intent(in) :: k
    call get_command_argument(k, a); read(a, *) rarg
  end function rarg

  subroutine do_eval()
    type(fv_bicop) :: b
    real(dp), allocatable :: u1(:), u2(:), out(:)
    real(dp) :: t0, tl, th, ti
    integer :: n, reps, r
    b%family = iarg(2); b%par = [rarg(3), rarg(4), rarg(5)]; n = iarg(6); reps = iarg(7)
    allocate(u1(n), u2(n), out(n))
    call random_seed(); call random_number(u1); call random_number(u2)
    u1 = 1e-6_dp + (1 - 2e-6_dp) * u1; u2 = 1e-6_dp + (1 - 2e-6_dp) * u2
    call b%logpdf(u1, u2, out)                    ! warm-up
    t0 = now(); do r = 1, reps; call b%logpdf(u1, u2, out); end do; tl = (now() - t0) / reps
    t0 = now(); do r = 1, reps; call b%hfunc1(u1, u2, out); end do; th = (now() - t0) / reps
    t0 = now(); do r = 1, reps; call b%hinv1(u1, u2, out); end do; ti = (now() - t0) / reps
    print '(a,es12.4,a,es12.4,a,es12.4)', 'logpdf_ns=', 1e9_dp * tl / n, ' hfunc1_ns=', 1e9_dp * th / n, &
      ' hinv1_ns=', 1e9_dp * ti / n
  end subroutine do_eval

  subroutine do_fit()
    type(fv_bicop) :: b
    type(fv_controls) :: ctl
    real(dp) :: t0, t
    integer :: reps, r
    call get_command_argument(2, a); call read_csv(trim(a), u)
    reps = iarg(4); ctl%nthreads = 1
    t0 = now()
    do r = 1, reps
      b%family = iarg(3); b%rotation = 0
      call b%fit(u(:, 1), u(:, 2), ctl)
    end do
    t = (now() - t0) / reps
    print '(a,es12.4,a,es16.8,a,3es14.6)', 'ms=', 1e3_dp * t, ' loglik=', b%loglik, ' par=', b%par
  end subroutine do_fit

  subroutine do_select()
    type(fv_bicop) :: b
    type(fv_controls) :: ctl
    real(dp) :: t0, t
    integer :: reps, r
    call get_command_argument(2, a); call read_csv(trim(a), u)
    reps = iarg(3)
    ctl%nthreads = 1; ctl%nfam = size(COMMON); ctl%family_set = 0; ctl%family_set(1:size(COMMON)) = COMMON
    t0 = now()
    do r = 1, reps
      call b%select(u(:, 1), u(:, 2), ctl)
    end do
    t = (now() - t0) / reps
    print '(a,es12.4,a,a,a,i0,a,es16.8)', 'ms=', 1e3_dp * t, ' family=', trim(fv_fam_name(b%family)), &
      ' rotation=', b%rotation, ' loglik=', b%loglik
  end subroutine do_select

  subroutine do_vine()
    type(fv_vinecop) :: vc
    type(fv_controls) :: ctl
    real(dp), allocatable :: dens(:), sim(:, :)
    real(dp) :: t0, ts, tp, tm, ll
    integer :: reps, r, n, d
    call get_command_argument(2, a); call read_csv(trim(a), u)
    n = size(u, 1); d = size(u, 2)
    ctl%nthreads = iarg(3); reps = iarg(4)
    ctl%nfam = size(COMMON); ctl%family_set = 0; ctl%family_set(1:size(COMMON)) = COMMON
    t0 = now()
    do r = 1, reps
      call vc%select(u, ctl)
    end do
    ts = (now() - t0) / reps
    ll = vc%loglik(u)
    allocate(dens(n), sim(n, d))
    t0 = now(); call vc%pdf(u, dens, ctl%nthreads); tp = now() - t0
    t0 = now(); call vc%simulate(n, sim, seed=1, nthreads=ctl%nthreads); tm = now() - t0
    print '(a,es12.4,a,es16.8,a,i0,a,es12.4,a,es12.4)', 'select_s=', ts, ' loglik=', ll, ' npars=', vc%npars(), &
      ' pdf_s=', tp, ' simulate_s=', tm
  end subroutine do_vine

  subroutine read_csv(fn, x)
    character(len=*), intent(in) :: fn
    real(dp), allocatable, intent(out) :: x(:, :)
    character(len=8192) :: ln
    integer :: uu, n, ios, j, q
    open(newunit=uu, file=fn, status='old', action='read')
    read(uu, '(a)') ln
    q = count([(ln(j:j) == ',', j = 1, len_trim(ln))]) + 1
    n = 1
    do
      read(uu, '(a)', iostat=ios) ln; if (ios /= 0) exit; n = n + 1
    end do
    rewind(uu); allocate(x(n, q))
    do j = 1, n
      read(uu, *) x(j, :)
    end do
    close(uu)
  end subroutine read_csv
end program fv_bench
