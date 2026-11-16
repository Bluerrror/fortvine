! Benchmark: bench <data.csv> [nreps]   (CSV of pseudo-observations, n rows x d cols)
! Times vine selection (full MLE, all parametric families), pdf, pdf_grid,
! rosenblatt and inverse_rosenblatt; prints loglik and npars.
program bench
  use fortvine_api
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  character(len=512) :: path, arg
  real(dp), allocatable :: u(:, :), p(:), w(:, :), b(:, :), ug(:, :, :), pg(:, :)
  type(fv_vinecop) :: vc
  type(fv_controls) :: ctl
  integer :: n, d, i, k, reps, r
  integer(8) :: t0, t1, rate
  real(dp) :: ts, tbest
  call get_command_argument(1, path)
  reps = 3
  if (command_argument_count() >= 2) then
    call get_command_argument(2, arg); read(arg, *) reps
  end if
  call read_csv(trim(path), u)
  n = size(u, 1); d = size(u, 2)
  print '(a,i0,a,i0,a,i0)', 'n = ', n, '  d = ', d, '  threads = ', fv_threads()
  tbest = huge(1.0_dp)
  do r = 1, reps
    call system_clock(t0, rate)
    call vc%init(d)
    call vc%select(u, ctl)
    call system_clock(t1)
    tbest = min(tbest, real(t1 - t0, dp) / rate)
  end do
  print '(a,f10.1,a,f12.2,a,i0)', 'select (MLE, default 23 set):    ', 1e3 * tbest, ' ms   loglik ', vc%loglik(u), &
    '  npars ', vc%npars()
  ctl%nfam = 12; ctl%family_set(1:12) = FV_PARAMETRIC
  tbest = huge(1.0_dp)
  do r = 1, reps
    call system_clock(t0)
    call vc%init(d)
    call vc%select(u, ctl)
    call system_clock(t1)
    tbest = min(tbest, real(t1 - t0, dp) / rate)
  end do
  print '(a,f10.1,a,f12.2,a,i0)', 'select (MLE, vinecopulib 12 set): ', 1e3 * tbest, ' ms   loglik ', vc%loglik(u), &
    '  npars ', vc%npars()
  ctl%fast_mle = .true.
  tbest = huge(1.0_dp)
  do r = 1, reps
    call system_clock(t0)
    call vc%init(d)
    call vc%select(u, ctl)
    call system_clock(t1)
    tbest = min(tbest, real(t1 - t0, dp) / rate)
  end do
  print '(a,f10.1,a,f12.2)', 'select (fast_mle, 12 set):         ', 1e3 * tbest, ' ms   loglik ', vc%loglik(u)
  allocate(p(n), w(n, d), b(n, d))
  call timeit('pdf', 1)
  call timeit('rosenblatt', 2)
  call timeit('inverse_rosenblatt', 3)
  allocate(ug(d, 300, 1000), pg(300, 1000))
  do i = 1, 1000
    do k = 1, 300
      ug(:, k, i) = u(i, :)
      ug(d, k, i) = (k - 0.5_dp) / 300
    end do
  end do
  call timeit('pdf_grid 1000 x 300 bins', 4)
contains
  subroutine timeit(name, what)
    character(*), intent(in) :: name
    integer, intent(in) :: what
    integer :: rr
    tbest = huge(1.0_dp)
    do rr = 1, max(reps, 5)
      call system_clock(t0)
      select case (what)
      case (1); call vc%pdf(u, p)
      case (2); call vc%rosenblatt(u, w)
      case (3); call vc%inverse_rosenblatt(w, b)
      case (4); call vc%pdf_grid(ug, d, pg)
      end select
      call system_clock(t1)
      tbest = min(tbest, real(t1 - t0, dp) / rate)
    end do
    print '(a,a,f10.2,a)', name, ': ', 1e3 * tbest, ' ms'
  end subroutine timeit

  subroutine read_csv(fn, x)
    character(*), intent(in) :: fn
    real(dp), allocatable, intent(out) :: x(:, :)
    integer :: un, ios, nr, nc, j
    character(len=4096) :: line
    open(newunit=un, file=fn, status='old', action='read')
    nr = 0
    do
      read(un, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) > 0) nr = nr + 1
    end do
    rewind(un)
    read(un, '(a)') line
    nc = 1
    do j = 1, len_trim(line)
      if (line(j:j) == ',') nc = nc + 1
    end do
    rewind(un)
    allocate(x(nr, nc))
    do j = 1, nr
      read(un, *) x(j, :)
    end do
    close(un)
  end subroutine read_csv
end program bench
