! vine_sample -- simulate the 4-dimensional D-vine pictured on the docs home page and write it as CSV
! (columns u1..u4, no header) to the file named on the command line.
program vine_sample
  use fortvine_api
  implicit none
  integer, parameter :: dp = kind(1.0d0), n = 1500
  type(fv_vinecop) :: v
  real(dp) :: u(n, 4)
  character(len=256) :: fn
  integer :: i, un
  call get_command_argument(1, fn)
  call v%set_dvine([1, 2, 3, 4])
  ! tree 1: Clayton (lower tail), Gumbel (upper tail), Tawn (asymmetric); tree 2: Frank, Gaussian; tree 3: Student t
  v%fam(1:3, 1) = [FV_CLAYTON, FV_GUMBEL, FV_TAWN]
  v%par(:, 1, 1) = [3.0_dp, 0.0_dp, 0.0_dp]
  v%par(:, 2, 1) = [2.5_dp, 0.0_dp, 0.0_dp]
  v%par(:, 3, 1) = [0.4_dp, 1.0_dp, 4.0_dp]
  v%fam(1:2, 2) = [FV_FRANK, FV_GAUSSIAN]
  v%par(:, 1, 2) = [-4.0_dp, 0.0_dp, 0.0_dp]
  v%par(:, 2, 2) = [0.4_dp, 0.0_dp, 0.0_dp]
  v%fam(1, 3) = FV_STUDENT; v%par(:, 1, 3) = [0.3_dp, 4.0_dp, 0.0_dp]
  call v%simulate(n, u, seed=7)
  open(newunit=un, file=trim(fn), status='replace')
  do i = 1, n
    write(un, '(3(es22.15,","),es22.15)') u(i, :)
  end do
  close(un)
end program vine_sample
