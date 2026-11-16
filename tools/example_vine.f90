! example_vine -- data for the figures of docs/examples/vine.md (tools/example_figures.py vine).
! Re-runs example/ex02_vine.f90 (same truth, seed 42, same controls) and writes, into the
! current directory:
!   vine_truth.txt, vine_select.txt, vine_fixed.txt   one line per edge of the true D-vine, the
!       vine selected by Dissmann's algorithm and the vine selected on the fixed D-vine 1-2-3-4-5:
!       tree edge a b fam rot p1 p2 p3 tau ncond cond(1..ncond)   (a, b conditioned; original labels)
!   vine_taumat.txt   the 5 x 5 empirical Kendall tau matrix of the sample
! stdout: log-likelihood and number of parameters of each vine.
program example_vine
  use fortvine_api
  implicit none
  integer, parameter :: dp = kind(1.0d0), n = 3000, d = 5
  type(fv_vinecop) :: truth, vc, vf
  type(fv_controls) :: ctl
  real(dp) :: u(n, d), tm(d, d)
  integer :: e, un

  call truth%set_dvine([1, 2, 3, 4, 5])
  do e = 1, d - 1
    truth%fam(e, 1) = FV_GUMBEL; truth%par(1, e, 1) = 2.0_dp
  end do
  do e = 1, d - 2
    truth%fam(e, 2) = FV_FRANK; truth%par(1, e, 2) = 3.0_dp
  end do
  call truth%simulate(n, u, seed=42)
  ctl%nthreads = 4
  call vc%select(u, ctl)
  call vf%set_dvine([1, 2, 3, 4, 5])
  call vf%select(u, ctl)

  call dump(truth, 'vine_truth.txt')
  call dump(vc, 'vine_select.txt')
  call dump(vf, 'vine_fixed.txt')
  print '(a,f12.3,1x,i0)', 'truth ', truth%loglik(u), truth%npars()
  print '(a,f12.3,1x,i0)', 'select', vc%loglik(u), vc%npars()
  print '(a,f12.3,1x,i0)', 'fixed ', vf%loglik(u), vf%npars()

  call fv_tau_mat(u, tm, 4)
  open(newunit=un, file='vine_taumat.txt', status='replace')
  do e = 1, d
    write(un, '(*(1x,f10.6))') tm(e, :)
  end do
  close(un)

contains

  subroutine dump(v, fn)
    type(fv_vinecop), intent(in) :: v
    character(*), intent(in) :: fn
    integer :: m(v%d, v%d), t, e, uo
    type(fv_bicop) :: b
    m = v%matrix()
    open(newunit=uo, file=fn, status='replace')
    do t = 1, v%d - 1
      do e = 1, v%d - t
        b%family = v%fam(e, t); b%rotation = v%rot(e, t); b%par = v%par(:, e, t)
        write(uo, '(i3,i3,i3,i3,1x,i6,i5,3es16.7,f10.5,i3,*(i3))') t, e, m(v%d - e + 1, e), m(t, e), &
          b%family, b%rotation, b%par, ktau(b, 300), t - 1, m(1:t - 1, e)
      end do
    end do
    close(uo)
  end subroutine dump

  ! Kendall's tau of a pair copula, 1 - 4 int h1 h2 du dv (midpoint rule)
  real(dp) function ktau(c, k)
    type(fv_bicop), intent(in) :: c
    integer, intent(in) :: k
    real(dp) :: a(k*k), bb(k*k), h1(k*k), h2(k*k)
    integer :: ii, jj
    if (c%family == FV_INDEP) then
      ktau = 0.0_dp; return
    end if
    do jj = 1, k
      do ii = 1, k
        a(ii + (jj - 1) * k) = (ii - 0.5_dp) / k; bb(ii + (jj - 1) * k) = (jj - 0.5_dp) / k
      end do
    end do
    call c%hfunc1(a, bb, h1); call c%hfunc2(a, bb, h2)
    ktau = 1.0_dp - 4.0_dp * sum(h1 * h2) / real(k, dp)**2
  end function ktau
end program example_vine
