! fortvine_r.f90 -- the two pieces of library metadata the R package needs and the C interface
! of src/fortvine.f90 does not export: the family table (codes, names, numbers of parameters)
! and the D-vine structure array. Both are taken from fortvine_api; nothing is computed here.
! This file belongs to the R package only (it is not a copy of src/).
module fortvine_r
  use iso_c_binding, only: c_int, c_char
  use fortvine_api, only: fv_bicop, fv_vinecop, fv_fam_name, fv_generic_families, fam_is_np, FV_NCLASSIC
  implicit none
  private
contains

  ! codes(1:n), npars(1:n) (-1 for the nonparametric tll / bernstein), names(48, 1:n);
  ! n = -(number of families) if nmax is too small
  subroutine fvr_families(nmax, n, codes, npars, names) bind(C, name="fvr_families")
    integer(c_int), value :: nmax
    integer(c_int), intent(out) :: n, codes(nmax), npars(nmax)
    character(kind=c_char), intent(out) :: names(48, nmax)
    type(fv_bicop) :: b
    character(len=48) :: nm
    integer :: i, k
    associate (g => fv_generic_families())
    n = FV_NCLASSIC + size(g)
    if (n > nmax) then
      n = -n; return
    end if
    do i = 1, n
      if (i <= FV_NCLASSIC) then
        codes(i) = i - 1
      else
        codes(i) = g(i - FV_NCLASSIC)
      end if
      b%family = codes(i)
      npars(i) = b%npars()
      if (fam_is_np(codes(i))) npars(i) = -1
      nm = fv_fam_name(codes(i))
      do k = 1, 48
        names(k, i) = nm(k:k)
      end do
    end do
    end associate
  end subroutine fvr_families

  ! natural-order structure array st(edge, tree) of the D-vine fv_vinecop%set_dvine builds
  subroutine fvr_set_dvine(d, order, trunc, st) bind(C, name="fvr_set_dvine")
    integer(c_int), value :: d, trunc
    integer(c_int), intent(in) :: order(d)
    integer(c_int), intent(out) :: st(d, d)
    type(fv_vinecop) :: vc
    call vc%set_dvine(order, trunc)
    st = vc%st
  end subroutine fvr_set_dvine
end module fortvine_r
