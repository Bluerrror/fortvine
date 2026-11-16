! fortvine_py.f90 -- bind(C) wrappers used by the Python package for the parts of the
! typed Fortran API (fortvine_api) that src/'s C interface does not export: family names,
! parameter counts and the family list, the pair-copula distribution function, D-/C-vine
! structures and the R-vine matrix.  Each routine only forwards to fortvine_api.
module fortvine_py
  use iso_c_binding, only: c_int, c_double, c_char
  use fortvine_api, only: fv_bicop, fv_vinecop, fv_fam_name, fv_generic_families, FV_NCLASSIC
  implicit none
contains

  ! fv_fam_name(fam), blank padded to 48 characters ('invalid' for an unknown code)
  subroutine fvpy_fam_name(fam, buf) bind(C, name="fvpy_fam_name")
    integer(c_int), value :: fam
    character(kind=c_char), intent(out) :: buf(48)
    character(len=48) :: nm
    integer :: i
    nm = fv_fam_name(fam)
    do i = 1, 48
      buf(i) = nm(i:i)
    end do
  end subroutine fvpy_fam_name

  ! fv_bicop%npars() of a family (0 for indep, tll, bernstein)
  integer(c_int) function fvpy_npars(fam) bind(C, name="fvpy_npars")
    integer(c_int), value :: fam
    type(fv_bicop) :: b
    b%family = fam
    fvpy_npars = b%npars()
  end function fvpy_npars

  ! the distinct families: codes 0..25, then fv_generic_families(); returns their number and
  ! writes at most nmax codes
  integer(c_int) function fvpy_families(nmax, codes) bind(C, name="fvpy_families")
    integer(c_int), value :: nmax
    integer(c_int), intent(out) :: codes(nmax)
    integer, allocatable :: g(:)
    integer :: j, m
    g = [(j, j = 0, FV_NCLASSIC - 1), fv_generic_families()]
    m = min(int(nmax), size(g))
    codes(1:m) = g(1:m)
    fvpy_families = size(g)
  end function fvpy_families

  ! fv_bicop%cdf on u(2, n)
  subroutine fvpy_bicop_cdf(fam, rot, par, n, u, out) bind(C, name="fvpy_bicop_cdf")
    integer(c_int), value :: fam, rot, n
    real(c_double), intent(in) :: par(3), u(2, n)
    real(c_double), intent(out) :: out(n)
    type(fv_bicop) :: b
    b%family = fam; b%rotation = rot; b%par = par
    call b%cdf(u(1, :), u(2, :), out)
  end subroutine fvpy_bicop_cdf

  ! fv_vinecop%set_dvine (kind 0) or %set_cvine (kind 1): struct array st(d, d) and the
  ! truncation level actually used
  subroutine fvpy_vine_structure(kind, d, order, trunc, st, trunc_out) bind(C, name="fvpy_vine_structure")
    integer(c_int), value :: kind, d, trunc
    integer(c_int), intent(in) :: order(d)
    integer(c_int), intent(out) :: st(d, d), trunc_out
    type(fv_vinecop) :: v
    if (kind == 0) then
      call v%set_dvine(order, trunc)
    else
      call v%set_cvine(order, trunc)
    end if
    st = v%st
    trunc_out = v%trunc
  end subroutine fvpy_vine_structure

  ! fv_vinecop%matrix(): the R-vine matrix m(d, d) in original labels
  subroutine fvpy_vine_matrix(d, order, trunc, st, m) bind(C, name="fvpy_vine_matrix")
    integer(c_int), value :: d, trunc
    integer(c_int), intent(in) :: order(d), st(d, d)
    integer(c_int), intent(out) :: m(d, d)
    type(fv_vinecop) :: v
    call v%init(d)
    v%order = order; v%st = st; v%trunc = trunc
    m = v%matrix()
  end subroutine fvpy_vine_matrix

end module fortvine_py
