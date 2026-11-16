! catalogue -- write the complete family catalogue of fortvine as Markdown (docs/families/catalogue.md)
! from the library itself, so the documentation always matches the code.
program catalogue
  use fortvine_api
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  character(len=24), parameter :: FORM(0:6) = [character(len=24) :: 'hand-coded', 'Archimedean', &
    'extreme-value', 'Archimax (1-par outer)', 'Archimax (2-par outer)', 'Khoudraji', 'nonparametric']
  integer, allocatable :: gen(:)
  integer :: k, u, f
  open(newunit=u, file='docs/families/catalogue.md', status='replace')
  write(u, '(a)') '# Complete family catalogue'
  write(u, '(a)') ''
  write(u, '(a)') 'Generated from the library (`tools/catalogue.f90`): every family code fortvine accepts, its'
  write(u, '(a)') 'name, number of parameters and parameter box. *Exch.* = exchangeable (C(u,v) = C(v,u));'
  write(u, '(a)') 'asymmetric families are transposed exactly when a vine edge is flipped. Every family except'
  write(u, '(a)') 'the rotationless ones also has 90, 180 and 270 degree rotations.'
  write(u, '(a)') ''
  write(u, '(a)') '## Hand-coded families (codes 0-25)'
  write(u, '(a)') ''
  write(u, '(a)') '| code | name | parameters | lower bounds | upper bounds | rotations |'
  write(u, '(a)') '|---:|---|---:|---|---|---|'
  do f = 0, 25
    call row_classic(f)
  end do
  gen = fv_generic_families()
  write(u, '(a)') ''
  write(u, '(a,i0,a)') '## Generic families (', size(gen), ' codes, `fv_generic_families()`)'
  write(u, '(a)') ''
  write(u, '(a)') 'Code = 10000 k + 100 og + ie (see [the generic engine](../manual/generic.md)).'
  do k = 1, 5
    write(u, '(a)') ''
    write(u, '(a,a)') '### ', trim(FORM(k))
    write(u, '(a)') ''
    write(u, '(a)') '| code | name | parameters | lower bounds | upper bounds | exch. | lower tail | upper tail |'
    write(u, '(a)') '|---:|---|---:|---|---|:-:|:-:|:-:|'
    do f = 1, size(gen)
      if (form_of(gen(f)) == k) call row_gen(gen(f))
    end do
  end do
  close(u)
  print '(a,i0,a)', 'wrote docs/families/catalogue.md (', 26 + size(gen), ' families)'
contains
  integer function form_of(c)
    integer, intent(in) :: c
    integer :: og, ie
    og = mod(c, 10000) / 100; ie = mod(c, 100)
    if (c >= 20000) then
      form_of = 5
    else if (ie == 0) then
      form_of = 1
    else if (og == 0) then
      form_of = 2
    else if (og <= 4) then
      form_of = 3
    else
      form_of = 4
    end if
  end function form_of
  character(len=60) function vec(x, np)
    real(dp), intent(in) :: x(3)
    integer, intent(in) :: np
    integer :: j
    vec = ''
    do j = 1, np
      vec = trim(vec) // ' ' // num(x(j))
    end do
    if (np == 0) vec = '-'
  end function vec
  ! a short, readable number: 0.01, 1, 2.5, 28, 1e-10
  function num(x) result(s)
    real(dp), intent(in) :: x
    character(len=16) :: s
    if (x == 0.0_dp) then
      s = '0'
    else if (abs(x) < 1.0e-3_dp .or. abs(x) >= 1.0e4_dp) then
      write(s, '(es8.1)') x
    else
      write(s, '(f12.4)') x
      s = adjustl(s)
      do while (s(len_trim(s):len_trim(s)) == '0')
        s(len_trim(s):len_trim(s)) = ' '
      end do
      if (s(len_trim(s):len_trim(s)) == '.') s(len_trim(s):len_trim(s)) = ' '
    end if
    s = adjustl(s)
  end function num
  ! escape | so a name like Archimax(PS|HuslerReiss) stays in one table cell
  function esc(s) result(r)
    character(len=*), intent(in) :: s
    character(len=:), allocatable :: r
    integer :: j
    r = ''
    do j = 1, len_trim(s)
      if (s(j:j) == '|') then
        r = r // '\|'
      else
        r = r // s(j:j)
      end if
    end do
  end function esc
  subroutine row_classic(c)
    integer, intent(in) :: c
    type(fv_bicop) :: b
    real(dp) :: lo(3), hi(3)
    integer :: np
    character(len=8) :: rot
    b%family = c; np = b%npars()
    call fv_par_bounds(c, lo, hi)
    rot = 'yes'
    if (any(c == [0, 1, 2, 5, 12, 13, 14, 24, 25])) rot = 'no'
    if (c == 24 .or. c == 25) then
      write(u, '(a,i0,a,a,a)') '| ', c, ' | ', esc(fv_fam_name(c)), ' | (smoothing) | - | - | no |'
      return
    end if
    write(u, '(a,i0,a,a,a,i0,a,a,a,a,a,a,a)') '| ', c, ' | ', esc(fv_fam_name(c)), ' | ', np, ' | ', &
      trim(adjustl(vec(lo, np))), ' | ', trim(adjustl(vec(hi, np))), ' | ', trim(rot), ' |'
  end subroutine row_classic
  subroutine row_gen(c)
    integer, intent(in) :: c
    real(dp) :: lo(3), hi(3)
    integer :: np
    np = gen_npar(c)
    call fv_par_bounds(c, lo, hi)
    write(u, '(a,i0,a,a,a,i0,a,a,a,a,a,a,a,a,a,a,a)') '| ', c, ' | ', esc(fv_fam_name(c)), ' | ', np, ' | ', &
      trim(adjustl(vec(lo, np))), ' | ', trim(adjustl(vec(hi, np))), ' | ', merge('yes', 'no ', gen_exch(c)), &
      ' | ', merge('yes', 'no ', gen_lt(c)), ' | ', merge('yes', 'no ', gen_ut(c)), ' |'
  end subroutine row_gen
end program catalogue
