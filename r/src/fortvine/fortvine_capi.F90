! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Shahab Aldin Shojaeezadeh. Part of fortvine (see LICENSE). If you use it, please cite:
! Shojaeezadeh, S. A. (2026). fortvine: a complete compendium to work on vine copula modeling.
! fortvine_capi.f90 -- C interface to the typed API (fortvine_api), used by the Python and R
! packages.  It adds what the fv_* routines of fortvine.f90 do not pass: observation weights,
! discrete margins, the nonparametric families (their estimate `aux`), every fv_controls field
! (select_trunc, select_threshold, fast_mle, ...), non-simplified vines (relax), pdf_grid, the
! fit cache, the abort flag, parameter bounds, family metadata and D-/C-vine structures.
! Each routine only converts arguments and forwards to fortvine_api; nothing is computed here.
!
! Conventions: arrays column-major, observations last (u(d, n), as in the fv_* routines);
! logicals as 0/1 ints.  Controls are passed as
!   ci(12) = selcrit, method, preselect, allow_rotations, fast_mle, tau_mode (= ci(6) of the fv_*
!            routines), tree_crit, tree_algo, trunc_lvl (= ci(7:9) of fv_vine_select), nthreads,
!            select_trunc, select_threshold
!   cr(2)  = threshold, psi0
! with the family set famset(nf).  A vine is passed as its arrays: d, trunc, order(d),
! st/fam/rot(d, d) (edge, tree), par(3, d, d), aux(naux, d, d) (naux = fvc_naux() when an edge is
! tll / bernstein, else 1), disc(d), and the relaxed edges nrel, erel(2, *), wa(d, *),
! wk(nk, *), pk(3, nk, *) (nk = fvc_nk(); arrays sized max(nrel, 1)).
module fortvine_capi
  use iso_c_binding, only: c_int, c_double, c_char, c_int64_t
  use iso_fortran_env, only: int64
  use fortvine, only: fv_abort
  use fortvine_api, only: fv_bicop, fv_vinecop, fv_controls, fv_fam_name, fv_generic_families, fam_is_np, &
    fv_par_bounds, gen_flip, gen_is, NAUX, FV_NCLASSIC, FV_NFAM, fv_fitcache_enable, fv_fitcache_clear, &
    fv_fitcache_disable, fv_fitcache_stats
  implicit none
  private
  integer, parameter :: dp = kind(1.0d0), NK = 5   ! NK: anchors of a relaxed edge (vc_relax)
contains

  ! ------------------------------------------------------------ metadata
  integer(c_int) function fvc_naux() bind(C, name="fvc_naux")
    fvc_naux = NAUX
  end function fvc_naux

  integer(c_int) function fvc_nk() bind(C, name="fvc_nk")
    fvc_nk = NK
  end function fvc_nk

  ! fv_fam_name(fam), blank padded to 48 characters ('invalid' for an unknown code)
  subroutine fvc_fam_name(fam, buf) bind(C, name="fvc_fam_name")
    integer(c_int), value :: fam
    character(kind=c_char), intent(out) :: buf(48)
    character(len=48) :: nm
    integer :: i
    nm = fv_fam_name(fam)
    do i = 1, 48
      buf(i) = nm(i:i)
    end do
  end subroutine fvc_fam_name

  ! the distinct families (codes 0..25, then fv_generic_families()) with their numbers of
  ! parameters (-1: tll / bernstein); returns their number and writes at most nmax of them
  integer(c_int) function fvc_families(nmax, codes, npars) bind(C, name="fvc_families")
    integer(c_int), value :: nmax
    integer(c_int), intent(out) :: codes(nmax), npars(nmax)
    integer, allocatable :: g(:)
    type(fv_bicop) :: b
    integer :: j
    g = [(j, j = 0, FV_NCLASSIC - 1), fv_generic_families()]
    do j = 1, min(int(nmax), size(g))
      codes(j) = g(j)
      b%family = g(j)
      npars(j) = b%npars()
      if (fam_is_np(g(j))) npars(j) = -1
    end do
    fvc_families = size(g)
  end function fvc_families

  ! fv_par_bounds: the parameter box of a family
  subroutine fvc_par_bounds(fam, lb, ub) bind(C, name="fvc_par_bounds")
    integer(c_int), value :: fam
    real(c_double), intent(out) :: lb(3), ub(3)
    call fv_par_bounds(fam, lb, ub)
  end subroutine fvc_par_bounds

  ! gen_flip: the generic family and parameters of the transposed copula C(u2, u1)
  ! (generic codes only; other codes are returned unchanged)
  subroutine fvc_gen_flip(fam, par) bind(C, name="fvc_gen_flip")
    integer(c_int), intent(inout) :: fam
    real(c_double), intent(inout) :: par(3)
    integer :: f
    f = fam
    if (gen_is(f)) call gen_flip(f, par)
    fam = f
  end subroutine fvc_gen_flip

  ! ------------------------------------------------------------ run-time switches
  subroutine fvc_set_abort(flag) bind(C, name="fvc_set_abort")
    integer(c_int), value :: flag
    fv_abort = flag /= 0
  end subroutine fvc_set_abort

  subroutine fvc_fitcache(op, capacity, hits, misses, used) bind(C, name="fvc_fitcache")
    integer(c_int), value :: op, capacity          ! op: 0 disable, 1 enable(capacity), 2 clear, 3 stats
    integer(c_int64_t), intent(out) :: hits, misses
    integer(c_int), intent(out) :: used
    integer(int64) :: h, m
    integer :: us
    select case (op)
    case (0); call fv_fitcache_disable()
    case (1); call fv_fitcache_enable(int(capacity))
    case (2); call fv_fitcache_clear()
    end select
    call fv_fitcache_stats(h, m, us)
    hits = h; misses = m; used = us
  end subroutine fvc_fitcache

  ! ------------------------------------------------------------ pair copula
  ! what: 0 log density, 1 density, 2 h1, 3 h2, 4 h1 inverse, 5 h2 inverse, 6 distribution function
  subroutine fvc_bicop_eval(fam, rot, par, naux, aux, n, u, what, out) bind(C, name="fvc_bicop_eval")
    integer(c_int), value :: fam, rot, naux, n, what
    real(c_double), intent(in) :: par(3), aux(naux), u(2, n)
    real(c_double), intent(out) :: out(n)
    type(fv_bicop) :: b
    b%family = fam; b%rotation = rot; b%par = par; b%aux = aux
    if (what == 6) then
      call b%cdf(u(1, :), u(2, :), out)
    else
      call bc_eval_c(b, u, what, out)
    end if
  end subroutine fvc_bicop_eval

  subroutine bc_eval_c(b, u, what, out)
    type(fv_bicop), intent(in) :: b
    real(dp), intent(in) :: u(:, :)
    integer, intent(in) :: what
    real(dp), intent(out) :: out(:)
    select case (what)
    case (0); call b%logpdf(u(1, :), u(2, :), out)
    case (1); call b%pdf(u(1, :), u(2, :), out)
    case (2); call b%hfunc1(u(1, :), u(2, :), out)
    case (3); call b%hfunc2(u(1, :), u(2, :), out)
    case (4); call b%hinv1(u(1, :), u(2, :), out)
    case default; call b%hinv2(u(1, :), u(2, :), out)
    end select
  end subroutine bc_eval_c

  ! fit family fam / rotation rot (sel = 0) or select among famset (sel = 1) on u(2, n), with
  ! weights w(n) when usew = 1; aux(naux) receives a tll / bernstein estimate (naux = fvc_naux()
  ! then), npe the number of parameters (effective degrees of freedom for tll / bernstein)
  subroutine fvc_bicop_fit(n, u, usew, w, sel, nf, famset, ci, cr, fam, rot, par, ll, naux, aux, npe) &
      bind(C, name="fvc_bicop_fit")
    integer(c_int), value :: n, usew, sel, nf, naux
    real(c_double), intent(in) :: u(2, n), w(n), cr(2)
    integer(c_int), intent(in) :: famset(nf), ci(12)
    integer(c_int), intent(inout) :: fam, rot
    real(c_double), intent(out) :: par(3), ll, aux(naux), npe
    type(fv_bicop) :: b
    type(fv_controls) :: c
    c = ctl_of(nf, famset, ci, cr)
    b%family = fam; b%rotation = rot
    if (sel == 0) then
      if (usew /= 0) then
        call b%fit(u(1, :), u(2, :), c, w)
      else
        call b%fit(u(1, :), u(2, :), c)
      end if
    else
      if (usew /= 0) then
        call b%select(u(1, :), u(2, :), c, w)
      else
        call b%select(u(1, :), u(2, :), c)
      end if
    end if
    fam = b%family; rot = b%rotation; par = b%par; ll = b%loglik
    aux = 0.0_dp
    aux(1:min(naux, size(b%aux))) = b%aux(1:min(naux, size(b%aux)))
    npe = b%npars_eff()
  end subroutine fvc_bicop_fit

  ! ------------------------------------------------------------ vine copula
  ! Dissmann selection (fv_vinecop%select) on u(d, n): weights w(n) if usew = 1, discrete
  ! variables disc(d) with um(d, n) = F(x-), a fixed structure (order, st, trunc on input) if
  ! fixed = 1.  Returns the vine arrays, its log-likelihood ll, number of parameters npe and the
  ! threshold used; info = 0 ok, -2 if naux is too small for a tll / bernstein family set.
  subroutine fvc_vine_select(n, d, u, um, disc, usew, w, fixed, nf, famset, ci, cr, trunc, order, st, fam, rot, &
                             par, edge_ll, naux, aux, ll, npe, thr, info) bind(C, name="fvc_vine_select")
    integer(c_int), value :: n, d, usew, fixed, nf, naux
    real(c_double), intent(in) :: u(d, n), um(d, n), w(n), cr(2)
    integer(c_int), intent(in) :: disc(d), famset(nf), ci(12)
    integer(c_int), intent(inout) :: trunc, order(d), st(d, d)
    integer(c_int), intent(out) :: fam(d, d), rot(d, d), info
    real(c_double), intent(out) :: par(3, d, d), edge_ll(d, d), aux(naux, d, d), ll, npe, thr
    type(fv_vinecop) :: v
    type(fv_controls) :: c
    real(dp), allocatable :: ut(:, :), umt(:, :)
    logical :: dsc(d)
    c = ctl_of(nf, famset, ci, cr)
    if (naux < NAUX .and. any(fam_is_np(c%family_set(1:c%nfam)))) then
      info = -2; return
    end if
    call v%init(d)
    if (fixed /= 0) then
      v%order = order; v%st = st; v%trunc = trunc
    end if
    dsc = disc /= 0
    ut = transpose(u)
    umt = ut
    if (any(dsc)) umt = transpose(um)
    if (usew /= 0) then
      call v%select(ut, c, fixed=fixed /= 0, info=info, weights=w, discrete=dsc, u_minus=umt)
    else
      call v%select(ut, c, fixed=fixed /= 0, info=info, discrete=dsc, u_minus=umt)
    end if
    if (info /= 0) return
    trunc = v%trunc; order = v%order; st = v%st; fam = v%fam; rot = v%rot; par = v%par; edge_ll = v%edge_ll
    aux = 0.0_dp
    if (size(v%aux, 1) > 1) aux(1:min(naux, size(v%aux, 1)), :, :) = v%aux(1:min(naux, size(v%aux, 1)), :, :)
    ll = sum(v%edge_ll(:, 1:v%trunc))
    npe = v%npars_eff()
    thr = v%threshold
  end subroutine fvc_vine_select

  ! what: 0 density out(1, :), 1 Rosenblatt transform (with discrete variables: randomised with
  ! seed, or the conditional CDF F(x|...) if seed = 0), 2 inverse Rosenblatt, 3 n draws (u
  ! unused; reproducible for seed > 0).  u, um(d, n) as in fvc_vine_select; out(d, n).
  subroutine fvc_vine_eval(d, trunc, order, st, fam, rot, par, naux, aux, disc, nrel, erel, wa, wk, pk, &
                           what, n, u, um, seed, nthreads, out) bind(C, name="fvc_vine_eval")
    integer(c_int), value :: d, trunc, naux, nrel, what, n, seed, nthreads
    integer(c_int), intent(in) :: order(d), st(d, d), fam(d, d), rot(d, d), disc(d), erel(2, max(nrel, 1))
    real(c_double), intent(in) :: par(3, d, d), aux(naux, d, d), wa(d, max(nrel, 1)), wk(NK, max(nrel, 1)), &
      pk(3, NK, max(nrel, 1)), u(d, n), um(d, n)
    real(c_double), intent(out) :: out(d, n)
    type(fv_vinecop) :: v
    real(dp), allocatable :: o(:, :), p(:)
    call load(v, d, trunc, order, st, fam, rot, par, naux, aux, disc, nrel, erel, wa, wk, pk)
    allocate(o(n, d))
    select case (what)
    case (0)
      allocate(p(n))
      if (any(v%disc)) then
        call v%pdf(transpose(u), p, nthreads, transpose(um))
      else
        call v%pdf(transpose(u), p, nthreads)
      end if
      out(1, :) = p
      return
    case (1)
      if (any(v%disc)) then
        call v%rosenblatt(transpose(u), o, nthreads, transpose(um), randomize=seed /= 0, seed=seed)
      else
        call v%rosenblatt(transpose(u), o, nthreads)
      end if
    case (2)
      call v%inverse_rosenblatt(transpose(u), o, nthreads)
    case default
      if (seed > 0) then
        call v%simulate(n, o, seed, nthreads)
      else
        call v%simulate(n, o, nthreads=nthreads)
      end if
    end select
    out = transpose(o)
  end subroutine fvc_vine_eval

  ! fv_vinecop%pdf_grid: ug(d, nb, n) varies only variable tgt (1-based) across the nb bins;
  ! ugm = F(x-) on the grid (discrete variables); out(nb, n)
  subroutine fvc_vine_pdf_grid(d, trunc, order, st, fam, rot, par, naux, aux, disc, nrel, erel, wa, wk, pk, &
                               nb, n, ug, ugm, tgt, nthreads, out) bind(C, name="fvc_vine_pdf_grid")
    integer(c_int), value :: d, trunc, naux, nrel, nb, n, tgt, nthreads
    integer(c_int), intent(in) :: order(d), st(d, d), fam(d, d), rot(d, d), disc(d), erel(2, max(nrel, 1))
    real(c_double), intent(in) :: par(3, d, d), aux(naux, d, d), wa(d, max(nrel, 1)), wk(NK, max(nrel, 1)), &
      pk(3, NK, max(nrel, 1)), ug(d, nb, n), ugm(d, nb, n)
    real(c_double), intent(out) :: out(nb, n)
    type(fv_vinecop) :: v
    call load(v, d, trunc, order, st, fam, rot, par, naux, aux, disc, nrel, erel, wa, wk, pk)
    call v%pdf_grid(ug, tgt, out, nthreads, ugm)
  end subroutine fvc_vine_pdf_grid

  ! fv_vinecop%relax (non-simplified edges) on u(d, n); the vine arrays are updated (later trees
  ! are refitted) and the relaxed edges written to erel, wa, wk, pk (room for nmax edges);
  ! nrel = their number (-1 if nmax is too small).  tgt, mode as in fv_vinecop%relax.
  subroutine fvc_vine_relax(n, d, u, usew, w, tgt, mode, nf, famset, ci, cr, trunc, order, st, fam, rot, par, &
                            edge_ll, naux, aux, disc, nmax, nrel, erel, wa, wk, pk, ll, npe) bind(C, name="fvc_vine_relax")
    integer(c_int), value :: n, d, usew, tgt, mode, nf, naux, nmax
    real(c_double), intent(in) :: u(d, n), w(n), cr(2)
    integer(c_int), intent(in) :: famset(nf), ci(12), trunc, order(d), st(d, d), disc(d)
    integer(c_int), intent(inout) :: fam(d, d), rot(d, d)
    real(c_double), intent(inout) :: par(3, d, d), edge_ll(d, d), aux(naux, d, d)
    integer(c_int), intent(out) :: nrel, erel(2, nmax)
    real(c_double), intent(out) :: wa(d, nmax), wk(NK, nmax), pk(3, NK, nmax), ll, npe
    type(fv_vinecop) :: v
    integer :: e1(2, 1), nr
    real(dp) :: w1(d, 1), k1(NK, 1), p1(3, NK, 1)
    e1 = 0; w1 = 0; k1 = 0; p1 = 0
    call load(v, d, trunc, order, st, fam, rot, par, naux, aux, disc, 0, e1, w1, k1, p1)
    v%edge_ll = edge_ll
    if (usew /= 0) then
      call v%relax(transpose(u), tgt, mode, ctl_of(nf, famset, ci, cr), w, nr)
    else
      call v%relax(transpose(u), tgt, mode, ctl_of(nf, famset, ci, cr), nrel_out=nr)
    end if
    fam = v%fam; rot = v%rot; par = v%par; edge_ll = v%edge_ll
    if (size(v%aux, 1) == naux) aux = v%aux
    ll = sum(v%edge_ll(:, 1:v%trunc)); npe = v%npars_eff()
    nrel = v%nrel
    if (v%nrel > nmax) then
      nrel = -1; return
    end if
    if (v%nrel > 0) then
      erel(:, 1:v%nrel) = v%erel; wa(:, 1:v%nrel) = v%wa; wk(:, 1:v%nrel) = v%wk; pk(:, :, 1:v%nrel) = v%pk
    end if
  end subroutine fvc_vine_relax

  ! fv_vinecop%rel_params: parameters prel(3, nrel, n) of the relaxed edges at the samples uc(d, n)
  ! (original order); the relaxed edges as returned by fvc_vine_relax
  subroutine fvc_vine_rel_params(d, nrel, wa, wk, pk, n, uc, prel) bind(C, name="fvc_vine_rel_params")
    integer(c_int), value :: d, nrel, n
    real(c_double), intent(in) :: wa(d, max(nrel, 1)), wk(NK, max(nrel, 1)), pk(3, NK, max(nrel, 1)), uc(d, n)
    real(c_double), intent(out) :: prel(3, max(nrel, 1), n)
    type(fv_vinecop) :: v
    if (nrel <= 0) return
    v%d = d; v%nrel = nrel
    v%wa = wa(:, 1:nrel); v%wk = wk(:, 1:nrel); v%pk = pk(:, :, 1:nrel)
    prel(:, 1:nrel, :) = v%rel_params(uc)
  end subroutine fvc_vine_rel_params

  ! fv_vinecop%set_dvine (kind 0) or %set_cvine (kind 1): st(d, d) and the truncation level used
  subroutine fvc_vine_structure(kind, d, order, trunc, st, trunc_out) bind(C, name="fvc_vine_structure")
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
  end subroutine fvc_vine_structure

  ! fv_vinecop%matrix(): the R-vine matrix m(d, d) in original labels
  subroutine fvc_vine_matrix(d, order, trunc, st, m) bind(C, name="fvc_vine_matrix")
    integer(c_int), value :: d, trunc
    integer(c_int), intent(in) :: order(d), st(d, d)
    integer(c_int), intent(out) :: m(d, d)
    type(fv_vinecop) :: v
    call v%init(d)
    v%order = order; v%st = st; v%trunc = trunc
    m = v%matrix()
  end subroutine fvc_vine_matrix

  ! ------------------------------------------------------------ helpers
  function ctl_of(nf, famset, ci, cr) result(c)
    integer, intent(in) :: nf, famset(:), ci(12)
    real(dp), intent(in) :: cr(2)
    type(fv_controls) :: c
    c%nfam = min(nf, FV_NFAM)
    c%family_set(1:c%nfam) = famset(1:c%nfam)
    c%selcrit = ci(1); c%method = ci(2); c%preselect = ci(3) /= 0; c%allow_rotations = ci(4) /= 0
    c%fast_mle = ci(5) /= 0; c%tau_mode = ci(6)
    c%tree_crit = ci(7); c%tree_algo = ci(8); c%trunc_lvl = ci(9); c%nthreads = ci(10)
    c%select_trunc = ci(11) /= 0; c%select_threshold = ci(12) /= 0
    c%threshold = cr(1); c%psi0 = cr(2)
  end function ctl_of

  subroutine load(v, d, trunc, order, st, fam, rot, par, naux, aux, disc, nrel, erel, wa, wk, pk)
    type(fv_vinecop), intent(out) :: v
    integer, intent(in) :: d, trunc, order(d), st(d, d), fam(d, d), rot(d, d), naux, disc(d), nrel, erel(:, :)
    real(dp), intent(in) :: par(3, d, d), aux(naux, d, d), wa(:, :), wk(:, :), pk(:, :, :)
    call v%init(d)
    v%trunc = trunc; v%order = order; v%st = st; v%fam = fam; v%rot = rot; v%par = par
    v%aux = aux                      ! (reallocated to naux rows)
    v%disc = disc /= 0
    if (nrel > 0) then               ! erel stays unallocated for a simplified vine
      v%nrel = nrel
      v%erel = erel(:, 1:nrel); v%wa = wa(:, 1:nrel); v%wk = wk(:, 1:nrel); v%pk = pk(:, :, 1:nrel)
    end if
  end subroutine load

end module fortvine_capi
