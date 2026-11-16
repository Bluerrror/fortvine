# `fv_vinecop`

A regular vine copula. Manual: [Vine copulas](../manual/vinecop.md).

```fortran
type, public :: fv_vinecop
  integer :: d = 0
  integer :: trunc = 0
  integer,  allocatable :: order(:)          ! structure order (1-based labels)
  integer,  allocatable :: st(:, :)          ! natural-order structure array (edge, tree)
  integer,  allocatable :: fam(:, :), rot(:, :)
  real(dp), allocatable :: par(:, :, :)      ! (3, edge, tree)
  real(dp), allocatable :: edge_ll(:, :)
  real(dp), allocatable :: aux(:, :, :)      ! (NAUX or 1, edge, tree)
  logical,  allocatable :: disc(:)           ! discrete variables (original labels)
  real(dp) :: threshold = 0.0_dp
  integer :: nrel = 0                        ! non-simplified (relaxed) edges
  integer,  allocatable :: erel(:, :)        ! (2, nrel): edge, tree
  real(dp), allocatable :: wa(:, :)          ! (d, nrel) loadings of w
  real(dp), allocatable :: wk(:, :)          ! (NK, nrel) anchor values of w
  real(dp), allocatable :: pk(:, :, :)       ! (3, NK, nrel) parameters at the anchors
contains
  procedure :: init, set_dvine, set_cvine, select
  procedure :: pdf, pdf_grid, loglik, rosenblatt, inverse_rosenblatt, simulate
  procedure :: npars, npars_eff, matrix, print
  procedure :: relax, unrelax, rel_params
end type fv_vinecop
```

## Fields

| name | type | description |
|---|---|---|
| `d` | `integer` | dimension |
| `trunc` | `integer` | number of trees stored (0 before a structure is set) |
| `order` | `integer(d)` | structure order |
| `st` | `integer(d, d)` | `st(e, t)`: edge `e` of tree `t` joins `order(e)` and `order(st(e, t))` given `order(st(e, 1:t-1))` |
| `fam`, `rot` | `integer(d, d)` | family and rotation of edge `(e, t)` |
| `par` | `real(dp)(3, d, d)` | parameters of edge `(e, t)` in `par(:, e, t)` |
| `edge_ll` | `real(dp)(d, d)` | edge log-likelihoods from the last `select` (or `relax`) |
| `aux` | `real(dp)(:, d, d)` | nonparametric estimates; first extent `NAUX` if TLL / Bernstein were candidates, else 1 |
| `disc` | `logical(d)` | discrete flags of the last `select` |
| `threshold` | `real(dp)` | threshold of the last `select` |
| `nrel`, `erel`, `wa`, `wk`, `pk` | | relaxed edges, see [`relax`](#relax) (`NK = 5` anchors) |

---

## init

```fortran
subroutine init(self, d)
  class(fv_vinecop), intent(inout) :: self
  integer, intent(in) :: d
```

Allocates all arrays for dimension `d` and resets them: `order = [1, ..., d]`, `st`, `fam`,
`rot`, `par`, `edge_ll`, `aux` (first extent 1) to zero, `disc = .false.`, `trunc = 0`, no
relaxed edges.

## set_dvine, set_cvine

```fortran
subroutine set_dvine(self, order, trunc)
subroutine set_cvine(self, order, trunc)
  class(fv_vinecop), intent(inout) :: self
  integer, intent(in) :: order(:)
  integer, intent(in), optional :: trunc
```

| name | type | intent | description |
|---|---|---|---|
| `order` | `integer(:)` | in | permutation of `1..d`; `d = size(order)` |
| `trunc` | `integer` | in, optional | number of trees, clipped to `[0, d-1]`; default `d - 1` |

Call `init(d)`, then set a D-vine (`st(e, t) = e + t`: path `order(1)-order(2)-...`) or a
C-vine (`st(e, t) = d - t + 1`: the root of tree `t` is `order(d - t + 1)`). All edges are
independence until families are assigned or `select` is called.

```fortran
call vc%set_cvine([4, 3, 2, 1])      ! roots 1, 2, 3
```

## select

```fortran
subroutine select(self, u, ctl, fixed, info, weights, discrete, u_minus)
  class(fv_vinecop), intent(inout) :: self
  real(dp), intent(in) :: u(:, :)
  type(fv_controls), intent(in), optional :: ctl
  logical, intent(in), optional :: fixed
  integer, intent(out), optional :: info
  real(dp), intent(in), optional :: weights(:), u_minus(:, :)
  logical, intent(in), optional :: discrete(:)
```

| name | type | intent | description |
|---|---|---|---|
| `u` | `real(dp)(n, d)` | in | pseudo-observations; for discrete columns $F(x)$ |
| `ctl` | `type(fv_controls)` | in, optional | all fields are used; default `fv_controls()` |
| `fixed` | `logical` | in, optional | keep the stored structure; default `.true.` if `self%trunc > 0`, else `.false.` |
| `info` | `integer` | out, optional | 0 on success; 2, 3, 4 internal failure |
| `weights` | `real(dp)(n)` | in, optional | observation weights |
| `discrete` | `logical(d)` | in, optional | discrete columns |
| `u_minus` | `real(dp)(n, d)` | in, optional | $F(x^-)$ of the discrete columns; required if any `discrete` is true |

Sets `order`, `st`, `trunc`, `fam`, `rot`, `par`, `edge_ll`, `aux`, `disc`, `threshold`;
removes relaxed edges. If `self%d` differs from `size(u, 2)` the object is re-initialised.
With `fixed`, the number of trees is the smaller of `self%trunc` and `ctl%trunc_lvl`
(`d - 1` when `trunc_lvl = -1`), possibly reduced further by `select_trunc`.

```fortran
call vc%select(u, ctl, fixed=.false., info=ierr)
```

## pdf

```fortran
subroutine pdf(self, u, out, nthreads, u_minus)
  class(fv_vinecop), intent(in) :: self
  real(dp), intent(in)  :: u(:, :)
  real(dp), intent(out) :: out(:)
  integer, intent(in), optional :: nthreads
  real(dp), intent(in), optional :: u_minus(:, :)
```

| name | type | intent | description |
|---|---|---|---|
| `u` | `real(dp)(n, d)` | in | points |
| `out` | `real(dp)(n)` | out | joint copula density |
| `nthreads` | `integer` | in, optional | threads (default `fv_threads()`) |
| `u_minus` | `real(dp)(n, d)` | in, optional | required if the vine has discrete variables |

Relaxed edges are evaluated with their sample-dependent parameters.

## pdf_grid

```fortran
subroutine pdf_grid(self, ug, tgt, out, nthreads, ugm)
  class(fv_vinecop), intent(in) :: self
  real(dp), intent(in)  :: ug(:, :, :)
  integer, intent(in)   :: tgt
  real(dp), intent(out) :: out(:, :)
  integer, intent(in), optional :: nthreads
  real(dp), intent(in), optional :: ugm(:, :, :)
```

| name | type | intent | description |
|---|---|---|---|
| `ug` | `real(dp)(d, bins, n)` | in | for each of `n` samples, `bins` points in original variable order |
| `tgt` | `integer` | in | 1-based label of the only variable that varies across bins |
| `out` | `real(dp)(bins, n)` | out | joint density at every grid point |
| `nthreads` | `integer` | in, optional | threads |
| `ugm` | `real(dp)(d, bins, n)` | in, optional | $F(x^-)$ grid, required if the vine has discrete variables |

Note the layout: variables first, observations last (unlike `pdf`).

## loglik

```fortran
real(dp) function loglik(self, u, u_minus)
  class(fv_vinecop), intent(in) :: self
  real(dp), intent(in) :: u(:, :)
  real(dp), intent(in), optional :: u_minus(:, :)
```

$\sum_i \log c(u_i)$ over the rows of `u(n, d)`.

## rosenblatt

```fortran
subroutine rosenblatt(self, u, out, nthreads, u_minus, out_minus, randomize, seed)
  class(fv_vinecop), intent(in) :: self
  real(dp), intent(in)  :: u(:, :)
  real(dp), intent(out) :: out(:, :)
  integer, intent(in), optional :: nthreads, seed
  real(dp), intent(in), optional :: u_minus(:, :)
  real(dp), intent(out), optional :: out_minus(:, :)
  logical, intent(in), optional :: randomize
```

| name | type | intent | description |
|---|---|---|---|
| `u` | `real(dp)(n, d)` | in | points |
| `out` | `real(dp)(n, d)` | out | Rosenblatt transform, original column labels; column `order(d)` equals the input |
| `nthreads` | `integer` | in, optional | threads |
| `u_minus` | `real(dp)(n, d)` | in, optional | required for discrete vines (`error stop` otherwise) |
| `out_minus` | `real(dp)(n, d)` | out, optional | discrete vines: $F(x_k^- \mid \cdot)$ |
| `randomize` | `logical` | in, optional | discrete vines: Brockwell randomisation, default `.true.` |
| `seed` | `integer` | in, optional | seed of the randomisation, default 1 |

Relaxed edges are evaluated with `par` (their middle-anchor parameters).

## inverse_rosenblatt

```fortran
subroutine inverse_rosenblatt(self, u, out, nthreads)
  class(fv_vinecop), intent(in) :: self
  real(dp), intent(in)  :: u(:, :)
  real(dp), intent(out) :: out(:, :)
  integer, intent(in), optional :: nthreads
```

Maps independent uniforms `u(n, d)` to `out(n, d)`; the inverse of the continuous
`rosenblatt`.

## simulate

```fortran
subroutine simulate(self, n, out, seed, nthreads)
  class(fv_vinecop), intent(in) :: self
  integer, intent(in) :: n
  real(dp), intent(out) :: out(:, :)
  integer, intent(in), optional :: seed, nthreads
```

| name | type | intent | description |
|---|---|---|---|
| `n` | `integer` | in | number of draws |
| `out` | `real(dp)(n, d)` | out | sample of the vine copula |
| `seed` | `integer` | in, optional | reseeds the intrinsic `random_number` generator (affects the caller's generator state) |
| `nthreads` | `integer` | in, optional | threads |

## npars, npars_eff

```fortran
integer  function npars(self)
real(dp) function npars_eff(self)
```

Number of parameters over the trees `1..trunc`: effective degrees of freedom for TLL /
Bernstein edges, and `(NK - 1) * npars` additional parameters per relaxed edge. `npars` is the
rounded value.

## matrix

```fortran
function matrix(self) result(m)
  class(fv_vinecop), intent(in) :: self
  integer :: m(self%d, self%d)
```

The vinecopulib R-vine matrix in original labels: `m(t, e) = order(st(e, t))` for
`t <= trunc`, `m(d - j + 1, j) = order(j)`, zeros elsewhere.

## print

```fortran
subroutine print(self, unit)
  class(fv_vinecop), intent(in) :: self
  integer, intent(in), optional :: unit
```

Writes a summary to `unit` (default 6): `d`, `trunc`, `order`, then one line per edge with
family name (first 12 characters), rotation and parameters, then the relaxed edges.

## relax

```fortran
subroutine relax(self, u, tgt, mode, ctl, weights, nrel_out)
  class(fv_vinecop), intent(inout) :: self
  real(dp), intent(in) :: u(:, :)
  integer, intent(in) :: tgt, mode
  type(fv_controls), intent(in) :: ctl
  real(dp), intent(in), optional :: weights(:)
  integer, intent(out), optional :: nrel_out
```

| name | type | intent | description |
|---|---|---|---|
| `u` | `real(dp)(n, d)` | in | the data the vine was fitted to |
| `tgt` | `integer` | in | relax only edges whose conditioned pair contains this label and whose conditioning set does not; 0 for any |
| `mode` | `integer` | in | 1: only edges rejected by the constant-conditional-correlation test at 5 %; 2: all eligible edges; $\le 0$: nothing |
| `ctl` | `type(fv_controls)` | in | controls for the refits (`nthreads`, `method`, ...) |
| `weights` | `real(dp)(n)` | in, optional | observation weights |
| `nrel_out` | `integer` | out, optional | number of relaxed edges |

Requires `trunc >= 2` and no discrete variables (otherwise returns without change). See
[Non-simplified edges](../manual/vinecop.md#non-simplified-edges-relax).

## unrelax

```fortran
subroutine unrelax(self)
```

Deletes the relaxed-edge records (`nrel = 0`). Does not restore the parameters `relax`
changed.

## rel_params

```fortran
function rel_params(self, uc) result(prel)
  class(fv_vinecop), intent(in) :: self
  real(dp), intent(in) :: uc(:, :)
  real(dp) :: prel(3, self%nrel, size(uc, 2))
```

Parameters of the relaxed edges at the samples `uc(d, n)` (variables first), linear in $w$
between the anchors and constant outside them.

## Example

```fortran
type(fv_vinecop) :: vc
real(dp) :: sim(1000, 3)
call vc%set_dvine([1, 2, 3])
vc%fam(1:2, 1) = FV_CLAYTON; vc%par(1, 1:2, 1) = 2.0_dp
call vc%simulate(1000, sim, seed=1)
```

Complete programs are in the [manual](../manual/vinecop.md#example).
