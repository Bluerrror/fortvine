# Functions

Module-level procedures of `fortvine_api`. Kernel routines of the module `fortvine` that are
not re-exported here (for example `bicop_eval_x`) are internal; their C counterparts are
documented in the [C interface](../manual/c-api.md).

## Utilities

### fv_threads

```fortran
integer function fv_threads(req)
  integer, intent(in), optional :: req
```

| name | type | intent | description |
|---|---|---|---|
| `req` | `integer` | in, optional | requested thread count |

Returns `req` if it is positive; otherwise the value of the environment variable
`FORTVINE_THREADS` if that is a positive integer; otherwise `max(1, omp_get_num_procs() / 2)`
(1 without OpenMP). See [Performance](../manual/performance.md#how-many-threads).

```fortran
print *, fv_threads()        ! what nthreads = 0 resolves to
```

### fv_pseudo_obs_nd

```fortran
subroutine fv_pseudo_obs_nd(x, u, nthreads)
  real(dp), intent(in)  :: x(:, :)
  real(dp), intent(out) :: u(:, :)
  integer, intent(in), optional :: nthreads
```

| name | type | intent | description |
|---|---|---|---|
| `x` | `real(dp)(n, d)` | in | data |
| `u` | `real(dp)(n, d)` | out | pseudo-observations: rank / $(n+1)$ per column, ties get their average rank |
| `nthreads` | `integer` | in, optional | threads (parallel over columns) |

### fv_kendall_tau

```fortran
real(dp) function fv_kendall_tau(x, y)
  real(dp), intent(in) :: x(:), y(:)
```

Kendall's $\tau_b$ (tie-corrected) of two vectors of equal length, $O(n \log n)$ (Knight 1966).

### fv_tau_mat

```fortran
subroutine fv_tau_mat(u, tm, nthreads)
  real(dp), intent(in)  :: u(:, :)
  real(dp), intent(out) :: tm(:, :)
  integer, intent(in), optional :: nthreads
```

| name | type | intent | description |
|---|---|---|---|
| `u` | `real(dp)(n, d)` | in | data |
| `tm` | `real(dp)(d, d)` | out | symmetric matrix of pairwise $\tau_b$, unit diagonal |
| `nthreads` | `integer` | in, optional | threads (parallel over pairs) |

### fv_dependence

```fortran
real(dp) function fv_dependence(x, y, kind, weights)
  real(dp), intent(in) :: x(:), y(:)
  integer, intent(in) :: kind
  real(dp), intent(in), optional :: weights(:)
```

| name | type | intent | description |
|---|---|---|---|
| `x`, `y` | `real(dp)(:)` | in | pseudo-observations of equal length $n$ |
| `kind` | `integer` | in | `FV_TREE_TAU`, `FV_TREE_RHO`, `FV_TREE_JOE`, `FV_TREE_HOEFFD`, `FV_TREE_XI`, `FV_TREE_MI`, `FV_TREE_FERREIRA` |
| `weights` | `real(dp)(n)` | in, optional | observation weights; used by $\tau$, $\rho$ and Joe only |

The absolute value of the dependence measure that vine selection uses as edge weight, see
[`tree_crit`](../manual/controls.md#tree_crit). Returns 0 for $n \le 10$.

```fortran
program funcs_demo
  use fortvine_api
  implicit none
  integer, parameter :: dp = kind(1.0d0), n = 1000
  real(dp) :: x(n, 3), u(n, 3), tm(3, 3)
  character(len=8) :: kinds(0:6) = [character(len=8) :: 'tau', 'rho', 'joe', 'hoeffd', 'xi', 'mi', 'ferreira']
  integer :: k, i

  call random_seed(put=[(7 * i, i = 1, 33)])
  call random_number(x)
  x(:, 2) = x(:, 1) + 0.3_dp * x(:, 2)              ! monotone dependence
  x(:, 3) = (x(:, 1) - 0.5_dp)**2 + 0.02_dp * x(:, 3)   ! non-monotone dependence
  call fv_pseudo_obs_nd(x, u)
  call fv_tau_mat(u, tm)
  print '(a,3f8.3)', 'tau matrix row 1:', tm(1, :)
  do k = 0, 6
    print '(a8,a,f7.3,a,f7.3)', kinds(k), '  (1,2):', fv_dependence(u(:, 1), u(:, 2), k), &
      '   (1,3):', fv_dependence(u(:, 1), u(:, 3), k)
  end do
end program funcs_demo
```

Output:

```text
tau matrix row 1:   1.000   0.818  -0.087
tau       (1,2):  0.818   (1,3):  0.087
rho       (1,2):  0.961   (1,3):  0.098
joe       (1,2):  1.005   (1,3):  0.003
hoeffd    (1,2):  0.724   (1,3):  0.214
xi        (1,2):  0.753   (1,3):  0.892
mi        (1,2):  1.321   (1,3):  2.389
ferreira  (1,2):  0.855   (1,3):  0.103
```

The pair (1, 3) is a noisy parabola: strongly dependent but almost unconcordant. Kendall's
$\tau$, Spearman's $\rho$ and the Gaussian measure `joe` barely see it; Chatterjee's $\xi$,
mutual information and Hoeffding's $D$ do.

## Families

### fv_fam_name

```fortran
pure function fv_fam_name(fam) result(nm)
  integer, intent(in) :: fam
  character(len=48) :: nm
```

Name of any family code: `FV_FAMILY_NAME(fam)` for 0 to 25 (`'gaussian'`, `'bb1'`, ...), the
generic name for codes $\ge 10000$ (see [Names](../manual/generic.md#names)), `'invalid'` for
invalid generic codes. Use `trim`.

### fv_par_bounds

```fortran
pure subroutine fv_par_bounds(fam, lb, ub)
  integer, intent(in) :: fam
  real(dp), intent(out) :: lb(3), ub(3)
```

The box searched by maximum likelihood for family `fam`; unused entries are 0. Independence,
TLL and Bernstein return zeros. Table in [Pair copulas](../manual/bicop.md#parameter-bounds).

### fv_generic_families

```fortran
pure function fv_generic_families(new_only) result(cs)
  logical, intent(in), optional :: new_only
  integer, allocatable :: cs(:)
```

| name | type | intent | description |
|---|---|---|---|
| `new_only` | `logical` | in, optional | `.true.` (default): leave out codes for which `gen_dup` is true |

Every valid generic code in increasing order: 85 codes with `new_only = .true.`, 114 with
`.false.`.

```fortran
integer, allocatable :: gen(:)
gen = fv_generic_families()
```

### fam_is_np

```fortran
pure elemental logical function fam_is_np(fam)
  integer, intent(in) :: fam
```

True for the nonparametric families 24 (TLL) and 25 (Bernstein).

### gen_is

```fortran
pure elemental logical function gen_is(fam)
```

True for generic codes (`fam >= 10000`); says nothing about validity.

### gen_valid

```fortran
pure elemental logical function gen_valid(fam)
```

True if `fam` is a well-defined, identifiable generic code with 1 to 3 parameters; the rules
are in [Valid codes](../manual/generic.md#valid-codes-and-duplicates).

### gen_dup

```fortran
pure elemental logical function gen_dup(fam)
```

True if the copula of generic code `fam` is also available as a hand-coded family or under
another generic code. Meaningful only for valid codes.

### gen_npar

```fortran
pure elemental integer function gen_npar(fam)
```

Number of parameters of a generic code: outer generator (0 to 2) + inner function (0 to 3) +
2 for the Khoudraji wrapper.

### gen_exch

```fortran
pure elemental logical function gen_exch(fam)
```

True if the generic copula is exchangeable, $C(u, v) = C(v, u)$: no Khoudraji wrapper and an
inner index of at most 4.

### gen_lt, gen_ut

```fortran
pure elemental logical function gen_lt(fam)
pure elemental logical function gen_ut(fam)
```

Lower- and upper-tail dependence of the unrotated generic copula, as used by the preselection
(see [Tail dependence](../manual/generic.md#tail-dependence-symmetry-and-transposition)).

### gen_flip

```fortran
pure subroutine gen_flip(fam, par)
  integer, intent(inout) :: fam
  real(dp), intent(inout) :: par(3)
```

| name | type | intent | description |
|---|---|---|---|
| `fam` | `integer` | inout | generic code; may change (asymmetric logistic types 1 and 2 swap) |
| `par` | `real(dp)(3)` | inout | parameters, replaced by those of the transpose |

Replaces `(fam, par)` by the transposed copula $C^\top(u, v) = C(v, u)$. No change for
exchangeable codes. See [Edge flips](../manual/rotations.md#asymmetric-families-and-edge-flips).

```fortran
f = 10006; p = [2.5_dp, 0.4_dp, 0.0_dp]
call gen_flip(f, p)            ! f = 10007 (Tawn type 2), p unchanged
```

## Fit cache

The cache memoises pair fits across calls; see [Fit cache](../manual/fitting.md#fit-cache).

### fv_fitcache_enable

```fortran
subroutine fv_fitcache_enable(capacity)
  integer, intent(in) :: capacity
```

Allocates an empty cache with room for `max(16, capacity)` fits and resets the counters.
Replaces any existing cache.

### fv_fitcache_capacity

```fortran
integer function fv_fitcache_capacity()
```

Size of the cache table; 0 when disabled.

### fv_fitcache_clear

```fortran
subroutine fv_fitcache_clear()
```

Empties the cache (keeps it enabled) and resets the hit and miss counters.

### fv_fitcache_disable

```fortran
subroutine fv_fitcache_disable()
```

Frees the cache; subsequent fits are not cached.

### fv_fitcache_stats

```fortran
subroutine fv_fitcache_stats(hits, misses, used)
  integer(int64), intent(out) :: hits, misses
  integer, intent(out), optional :: used
```

| name | type | intent | description |
|---|---|---|---|
| `hits` | `integer(int64)` | out | lookups answered from the cache since `enable` / `clear` |
| `misses` | `integer(int64)` | out | lookups that had to fit |
| `used` | `integer` | out, optional | occupied slots |

`int64` is the kind from `iso_fortran_env`.
