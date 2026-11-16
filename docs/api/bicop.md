# `fv_bicop`

A bivariate (pair) copula. Manual: [Pair copulas](../manual/bicop.md).

```fortran
type, public :: fv_bicop
  integer  :: family   = FV_INDEP
  integer  :: rotation = 0
  real(dp) :: par(3)   = 0.0_dp
  real(dp) :: loglik   = 0.0_dp
  real(dp), allocatable :: aux(:)        ! nonparametric estimate (tll, bernstein)
contains
  procedure :: pdf, logpdf, cdf, hfunc1, hfunc2, hinv1, hinv2
  procedure :: fit, select, npars, npars_eff
end type fv_bicop
```

## Fields

| name | type | default | description |
|---|---|---|---|
| `family` | `integer` | `FV_INDEP` | family code, see [constants](constants.md) and the [catalogue](../families/catalogue.md) |
| `rotation` | `integer` | 0 | 0, 90, 180 or 270; ignored by rotationless families |
| `par` | `real(dp)(3)` | 0 | parameters, first `npars()` used |
| `loglik` | `real(dp)` | 0 | log-likelihood of the last `fit` / `select` (weighted if weights were given) |
| `aux` | `real(dp), allocatable (:)` | unallocated | TLL / Bernstein estimate. Allocated by `fit` (length `NAUX` for TLL / Bernstein, else 1) and `select` (length `NAUX` if the family set contains 24 or 25, else 1) |

---

## pdf

```fortran
subroutine pdf(self, u1, u2, out)
  class(fv_bicop), intent(in) :: self
  real(dp), intent(in)  :: u1(:), u2(:)
  real(dp), intent(out) :: out(:)
```

| name | type | intent | description |
|---|---|---|---|
| `u1`, `u2` | `real(dp)(:)` | in | evaluation points, same length $m$ |
| `out` | `real(dp)(:)` | out | copula density $c(u_1, u_2)$, length $m$; never below `tiny(1.0_dp)` |

```fortran
call cop%pdf([0.2_dp, 0.8_dp], [0.3_dp, 0.9_dp], dens)
```

## logpdf

```fortran
subroutine logpdf(self, u1, u2, out)
```

Same arguments as `pdf`; `out` $= \log c(u_1, u_2)$, computed directly in log space (no
floor).

## cdf

```fortran
subroutine cdf(self, u1, u2, out)
```

Same arguments as `pdf`; `out` $= C(u_1, u_2)$. Closed form for every family except the
Student copula (Gauss-Legendre quadrature of $h_1$, accuracy about $10^{-6}$). Runs on the
OpenMP default team.

## hfunc1, hfunc2

```fortran
subroutine hfunc1(self, u1, u2, out)
subroutine hfunc2(self, u1, u2, out)
```

| name | type | intent | description |
|---|---|---|---|
| `u1`, `u2` | `real(dp)(:)` | in | points |
| `out` | `real(dp)(:)` | out | `hfunc1`: $\partial C/\partial u_1 = P(U_2 \le u_2 \mid U_1 = u_1)$; `hfunc2`: $\partial C/\partial u_2 = P(U_1 \le u_1 \mid U_2 = u_2)$ |

## hinv1, hinv2

```fortran
subroutine hinv1(self, u1, u2, out)
subroutine hinv2(self, u1, u2, out)
```

| name | type | intent | description |
|---|---|---|---|
| `u1` | `real(dp)(:)` | in | `hinv1`: conditioning value $u_1$; `hinv2`: level $w$ |
| `u2` | `real(dp)(:)` | in | `hinv1`: level $w$; `hinv2`: conditioning value $u_2$ |
| `out` | `real(dp)(:)` | out | `hinv1`: $v$ with $h_1(u_1, v) = w$; `hinv2`: $v$ with $h_2(v, u_2) = w$ |

```fortran
call random_number(u1); call random_number(w)
call cop%hinv1(u1, w, u2)          ! (u1, u2) is a sample of the copula
```

## fit

```fortran
subroutine fit(self, u1, u2, ctl, weights)
  class(fv_bicop), intent(inout) :: self
  real(dp), intent(in) :: u1(:), u2(:)
  type(fv_controls), intent(in), optional :: ctl
  real(dp), intent(in), optional :: weights(:)
```

Estimates `par` (and `aux` for TLL / Bernstein) for the current `family` and `rotation`.

| name | type | intent | description |
|---|---|---|---|
| `u1`, `u2` | `real(dp)(:)` | in | pseudo-observations, length $n$ |
| `ctl` | `type(fv_controls)` | in, optional | uses `method`, `fast_mle`, `tau_mode`, `nthreads`; default `fv_controls()` |
| `weights` | `real(dp)(:)` | in, optional | observation weights, length $n$ |

Sets `par`, `loglik`, `aux`. `FV_ITAU` is only meaningful for families 0 to 6; see the
[warning in the manual](../manual/bicop.md#fitting-a-given-family).

```fortran
cop%family = FV_BB1; cop%rotation = 0
call cop%fit(u1, u2)
```

## select

```fortran
subroutine select(self, u1, u2, ctl, weights)
  class(fv_bicop), intent(inout) :: self
  real(dp), intent(in) :: u1(:), u2(:)
  type(fv_controls), intent(in), optional :: ctl
  real(dp), intent(in), optional :: weights(:)
```

Chooses `family`, `rotation` and `par` among `ctl%family_set(1:ctl%nfam)`.

| name | type | intent | description |
|---|---|---|---|
| `u1`, `u2` | `real(dp)(:)` | in | pseudo-observations, length $n$ |
| `ctl` | `type(fv_controls)` | in, optional | uses `family_set`, `nfam`, `selcrit`, `psi0`, `method`, `fast_mle`, `preselect`, `allow_rotations`, `tau_mode`, `nthreads` |
| `weights` | `real(dp)(:)` | in, optional | observation weights, length $n$ |

Sets `family`, `rotation`, `par`, `loglik` and `aux`.

```fortran
type(fv_controls) :: ctl
ctl%selcrit = FV_AIC
call cop%select(u1, u2, ctl)
```

## npars

```fortran
integer function npars(self)
```

Number of parameters of `family`: 0 (independence), 1, 2 or 3 (for generic codes
[`gen_npar`](functions.md#gen_npar)); for TLL / Bernstein the rounded effective degrees of
freedom.

## npars_eff

```fortran
real(dp) function npars_eff(self)
```

As `npars`, as a real number; for TLL / Bernstein the unrounded effective degrees of freedom
`aux(2)` (0 before a fit).

## Example

```fortran
type(fv_bicop) :: cop
real(dp) :: d(1)
cop%family = FV_JOE; cop%rotation = 90; cop%par(1) = 2.5_dp
call cop%pdf([0.9_dp], [0.1_dp], d)
```

A complete program is in the [manual](../manual/bicop.md#a-complete-example).
