# User manual

This manual explains how fortvine is organised and how its pieces fit together. The
[API reference](../api/index.md) lists every routine and argument; the manual is where the
behaviour behind those signatures is described: which defaults are used, what the fitting
code actually does, and which options interact.

## Two levels of interface

fortvine ships two Fortran modules:

| module | file | what it is for |
|---|---|---|
| `fortvine_api` | `src/fortvine_api.f90` | the **typed interface**: derived types with type-bound procedures, Fortran-natural array layout `u(n, d)`. This is what Fortran users should `use`. |
| `fortvine` | `src/fortvine.f90` | the **kernels**: evaluation, fitting, selection, Rosenblatt transforms, the generic family engine and the `bind(C)` entry points. Arrays have observations in the *last* dimension (`u(d, n)`). |

Everything a Fortran program needs is re-exported by `fortvine_api`, including the family
constants, the fit cache controls and the `gen_*` helpers of the generic families. The kernel
module is only of direct interest if you call the library from another language (see the
[C interface](c-api.md)).

## The three types

`fv_bicop`
:   One bivariate (pair) copula: a family code, a rotation, up to three parameters, the
    log-likelihood of the last fit and, for the two nonparametric families, an estimate vector
    `aux`. Methods evaluate the density, distribution function, both h-functions and their
    inverses, fit parameters for a given family, or select the family. See
    [Pair copulas](bicop.md).

`fv_vinecop`
:   A regular vine copula in `d` dimensions, stored as a structure (an order vector plus a
    natural-order structure array) and one pair copula per edge. Methods select structure,
    families and parameters (Dissmann et al. 2013), and evaluate the density, the Rosenblatt
    transform and its inverse, and simulate. See [Vine copulas](vinecop.md).

`fv_controls`
:   A plain record of options shared by every fit and selection: family set, selection
    criterion, estimation method, tree criterion, truncation, thresholding, threads. Every
    field has a default, so `type(fv_controls) :: ctl` is a valid configuration. See
    [Controls](controls.md).

## A typical session

The usual sequence is: transform the data to the copula scale, select a model, then use it.

```fortran
program workflow
  use fortvine_api
  implicit none
  integer, parameter :: dp = kind(1.0d0), n = 1500, d = 4
  real(dp) :: z(n, d), x(n, d), u(n, d), sim(5, d)
  type(fv_vinecop) :: vc
  type(fv_controls) :: ctl
  integer :: i

  ! 1. raw data with arbitrary margins (here: a skewed, dependent toy sample)
  call random_seed(put=[(17 * i, i = 1, 33)])
  call random_number(z)
  z = -log(z)                                   ! exponential noise
  x(:, 1) = z(:, 1)
  x(:, 2) = z(:, 1) + z(:, 2)
  x(:, 3) = sqrt(x(:, 2)) + 0.5_dp * z(:, 3)
  x(:, 4) = exp(-z(:, 4)) * x(:, 1)

  ! 2. pseudo-observations: ranks / (n + 1), column by column
  call fv_pseudo_obs_nd(x, u)

  ! 3. structure, families and parameters in one call
  ctl%nthreads = 4
  call vc%select(u, ctl)
  call vc%print()

  ! 4. use the fitted model
  print '(a,f10.2,a,i0)', 'log-likelihood', vc%loglik(u), '   parameters ', vc%npars()
  call vc%simulate(5, sim, seed=1)
  print '(a,4f8.4)', 'first simulated row:', sim(1, :)
end program workflow
```

Output (gfortran 16.1, 4 threads):

```text
fortvine vine copula: d = 4, trunc = 3
order: 3 2 1 4
  tree 1 edge 1:         tawn rot 180 par  1.0000E+00  5.9113E-01  4.1699E+00
  tree 1 edge 2:         tawn rot   0 par  5.1433E-01  1.0000E+00  3.8241E+00
  tree 1 edge 3:         tawn rot   0 par  6.7211E-01  1.0000E+00  3.7311E+00
  tree 2 edge 1:       gumbel rot 180 par  1.0480E+00  0.0000E+00  0.0000E+00
  tree 2 edge 2:       gumbel rot 180 par  1.0573E+00  0.0000E+00  0.0000E+00
  tree 3 edge 1:        indep rot   0 par  0.0000E+00  0.0000E+00  0.0000E+00
log-likelihood   2242.61   parameters 11
first simulated row:  0.2909  0.8949  0.9592  0.3023
```

Compile against the built library (see [Installation](../getting-started/installation.md)):

```bash
gfortran -O2 -ffree-line-length-none -fopenmp -I build workflow.f90 build/libfortvine.a -o workflow
```

!!! note "Copula scale"
    Every routine expects data on the unit square or cube. fortvine does not model margins:
    use `fv_pseudo_obs_nd` (ranks) or your own fitted marginal distributions to map the data to
    $(0,1)$ first. Values are clamped internally to $[10^{-10}, 1-10^{-10}]$, so exact 0 or 1
    do not crash, but they carry no information.

## Conventions used throughout

* **Array layout (Fortran API).** Data matrices are `u(n, d)`: rows are observations,
  columns are variables. Pair data are two vectors `u1(n)`, `u2(n)`.
* **h-functions.** `hfunc1(u1, u2)` $= P(U_2 \le u_2 \mid U_1 = u_1) = \partial C/\partial u_1$,
  `hfunc2(u1, u2)` $= P(U_1 \le u_1 \mid U_2 = u_2) = \partial C/\partial u_2$. Details in
  [Pair copulas](bicop.md#h-functions-and-their-inverses).
* **Rotations** are counter-clockwise rotations of the data by 0, 90, 180 or 270 degrees,
  following vinecopulib. See [Rotations and asymmetry](rotations.md).
* **Family codes** are integers: 0 to 25 for the hand-coded families, $\ge 10000$ for the
  generic engine. The constants are listed in [Family codes and constants](../api/constants.md)
  and the full list in the [catalogue](../families/catalogue.md).
* **Variable labels** in a vine are 1-based and refer to the columns of `u`.
* **Threads.** Every heavy routine runs on OpenMP threads; the count comes from
  `fv_controls%nthreads`, an optional `nthreads` argument, the `FORTVINE_THREADS`
  environment variable or half the logical CPUs, in that order. See
  [Performance](performance.md).

## Where to go next

| task | page |
|---|---|
| evaluate, fit or choose a single pair copula | [Pair copulas](bicop.md) |
| select, evaluate, transform and simulate a vine | [Vine copulas](vinecop.md) |
| understand what the optimiser does and how to make it faster | [Fitting and selection](fitting.md) |
| change any option | [Controls](controls.md) |
| negative dependence, tails, asymmetric families | [Rotations and asymmetry](rotations.md) |
| count data or mixed continuous/discrete data | [Discrete and mixed data](discrete.md) |
| TLL and Bernstein estimates | [Nonparametric families](nonparametric.md) |
| the 85 extra families built from generators and Pickands functions | [The generic engine](generic.md) |
| threads, build flags, cost of options | [Performance](performance.md) |
| calling from C, Python, Julia or R | [C interface](c-api.md) |
