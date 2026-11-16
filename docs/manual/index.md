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

=== "Fortran"

    ```fortran title="workflow.f90"
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

=== "Python"

    ```python title="workflow.py"
    import numpy as np
    import fortvine as fv

    n, d = 1500, 4

    # 1. raw data with arbitrary margins (here: a skewed, dependent toy sample)
    z = -np.log(np.random.default_rng(17).random((n, d)))      # exponential noise
    x = np.empty((n, d))
    x[:, 0] = z[:, 0]
    x[:, 1] = z[:, 0] + z[:, 1]
    x[:, 2] = np.sqrt(x[:, 1]) + 0.5 * z[:, 2]
    x[:, 3] = np.exp(-z[:, 3]) * x[:, 0]

    # 2. pseudo-observations: ranks / (n + 1), column by column
    u = fv.to_pseudo_obs(x)

    # 3. structure, families and parameters in one call
    vc = fv.Vinecop.from_data(u, num_threads=4)
    print(vc)

    # 4. use the fitted model
    print(f"log-likelihood{vc.loglik(u):10.2f}   parameters {vc.npars}")
    sim = vc.simulate(5, seed=1)
    print("first simulated row:" + "".join(f"{s:8.4f}" for s in sim[0]))
    ```

=== "R"

    ```r title="workflow.R"
    library(fortvine)
    n <- 1500; d <- 4

    # 1. raw data with arbitrary margins (here: a skewed, dependent toy sample)
    set.seed(17)
    z <- matrix(rexp(n * d), n, d)                             # exponential noise
    x <- matrix(0, n, d)
    x[, 1] <- z[, 1]
    x[, 2] <- z[, 1] + z[, 2]
    x[, 3] <- sqrt(x[, 2]) + 0.5 * z[, 3]
    x[, 4] <- exp(-z[, 4]) * x[, 1]

    # 2. pseudo-observations: ranks / (n + 1), column by column
    u <- pseudo_obs(x)

    # 3. structure, families and parameters in one call
    vc <- vinecop(u, cores = 4)
    vc

    # 4. use the fitted model
    cat(sprintf("log-likelihood%10.2f   parameters %d\n", logLik(vc), vc$npars))
    set.seed(1)
    sim <- rvinecop(5, vc)
    cat("first simulated row:", sprintf("%8.4f", sim[1, ]), "\n", sep = "")
    ```

Output (Fortran: gfortran 16.1, 4 threads):

=== "Fortran"

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

=== "Python"

    ```text
    fortvine vine copula: d = 4, trunc = 3
    order: 3 2 1 4
      tree 1 edge 1: tawn         rot 180 par 1.0000 0.5916 4.0834
      tree 1 edge 2: tawn         rot   0 par 0.5120 1.0000 3.9125
      tree 1 edge 3: tawn         rot   0 par 0.6868 1.0000 3.9937
      tree 2 edge 1: indep        rot   0 par
      tree 2 edge 2: indep        rot   0 par
      tree 3 edge 1: indep        rot   0 par
    log-likelihood   2278.20   parameters 9
    first simulated row:  0.9210  0.9385  0.6721  0.9486
    ```

=== "R"

    ```text
    Vine copula ('vinecop_dist'): d = 4, trunc_lvl = 3, order = 3 2 1 4
      npars = 10, loglik = 2392.89 (n = 1500)
     tree edge conditioned conditioning family rotation           parameters
        1    1        3, 2                tawn      180     1, 0.6212, 4.164
        1    2        2, 1                tawn        0 0.5333,     1, 4.124
        1    3        1, 4                tawn        0 0.6527,     1,  4.23
        2    1        3, 1            2  indep        0                 none
        2    2        2, 4            1 gumbel        0                1.042
        3    1        3, 4         2, 1  indep        0                 none
    log-likelihood   2392.89   parameters 10
    first simulated row:  0.3022  0.7978  0.5475  0.4977
    ```

The raw data are random draws, so the three programs work on different samples and select
slightly different vines. The Fortran program is compiled against the built library (see
[Installation](../getting-started/installation.md)):

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
