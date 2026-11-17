---
hide:
  - navigation
  - toc
---

<div class="fv-hero" markdown>

<img class="fv-logo" src="assets/logo.svg" alt="">

<div markdown>

# fortvine

<p class="fv-tagline">Vine copulas in modern Fortran: 111 pair-copula families with exact
parameter scores, structure selection, densities, Rosenblatt transforms and simulation.</p>

[Get started](getting-started/installation.md){ .md-button .md-button--primary }
[Browse the families](families/index.md){ .md-button }

</div>
</div>

<figure class="fv-strip" markdown>
![Density contours of eight fortvine copula families](assets/hero.png)
<figcaption>Eight of the 111 families, drawn by fortvine in normal scores: Gaussian, Student t,
Clayton, Gumbel, Frank, Tawn's asymmetric logistic, an Archimax copula (gamma generator with an
asymmetric logistic Pickands function) and a Khoudraji-Joe copula.</figcaption>
</figure>

fortvine fits, evaluates and simulates bivariate (pair) copulas and regular-vine copulas. It is
standard Fortran 2008 with OpenMP, needs nothing but the compiler, and can be called from
Fortran, C, Python, R or any language that speaks the C ABI.

## Thirty seconds of fortvine

=== "Fortran"

    ```fortran title="hello_fortvine.f90"
    program hello_fortvine
      use fortvine_api
      implicit none
      integer, parameter :: dp = kind(1.0d0), n = 2000
      type(fv_bicop) :: cop
      real(dp) :: u1(n), w(n), u2(n)

      cop%family = FV_CLAYTON; cop%par(1) = 2.0_dp       ! Clayton copula, Kendall's tau = 0.5
      call random_number(u1); call random_number(w)
      call cop%hinv1(u1, w, u2)                          ! simulate (u1, u2) from it

      call cop%select(u1, u2)                            ! choose family, rotation and parameters by BIC
      print *, trim(fv_fam_name(cop%family)), cop%rotation, cop%par(1)
    end program hello_fortvine
    ```

    ```text
     clayton           0   2.1450855929878818
    ```

    The sample is not seeded, so the last digits differ from run to run.

=== "Python"

    ```python title="hello_fortvine.py"
    import fortvine as fv

    cop = fv.Bicop("clayton", 0, [2.0])              # Clayton copula, Kendall's tau = 0.5
    u = cop.simulate(2000, seed=1)                   # simulate (u1, u2) from it

    best = fv.Bicop().select(u)                      # choose family, rotation and parameters by BIC
    print(best.family, best.rotation, best.parameters[0])
    ```

    ```text
    clayton 0 1.956774237371629
    ```

=== "R"

    ```r title="hello_fortvine.R"
    library(fortvine)
    cop <- bicop_dist("clayton", 0, 2)               # Clayton copula, Kendall's tau = 0.5
    set.seed(1)
    u <- rbicop(2000, cop)                           # simulate (u1, u2) from it

    best <- bicop(u)                                 # choose family, rotation and parameters by BIC
    best
    ```

    ```text
    Bivariate copula ('bicop_dist'): family = clayton (3), rotation = 0, parameters = 1.879
    ```

The [Python](python.md) and [R](r.md) packages are bindings: every number above is computed by
the same compiled Fortran library. The three samples come from different random number
generators, so the estimates differ. `bash build.sh` builds `build/libfortvine.a`, the module
files and a shared library; see [Installation](getting-started/installation.md) and the
[Quick start](getting-started/quickstart.md).

## What is in the box

| Area | What fortvine provides |
|---|---|
| Families | elliptical (Gaussian, t); Archimedean (Clayton, Gumbel, Frank, Joe, AMH, BB1-BB3, BB6-BB10 and 6 new Laplace-transform compositions); extreme-value (Tawn, Galambos, Hüsler-Reiss, t-EV, BB5 and 9 new asymmetric and mixed Pickands families); Archimax (BB4 and 67 new); Khoudraji (3 new); Plackett, FGM; nonparametric TLL and Bernstein. All with 90/180/270 degree rotations where they apply. See [Families](families/index.md). |
| Pair copulas | density, log-density, CDF, both h-functions and their inverses, maximum-likelihood fit with exact scores, family selection (log-likelihood, AIC, BIC, mBIC), Kendall's tau inversion, weights, discrete margins. See [Pair copulas](manual/bicop.md). |
| Vines | R-, C- and D-vine structures; Dissmann selection with Kendall's tau, Spearman's rho, Hoeffding's D, Chatterjee's xi, mutual information and other tree criteria; truncation and thresholding by mBICv; fixed structures. See [Vine copulas](manual/vinecop.md). |
| Vine functions | density on samples and on prediction grids, log-likelihood, Rosenblatt and inverse Rosenblatt, simulation, parameter counts, structure matrix, printing. |
| Interfaces | a typed Fortran API (`fortvine_api`), a C ABI (`bind(C)` kernels and `fortvine_capi`), and [Python](python.md) and [R](r.md) packages with every feature of the Fortran API. See [C interface](manual/c-api.md). |
| Quality | test suites over every family and rotation (density / h-function / CDF consistency, inverses, MLE recovery); the generic engine reproduces the hand-coded families to 1e-15-1e-9. See [Validation](theory/validation.md). |

## How it compares

fortvine follows the conventions of [vinecopulib](https://github.com/vinecopulib/vinecopulib)
(rotations, structure arrays, selection algorithm), so the two can be compared one to one: on
the benchmark data they reach the same maximum log-likelihoods, and fortvine gets there
7-12 times faster, from Fortran or from Python (2-3 times slower from R on Windows, still
7-9 times faster than rvinecopulib).

![Vine selection time against dimension for fortvine, pyvinecopulib and torchvine](assets/bench_vine.png)

| | fortvine | vinecopulib / pyvinecopulib | torchvine |
|---|---|---|---|
| language | Fortran 2008 + OpenMP | C++ (Python bindings) | Python / PyTorch |
| parametric families | 108 | 11 | 10 |
| asymmetric families | 32 | 1 (Tawn) | — |
| nonparametric | TLL, Bernstein | TLL | TLL |
| exact parameter gradients | yes (dual numbers) | yes (analytic, 1.0) | autograd |
| GPU | no | no | yes |

All timings and the method are on the [Benchmarks](benchmarks.md) page.

## Citing

If you use fortvine — from Fortran, C, Python or R — in published work, please cite:

> Shojaeezadeh, S. A. (2026). *fortvine: a complete compendium to work on vine copula modeling.*
> In preparation.

```bibtex
@unpublished{shojaeezadeh2026fortvine,
  author = {Shojaeezadeh, Shahab Aldin},
  title  = {fortvine: a complete compendium to work on vine copula modeling},
  year   = {2026},
  note   = {In preparation. Software: https://github.com/Bluerrror/fortvine}
}
```

The same reference is returned by `citation("fortvine")` in R, by `fortvine.cite()` in Python and by
GitHub's "Cite this repository" button.

## Licence

fortvine is free software, released under the
[GNU General Public License, version 3 or later](https://github.com/Bluerrror/fortvine/blob/main/LICENSE).
You may use, study, modify and redistribute it; distributed copies and derived works must remain
under the GPL and keep the copyright notice of Shahab Aldin Shojaeezadeh.
