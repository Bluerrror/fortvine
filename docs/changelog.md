# Changelog

All notable changes to fortvine are recorded here. Versions follow
[semantic versioning](https://semver.org).

## 1.0.0 (2026-11)

First public release.

### Pair copulas

* **111 families.** 26 hand-coded families (codes 0-25): independence, Gaussian, Student t,
  Clayton, Gumbel, Frank, Joe, BB1-BB10, Tawn, Plackett, FGM, AMH, Galambos, Hüsler-Reiss,
  t-EV, and the nonparametric transformation local-likelihood (TLL) and Bernstein
  estimators; plus 85 generic families (codes 10000 and up) built by one engine,
  $C = \Psi(\ell(\Phi(u), \Phi(v)))$: two-parameter Archimedean compositions, extreme-value
  copulas (Tawn mixed, asymmetric logistic and negative logistic models and their types,
  Khoudraji-Hüsler-Reiss, ...), Archimax copulas with one- and two-parameter outer
  generators, and Khoudraji asymmetrisations of Clayton, Joe and Frank.
  See the [catalogue](families/catalogue.md).
* Rotations by 90, 180 and 270 degrees (vinecopulib convention); exact transposition of
  asymmetric families when a vine edge is flipped.
* Log-density, density, distribution function, both h-functions and their inverses for
  every family, vectorised and parallel; log-space evaluation with cancellation-free and
  asymptotic forms throughout ([Numerics](theory/numerics.md)).
* Maximum likelihood with **exact parameter scores** from dual numbers (each family written
  once and compiled for real and dual types), projected Fisher-scoring (BHHH) steps with
  Levenberg-Marquardt damping, subsample-then-polish fitting, and analytic second
  derivatives for the one-parameter families ([Exact parameter scores](theory/scores.md));
  estimation by inversion of Kendall's tau as an option.
* Family selection by log-likelihood, AIC, BIC or mBIC, with vinecopulib-style
  preselection and racing of candidate fits.

### Vine copulas

* Regular vines with the natural-order structure array; fixed C- and D-vine constructors;
  R-vine matrix export in original labels.
* Structure selection by Dissmann's algorithm with Prim's or Kruskal's maximum spanning
  tree and seven edge criteria (Kendall's tau, Spearman's rho, Joe's criterion, Hoeffding's
  D, Chatterjee's xi, KSG mutual information, Ferreira's criterion).
* Truncation and thresholding, fixed or chosen by mBICv (vinecopulib's
  `select_trunc_lvl` and `select_threshold`).
* Vine density and log-likelihood, a fast density on prediction grids (`pdf_grid`),
  Rosenblatt and inverse Rosenblatt transforms, reproducible simulation.
* **Observation weights** for every fit and selection.
* **Discrete and mixed margins** (vinecopulib convention, $u = F(x)$, $u^- = F(x^-)$),
  including Brockwell's randomised Rosenblatt transform.
* Non-simplified edges (`relax`): pair-copula parameters varying with the conditioning
  variables, with a constant-conditional-correlation test.
* A pair-fit memo cache for repeated selections (tuning loops).

### Infrastructure

* Pure Fortran 2008 with OpenMP; no external dependencies. Likelihood sums are independent
  of the thread count.
* Typed Fortran API (`fortvine_api`: `fv_bicop`, `fv_vinecop`, `fv_controls`) and a
  `bind(C)` interface exported from a shared library for C, C++, Python and other languages.
* Build with `build.sh` (static and shared library, tests, examples, benchmark driver) or
  fpm.
* Two property test suites (`test/test_fortvine.f90`, `test/test_genfam.f90`) covering
  every family and rotation, vines, weights, discrete data, criteria and the generic engine
  ([Validation](theory/validation.md)).
* Five documented example programs, a Python (ctypes) and a C example, and a cross-library
  benchmark against pyvinecopulib and torchvine.
