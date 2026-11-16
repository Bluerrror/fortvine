# Pair copula basics

**Program:** `example/ex01_pair_copula.f90` &nbsp;·&nbsp; **Topics:** `fv_bicop`, density,
inverse h-function, simulation, maximum likelihood, BIC family selection.

The program works with one bivariate copula from start to finish: it evaluates a Clayton
copula, draws a sample from it, refits the parameter, and finally lets fortvine choose the
family without being told.

## Code

```fortran
--8<-- "example/ex01_pair_copula.f90"
```

## Output

```text
Clayton(2) density on the diagonal:    5.3702    1.4810    2.1578
empirical Kendall tau:  0.514
fitted theta:  2.0422   log-likelihood:    892.73
selected: clayton  rotation 0  par   2.042   0.000   0.000
```

The first line is deterministic. The other three depend on the random sample: the program
calls `random_seed()` without arguments, which gfortran seeds from the operating system's
random source, so each run prints slightly different values around the same truth.

## Walk-through

### Defining a copula

```fortran
truth%family = FV_CLAYTON; truth%par = [2.0_dp, 0.0_dp, 0.0_dp]
```

An `fv_bicop` is a plain derived type: `family` is an integer code (named constants
`FV_GAUSSIAN`, `FV_CLAYTON`, ... in `fortvine_api`), `rotation` is 0, 90, 180 or 270, and
`par(3)` holds up to three parameters with unused slots left at 0. Nothing has to be
initialised or allocated. The Clayton copula has the distribution function and density

$$
C(u, v) = \bigl(u^{-\theta} + v^{-\theta} - 1\bigr)^{-1/\theta},
\qquad
c(u, v) = (1 + \theta)\,(uv)^{-1-\theta} \bigl(u^{-\theta} + v^{-\theta} - 1\bigr)^{-2 - 1/\theta},
$$

with Kendall's $\tau = \theta / (\theta + 2)$, i.e. $\tau = 0.5$ for $\theta = 2$.

### Evaluating the density

```fortran
call truth%pdf([0.1_dp, 0.5_dp, 0.9_dp], [0.1_dp, 0.5_dp, 0.9_dp], dens)
```

All evaluation routines take arrays of points: here the three diagonal points
$(0.1, 0.1)$, $(0.5, 0.5)$, $(0.9, 0.9)$. At the centre the formula gives
$3 \cdot 0.25^{-3} \cdot 7^{-5/2} = 1.4810$, as printed. The density is largest near the
origin (5.37 at $(0.1, 0.1)$) because Clayton has lower tail dependence,
$\lambda_L = 2^{-1/\theta} \approx 0.71$; at $(0.9, 0.9)$ it is smaller (2.16) since
$\lambda_U = 0$. The same calling pattern serves `logpdf`, `cdf`, `hfunc1`, `hfunc2`,
`hinv1` and `hinv2`. Long arrays (more than 4096 points) are evaluated in parallel.

### Simulating

```fortran
call random_number(u1); call random_number(w)
call truth%hinv1(u1, w, u2)
```

`hinv1(u1, w, u2)` returns the $v$ that solves $h_1(v \mid u_1) = w$, where
$h_1(v \mid u) = \partial C(u, v)/\partial u$ is the conditional distribution of $V$ given
$U = u$. If $u_1$ and $w$ are independent uniforms, $(u_1, u_2)$ is a draw from the copula:
the inverse Rosenblatt transform in two dimensions. For Clayton, fortvine uses the closed
form

$$
h_1^{-1}(w \mid u) = \Bigl(1 + u^{-\theta}\bigl(w^{-\theta/(1+\theta)} - 1\bigr)\Bigr)^{-1/\theta},
$$

written so that the factor $w^{-\theta/(1+\theta)} - 1$ is computed with `expm1` and nothing
cancels for small $u$. The empirical Kendall tau of the 2000 draws (0.514 with the fixed seed)
scatters around the theoretical 0.5 with a standard error of roughly 0.01.

### Fitting a known family

```fortran
fitted%family = FV_CLAYTON
call fitted%fit(u1, u2)
```

`fit` keeps `family` and `rotation` and estimates `par` by maximum likelihood; the
maximised log-likelihood goes to `fitted%loglik`. Internally the start value comes from
inverting Kendall's tau, and the likelihood is maximised by Fisher-scoring steps with exact
parameter scores (see [Exact parameter scores](../theory/scores.md)). Passing an
`fv_controls` with `method = FV_ITAU` instead returns the tau-inversion estimate itself, and
`weights=` gives observation weights ([example](weights-discrete.md)).

### Selecting a family

```fortran
call chosen%select(u1, u2)
```

`select` with no controls tries the default candidate set: the 23 hand-coded parametric
families `0..22` (every parametric family except the slow-to-fit t-EV). A rotatable
family is tried in the two rotations that match the sign of the sample's Kendall tau
(0 and 180 degrees for positive dependence, 90 and 270 for negative). With `preselect`
on (the default), the correlations of the normal scores in the lower and the upper part of
the sample are compared, as in vinecopulib, and a rotation whose tail sits on the weaker
side is dropped before fitting. Every remaining candidate is fitted and the one with the
smallest BIC,
$-2\ell + k \log n$, wins. The data came from Clayton, and Clayton without rotation is
recovered with the same parameter as the dedicated fit.

## Variations

* Restrict or extend the candidates with a control object:
  ```fortran
  type(fv_controls) :: ctl
  ctl%nfam = 5; ctl%family_set(1:5) = FV_ONE_PAR     ! gaussian, clayton, gumbel, frank, joe
  ctl%selcrit = FV_AIC
  call chosen%select(u1, u2, ctl)
  ```
  See [Controls](../manual/controls.md) for every field.
* Evaluate a rotated copula: set `truth%rotation = 90` before calling `pdf`/`hinv1`; a 90 or
  270 degree rotation of Clayton models negative dependence with a tail in an off-diagonal
  corner ([Rotations and asymmetry](../manual/rotations.md)).
* `cop%cdf(u1, u2, out)` returns $C$ itself, closed form for every family except the
  Gaussian and Student ones, which are integrated by Gauss-Legendre quadrature over their
  exact h-function.
