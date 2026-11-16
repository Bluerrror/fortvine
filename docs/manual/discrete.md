# Discrete and mixed data

Counts, ordinal scores, rounded measurements and zero-inflated variables have marginal
distributions with jumps. A copula still links such margins (Sklar's theorem holds), but the
probability integral transform $F(X)$ is no longer uniform, and the likelihood of an
observation is a probability mass rather than a density. fortvine handles this with the
convention of vinecopulib's `var_types = "d"` (Panagiotelis, Czado & Joe 2012; Panagiotelis,
Czado, Joe & Stöber 2017): each discrete variable is described by two numbers per
observation.

!!! note "Fortran only"
    Discrete margins are available from Fortran only; the [Python](../python.md) and
    [R](../r.md) packages do not expose `discrete` and `u_minus`.

## Input: `u` and `u_minus`

For a discrete variable $X_j$ with distribution function $F_j$:

$$
u_{ij} = F_j(x_{ij}), \qquad u^-_{ij} = F_j(x_{ij}^-) = P(X_j < x_{ij}),
$$

so that $u_{ij} - u^-_{ij} = P(X_j = x_{ij})$. For integer data $F_j(x^-) = F_j(x - 1)$.
Pass

* `u(n, d)`: $F(x)$ for every column (continuous columns as usual),
* `u_minus(n, d)`: $F(x^-)$ for the discrete columns; entries of continuous columns are not
  read (copying `u` there is fine),
* `discrete(d)`: `.true.` for the discrete columns.

```fortran
call vc%select(u, ctl, discrete=[.false., .false., .true.], u_minus=um)
ll = vc%loglik(u, u_minus=um)
```

$F_j$ can be a fitted parametric margin or the empirical distribution function; with the
latter, using the denominator $n + 1$ (as `fv_pseudo_obs_nd` does for continuous columns)
keeps $u < 1$.

The discreteness flags are stored in `vc%disc`. Every evaluation of a vine with any discrete
variable then requires `u_minus` (`pdf`, `loglik`, `rosenblatt`) or `ugm` (`pdf_grid`).
`fv_bicop` methods have no `u_minus` argument; to fit or select a single pair with discrete
margins, use a two-dimensional `fv_vinecop`.

## The likelihood of a pair

Let $a, b$ be the two (conditional) margins of an edge, with $a^-, b^-$ their left limits.
The pair-copula likelihood used by fortvine is the joint probability divided by the product
of the marginal masses, so that it is a density with respect to the product of the margins
and directly comparable between candidates:

| margins of the observation | contribution |
|---|---|
| both continuous | $c(a, b)$ |
| $a$ discrete, $b$ continuous | $\dfrac{h_2(a, b) - h_2(a^-, b)}{a - a^-}$ |
| $a$ continuous, $b$ discrete | $\dfrac{h_1(a, b) - h_1(a, b^-)}{b - b^-}$ |
| both discrete | $\dfrac{C(a,b) - C(a^-,b) - C(a,b^-) + C(a^-,b^-)}{(a - a^-)(b - b^-)}$ |

A variable declared discrete whose observation has $a - a^- \le 10^{-12}$ is treated as
continuous for that observation, so mixed (partly continuous, partly atomic) margins such as
zero-inflated data are handled row by row.

## Conditional margins in higher trees

The next tree needs $F(b \mid a)$ and its left limit. If the conditioning variable $a$ is
discrete, the conditional distribution averages over its jump:

$$
F(b \mid a) = \frac{C(a, b) - C(a^-, b)}{a - a^-} = \frac{1}{a - a^-}\int_{a^-}^{a} h_1(s, b)\, ds .
$$

fortvine computes this from the closed-form CDF when the jump is at least $10^{-3}$ and the
family has one, and by Gauss-Legendre quadrature of the h-function otherwise (the Student
copula always uses quadrature). Both $F(b \mid a)$ and $F(b^- \mid a)$ are propagated, and a
conditional margin is discrete exactly when the variable it describes is.

## Fitting with discrete margins

* **Gaussian, Clayton, Gumbel, Frank, Joe:** the exact discrete likelihood with its score and
  curvature in closed form, maximised by Newton's method (as on continuous data).
* **Other one-parameter families:** Brent's method on the exact discrete likelihood.
* **Multi-parameter families:** the start grid followed by Nelder-Mead in logit coordinates
  on the exact discrete likelihood (value only: no dual-number scores on this path). The
  two-stage subsample scheme applies for $n > 1800$.
* **TLL, Bernstein:** the estimate is computed from the mid-points
  $\big(F(x) + F(x^-)\big)/2$ of the discrete margins, and its log-likelihood is then
  evaluated with the discrete formulas above.

Discrete fits are slower than continuous ones, mainly because the multi-parameter families
lose the exact scores and because every observation needs two to four copula evaluations. The
generic-family pre-screen is skipped on discrete edges, and [`relax`](vinecop.md#non-simplified-edges-relax)
is not available for vines with discrete variables.

## Rosenblatt transform

For a discrete variable the conditional distribution function $F(x_k \mid \cdot)$ evaluated
at the data is not uniform. `rosenblatt` therefore returns, by default, Brockwell's (2007)
randomised transform for the discrete columns,

$$
w_{ik} = W_{ik}\, F(x_{ik} \mid \dots) + (1 - W_{ik})\, F(x_{ik}^- \mid \dots),
$$

with independent uniforms $W_{ik}$; under the model the $w_{ik}$ are iid uniform, so the usual
goodness-of-fit diagnostics apply. Continuous columns are transformed as usual.

```fortran
call vc%rosenblatt(u, w, u_minus=um)                          ! randomised, seed 1
call vc%rosenblatt(u, w, u_minus=um, seed=7)                  ! another randomisation
call vc%rosenblatt(u, w, u_minus=um, randomize=.false., out_minus=wm)   ! F(x|.), F(x-|.)
```

The uniforms come from a counter-based hash of `(seed, row, column)`, so the result is
reproducible and independent of the number of threads and of the Fortran random generator.
`randomize = .false.` returns $F(x_k \mid \dots)$ itself, and `out_minus(n, d)` (optional)
receives $F(x_k^- \mid \dots)$.

!!! note "Simulation"
    `inverse_rosenblatt` and `simulate` do not use the discreteness flags: they return a sample
    $U$ of the copula on the continuous scale. Discrete observations are obtained by applying
    the marginal quantile functions, $x_j = F_j^{-1}(u_j)$, which is the standard latent-variable
    construction.

## Example

A latent Gaussian D-vine whose third variable is observed only as a count with four levels:

```fortran
program discrete_demo
  use fortvine_api
  implicit none
  integer, parameter :: dp = kind(1.0d0), n = 2000, d = 3
  type(fv_vinecop) :: gen, vc
  type(fv_controls) :: ctl
  real(dp) :: v(n, d), u(n, d), um(n, d), r(n, d)
  integer :: k(n), j, i

  ! latent Gaussian-copula sample; variable 3 is then made a count with 4 levels
  call gen%set_dvine([1, 2, 3])
  gen%fam(1:2, 1) = FV_GAUSSIAN; gen%par(1, 1:2, 1) = 0.7_dp
  call gen%simulate(n, v, seed=2)
  k = min(int(4.0_dp * v(:, 3)**1.5_dp), 3)            ! counts 0..3

  ! continuous columns: u_minus = u (not used); discrete column: empirical F(x) and F(x-)
  u = v; um = v
  do i = 1, n
    u(i, 3) = count(k <= k(i)) / (n + 1.0_dp)
    um(i, 3) = count(k < k(i)) / (n + 1.0_dp)
  end do

  ctl%nthreads = 4
  call vc%select(u, ctl, discrete=[.false., .false., .true.], u_minus=um)
  call vc%print()
  print '(a,f10.2)', 'discrete log-likelihood:', vc%loglik(u, u_minus=um)

  ! randomised (Brockwell) Rosenblatt transform: iid uniform under the model
  call vc%rosenblatt(u, r, u_minus=um, seed=3)
  do j = 1, d
    print '(a,i0,a,f7.4,a,f7.4)', 'column ', j, ': mean', sum(r(:, j)) / n, &
      '  variance', sum((r(:, j) - 0.5_dp)**2) / n
  end do
  print '(a,f8.4)', 'tau(r1, r3):', fv_kendall_tau(r(:, 1), r(:, 3))
end program discrete_demo
```

Output:

```text
fortvine vine copula: d = 3, trunc = 2
order: 1 2 3
  tree 1 edge 1:     gaussian rot   0 par  6.9863E-01  0.0000E+00  0.0000E+00
  tree 1 edge 2:     gaussian rot   0 par  6.9748E-01  0.0000E+00  0.0000E+00
  tree 2 edge 1:        indep rot   0 par  0.0000E+00  0.0000E+00  0.0000E+00
discrete log-likelihood:   1112.53
column 1: mean 0.4918  variance 0.0830
column 2: mean 0.4974  variance 0.0798
column 3: mean 0.4993  variance 0.0831
tau(r1, r3): -0.0190
```

Both latent correlations (true value 0.7) are recovered although one of the variables has
only four distinct values, and the transformed columns have the mean $1/2$ and variance
$1/12 \approx 0.0833$ of independent uniforms. Treating the count as continuous instead (its
mid-ranks $(u + u^-)/2$ passed as ordinary pseudo-observations) gives 0.679 for the second
parameter on the same data: the ties attenuate the dependence estimate, and the effect grows
as the number of distinct values shrinks.
