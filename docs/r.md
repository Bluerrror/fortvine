# R package

The R package `fortvine` lets you use the library from R. It is a binding, not a port: every
density, h-function, fit, family selection, vine selection, Rosenblatt transform and
pseudo-observation is computed by the Fortran code through its [C interface](manual/c-api.md).
The R functions only check and convert their arguments (an `n x d` R matrix becomes the
`u(d, n)` array of the C routines, family names become integer codes) and wrap the results in
R objects with `print`, `summary` and `logLik` methods. The function names follow
[rvinecopulib](https://github.com/vinecopulib/rvinecopulib), so code written for it needs few
changes.

## Installation

The package compiles the library from source with R's own toolchain, so you need a Fortran
compiler that R can use:

| platform | what to install |
|---|---|
| Windows | [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version (Rtools43 for R 4.3, Rtools44 for R 4.4, ...). It contains gfortran. |
| Linux | gfortran, e.g. `sudo apt-get install gfortran` (Debian, Ubuntu) or `sudo dnf install gcc-gfortran` (Fedora) |
| macOS | the gfortran build that CRAN recommends for your R version ([mac.r-project.org/tools](https://mac.r-project.org/tools/)) |

Then, from R:

```r
install.packages("remotes")
remotes::install_github("Bluerrror/fortvine", subdir = "r")
```

or, from a clone of the repository:

```bash
R CMD INSTALL r
```

OpenMP is used when R's compiler supports it; otherwise the library runs on one thread and
gives the same results. The `cores` argument of each function sets the number of threads
(default 1, as CRAN expects of R packages). Selection over large family sets is several times
faster with more threads, e.g. `cores = parallel::detectCores() %/% 2`.

## Quick start

A pair copula: build a model, evaluate the density and the h-function, simulate, fit and
select.

```r
library(fortvine)
cop <- bicop_dist("clayton", 0, 2)
u <- rbind(c(0.1, 0.1), c(0.5, 0.5), c(0.9, 0.9))
dbicop(u, cop)
#> [1] 5.370177 1.481004 2.157801
hbicop(u, 1, cop)
#> [1] 0.3562217 0.4319594 0.7703349

set.seed(42)
x <- rbicop(2000, cop)
bicop(x, family_set = "clayton", rotation = 0)$parameters   # fit one family and rotation
#> [1] 1.974684
bicop(x, family_set = "all")                                # select among families 0-22
#> Bivariate copula ('bicop_dist'): family = clayton (3), rotation = 0, parameters = 1.975
```

Families are given by their fortvine name (any case; `"t"` also means `"student"`) or code.
`fortvine_families()` lists the 111 families with their codes and numbers of parameters,
read from the library. The generic families (codes $\ge 10000$) work everywhere a family is
accepted:

```r
fam <- fortvine_families()
table(generic = fam$code >= 10000)
#> generic
#> FALSE  TRUE
#>    26    85
bicop(x, family_set = "generic", cores = 2)
#> Bivariate copula ('bicop_dist'): family = Arch(gamma.Sibuya) (11700), rotation = 0, parameters = 1.936, 1.064
```

## Functions

| R function | what it does | computed by |
|---|---|---|
| `bicop_dist(family, rotation, parameters)` | pair-copula model object | (checks only) |
| `dbicop(u, model, log = FALSE)` | density $c(u_1, u_2)$ or $\log c$ | `fv_bicop_eval` |
| `hbicop(u, cond_var, model, inverse = FALSE)` | $h_1 = \partial C / \partial u_1$, $h_2 = \partial C / \partial u_2$ and their inverses | `fv_bicop_eval` |
| `rbicop(n, model)` | simulation $\big(u_1, h_1^{-1}(u_1, w)\big)$ | `fv_bicop_eval` |
| `bicop(data, family_set, selcrit, rotation = NULL)` | select family, rotation and parameters; or, with `rotation`, fit one family | `fv_bicop_select`, `fv_bicop_fit` |
| `dvine_structure(order, trunc_lvl)` | D-vine structure | `fv_vinecop%set_dvine` |
| `vinecop_dist(pair_copulas, structure)` | vine model object from pair copulas | (checks only) |
| `vinecop(data, family_set, structure, selcrit, trunc_lvl, threshold, ...)` | Dissmann selection, or families on a fixed structure | `fv_vine_select` |
| `dvinecop(u, vinecop)` | vine density | `fv_vine_pdf` |
| `rosenblatt(x, model)`, `inverse_rosenblatt(u, model)` | Rosenblatt transform and its inverse | `fv_vine_rosenblatt`, `fv_vine_inverse_rosenblatt` |
| `rvinecop(n, vinecop)` | simulation: inverse Rosenblatt of uniforms | `fv_vine_inverse_rosenblatt` |
| `pseudo_obs(x)` | ranks / $(n + 1)$, ties averaged | `fv_pseudo_obs` |
| `fortvine_families()` | the 111 families: code, name, number of parameters | `fv_fam_name`, `fv_generic_families` |
| `logLik()`, `AIC()`, `BIC()`, `print()`, `summary()` | methods for fitted `bicop_dist` and `vinecop_dist` objects | (stored results) |

`family_set` takes names, codes and three groups: `"all"` (codes 0-22, the default set of the
library), `"parametric"` (0-11, the parametric families of vinecopulib) and `"generic"` (the 85
generic families). `selcrit` is one of `"loglik"`, `"aic"`, `"bic"`, `"mbic"` and `"mbicv"`.

**Random numbers.** `rbicop()` and `rvinecop()` draw their uniforms with R's `runif()`, so
`set.seed()` makes them reproducible. They do not use the counter-based generator behind the
`seed` argument of `fv_vinecop%simulate`, so the same seed gives different samples in R and
in Fortran.

**Not available from R.** The C interface has no routine for these, so the package leaves
them out: the distribution function $C(u_1, u_2)$ (`pbicop` in rvinecopulib), the
nonparametric families `tll` (24) and `bernstein` (25) (requesting them stops with an error),
observation weights, discrete margins, automatic truncation and threshold selection, and
non-simplified (relaxed) vines.

## Worked example

Daily log returns of four European stock indices, 1991-1998 (`EuStockMarkets`, shipped with
R): first a pair copula for DAX and CAC, then a vine for all four indices.

```r
library(fortvine)
x <- diff(log(EuStockMarkets))          # daily log returns of DAX, SMI, CAC, FTSE
u <- pseudo_obs(x)
dim(u)
#> [1] 1859    4

pair <- bicop(u[, c("DAX", "CAC")], family_set = "all", cores = 2)
summary(pair)
#>   family code rotation    parameters npars   loglik
#> 1    bb1    7      180 0.3035, 1.771     2 709.9664
c(AIC = AIC(pair), BIC = BIC(pair))
#>       AIC       BIC
#> -1415.933 -1404.877
```

The DAX-CAC pair is a survival BB1 copula (rotation 180), which has a stronger lower tail
than upper tail: large joint losses are more common than large joint gains. The vine:

```r
fit <- vinecop(u, family_set = "all", cores = 2)
fit
#> Vine copula ('vinecop_dist'): d = 4, trunc_lvl = 3, order = 2 1 3 4
#>   npars = 11, loglik = 2038.51 (n = 1859)
#>  tree edge conditioned conditioning   family rotation    parameters
#>     1    1    SMI, DAX                   bb3        0 1.577, 0.2891
#>     1    2    DAX, CAC                   bb1      180 0.3035, 1.771
#>     1    3   CAC, FTSE                   bb1      180 0.2627, 1.607
#>     2    1    SMI, CAC          DAX  student        0 0.2101,  11.9
#>     2    2   DAX, FTSE          CAC  student        0 0.3204, 10.88
#>     3    1   SMI, FTSE     DAX, CAC plackett        0          1.87
```

Tree 1 is the path SMI - DAX - CAC - FTSE. `summary(fit)` returns the same table as a data
frame, with the log-likelihood of each edge. The vine log-likelihood is the sum of the edge
log-likelihoods, and the density routine reproduces it:

```r
c(loglik = fit$loglik, sum_log_density = sum(log(dvinecop(u, fit))), npars = fit$npars)
#>          loglik sum_log_density           npars
#>        2038.507        2038.507          11.000
```

If the model fits, the Rosenblatt transform $w = T(u)$ has independent uniform columns, and the
inverse transform recovers $u$:

```r
w <- rosenblatt(u, fit)
round(cor(w, method = "kendall"), 3)
#>       [,1]   [,2]  [,3]   [,4]
#> [1,] 1.000  0.008 0.011  0.003
#> [2,] 0.008  1.000 0.007 -0.004
#> [3,] 0.011  0.007 1.000  0.002
#> [4,] 0.003 -0.004 0.002  1.000
max(abs(inverse_rosenblatt(w, fit) - u))
#> [1] 6.448175e-13

set.seed(1)
round(rvinecop(3, fit), 4)
#>        [,1]   [,2]   [,3]   [,4]
#> [1,] 0.2132 0.6958 0.6235 0.0618
#> [2,] 0.2926 0.1680 0.4114 0.2060
#> [3,] 0.3630 0.6858 0.3562 0.1766
```

The figure is drawn by `tools/r_example.R` with the package. Panel (a) is the density of the
selected DAX-CAC pair copula in normal scores, $c\big(\Phi(z_1), \Phi(z_2)\big)\,
\varphi(z_1)\varphi(z_2)$; panel (b) compares the DAX-FTSE returns with a sample of the same
size simulated from the vine (DAX and FTSE are not neighbours in tree 1, so this pair is
modelled through the higher trees); panel (c) shows two columns of the Rosenblatt transform.

<figure markdown>
![Density of the DAX-CAC pair copula, observed and simulated DAX-FTSE pairs, and the Rosenblatt transform](assets/r_example.png){ width="700" }
<figcaption>The EuStockMarkets vine fitted from R: (a) the selected DAX-CAC pair copula, (b) data
and vine simulation, (c) Rosenblatt-transformed data.</figcaption>
</figure>

## From rvinecopulib

The names and the main arguments are those of rvinecopulib. The differences:

* families use fortvine names and codes (`fortvine_families()`), and `family_set = "all"`
  means the hand-coded parametric families 0-22;
* vine structures are stored in the natural-order arrays of fortvine (`order`,
  `struct_array[edge, tree]`, see [Vine copulas](manual/vinecop.md#how-the-structure-is-stored)),
  not as an R-vine matrix;
* `bicop()` accepts `rotation` to fit a single family with a fixed rotation;
* there is no `pbicop()`, no `var_types`, no `weights` and no nonparametric family (see above).

## Citing

If you use fortvine, from R or otherwise, please cite:

> Shojaeezadeh, S. A. (2026). fortvine: a complete compendium to work on vine copula modeling.
> In preparation. Software: <https://github.com/Bluerrror/fortvine>

In R, `citation("fortvine")` prints this reference and a BibTeX entry.
