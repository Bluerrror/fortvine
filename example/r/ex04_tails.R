# ex04_tails -- tail dependence of five pair copulas, three ways: the limit from the family
# formula, the exact coefficient at a finite level, and the empirical coefficient of a simulated
# sample.  Then: which family does BIC select from each sample?
# The R version of example/ex04_tails.f90 (R CMD INSTALL r).
#
# The R package has no copula distribution function C (the C interface does not export it), so
# the finite-level coefficients are computed from the h-function h1(s, v) = dC(s, v)/ds:
#   lower tail  lambda_L(q) = C(q, q) / q = (1/q) * integral_0^q h1(s, q) ds
#   upper tail  lambda_U(q) = P(U1 > 1-q, U2 > 1-q) / q = (1/q) * integral_{1-q}^1 (1 - h1(s, 1-q)) ds
library(fortvine)
n <- 200000; m <- 2000

# the five copulas (parameters as in the Fortran program)
cops <- list(bicop_dist("gumbel", 0, 2),
             bicop_dist("clayton", 0, 2),
             bicop_dist("clayton", 180, 2),
             bicop_dist("tawn", 0, c(0.4, 0.9, 3)),       # psi1, psi2, theta
             bicop_dist(20300, 0, c(3, 0.5, 0.9)))        # Khoudraji-Joe: theta, a, b

# limits (lambda_L, lambda_U) from the family formulas
th <- 2
lim <- list(c(0, 2 - 2^(1 / th)),                         # Gumbel
            c(2^(-1 / th), 0),                            # Clayton
            c(0, 2^(-1 / th)),                            # survival Clayton
            c(0, 0.4 + 0.9 - (0.4^3 + 0.9^3)^(1 / 3)),    # Tawn
            c(0, 0.5 + 0.9 - (0.5^3 + 0.9^3)^(1 / 3)))    # Khoudraji-Joe

# n draws of (u1, u2) from cop by the inverse Rosenblatt transform, reproducible from seed
simulate <- function(cop, seed, size = n) {
  set.seed(seed)
  u1 <- runif(size); w <- runif(size)
  cbind(u1, hbicop(cbind(u1, w), 1, cop, inverse = TRUE))  # u2 = h1^-1(w | u1)
}

h1 <- function(s, v, cop) hbicop(cbind(s, v), 1, cop)
lam_lower <- function(cop, q)
  integrate(h1, 0, q, v = q, cop = cop, rel.tol = 1e-10)$value / q
lam_upper <- function(cop, q)                              # q is the distance from 1
  integrate(function(s) 1 - h1(s, 1 - q, cop), 1 - q, 1, rel.tol = 1e-10)$value / q

f8 <- function(x) paste(sprintf("%8.3f", x), collapse = "")
cat("                        |      lower tail lambda_L(q)         |      upper tail lambda_U(q)\n")
cat("copula              rot |      limit  q=1e-2  q=1e-4  sample |      limit  q=1e-2  q=1e-4  sample\n")
for (k in seq_along(cops)) {
  cop <- cops[[k]]
  u <- simulate(cop, 1000 + k)
  low <- mean(u[, 1] < 0.01 & u[, 2] < 0.01) / 0.01
  upp <- mean(u[, 1] > 0.99 & u[, 2] > 0.99) / 0.01
  cat(sprintf("%-20s%3d |   ", cop$family, cop$rotation),
      f8(c(lim[[k]][1], lam_lower(cop, 1e-2), lam_lower(cop, 1e-4), low)), " |   ",
      f8(c(lim[[k]][2], lam_upper(cop, 1e-2), lam_upper(cop, 1e-4), upp)), "\n", sep = "")
}

# BIC selection from m observations of each copula: default family set (hand-coded parametric
# families with rotations), then the same plus every generic family
cat("\nBIC selection from 2000 observations:\n")
for (k in seq_along(cops)) {
  cop <- cops[[k]]
  u <- simulate(cop, 1000 + k)[1:m, ]
  d <- bicop(u, family_set = "all", cores = 4)
  g <- bicop(u, family_set = list(0:23, "generic"), cores = 4)
  cat(sprintf("  %-14s default set: %-10s rot%3d   + generic families: %-14s rot%3d\n",
              cop$family, d$family, d$rotation, g$family, g$rotation))
}
