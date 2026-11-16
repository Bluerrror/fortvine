# ex01_pair_copula -- one bivariate copula: evaluate, simulate, fit a known family, select a family.
# The R version of example/ex01_pair_copula.f90 (R CMD INSTALL r).
library(fortvine)
n <- 2000

# a Clayton copula with theta = 2 (Kendall's tau = theta / (theta + 2) = 0.5)
truth <- bicop_dist("clayton", 0, 2)

# density at three points
g <- c(0.1, 0.5, 0.9)
dens <- dbicop(cbind(g, g), truth)
cat("Clayton(2) density on the diagonal:", sprintf("%10.4f", dens), "\n", sep = "")

# simulate by inverse Rosenblatt: u2 = h^-1(w | u1) with independent uniforms u1, w
set.seed(20261116)                                   # reproducible sample
u1 <- runif(n); w <- runif(n)
u2 <- hbicop(cbind(u1, w), 1, truth, inverse = TRUE)
u <- cbind(u1, u2)
cat(sprintf("empirical Kendall tau:%7.3f\n", cor(u1, u2, method = "kendall")))

# maximum likelihood fit of a given family
fitted <- bicop(u, family_set = "clayton", rotation = 0)
cat(sprintf("fitted theta:%8.4f   log-likelihood:%10.2f\n", fitted$parameters, fitted$loglik))

# family selection by BIC over the default family set (with rotations)
chosen <- bicop(u, family_set = "all", selcrit = "bic")
par <- c(chosen$parameters, 0, 0, 0)[1:3]
cat("selected: ", chosen$family, "  rotation ", chosen$rotation, "  par",
    sprintf("%8.3f", par), "\n", sep = "")
