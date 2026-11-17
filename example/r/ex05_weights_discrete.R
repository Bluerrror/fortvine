# ex05_weights_discrete.R -- the R version of example/ex05_weights_discrete.f90:
# (1) observation weights: a time-decay weighted fit that follows a change in dependence, and
# integer weights = duplicated rows; (2) a vine with a discrete (Poisson count) variable.
library(fortvine)
n <- 3000; LAM <- 2
set.seed(11)

# ---------------------------------------------------------------- (1) weights
u1 <- runif(n); w <- runif(n); h <- 1:(n / 2)
u2 <- c(hbicop(cbind(u1, w)[h, ], 1, "clayton", 0, 1, inverse = TRUE),      # tau 1/3
        hbicop(cbind(u1, w)[-h, ], 1, "clayton", 0, 4, inverse = TRUE))     # tau 2/3
u <- cbind(u1, u2)
fit <- bicop(u, family_set = "clayton", rotation = 0)
cat(sprintf("Clayton fit, equal weights:          theta = %6.3f\n", fit$parameters))
wt <- 0.998^(n - seq_len(n))                     # exponential time decay, half-life ~350 rows
fit <- bicop(u, family_set = "clayton", rotation = 0, weights = wt)
cat(sprintf("Clayton fit, time-decay weights:     theta = %6.3f   (effective n = %6.1f)\n",
            fit$parameters, sum(wt)^2 / sum(wt^2)))

wt <- ifelse(runif(n) > 0.5, 2, 1)               # integer weights = repeated rows
fit <- bicop(u, family_set = "clayton", rotation = 0, weights = wt)
cat(sprintf("weights 1/2:     theta = %.8f  loglik = %.4f\n", fit$parameters, fit$loglik))
fit <- bicop(u[rep(seq_len(n), wt), ], family_set = "clayton", rotation = 0)
cat(sprintf("duplicated rows: theta = %.8f  loglik = %.4f\n", fit$parameters, fit$loglik))

# ---------------------------------------------------------------- (2) a discrete variable
# D-vine 1-2-3: gaussian(0.7) and clayton(2) in tree 1, frank(2) in tree 2
pcs <- list(list(bicop_dist("gaussian", 0, 0.7), bicop_dist("clayton", 0, 2)),
            list(bicop_dist("frank", 0, 2)))
truth <- vinecop_dist(pcs, dvine_structure(1:3))
set.seed(5)
z <- rvinecop(n, truth)

# variable 2 becomes a Poisson(2) count; the data are F(x) of all three, then F(x-) of the count
k <- qpois(z[, 2], LAM)
u <- cbind(z[, 1], ppois(k, LAM), z[, 3], ppois(k - 1, LAM))
cat(sprintf("\ncount variable: %d zeros, maximum %d, mean %.3f\n", sum(k == 0), max(k), mean(k)))

vd <- vinecop(u, structure = dvine_structure(1:3), var_types = c("c", "d", "c"))   # discrete-aware
umid <- u[, 1:3]; umid[, 2] <- (u[, 2] + u[, 4]) / 2                               # naive: mid-point
vn <- vinecop(umid, structure = dvine_structure(1:3))
cat("discrete-aware fit:\n"); print(summary(vd)[, c("tree", "conditioned", "family", "parameters")])
cat("naive fit:\n"); print(summary(vn)[, c("tree", "conditioned", "family", "parameters")])

truth <- vinecop_dist(pcs, dvine_structure(1:3), var_types = c("c", "d", "c"))
cat(sprintf("loglik of the mixed data: truth %9.2f   discrete-aware fit %9.2f\n",
            sum(log(dvinecop(u, truth))), vd$loglik))
cat("first three mixed densities:", signif(dvinecop(u[1:3, ], vd), 4), "\n")

set.seed(1)
r <- rosenblatt(u, vd)                            # randomised Rosenblatt transform (Brockwell 2007)
cat("Rosenblatt column means (1/2):      ", round(colMeans(r), 4), "\n")
cat("Rosenblatt column variances (1/12): ", round(apply(r, 2, var), 4), "\n")
tc <- cor(r, method = "kendall")
cat("Kendall tau (1,2) (1,3) (2,3), ~0:  ", round(tc[cbind(c(1, 1, 2), c(2, 3, 3))], 4), "\n")
