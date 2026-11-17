# Basic checks of the R binding; base R only. Every number is computed by the Fortran library.
library(fortvine)
set.seed(20261116)

# the family table of the library: 26 hand-coded + 85 generic families
fams <- fortvine_families()
stopifnot(nrow(fams) == 111, !anyDuplicated(fams$code), !anyDuplicated(tolower(fams$name)),
          identical(fams$name[1:4], c("indep", "gaussian", "student", "clayton")))
stopifnot(inherits(try(bicop_dist("tll"), silent = TRUE), "try-error"))

# h-function and inverse h-function round trip, both conditioning variables, rotated family
cop <- bicop_dist("bb1", 90, c(0.8, 1.6))
u <- matrix(runif(400, 0.01, 0.99), ncol = 2)
for (cv in 1:2) {
  h <- hbicop(u, cv, cop)
  back <- if (cv == 1) hbicop(cbind(u[, 1], h), 1, cop, inverse = TRUE) else
    hbicop(cbind(h, u[, 2]), 2, cop, inverse = TRUE)
  stopifnot(max(abs(back - u[, 3 - cv])) < 1e-7)
}

# the density integrates to about 1 (midpoint rule on a 200 x 200 grid)
g <- (seq_len(200) - 0.5) / 200
grid <- as.matrix(expand.grid(g, g))
stopifnot(abs(mean(dbicop(grid, "frank", 0, 4)) - 1) < 1e-3,
          abs(mean(dbicop(grid, "gumbel", 180, 1.5)) - 1) < 1e-2)
stopifnot(all.equal(log(dbicop(grid[1:50, ], "joe", 0, 2)),
                    dbicop(grid[1:50, ], "joe", 0, 2, log = TRUE)))

# fitting recovers the parameter; selection finds the true family and rotation
x <- rbicop(3000, "clayton", 0, 3)
fit <- bicop(x, family_set = "clayton", rotation = 0)
stopifnot(abs(fit$parameters - 3) < 0.25, fit$loglik > 0)
sel <- bicop(x, family_set = "parametric", cores = 2)
stopifnot(sel$family == "clayton", sel$rotation == 0,
          inherits(logLik(sel), "logLik"), attr(logLik(sel), "df") == 1)
x <- rbicop(3000, "gumbel", 270, 2)
sel <- bicop(x, family_set = c("gaussian", "clayton", "gumbel", "frank", "joe"))
stopifnot(sel$family == "gumbel", sel$rotation == 270, abs(sel$parameters - 2) < 0.2)

# a vine: simulate a D-vine, select a vine, loglik = sum of log density, Rosenblatt round trip
pcs <- list(list(bicop_dist("gumbel", 0, 2), bicop_dist("clayton", 0, 2),
                 bicop_dist("gaussian", 0, 0.5)),
            list(bicop_dist("frank", 0, 3), bicop_dist("frank", 0, -2)),
            list(bicop_dist("indep")))
truth <- vinecop_dist(pcs, dvine_structure(c(2, 1, 4, 3)))
u <- rvinecop(1500, truth)
vc <- vinecop(u, family_set = "parametric", cores = 2)
stopifnot(abs(vc$loglik - sum(log(dvinecop(u, vc)))) < 1e-6 * abs(vc$loglik),
          as.numeric(logLik(vc)) == vc$loglik, nrow(summary(vc)) == 6)
w <- rosenblatt(u, vc)
stopifnot(max(abs(inverse_rosenblatt(w, vc) - u)) < 1e-8,
          max(abs(cor(w, method = "kendall")[upper.tri(diag(4))])) < 0.06)
# a fixed structure: only families and parameters are selected. The library re-encodes the
# structure, so compare the edges as unordered sets (a D-vine may come back in reverse order).
edges <- function(v) sort(vapply(strsplit(paste(summary(v)$conditioned, summary(v)$conditioning),
                                          "[, ]+"), function(s) paste(sort(s), collapse = ""), ""))
vd <- vinecop(u, family_set = "parametric", structure = truth$structure)
stopifnot(identical(edges(vd), edges(truth)), is.finite(vd$loglik))

# one generic family (code >= 10000): density, h round trip, fit
gen <- fams[fams$code >= 10000 & fams$npars == 2, ][1, ]
lb <- bicop(rbicop(2000, "frank", 0, 5), family_set = gen$code, rotation = 0)
stopifnot(lb$code == gen$code, length(lb$parameters) == 2, is.finite(lb$loglik))
gcop <- bicop_dist(gen$name, 0, lb$parameters)
stopifnot(abs(mean(dbicop(grid, gcop)) - 1) < 1e-2)
h <- hbicop(u[, 1:2], 1, gcop)
stopifnot(max(abs(hbicop(cbind(u[, 1], h), 1, gcop, inverse = TRUE) - u[, 2])) < 1e-6)

# pseudo-observations: average ranks / (n + 1)
stopifnot(all.equal(pseudo_obs(c(3.1, 0.2, 5, 0.2)), c(3, 1.5, 4, 1.5) / 5))

# non-simplified vine: the (1, 3 | 2) correlation grows with u2
set.seed(6)
ur <- matrix(runif(1500), 500, 3)
rr <- 0.9 * (2 * ur[, 2] - 1)
ur[, 3] <- pnorm(rr * qnorm(ur[, 1]) + sqrt(1 - rr^2) * qnorm(ur[, 3]))
fr <- vinecop(ur, family_set = "gaussian", structure = dvine_structure(1:3), cores = 2)
rel <- vinecop_relax(fr, ur, mode = 2, cores = 2)
p <- vinecop_rel_params(rel, ur)
stopifnot(!is.null(rel$relaxed), dim(p) == c(500, ncol(rel$relaxed$erel), 3),
          cor(p[, 1, 1], ur[, 2]) > 0.9, rel$loglik > fr$loglik,
          abs(rel$loglik - sum(log(dvinecop(ur, rel)))) < 1e-6)
un <- vinecop_unrelax(rel, ur)
stopifnot(is.null(un$relaxed), un$npars < rel$npars, is.finite(un$loglik))
cat("fortvine tests passed\n")
