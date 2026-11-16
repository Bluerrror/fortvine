# ex02_vine -- a 5-dimensional vine copula: simulate from a known D-vine, select a new vine from
# the sample (structure, families, parameters), then density, Rosenblatt transform and simulation.
# The R version of example/ex02_vine.f90 (R CMD INSTALL r).
library(fortvine)
n <- 3000; d <- 5

# a D-vine 1-2-3-4-5 with Gumbel pairs in tree 1 and Frank pairs in tree 2, independence above
pcs <- list(rep(list(bicop_dist("gumbel", 0, 2)), d - 1),
            rep(list(bicop_dist("frank", 0, 3)), d - 2))
truth <- vinecop_dist(pcs, dvine_structure(1:d))
set.seed(42)
u <- rvinecop(n, truth)

# select a regular vine: Dissmann's algorithm, all parametric families, BIC
vc <- vinecop(u, family_set = "all", selcrit = "bic", cores = 4)
print(vc)
cat(sprintf("log-likelihood:%12.2f   parameters:%d\n", sum(log(dvinecop(u, vc))), vc$npars))

# density, Rosenblatt transform (independent uniforms if the model is right) and new samples
dens <- dvinecop(u, vc)
w <- rosenblatt(u, vc)
cat(sprintf("Kendall tau of the first two Rosenblatt columns (should be ~0):%8.4f\n",
            cor(w[, 1], w[, 2], method = "kendall")))
set.seed(7)
sim <- rvinecop(10, vc)
cat("first simulated row:", sprintf("%8.4f", sim[1, ]), "\n", sep = "")
for (t in 1:2)
  cat("tree ", t, " families: ",
      paste(vapply(vc$pair_copulas[[t]], function(b) b$family, ""), collapse = " "), "\n", sep = "")
