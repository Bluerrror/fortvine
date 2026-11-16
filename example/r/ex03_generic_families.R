# ex03_generic_families -- the generic engine (codes 10000+): list the families, use an Archimax
# copula and an asymmetric Khoudraji copula, and select among classic + generic families.
# The R version of example/ex03_generic_families.f90 (R CMD INSTALL r).
library(fortvine)
n <- 3000

fam <- fortvine_families()
gen <- fam[fam$code >= 10000, ]                                # the 85 new families
cat("generic families: ", nrow(gen), "; the first five:\n", sep = "")
for (k in 1:5)
  cat(sprintf("%7d  %-40s  (%d parameters)\n", gen$code[k], gen$name[k], gen$npars[k]))

# Archimax(gamma | Husler-Reiss) = code 10203: gamma (Clayton) generator, Husler-Reiss Pickands
set.seed(20261116)                                             # reproducible sample
u1 <- runif(n); w <- runif(n)
truth <- bicop_dist(10203, 0, c(1.0, 1.5))
u <- cbind(u1, hbicop(cbind(u1, w), 1, truth, inverse = TRUE))
fitted <- bicop(u, family_set = 10203, rotation = 0)
cat("fitted ", fitted$family, ":", sprintf("%8.3f", fitted$parameters), "\n", sep = "")

# Khoudraji-Clayton = code 20200: asymmetric, C(u,v) = u^(1-a) v^(1-b) Clayton(u^a, v^b)
truth <- bicop_dist(20200, 0, c(2.0, 0.6, 0.95))
u <- cbind(u1, hbicop(cbind(u1, w), 1, truth, inverse = TRUE))
fitted <- bicop(u, family_set = 20200, rotation = 0)
cat("fitted Khoudraji-Clayton (theta, a, b):", sprintf("%8.3f", fitted$parameters), "\n", sep = "")

# selection among every parametric family + every generic family
chosen <- bicop(u, family_set = c(0:23, gen$code), cores = 4)
cat("selected among all families: ", chosen$family, "  rotation ", chosen$rotation, "\n", sep = "")
