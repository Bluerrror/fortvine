# R part of bench/run_bench.py: times the fortvine R package and rvinecopulib on the data files
# run_bench.py writes (pair_<family>.csv, vine_<d>.csv) and writes the results to <dir>/r.json.
#   Rscript bench/bench_r.R <dir> <n_eval> <quick 0/1>
# Both packages export the same function names, so every call is namespace-qualified.
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]; n_eval <- as.integer(args[2]); quick <- args[3] == "1"
spec <- jsonlite::fromJSON(file.path(dir, "spec.json"), simplifyVector = FALSE)
have_rv <- requireNamespace("rvinecopulib", quietly = TRUE)
fv <- asNamespace("fortvine")
if (have_rv) rv <- asNamespace("rvinecopulib")

# mean time of f() in seconds: one warm-up call, then repeats until 0.3 s have passed (at least reps);
# a warm-up over 2 s is the timing itself (as in run_bench.py)
timed <- function(f, reps = 1) {
  t0 <- proc.time()[["elapsed"]]
  r <- f()
  el <- proc.time()[["elapsed"]] - t0
  if (el > 2) return(list(t = el, r = r))
  k <- 0; t0 <- proc.time()[["elapsed"]]
  repeat {
    r <- f(); k <- k + 1
    el <- proc.time()[["elapsed"]] - t0
    if (k >= reps && el > 0.3) break
  }
  list(t = el / k, r = r)
}
rvname <- function(f) if (f == "student") "t" else f
res <- list(eval = list(), fit = list(), select = list(), vine = list(),
            versions = list(R = R.version.string, fortvine = as.character(packageVersion("fortvine")),
                            rvinecopulib = if (have_rv) as.character(packageVersion("rvinecopulib")) else NA))

set.seed(1)
u <- matrix(runif(2 * n_eval, 1e-6, 1 - 1e-6), n_eval, 2)
for (f in spec$fams) {
  p <- unlist(f$par)
  ns <- function(g) 1e9 * timed(g)$t / n_eval
  row <- list(family = f$name)
  row$fortvine_r <- list(logpdf = ns(function() fv$dbicop(u, f$name, 0, p, log = TRUE)),
                         hfunc1 = ns(function() fv$hbicop(u, 1, f$name, 0, p)),
                         hinv1 = ns(function() fv$hbicop(u, 1, f$name, 0, p, inverse = TRUE)))
  if (have_rv)
    row$rvinecopulib <- list(logpdf = ns(function() rv$dbicop(u, rvname(f$name), 0, p)),
                             hfunc1 = ns(function() rv$hbicop(u, 1, rvname(f$name), 0, p)),
                             hinv1 = ns(function() rv$hbicop(u, 1, rvname(f$name), 0, p, inverse = TRUE)))
  res$eval[[length(res$eval) + 1]] <- row
  message("eval ", f$name)
}

common <- unlist(spec$common)
for (f in spec$fams) {
  x <- as.matrix(read.csv(file.path(dir, paste0("pair_", f$name, ".csv")), header = FALSE))
  reps <- if (quick) 3 else 10
  a <- timed(function() fv$bicop(x, family_set = f$name, rotation = 0), reps)
  rowf <- list(family = f$name, fortvine_r = list(ms = 1e3 * a$t, loglik = a$r$loglik))
  a <- timed(function() fv$bicop(x, family_set = common, selcrit = "bic"))
  rows <- list(truth = f$name, fortvine_r = list(ms = 1e3 * a$t, family = a$r$family, loglik = a$r$loglik))
  if (have_rv) {
    b <- timed(function() rv$bicop(x, family_set = rvname(f$name), allow_rotations = FALSE, presel = FALSE), reps)
    rowf$rvinecopulib <- list(ms = 1e3 * b$t, loglik = b$r$loglik)
    b <- timed(function() rv$bicop(x, family_set = vapply(common, rvname, ""), selcrit = "bic"))
    fam <- if (b$r$family == "t") "student" else b$r$family
    rows$rvinecopulib <- list(ms = 1e3 * b$t, family = fam, loglik = b$r$loglik)
  }
  res$fit[[length(res$fit) + 1]] <- rowf
  res$select[[length(res$select) + 1]] <- rows
  message("fit/select ", f$name)
}

for (d in unlist(spec$dims)) {
  x <- as.matrix(read.csv(file.path(dir, paste0("vine_", d, ".csv")), header = FALSE))
  for (thr in c(1L, 8L)) {
    row <- list(d = d, threads = thr)
    a <- timed(function() fv$vinecop(x, family_set = common, selcrit = "bic", cores = thr))
    v <- a$r
    row$fortvine_r <- list(select_s = a$t, loglik = v$loglik, npars = v$npars,
                           pdf_s = timed(function() fv$dvinecop(x, v, cores = thr))$t,
                           simulate_s = timed(function() fv$rvinecop(nrow(x), v, cores = thr))$t)
    if (have_rv) {
      b <- timed(function() rv$vinecop(x, family_set = vapply(common, rvname, ""), selcrit = "bic",
                                       cores = thr))
      w <- b$r
      row$rvinecopulib <- list(select_s = b$t, loglik = w$loglik, npars = w$npars,
                               pdf_s = timed(function() rv$dvinecop(x, w, cores = thr))$t,
                               simulate_s = timed(function() rv$rvinecop(nrow(x), w, cores = thr))$t)
    }
    res$vine[[length(res$vine) + 1]] <- row
    message("vine d=", d, " threads=", thr)
  }
}
jsonlite::write_json(res, file.path(dir, "r.json"), auto_unbox = TRUE, digits = NA)
