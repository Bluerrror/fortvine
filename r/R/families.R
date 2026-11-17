# Family table and the conversion of family names / codes / sets to the integer codes of the
# C interface. The table itself comes from the Fortran library (fvc_families, fvc_fam_name).

.fv <- new.env(parent = emptyenv())

#' List the pair-copula families of fortvine
#'
#' Returns the 111 distinct families of the library: the 26 hand-coded families (codes 0-25)
#' and the 85 generic families (codes 10000 and above). The names and parameter counts are
#' read from the Fortran library.
#'
#' @return A data frame with columns `code`, `name` and `npars` (`NA` for the nonparametric
#'   families `tll` (24) and `bernstein` (25), whose number of parameters is the effective
#'   degrees of freedom of each fit).
#' @examples
#' fams <- fortvine_families()
#' nrow(fams)
#' head(fams[fams$code >= 10000, ])
#' @export
fortvine_families <- function() {
  if (is.null(.fv$families)) {
    r <- .Call(C_families)
    .fv$families <- data.frame(code = r$code, name = r$name, npars = r$npars,
                               stringsAsFactors = FALSE)
    .fv$naux <- r$naux
    .fv$nk <- r$nk
  }
  .fv$families
}

# family names or codes -> integer codes (vector)
.fam_code <- function(family) {
  tab <- fortvine_families()
  if (is.numeric(family)) {
    i <- match(as.integer(family), tab$code)
  } else {
    f <- tolower(as.character(family))
    f[f == "t"] <- "student"
    num <- grepl("^[0-9]+$", f)
    i <- match(f, tolower(tab$name))
    i[num] <- match(as.integer(f[num]), tab$code)
  }
  if (anyNA(i))
    stop("unknown fortvine family: ", paste(family[is.na(i)], collapse = ", "),
         " (see fortvine_families())", call. = FALSE)
  tab$code[i]
}

# length of the estimate of a tll / bernstein pair copula, anchors of a relaxed edge
.naux <- function() {
  fortvine_families()
  .fv$naux
}
.nk <- function() {
  fortvine_families()
  .fv$nk
}
.is_np <- function(code) code %in% 24:25

# family_set: names, codes and the groups "all" (0-22, fortvine's default set), "parametric"
# (0-11, vinecopulib's parametric set), "nonparametric" (tll and bernstein) and "generic" (the
# 85 generic families)
.family_set <- function(family_set) {
  tab <- fortvine_families()
  groups <- list(all = 0:22, parametric = 0:11, nonparametric = 24:25,
                 generic = tab$code[tab$code >= 10000])
  out <- integer()
  for (f in as.list(family_set)) {
    g <- if (is.character(f)) groups[[tolower(f)]] else NULL
    out <- c(out, if (is.null(g)) .fam_code(f) else g)
  }
  unique(as.integer(out))
}

.fam_name <- function(code) {
  tab <- fortvine_families()
  tab$name[match(code, tab$code)]
}

.fam_npars <- function(code) {
  tab <- fortvine_families()
  tab$npars[match(code, tab$code)]
}

.match_code <- function(x, choices, what) {
  i <- match(tolower(x), choices)
  if (length(i) != 1L || is.na(i))
    stop(what, " must be one of ", paste(dQuote(choices, FALSE), collapse = ", "), call. = FALSE)
  i - 1L
}

.cores <- function(cores) as.integer(max(1, min(cores, 30)))

# the routines without a thread argument (pair-copula evaluation) read FORTVINE_THREADS;
# set it for the call so that they also use `cores` threads (and not half the CPUs)
.with_threads <- function(cores, expr) {
  old <- Sys.getenv("FORTVINE_THREADS", unset = NA)
  Sys.setenv(FORTVINE_THREADS = .cores(cores))
  on.exit(if (is.na(old)) Sys.unsetenv("FORTVINE_THREADS") else
    Sys.setenv(FORTVINE_THREADS = old))
  expr
}

# ci(12), cr(2) of the C interface (fv_controls); ci[10] is the number of threads
.controls <- function(selcrit = "bic", par_method = "mle", tree_crit = "tau",
                      tree_algorithm = "mst_prim", tau_mode = "exact", trunc = -1L, cores = 1,
                      presel = TRUE, allow_rotations = TRUE, fast_mle = FALSE,
                      select_trunc = FALSE, select_threshold = FALSE, threshold = 0,
                      psi0 = 0.9) {
  .flag <- function(x, what) {
    if (!is.logical(x) || length(x) != 1L || is.na(x))
      stop(what, " must be TRUE or FALSE", call. = FALSE)
    as.integer(x)
  }
  if (!is.numeric(psi0) || length(psi0) != 1L || !(psi0 > 0 && psi0 < 1))
    stop("psi0 must be in (0, 1)", call. = FALSE)
  if (!is.numeric(threshold) || length(threshold) != 1L || !(threshold >= 0 && threshold <= 1))
    stop("threshold must be in [0, 1] (or NA to select it)", call. = FALSE)
  # order of fortvine_capi: selcrit, method, preselect, allow_rotations, fast_mle, tau_mode,
  # tree_crit, tree_algo, trunc_lvl, nthreads, select_trunc, select_threshold
  list(ci = c(.match_code(selcrit, c("loglik", "aic", "bic", "mbic", "mbicv"), "selcrit"),
              .match_code(par_method, c("mle", "itau"), "par_method"),
              .flag(presel, "presel"), .flag(allow_rotations, "allow_rotations"),
              .flag(fast_mle, "fast_mle"),
              .match_code(tau_mode, c("auto", "exact", "approx"), "tau_mode"),
              .match_code(tree_crit, c("tau", "rho", "joe", "hoeffd", "xi", "mi", "ferreira"),
                          "tree_crit"),
              .match_code(tree_algorithm, c("mst_prim", "mst_kruskal"), "tree_algorithm"),
              as.integer(trunc), .cores(cores), as.integer(select_trunc),
              as.integer(select_threshold)),
       cr = as.double(c(threshold, psi0)))
}

.weights <- function(weights, n) {
  if (is.null(weights) || !length(weights)) return(double())
  w <- as.double(weights)
  if (length(w) != n || anyNA(w) || any(w < 0))
    stop("weights must be ", n, " non-negative numbers (one per observation)", call. = FALSE)
  w
}
