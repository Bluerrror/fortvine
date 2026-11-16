# Pair copulas: model objects, evaluation, simulation, fitting and selection. Every number is
# computed by fv_bicop_eval / fv_bicop_fit / fv_bicop_select of the Fortran library.

#' Pair-copula model
#'
#' Builds a pair-copula model from a family, a rotation and its parameters.
#'
#' @param family a fortvine family name (case-insensitive, `"t"` is accepted for `"student"`)
#'   or integer code, generic codes 10000+ included; see [fortvine_families()].
#' @param rotation 0, 90, 180 or 270 (counter-clockwise, as in vinecopulib).
#' @param parameters numeric vector with as many entries as the family has parameters.
#' @return An object of class `bicop_dist`: a list with `family` (name), `code`, `rotation`,
#'   `parameters` and `npars`; fitted models also hold `loglik`, `nobs` and `selcrit`.
#' @examples
#' cop <- bicop_dist("clayton", 0, 2)
#' cop
#' @export
bicop_dist <- function(family = "indep", rotation = 0, parameters = numeric()) {
  code <- .fam_code(family)
  if (length(code) != 1L) stop("give one family", call. = FALSE)
  np <- .fam_npars(code)
  if (length(parameters) != np)
    stop(sprintf("family '%s' takes %d parameter(s), %d given", .fam_name(code), np,
                 length(parameters)), call. = FALSE)
  if (length(rotation) != 1L || !rotation %in% c(0, 90, 180, 270))
    stop("rotation must be 0, 90, 180 or 270", call. = FALSE)
  structure(list(family = .fam_name(code), code = code, rotation = as.integer(rotation),
                 parameters = as.numeric(parameters), npars = np), class = "bicop_dist")
}

.as_bicop <- function(family, rotation, parameters)
  if (inherits(family, "bicop_dist")) family else bicop_dist(family, rotation, parameters)

.par3 <- function(b) c(b$parameters, 0, 0, 0)[1:3]

# n x 2 matrix (or a length-2 vector) -> 2 x n double matrix of the C interface
.pairs <- function(u) {
  u <- if (is.null(dim(u))) matrix(as.double(u), ncol = 2) else as.matrix(u)
  if (ncol(u) != 2L) stop("u must have two columns", call. = FALSE)
  storage.mode(u) <- "double"
  t(u)
}

.bicop_eval <- function(u, b, what, cores)
  .Call(C_bicop_eval, b$code, b$rotation, .par3(b), .pairs(u), as.integer(what), .cores(cores))

#' Pair-copula density, h-functions and simulation
#'
#' `dbicop()` evaluates the density, `hbicop()` the h-functions and their inverses, and
#' `rbicop()` simulates. The model is given either by `family`, `rotation` and `parameters`
#' or by a [bicop_dist()] object in place of `family`.
#'
#' `hbicop(u, 1)` is \eqn{h_1(u_1, u_2) = P(U_2 \le u_2 \mid U_1 = u_1)} and `hbicop(u, 2)` is
#' \eqn{h_2(u_1, u_2) = P(U_1 \le u_1 \mid U_2 = u_2)}. With `inverse = TRUE`, `hbicop(u, 1)`
#' returns \eqn{v} with \eqn{h_1(u_1, v) = u_2} and `hbicop(u, 2)` returns \eqn{v} with
#' \eqn{h_2(v, u_2) = u_1}.
#'
#' `rbicop()` draws independent uniforms with [stats::runif()], so it follows [set.seed()],
#' and maps them through the inverse h-function of the library:
#' \eqn{(u_1, h_1^{-1}(u_1, w))}.
#'
#' The distribution function (`pbicop` in rvinecopulib) is not available: the fortvine C
#' interface does not export it.
#'
#' @param u an `n x 2` matrix (or a vector of length 2) of points in \eqn{(0, 1)^2}.
#' @param family,rotation,parameters the model, see [bicop_dist()]; or a `bicop_dist` object
#'   as `family`.
#' @param log return the log-density.
#' @param cond_var 1 or 2, the conditioning variable.
#' @param inverse compute the inverse h-function.
#' @param n number of samples.
#' @param cores number of OpenMP threads (used for more than 4096 points).
#' @return `dbicop()` and `hbicop()`: a numeric vector of length `n`; `rbicop()`: an `n x 2`
#'   matrix.
#' @examples
#' cop <- bicop_dist("gumbel", 0, 2)
#' u <- rbind(c(0.2, 0.3), c(0.7, 0.9))
#' dbicop(u, cop)
#' dbicop(u, "gumbel", 0, 2, log = TRUE)
#' w <- hbicop(u, 1, cop)
#' hbicop(cbind(u[, 1], w), 1, cop, inverse = TRUE)   # gives back u[, 2]
#' set.seed(1)
#' head(rbicop(5, cop))
#' @name dbicop
#' @export
dbicop <- function(u, family, rotation = 0, parameters = numeric(), log = FALSE, cores = 1)
  .bicop_eval(u, .as_bicop(family, rotation, parameters), if (log) 0L else 1L, cores)

#' @rdname dbicop
#' @export
hbicop <- function(u, cond_var, family, rotation = 0, parameters = numeric(), inverse = FALSE,
                   cores = 1) {
  if (!cond_var %in% 1:2) stop("cond_var must be 1 or 2", call. = FALSE)
  .bicop_eval(u, .as_bicop(family, rotation, parameters), 1L + cond_var + 2L * inverse, cores)
}

#' @rdname dbicop
#' @export
rbicop <- function(n, family, rotation = 0, parameters = numeric(), cores = 1) {
  b <- .as_bicop(family, rotation, parameters)
  u <- matrix(runif(2 * n), n, 2)
  u[, 2] <- .bicop_eval(u, b, 4L, cores)
  u
}

# ci(6) of the pair routines
.ci6 <- function(selcrit, par_method, presel, allow_rotations) {
  c(.match_code(selcrit, c("loglik", "aic", "bic", "mbic", "mbicv"), "selcrit"),
    .match_code(par_method, c("mle", "itau"), "par_method"),
    as.integer(presel), as.integer(allow_rotations), 0L, 1L)
}

.pseudo_check <- function(u) {
  if (anyNA(u) || any(u <= 0 | u >= 1))
    stop("data must be pseudo-observations in (0, 1) without NA; see pseudo_obs()",
         call. = FALSE)
}

#' Fit or select a pair copula
#'
#' Selects the family, rotation and parameters among `family_set` (as `Bicop::select` of
#' vinecopulib), or, with `rotation` given and a single family, fits that family and rotation
#' only. Both are done by the Fortran library.
#'
#' @param data an `n x 2` matrix of pseudo-observations.
#' @param family_set family names or codes, and the groups `"all"` (codes 0-22, fortvine's
#'   default set), `"parametric"` (0-11, vinecopulib's parametric families) and `"generic"`
#'   (the 85 generic families).
#' @param par_method `"mle"` or `"itau"` (Kendall's tau inversion, one-parameter families).
#' @param selcrit `"loglik"`, `"aic"`, `"bic"`, `"mbic"` or `"mbicv"`.
#' @param psi0 prior probability of a non-independence copula (mBIC).
#' @param presel preselect families from the empirical tail asymmetry.
#' @param allow_rotations also consider rotated families.
#' @param rotation `NULL` to select; 0, 90, 180 or 270 to fit the single family of
#'   `family_set` with that rotation.
#' @param cores number of OpenMP threads.
#' @return A [bicop_dist()] object with `loglik`, `nobs` and `selcrit`.
#' @examples
#' set.seed(2)
#' u <- rbicop(500, "clayton", 0, 3)
#' fit <- bicop(u, family_set = c("gaussian", "clayton", "gumbel", "frank"))
#' summary(fit)
#' logLik(fit)
#' bicop(u, family_set = "clayton", rotation = 0)$parameters
#' @export
bicop <- function(data, family_set = "all", par_method = "mle", selcrit = "bic", psi0 = 0.9,
                  presel = TRUE, allow_rotations = TRUE, rotation = NULL, cores = 1) {
  u <- .pairs(data)
  .pseudo_check(u)
  ci <- .ci6(selcrit, par_method, presel, allow_rotations)
  fs <- .family_set(family_set)
  if (is.null(rotation)) {
    r <- .Call(C_bicop_select, u, fs, ci, as.double(psi0), .cores(cores))
  } else {
    if (length(fs) != 1L) stop("with rotation given, family_set must hold one family",
                               call. = FALSE)
    bicop_dist(fs, rotation, numeric(.fam_npars(fs)))     # checks the rotation
    r <- .Call(C_bicop_fit, u, fs, as.integer(rotation), ci, .cores(cores))
    r$family <- fs
    r$rotation <- as.integer(rotation)
  }
  b <- bicop_dist(r$family, r$rotation, r$par[seq_len(.fam_npars(r$family))])
  b$loglik <- r$loglik
  b$nobs <- ncol(u)
  b$selcrit <- selcrit
  b
}

.fmt_par <- function(p) if (length(p)) paste(formatC(p, digits = 4, format = "g"),
                                             collapse = ", ") else "none"

#' @export
print.bicop_dist <- function(x, ...) {
  cat("Bivariate copula ('bicop_dist'): family = ", x$family, " (", x$code, "), rotation = ",
      x$rotation, ", parameters = ", .fmt_par(x$parameters), "\n", sep = "")
  invisible(x)
}

#' Summaries and log-likelihoods of fortvine models
#'
#' `summary()` returns a data frame: one row for a pair copula, one row per edge for a vine
#' (tree, edge, conditioned and conditioning variables, family, rotation, parameters,
#' number of parameters and edge log-likelihood). `logLik()` returns the log-likelihood of a
#' fitted model as a `logLik` object, so [stats::AIC()] and [stats::BIC()] work.
#'
#' @param object a `bicop_dist` or `vinecop_dist` object.
#' @param ... unused.
#' @return `summary()`: a data frame; `logLik()`: an object of class `logLik`.
#' @examples
#' set.seed(3)
#' fit <- bicop(rbicop(300, "frank", 0, 4), family_set = "parametric")
#' summary(fit)
#' AIC(fit)
#' @name summary.bicop_dist
#' @export
summary.bicop_dist <- function(object, ...) {
  data.frame(family = object$family, code = object$code, rotation = object$rotation,
             parameters = .fmt_par(object$parameters), npars = object$npars,
             loglik = if (is.null(object$loglik)) NA_real_ else object$loglik,
             stringsAsFactors = FALSE)
}

#' @rdname summary.bicop_dist
#' @export
logLik.bicop_dist <- function(object, ...) {
  if (is.null(object$loglik)) stop("the model was not fitted to data", call. = FALSE)
  structure(object$loglik, df = object$npars, nobs = object$nobs, class = "logLik")
}
