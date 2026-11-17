# Pair copulas: model objects, evaluation, simulation, fitting and selection. Every number is
# computed by fvc_bicop_eval / fvc_bicop_fit of the Fortran library (src/fortvine_capi.f90).

#' Pair-copula model
#'
#' Builds a pair-copula model from a family, a rotation and its parameters.
#'
#' The nonparametric families `"tll"` (code 24) and `"bernstein"` (25) are described by an
#' estimate instead of parameters: get them from [bicop()] (or from a fitted vine), which
#' stores the estimate in `aux`. Their `npars` is the effective degrees of freedom of the fit.
#'
#' @param family a fortvine family name (case-insensitive, `"t"` is accepted for `"student"`)
#'   or integer code, generic codes 10000+ included; see [fortvine_families()].
#' @param rotation 0, 90, 180 or 270 (counter-clockwise, as in vinecopulib). Ignored (set to 0)
#'   for `tll` and `bernstein`.
#' @param parameters numeric vector with as many entries as the family has parameters (none for
#'   `tll` and `bernstein`).
#' @param aux the estimate of a `tll` / `bernstein` copula (a numeric vector of length 933, as
#'   stored by [bicop()] in `$aux`); `NULL` for the parametric families.
#' @return An object of class `bicop_dist`: a list with `family` (name), `code`, `rotation`,
#'   `parameters` and `npars` (and `aux` for `tll` / `bernstein`); fitted models also hold
#'   `loglik`, `nobs` and `selcrit`.
#' @examples
#' cop <- bicop_dist("clayton", 0, 2)
#' cop
#' @export
bicop_dist <- function(family = "indep", rotation = 0, parameters = numeric(), aux = NULL) {
  code <- .fam_code(family)
  if (length(code) != 1L) stop("give one family", call. = FALSE)
  if (length(rotation) != 1L || !rotation %in% c(0, 90, 180, 270))
    stop("rotation must be 0, 90, 180 or 270", call. = FALSE)
  if (.is_np(code)) {
    if (is.null(aux) || !is.numeric(aux) || length(aux) != .naux() || aux[1] != code)
      stop(sprintf("family '%s' is nonparametric: its estimate comes from bicop() ",
                   .fam_name(code)), "(or a fitted vine), not from parameters", call. = FALSE)
    if (length(parameters)) stop("tll and bernstein take no parameters", call. = FALSE)
    return(structure(list(family = .fam_name(code), code = code, rotation = 0L,
                          parameters = numeric(), npars = aux[2], aux = as.double(aux)),
                     class = "bicop_dist"))
  }
  if (!is.null(aux)) stop("aux is only used by tll and bernstein", call. = FALSE)
  np <- .fam_npars(code)
  if (length(parameters) != np)
    stop(sprintf("family '%s' takes %d parameter(s), %d given", .fam_name(code), np,
                 length(parameters)), call. = FALSE)
  structure(list(family = .fam_name(code), code = code, rotation = as.integer(rotation),
                 parameters = as.numeric(parameters), npars = np), class = "bicop_dist")
}

.as_bicop <- function(family, rotation, parameters)
  if (inherits(family, "bicop_dist")) family else bicop_dist(family, rotation, parameters)

.par3 <- function(b) c(b$parameters, 0, 0, 0)[1:3]
.aux_of <- function(b) if (is.null(b$aux)) 0 else b$aux

# n x 2 matrix (or a length-2 vector) -> 2 x n double matrix of the C interface
.pairs <- function(u) {
  u <- if (is.null(dim(u))) matrix(as.double(u), ncol = 2) else as.matrix(u)
  if (ncol(u) != 2L) stop("u must have two columns", call. = FALSE)
  storage.mode(u) <- "double"
  t(u)
}

.bicop_eval <- function(u, b, what, cores)
  .with_threads(cores, .Call(C_bicop_eval, as.integer(b$code), as.integer(b$rotation),
                             .par3(b), .aux_of(b), .pairs(u), as.integer(what)))

#' Pair-copula density, distribution, h-functions and simulation
#'
#' `dbicop()` evaluates the density, `pbicop()` the distribution function, `hbicop()` the
#' h-functions and their inverses, and `rbicop()` simulates. The model is given either by
#' `family`, `rotation` and `parameters` or by a [bicop_dist()] object in place of `family`.
#'
#' `hbicop(u, 1)` is \eqn{h_1(u_1, u_2) = P(U_2 \le u_2 \mid U_1 = u_1)} and `hbicop(u, 2)` is
#' \eqn{h_2(u_1, u_2) = P(U_1 \le u_1 \mid U_2 = u_2)}. With `inverse = TRUE`, `hbicop(u, 1)`
#' returns \eqn{v} with \eqn{h_1(u_1, v) = u_2} and `hbicop(u, 2)` returns \eqn{v} with
#' \eqn{h_2(v, u_2) = u_1}. `pbicop()` is \eqn{C(u_1, u_2)} as computed by the library (closed
#' form, or quadrature of the h-function for the Gaussian and Student copulas).
#'
#' `rbicop()` draws independent uniforms with [stats::runif()], so it follows [set.seed()],
#' and maps them through the inverse h-function of the library:
#' \eqn{(u_1, h_1^{-1}(u_1, w))}.
#'
#' @param u an `n x 2` matrix (or a vector of length 2) of points in \eqn{(0, 1)^2}.
#' @param family,rotation,parameters the model, see [bicop_dist()]; or a `bicop_dist` object
#'   as `family`.
#' @param log return the log-density.
#' @param cond_var 1 or 2, the conditioning variable.
#' @param inverse compute the inverse h-function.
#' @param n number of samples.
#' @param cores number of OpenMP threads (used for more than 4096 points).
#' @return `dbicop()`, `pbicop()` and `hbicop()`: a numeric vector of length `n`; `rbicop()`:
#'   an `n x 2` matrix.
#' @examples
#' cop <- bicop_dist("gumbel", 0, 2)
#' u <- rbind(c(0.2, 0.3), c(0.7, 0.9))
#' dbicop(u, cop)
#' dbicop(u, "gumbel", 0, 2, log = TRUE)
#' pbicop(u, cop)
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
pbicop <- function(u, family, rotation = 0, parameters = numeric(), cores = 1)
  .bicop_eval(u, .as_bicop(family, rotation, parameters), 6L, cores)

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
#'   default set), `"parametric"` (0-11, vinecopulib's parametric families),
#'   `"nonparametric"` (`"tll"` and `"bernstein"`, codes 24 and 25) and `"generic"` (the 85
#'   generic families).
#' @param par_method `"mle"` or `"itau"` (Kendall's tau inversion, one-parameter families).
#' @param selcrit `"loglik"`, `"aic"`, `"bic"`, `"mbic"` or `"mbicv"`.
#' @param psi0 prior probability of a non-independence copula (mBIC).
#' @param presel preselect families from the empirical tail asymmetry.
#' @param allow_rotations also consider rotated families.
#' @param rotation `NULL` to select; 0, 90, 180 or 270 to fit the single family of
#'   `family_set` with that rotation.
#' @param cores number of OpenMP threads.
#' @param weights optional vector of `n` non-negative observation weights.
#' @param fast_mle with `par_method = "mle"`, fit the one-parameter Gaussian, Clayton, Gumbel,
#'   Frank and Joe copulas by Kendall's tau inversion only.
#' @param tau_mode how Kendall's tau is computed for start values and preselection:
#'   `"exact"`, `"approx"` (normal-score approximation) or `"auto"` (exact below 800
#'   observations).
#' @return A [bicop_dist()] object with `loglik`, `nobs` and `selcrit`. For `tll` and
#'   `bernstein`, `aux` holds the estimate and `npars` the effective degrees of freedom, which
#'   `logLik()`, `AIC()` and `BIC()` use.
#' @examples
#' set.seed(2)
#' u <- rbicop(500, "clayton", 0, 3)
#' fit <- bicop(u, family_set = c("gaussian", "clayton", "gumbel", "frank"))
#' summary(fit)
#' logLik(fit)
#' bicop(u, family_set = "clayton", rotation = 0)$parameters
#' np <- bicop(u, family_set = "tll")                # nonparametric
#' c(edf = np$npars, loglik = np$loglik)
#' w <- rep(1:2, length.out = 500)                   # weighted fit
#' bicop(u, family_set = "clayton", rotation = 0, weights = w)$parameters
#' @export
bicop <- function(data, family_set = "all", par_method = "mle", selcrit = "bic", psi0 = 0.9,
                  presel = TRUE, allow_rotations = TRUE, rotation = NULL, cores = 1,
                  weights = numeric(), fast_mle = FALSE, tau_mode = "exact") {
  u <- .pairs(data)
  .pseudo_check(u)
  ctl <- .controls(selcrit = selcrit, par_method = par_method, tau_mode = tau_mode,
                   cores = cores, presel = presel, allow_rotations = allow_rotations,
                   fast_mle = fast_mle, psi0 = psi0)
  w <- .weights(weights, ncol(u))
  fs <- .family_set(family_set)
  naux <- if (any(.is_np(fs))) .naux() else 1L
  if (is.null(rotation)) {
    r <- .Call(C_bicop_fit, u, w, 1L, fs, ctl$ci, ctl$cr, 0L, 0L, naux)
  } else {
    if (length(fs) != 1L) stop("with rotation given, family_set must hold one family",
                               call. = FALSE)
    if (length(rotation) != 1L || !rotation %in% c(0, 90, 180, 270))
      stop("rotation must be 0, 90, 180 or 270", call. = FALSE)
    r <- .Call(C_bicop_fit, u, w, 0L, fs, ctl$ci, ctl$cr, fs, as.integer(rotation), naux)
  }
  b <- if (.is_np(r$family)) bicop_dist(r$family, 0, aux = r$aux) else
    bicop_dist(r$family, r$rotation, r$par[seq_len(.fam_npars(r$family))])
  b$npars <- r$npars
  b$loglik <- r$loglik
  b$nobs <- ncol(u)
  b$selcrit <- selcrit
  b
}

#' Parameter bounds and transposed generic families
#'
#' `bicop_par_bounds()` returns the parameter box of a family as used by the fitting routines
#' of the library (`fv_par_bounds`). `gen_flip()` returns the generic family and parameters of
#' the transposed copula \eqn{C(u_2, u_1)} (`gen_flip` of the library); other families are
#' returned unchanged. The rotation is relabelled (90 and 270 swap) so that
#' `dbicop(u[, 2:1], gen_flip(cop))` equals `dbicop(u, cop)`.
#'
#' @param family a family name or code, or a [bicop_dist()] object.
#' @param rotation,parameters the model when `family` is not a `bicop_dist` object.
#' @return `bicop_par_bounds()`: a 2-row matrix (`lower`, `upper`) with one column per
#'   parameter; `gen_flip()`: a [bicop_dist()] object.
#' @examples
#' bicop_par_bounds("bb1")
#' fams <- fortvine_families()
#' kh <- fams[fams$code >= 10000 & grepl("^Khoudraji", fams$name), ][1, ]
#' kh$name
#' @export
bicop_par_bounds <- function(family) {
  code <- if (inherits(family, "bicop_dist")) family$code else .fam_code(family)
  if (length(code) != 1L) stop("give one family", call. = FALSE)
  np <- if (.is_np(code)) 0L else .fam_npars(code)
  r <- .Call(C_par_bounds, as.integer(code))
  rbind(lower = r$lower[seq_len(np)], upper = r$upper[seq_len(np)])
}

#' @rdname bicop_par_bounds
#' @export
gen_flip <- function(family, rotation = 0, parameters = numeric()) {
  b <- .as_bicop(family, rotation, parameters)
  if (b$code < 10000) return(b)
  r <- .Call(C_gen_flip, as.integer(b$code), .par3(b))
  rot <- c(0L, 270L, 180L, 90L)[match(b$rotation, c(0L, 90L, 180L, 270L))]
  bicop_dist(r$family, rot, r$par[seq_len(.fam_npars(r$family))])
}

.fmt_par <- function(p) if (length(p)) paste(formatC(p, digits = 4, format = "g"),
                                             collapse = ", ") else "none"
.fmt_np <- function(b) if (.is_np(b$code)) sprintf("none (nonparametric, edf %.3g)", b$npars) else
  .fmt_par(b$parameters)

#' @export
print.bicop_dist <- function(x, ...) {
  cat("Bivariate copula ('bicop_dist'): family = ", x$family, " (", x$code, "), rotation = ",
      x$rotation, ", parameters = ", .fmt_np(x), "\n", sep = "")
  invisible(x)
}

#' Summaries and log-likelihoods of fortvine models
#'
#' `summary()` returns a data frame: one row for a pair copula, one row per edge for a vine
#' (tree, edge, conditioned and conditioning variables, family, rotation, parameters,
#' number of parameters and edge log-likelihood). `logLik()` returns the log-likelihood of a
#' fitted model as a `logLik` object, so [stats::AIC()] and [stats::BIC()] work; its `df` is
#' the number of parameters, with the effective degrees of freedom for `tll` and `bernstein`
#' and the anchor parameters of relaxed edges.
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
