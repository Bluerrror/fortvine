# Vine copulas: structure, model objects, selection, density, Rosenblatt transforms and
# simulation, through fv_vine_select / fv_vine_pdf / fv_vine_(inverse_)rosenblatt.
#
# A vine is stored in the arrays of the C interface (natural-order convention):
#   order (d), st, fam, rot (d x d, [edge, tree]), par (3 x d x d, [k, edge, tree]).
# Edge e of tree t joins order[e] and order[st[e, t]] given order[st[e, 1:(t - 1)]].

.trunc <- function(trunc_lvl, d) as.integer(if (is.infinite(trunc_lvl)) d - 1 else
                                             max(0, min(d - 1, trunc_lvl)))

#' D-vine structure
#'
#' The structure of a D-vine (path) in the given order, as built by `set_dvine` of the
#' Fortran library: tree 1 links `order[1]-order[2]`, `order[2]-order[3]`, ...
#'
#' @param order a permutation of `1:d`.
#' @param trunc_lvl number of trees (`Inf` for \eqn{d - 1}).
#' @return An object of class `fv_structure` with `order`, `trunc_lvl` and `struct_array`
#'   (`d x d`, indexed `[edge, tree]`).
#' @examples
#' dvine_structure(1:4)
#' @export
dvine_structure <- function(order, trunc_lvl = Inf) {
  order <- as.integer(order)
  d <- length(order)
  if (d < 2L || !identical(sort(order), seq_len(d)))
    stop("order must be a permutation of 1:d with d >= 2", call. = FALSE)
  tr <- .trunc(trunc_lvl, d)
  structure(list(order = order, trunc_lvl = tr, struct_array = .Call(C_set_dvine, order, tr)),
            class = "fv_structure")
}

#' @export
print.fv_structure <- function(x, ...) {
  cat("Vine structure ('fv_structure'): d = ", length(x$order), ", trunc_lvl = ", x$trunc_lvl,
      ", order = ", paste(x$order, collapse = " "), "\n", sep = "")
  invisible(x)
}

# build the R object from the C-interface arrays
.new_vine <- function(order, trunc, st, fam, rot, par, edge_ll = NULL) {
  d <- length(order)
  storage.mode(st) <- storage.mode(fam) <- storage.mode(rot) <- "integer"
  storage.mode(par) <- "double"
  pcs <- lapply(seq_len(trunc), function(t) lapply(seq_len(d - t), function(e) {
    b <- bicop_dist(fam[e, t], rot[e, t], par[seq_len(.fam_npars(fam[e, t])), e, t])
    if (!is.null(edge_ll)) b$loglik <- edge_ll[e, t]
    b
  }))
  npars <- sum(vapply(unlist(pcs, recursive = FALSE), function(b) b$npars, 0))
  str <- structure(list(order = as.integer(order), trunc_lvl = as.integer(trunc),
                        struct_array = st), class = "fv_structure")
  structure(list(order = as.integer(order), trunc = as.integer(trunc), struct_array = st,
                 fam = fam, rot = rot, par = par, pair_copulas = pcs, npars = npars, d = d,
                 structure = str), class = "vinecop_dist")
}

#' Vine copula model
#'
#' Builds a vine copula from a structure and its pair copulas.
#'
#' @param pair_copulas a list of trees; tree `t` is a list of `d - t` [bicop_dist()] objects,
#'   one per edge. Trees beyond the length of the list are independence.
#' @param structure a structure from [dvine_structure()] or the `structure` of a fitted vine.
#' @return An object of class `vinecop_dist`.
#' @examples
#' pcs <- list(list(bicop_dist("gumbel", 0, 2), bicop_dist("gumbel", 0, 2)),
#'             list(bicop_dist("frank", 0, 3)))
#' vc <- vinecop_dist(pcs, dvine_structure(1:3))
#' summary(vc)
#' @export
vinecop_dist <- function(pair_copulas, structure) {
  if (!inherits(structure, "fv_structure")) stop("structure must come from dvine_structure()",
                                                 call. = FALSE)
  d <- length(structure$order)
  trunc <- min(structure$trunc_lvl, length(pair_copulas))
  fam <- rot <- matrix(0L, d, d)
  par <- array(0, c(3, d, d))
  for (t in seq_len(trunc)) {
    if (length(pair_copulas[[t]]) != d - t)
      stop(sprintf("tree %d needs %d pair copulas", t, d - t), call. = FALSE)
    for (e in seq_len(d - t)) {
      b <- pair_copulas[[t]][[e]]
      if (!inherits(b, "bicop_dist")) stop("pair copulas must be bicop_dist objects",
                                           call. = FALSE)
      fam[e, t] <- b$code
      rot[e, t] <- b$rotation
      par[, e, t] <- .par3(b)
    }
  }
  .new_vine(structure$order, trunc, structure$struct_array, fam, rot, par)
}

# n x d matrix -> d x n double matrix
.cols <- function(u, d = NULL) {
  u <- as.matrix(u)
  if (!is.null(d) && ncol(u) != d) stop("data must have ", d, " columns", call. = FALSE)
  storage.mode(u) <- "double"
  t(u)
}

#' Select a vine copula
#'
#' Dissmann's algorithm as implemented by `fv_vinecop%select`: the structure (maximum
#' spanning trees of `tree_crit`), the pair-copula families among `family_set` and their
#' parameters, all chosen by the Fortran library. With `structure` given only the families
#' and parameters are selected.
#'
#' @param data an `n x d` matrix of pseudo-observations.
#' @param family_set,par_method,selcrit,psi0,presel,allow_rotations see [bicop()]; `selcrit`
#'   `"mbicv"` is the vine criterion of Nagler et al. (2019).
#' @param structure `NA` to select the structure, or an [dvine_structure()] object.
#' @param trunc_lvl number of trees (`Inf` for all).
#' @param tree_crit `"tau"`, `"rho"`, `"joe"`, `"hoeffd"`, `"xi"`, `"mi"` or `"ferreira"`.
#' @param tree_algorithm `"mst_prim"` or `"mst_kruskal"`.
#' @param threshold edges with \eqn{|\tau|} below it are set to independence.
#' @param cores number of OpenMP threads.
#' @return A [vinecop_dist()] object with `loglik`, `nobs` and `names` (the column names of
#'   `data`, used by `print()` and `summary()`); each pair copula holds its edge log-likelihood.
#' @examples
#' set.seed(4)
#' z <- matrix(rnorm(1500), 500, 3)
#' z[, 2] <- z[, 2] + z[, 1]
#' z[, 3] <- z[, 3] + z[, 2]
#' fit <- vinecop(pseudo_obs(z), family_set = "parametric", cores = 2)
#' fit
#' summary(fit)
#' @export
vinecop <- function(data, family_set = "all", structure = NA, par_method = "mle",
                    selcrit = "bic", psi0 = 0.9, presel = TRUE, allow_rotations = TRUE,
                    trunc_lvl = Inf, tree_crit = "tau", tree_algorithm = "mst_prim",
                    threshold = 0, cores = 1) {
  u <- .cols(data)
  .pseudo_check(u)
  d <- nrow(u)
  if (d < 2L) stop("data must have at least two columns", call. = FALSE)
  fixed <- inherits(structure, "fv_structure")
  if (fixed) {
    if (length(structure$order) != d) stop("structure dimension differs from the data",
                                           call. = FALSE)
    trunc_lvl <- min(trunc_lvl, structure$trunc_lvl)
    order <- structure$order
    st_in <- structure$struct_array
  } else {
    order <- integer(d)
    st_in <- matrix(0L, d, d)
  }
  tr <- .trunc(trunc_lvl, d)
  ci <- c(.ci6(selcrit, par_method, presel, allow_rotations),
          .match_code(tree_crit, c("tau", "rho", "joe", "hoeffd", "xi", "mi", "ferreira"),
                      "tree_crit"),
          .match_code(tree_algorithm, c("mst_prim", "mst_kruskal"), "tree_algorithm"),
          tr, as.integer(fixed))
  r <- .Call(C_vine_select, u, .family_set(family_set), ci, as.double(c(threshold, psi0)),
             as.integer(order), st_in, .cores(cores))
  if (r$info != 0L) stop("fv_vine_select failed (info = ", r$info, ")", call. = FALSE)
  vc <- .new_vine(r$order, tr, r$st, r$fam, r$rot, r$par, r$edge_ll)
  vc$loglik <- sum(r$edge_ll[, seq_len(tr)])
  vc$nobs <- ncol(u)
  vc$selcrit <- selcrit
  vc$names <- colnames(data)
  vc
}

.vine_eval <- function(u, vc, what, cores) {
  if (!inherits(vc, "vinecop_dist")) stop("model must be a vinecop_dist object", call. = FALSE)
  .Call(C_vine_eval, .cols(u, vc$d), vc$order, vc$trunc, vc$struct_array, vc$fam, vc$rot,
        vc$par, as.integer(what), .cores(cores))
}

#' Vine copula density, simulation and Rosenblatt transforms
#'
#' `dvinecop()` evaluates the density, `rosenblatt()` the Rosenblatt transform (which maps a
#' sample of the model to independent uniforms) and `inverse_rosenblatt()` its inverse.
#' `rvinecop()` simulates by applying the inverse Rosenblatt transform of the library to
#' independent uniforms from [stats::runif()], so it follows [set.seed()].
#'
#' @param u,x an `n x d` matrix in \eqn{(0, 1)^d}.
#' @param vinecop,model a [vinecop_dist()] object.
#' @param n number of samples.
#' @param cores number of OpenMP threads.
#' @return `dvinecop()`: a vector of length `n`; the others: an `n x d` matrix.
#' @examples
#' vc <- vinecop_dist(list(list(bicop_dist("clayton", 0, 2), bicop_dist("gumbel", 0, 1.5))),
#'                    dvine_structure(1:3))
#' set.seed(5)
#' u <- rvinecop(200, vc)
#' head(dvinecop(u, vc))
#' w <- rosenblatt(u, vc)
#' max(abs(inverse_rosenblatt(w, vc) - u))
#' @name dvinecop
#' @export
dvinecop <- function(u, vinecop, cores = 1) .vine_eval(u, vinecop, 0L, cores)

#' @rdname dvinecop
#' @export
rvinecop <- function(n, vinecop, cores = 1)
  inverse_rosenblatt(matrix(runif(n * vinecop$d), n, vinecop$d), vinecop, cores)

#' @rdname dvinecop
#' @export
rosenblatt <- function(x, model, cores = 1) t(.vine_eval(x, model, 1L, cores))

#' @rdname dvinecop
#' @export
inverse_rosenblatt <- function(u, model, cores = 1) t(.vine_eval(u, model, 2L, cores))

.edge_table <- function(x) {
  lab <- if (is.null(x$names)) as.character(seq_len(x$d)) else x$names
  rows <- list()
  for (t in seq_len(x$trunc)) for (e in seq_len(x$d - t)) {
    b <- x$pair_copulas[[t]][[e]]
    st <- x$struct_array
    rows[[length(rows) + 1L]] <- data.frame(
      tree = t, edge = e,
      conditioned = paste(lab[x$order[e]], lab[x$order[st[e, t]]], sep = ", "),
      conditioning = paste(lab[x$order[st[e, seq_len(t - 1)]]], collapse = ", "),
      family = b$family, rotation = b$rotation, parameters = .fmt_par(b$parameters),
      npars = b$npars, loglik = if (is.null(b$loglik)) NA_real_ else b$loglik,
      stringsAsFactors = FALSE)
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

#' @export
print.vinecop_dist <- function(x, ...) {
  cat("Vine copula ('vinecop_dist'): d = ", x$d, ", trunc_lvl = ", x$trunc, ", order = ",
      paste(x$order, collapse = " "), "\n", sep = "")
  cat("  npars = ", x$npars, if (!is.null(x$loglik))
        sprintf(", loglik = %.2f (n = %d)", x$loglik, x$nobs), "\n", sep = "")
  tab <- .edge_table(x)
  if (nrow(tab)) print(tab[, c("tree", "edge", "conditioned", "conditioning", "family",
                               "rotation", "parameters")], row.names = FALSE)
  invisible(x)
}

#' @rdname summary.bicop_dist
#' @export
summary.vinecop_dist <- function(object, ...) .edge_table(object)

#' @rdname summary.bicop_dist
#' @export
logLik.vinecop_dist <- function(object, ...) {
  if (is.null(object$loglik)) stop("the model was not fitted to data", call. = FALSE)
  structure(object$loglik, df = object$npars, nobs = object$nobs, class = "logLik")
}

#' Pseudo-observations
#'
#' Column-wise ranks divided by \eqn{n + 1}, ties sharing their mean rank, computed by
#' `fv_pseudo_obs` of the Fortran library.
#'
#' @param x a numeric matrix or data frame (`n x d`), or a vector.
#' @param cores number of OpenMP threads.
#' @return A matrix of the same shape as `x` with values in \eqn{(0, 1)}.
#' @examples
#' pseudo_obs(c(3.1, 0.2, 5.0, 0.2))
#' @export
pseudo_obs <- function(x, cores = 1) {
  vec <- is.null(dim(x))
  u <- t(.Call(C_pseudo_obs, .cols(if (vec) matrix(x) else x), .cores(cores)))
  if (vec) drop(u) else {
    dimnames(u) <- dimnames(as.matrix(x))
    u
  }
}
