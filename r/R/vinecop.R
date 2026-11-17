# Vine copulas: structure, model objects, selection, density, Rosenblatt transforms,
# simulation, non-simplified edges and the density on a grid, through fvc_vine_select,
# fvc_vine_eval, fvc_vine_pdf_grid, fvc_vine_relax, fvc_vine_structure and fvc_vine_matrix.
#
# A vine is stored in the arrays of the C interface (natural-order convention):
#   order (d), st, fam, rot (d x d, [edge, tree]), par (3 x d x d, [k, edge, tree]),
#   aux (933 x d x d, the estimates of tll / bernstein edges; NULL without them),
#   var_types (d, "c" or "d"), relaxed (NULL or the non-simplified edges of relax:
#   erel (2 x nrel: edge, tree), wa (d x nrel), wk (nk x nrel), pk (3 x nk x nrel)).
# Edge e of tree t joins order[e] and order[st[e, t]] given order[st[e, 1:(t - 1)]].

.trunc <- function(trunc_lvl, d) {
  if (!is.numeric(trunc_lvl) || length(trunc_lvl) != 1L || is.na(trunc_lvl))
    stop("trunc_lvl must be a number (Inf for all trees)", call. = FALSE)
  as.integer(if (is.infinite(trunc_lvl)) d - 1 else max(0, min(d - 1, trunc_lvl)))
}

.structure <- function(kind, order, trunc_lvl) {
  order <- as.integer(order)
  d <- length(order)
  if (d < 2L || !identical(sort(order), seq_len(d)))
    stop("order must be a permutation of 1:d with d >= 2", call. = FALSE)
  r <- .Call(C_vine_structure, as.integer(kind), order, .trunc(trunc_lvl, d))
  structure(list(order = order, trunc_lvl = r$trunc, struct_array = r$st),
            class = "fv_structure")
}

#' D-vine and C-vine structures, R-vine matrix
#'
#' `dvine_structure()` builds a D-vine (path) in the given order, as `set_dvine` of the
#' Fortran library: tree 1 links `order[1]-order[2]`, `order[2]-order[3]`, ...
#' `cvine_structure()` builds a C-vine (a star in every tree), as `set_cvine`: the root of tree
#' `t` is `order[d - t + 1]`, so the **last** element of `order` is the root of the first tree.
#' `as_rvine_matrix()` returns the R-vine matrix of vinecopulib in the original variable labels
#' (`fv_vinecop%matrix`): column `e` holds in row `t` the partner of `order[e]` in tree `t`, and
#' the anti-diagonal holds `order`.
#'
#' @param order a permutation of `1:d`.
#' @param trunc_lvl number of trees (`Inf` for \eqn{d - 1}).
#' @param x an `fv_structure` or a [vinecop_dist()] object.
#' @return `dvine_structure()`, `cvine_structure()`: an object of class `fv_structure` with
#'   `order`, `trunc_lvl` and `struct_array` (`d x d`, indexed `[edge, tree]`).
#'   `as_rvine_matrix()`: a `d x d` integer matrix (zero above the truncation level).
#' @examples
#' dvine_structure(1:4)
#' cv <- cvine_structure(4:1)       # variable 1 is the root of tree 1
#' as_rvine_matrix(cv)
#' @export
dvine_structure <- function(order, trunc_lvl = Inf) .structure(0L, order, trunc_lvl)

#' @rdname dvine_structure
#' @export
cvine_structure <- function(order, trunc_lvl = Inf) .structure(1L, order, trunc_lvl)

#' @rdname dvine_structure
#' @export
as_rvine_matrix <- function(x) {
  if (inherits(x, "vinecop_dist")) x <- x$structure
  if (!inherits(x, "fv_structure"))
    stop("x must be an fv_structure or a vinecop_dist object", call. = FALSE)
  .Call(C_vine_matrix, as.integer(x$order), as.integer(x$trunc_lvl), x$struct_array)
}

#' @export
print.fv_structure <- function(x, ...) {
  cat("Vine structure ('fv_structure'): d = ", length(x$order), ", trunc_lvl = ", x$trunc_lvl,
      ", order = ", paste(x$order, collapse = " "), "\n", sep = "")
  invisible(x)
}

.var_types <- function(var_types, d = length(var_types)) {
  vt <- as.character(var_types)
  if (length(vt) != d || anyNA(vt) || !all(vt %in% c("c", "d")))
    stop("var_types must hold ", d, " entries \"c\" (continuous) or \"d\" (discrete)",
         call. = FALSE)
  vt
}

# build the R object from the C-interface arrays
.new_vine <- function(order, trunc, st, fam, rot, par, edge_ll = NULL, aux = NULL,
                      var_types = NULL, relaxed = NULL, npars = NULL) {
  d <- length(order)
  storage.mode(st) <- storage.mode(fam) <- storage.mode(rot) <- "integer"
  storage.mode(par) <- "double"
  if (!is.null(aux) && dim(aux)[1] < 2L) aux <- NULL
  if (is.null(var_types)) var_types <- rep("c", d)
  pcs <- lapply(seq_len(trunc), function(t) lapply(seq_len(d - t), function(e) {
    b <- if (.is_np(fam[e, t])) bicop_dist(fam[e, t], 0, aux = aux[, e, t]) else
      bicop_dist(fam[e, t], rot[e, t], par[seq_len(.fam_npars(fam[e, t])), e, t])
    if (!is.null(edge_ll)) b$loglik <- edge_ll[e, t]
    b
  }))
  if (is.null(npars)) {
    npars <- sum(vapply(unlist(pcs, recursive = FALSE), function(b) as.double(b$npars), 0))
    for (k in seq_len(if (is.null(relaxed)) 0L else ncol(relaxed$erel)))
      npars <- npars + (.nk() - 1) * pcs[[relaxed$erel[2, k]]][[relaxed$erel[1, k]]]$npars
  }
  str <- structure(list(order = as.integer(order), trunc_lvl = as.integer(trunc),
                        struct_array = st), class = "fv_structure")
  structure(list(order = as.integer(order), trunc = as.integer(trunc), struct_array = st,
                 fam = fam, rot = rot, par = par, aux = aux, var_types = var_types,
                 relaxed = relaxed, edge_ll = edge_ll, pair_copulas = pcs, npars = npars,
                 d = d, structure = str), class = "vinecop_dist")
}

#' Vine copula model
#'
#' Builds a vine copula from a structure and its pair copulas.
#'
#' @param pair_copulas a list of trees; tree `t` is a list of `d - t` [bicop_dist()] objects,
#'   one per edge. Trees beyond the length of the list are independence.
#' @param structure a structure from [dvine_structure()] or [cvine_structure()], or the
#'   `structure` of a fitted vine.
#' @param var_types `"c"` (continuous) or `"d"` (discrete) for each variable; see [vinecop()]
#'   for the data layout with discrete variables.
#' @return An object of class `vinecop_dist`.
#' @examples
#' pcs <- list(list(bicop_dist("gumbel", 0, 2), bicop_dist("gumbel", 0, 2)),
#'             list(bicop_dist("frank", 0, 3)))
#' vc <- vinecop_dist(pcs, dvine_structure(1:3))
#' summary(vc)
#' @export
vinecop_dist <- function(pair_copulas, structure, var_types = rep("c", length(structure$order))) {
  if (!inherits(structure, "fv_structure"))
    stop("structure must come from dvine_structure() or cvine_structure()", call. = FALSE)
  d <- length(structure$order)
  var_types <- .var_types(var_types, d)
  trunc <- min(structure$trunc_lvl, length(pair_copulas))
  fam <- rot <- matrix(0L, d, d)
  par <- array(0, c(3, d, d))
  aux <- NULL
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
      if (.is_np(b$code)) {
        if (is.null(aux)) aux <- array(0, c(.naux(), d, d))
        aux[, e, t] <- b$aux
      }
    }
  }
  .new_vine(structure$order, trunc, structure$struct_array, fam, rot, par, aux = aux,
            var_types = var_types)
}

# data (n x (d + k), the F(x-) columns of the k discrete variables appended) -> u, um (d x n)
.split <- function(data, var_types, check = FALSE) {
  x <- as.matrix(data)
  storage.mode(x) <- "double"
  d <- length(var_types)
  disc <- var_types == "d"
  k <- sum(disc)
  if (ncol(x) != d + k)
    stop(if (k) sprintf(paste("data must have %d columns: the %d variables, then F(x-) of the",
                              "%d discrete ones"), d + k, d, k) else
           sprintf("data must have %d columns", d), call. = FALSE)
  u <- x[, seq_len(d), drop = FALSE]
  um <- u
  if (k) um[, disc] <- x[, d + seq_len(k)]
  if (check) {
    .pseudo_check(u[, !disc])
    if (k && (anyNA(x) || any(u[, disc] < 0 | u[, disc] > 1 | um[, disc] < 0 |
                              um[, disc] > u[, disc])))
      stop("discrete columns need 0 <= F(x-) <= F(x) <= 1", call. = FALSE)
  }
  list(u = t(u), um = t(um), disc = as.integer(disc))
}

#' Select a vine copula
#'
#' Dissmann's algorithm as implemented by `fv_vinecop%select`: the structure (maximum
#' spanning trees of `tree_crit`), the pair-copula families among `family_set` and their
#' parameters, all chosen by the Fortran library. With `structure` given only the families
#' and parameters are selected.
#'
#' **Discrete variables** follow rvinecopulib: with `k` entries `"d"` in `var_types`, `data`
#' has `d + k` columns, the `d` columns \eqn{F_j(x_j)} followed by the columns
#' \eqn{F_j(x_j^-)} of the `k` discrete variables (in the order of `var_types`). The same
#' layout is used by [dvinecop()] and [rosenblatt()] of the fitted model.
#'
#' @param data an `n x d` matrix of pseudo-observations (`n x (d + k)` with discrete
#'   variables).
#' @param family_set,par_method,selcrit,psi0,presel,allow_rotations,weights,fast_mle,tau_mode
#'   see [bicop()]; `selcrit` `"mbicv"` is the vine criterion of Nagler et al. (2019).
#' @param structure `NA` to select the structure, or an `fv_structure` from
#'   [dvine_structure()] / [cvine_structure()].
#' @param trunc_lvl number of trees (`Inf` for all); `NA` selects the truncation level by
#'   mBICv (vinecopulib's `select_trunc_lvl`).
#' @param tree_crit `"tau"`, `"rho"`, `"joe"`, `"hoeffd"`, `"xi"`, `"mi"` or `"ferreira"`.
#' @param tree_algorithm `"mst_prim"` or `"mst_kruskal"`.
#' @param threshold edges with a tree criterion below it are set to independence; `NA`
#'   selects the threshold by mBICv (vinecopulib's `select_threshold`).
#' @param cores number of OpenMP threads.
#' @param var_types `"c"` (continuous) or `"d"` (discrete) for each variable; `NULL` for all
#'   continuous.
#' @return A [vinecop_dist()] object with `loglik`, `nobs`, `threshold` (the value used),
#'   `names` (the column names of `data`, used by `print()` and `summary()`) and `npars` (with
#'   the effective degrees of freedom of `tll` / `bernstein` edges); each pair copula holds its
#'   edge log-likelihood.
#' @examples
#' set.seed(4)
#' z <- matrix(rnorm(1500), 500, 3)
#' z[, 2] <- z[, 2] + z[, 1]
#' z[, 3] <- z[, 3] + z[, 2]
#' fit <- vinecop(pseudo_obs(z), family_set = "parametric", cores = 2)
#' fit
#' summary(fit)
#' # automatic truncation and threshold
#' vinecop(pseudo_obs(z), family_set = "parametric", trunc_lvl = NA, threshold = NA)$trunc
#' # variable 3 observed as a count: F(x) and F(x-) of the count, appended as column 4
#' k <- cut(z[, 3], c(-Inf, -1, 0, 1, Inf), labels = FALSE)
#' Fx <- ecdf(k)
#' x <- cbind(pseudo_obs(z[, 1:2]), Fx(k) * 500 / 501, Fx(k - 1) * 500 / 501)
#' vd <- vinecop(x, family_set = "parametric", var_types = c("c", "c", "d"))
#' summary(vd)[, c("tree", "family", "parameters")]
#' @export
vinecop <- function(data, family_set = "all", structure = NA, par_method = "mle",
                    selcrit = "bic", psi0 = 0.9, presel = TRUE, allow_rotations = TRUE,
                    trunc_lvl = Inf, tree_crit = "tau", tree_algorithm = "mst_prim",
                    threshold = 0, cores = 1, var_types = NULL, weights = numeric(),
                    fast_mle = FALSE, tau_mode = "exact") {
  if (is.null(var_types)) var_types <- rep("c", NCOL(data))
  var_types <- .var_types(var_types)
  d <- length(var_types)
  if (d < 2L) stop("data must have at least two variables", call. = FALSE)
  x <- .split(data, var_types, check = TRUE)
  n <- ncol(x$u)
  if (inherits(structure, "vinecop_dist")) structure <- structure$structure
  fixed <- inherits(structure, "fv_structure")
  sel_tr <- length(trunc_lvl) == 1L && is.na(trunc_lvl)
  tr <- if (sel_tr) -1L else .trunc(trunc_lvl, d)
  if (fixed) {
    if (length(structure$order) != d) stop("structure dimension differs from the data",
                                           call. = FALSE)
    if (!sel_tr) tr <- min(tr, structure$trunc_lvl)
    order <- structure$order
    st_in <- structure$struct_array
    trunc_in <- structure$trunc_lvl
  } else {
    order <- seq_len(d)
    st_in <- matrix(0L, d, d)
    trunc_in <- 0L
  }
  sel_thr <- length(threshold) == 1L && is.na(threshold)
  ctl <- .controls(selcrit = selcrit, par_method = par_method, tree_crit = tree_crit,
                   tree_algorithm = tree_algorithm, tau_mode = tau_mode, trunc = tr,
                   cores = cores, presel = presel, allow_rotations = allow_rotations,
                   fast_mle = fast_mle, select_trunc = sel_tr, select_threshold = sel_thr,
                   threshold = if (sel_thr) 0 else threshold, psi0 = psi0)
  fs <- .family_set(family_set)
  naux <- if (any(.is_np(fs))) .naux() else 1L
  storage.mode(st_in) <- "integer"
  r <- .Call(C_vine_select, x$u, x$um, x$disc, .weights(weights, n), as.integer(fixed), fs,
             ctl$ci, ctl$cr, as.integer(trunc_in), as.integer(order), st_in, naux)
  if (r$info != 0L) stop("fvc_vine_select failed (info = ", r$info, ")", call. = FALSE)
  vc <- .new_vine(r$order, r$trunc, r$st, r$fam, r$rot, r$par, r$edge_ll, r$aux, var_types,
                  npars = r$npars)
  vc$loglik <- r$loglik
  vc$nobs <- n
  vc$selcrit <- selcrit
  vc$threshold <- r$threshold
  vc$names <- colnames(data)[seq_len(d)]
  vc$controls <- ctl
  vc
}

# the vine as the list of arrays that init.c reads
.vine_c <- function(vc, relaxed = vc$relaxed) {
  if (!inherits(vc, "vinecop_dist")) stop("model must be a vinecop_dist object", call. = FALSE)
  d <- vc$d
  nk <- .nk()
  aux <- vc$aux
  if (is.null(aux)) aux <- array(0, c(1L, d, d))
  vt <- if (is.null(vc$var_types)) rep("c", d) else vc$var_types
  v <- list(order = as.integer(vc$order), trunc = as.integer(vc$trunc),
            st = vc$struct_array, fam = vc$fam, rot = vc$rot, par = vc$par, aux = aux,
            disc = as.integer(vt == "d"), nrel = 0L, erel = matrix(0L, 2, 1),
            wa = matrix(0, d, 1), wk = matrix(0, nk, 1), pk = array(0, c(3, nk, 1)))
  if (!is.null(relaxed)) {
    v$nrel <- ncol(relaxed$erel)
    v$erel <- relaxed$erel
    v$wa <- relaxed$wa
    v$wk <- relaxed$wk
    v$pk <- relaxed$pk
  }
  storage.mode(v$st) <- storage.mode(v$fam) <- storage.mode(v$rot) <- "integer"
  storage.mode(v$erel) <- "integer"
  storage.mode(v$par) <- storage.mode(v$aux) <- storage.mode(v$wa) <- "double"
  storage.mode(v$wk) <- storage.mode(v$pk) <- "double"
  v
}

.vine_eval <- function(u, vc, what, cores, seed = 0L) {
  v <- .vine_c(vc)
  s <- .split(u, if (what == 2L) rep("c", vc$d) else vc$var_types %||% rep("c", vc$d))
  .Call(C_vine_eval, v, s$u, s$um, as.integer(what), as.integer(seed), .cores(cores))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Vine copula density, simulation and Rosenblatt transforms
#'
#' `dvinecop()` evaluates the density, `rosenblatt()` the Rosenblatt transform (which maps a
#' sample of the model to independent uniforms) and `inverse_rosenblatt()` its inverse.
#' `rvinecop()` simulates by applying the inverse Rosenblatt transform of the library to
#' independent uniforms from [stats::runif()], so it follows [set.seed()].
#'
#' For a vine with discrete variables (`var_types` of [vinecop()]), `u` / `x` have the
#' `d + k` columns described there, `dvinecop()` returns the mixed density (a probability
#' mass in the discrete directions) and `rosenblatt()` returns, for the discrete columns,
#' Brockwell's randomised transform \eqn{W F(x \mid \cdot) + (1 - W) F(x^- \mid \cdot)}, whose
#' uniforms \eqn{W} come from a seed drawn with R's generator (so [set.seed()] makes it
#' reproducible); `randomize_discrete = FALSE` returns \eqn{F(x \mid \cdot)} itself.
#' `inverse_rosenblatt()` and `rvinecop()` work on the continuous (copula) scale.
#'
#' Relaxed (non-simplified) edges from [vinecop_relax()] are used by `dvinecop()` and
#' [dvinecop_grid()]; the Rosenblatt transforms and `rvinecop()` use the simplified parameters
#' (the middle anchor), as the library does.
#'
#' @param u,x an `n x d` matrix in \eqn{(0, 1)^d} (`n x (d + k)` with discrete variables for
#'   `dvinecop()` and `rosenblatt()`).
#' @param vinecop,model a [vinecop_dist()] object.
#' @param n number of samples.
#' @param cores number of OpenMP threads.
#' @param randomize_discrete use the randomised transform for discrete variables.
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
dvinecop <- function(u, vinecop, cores = 1) .vine_eval(u, vinecop, 0L, cores)[1, ]

#' @rdname dvinecop
#' @export
rvinecop <- function(n, vinecop, cores = 1)
  inverse_rosenblatt(matrix(runif(n * vinecop$d), n, vinecop$d), vinecop, cores)

#' @rdname dvinecop
#' @export
rosenblatt <- function(x, model, cores = 1, randomize_discrete = TRUE) {
  disc <- any(model$var_types == "d")
  seed <- if (disc && isTRUE(randomize_discrete)) sample.int(.Machine$integer.max, 1L) else 0L
  t(.vine_eval(x, model, 1L, cores, seed))
}

#' @rdname dvinecop
#' @export
inverse_rosenblatt <- function(u, model, cores = 1) t(.vine_eval(u, model, 2L, cores))

#' Vine density on a grid of one variable
#'
#' Evaluates the vine density at points that differ only in variable `tgt`
#' (`fv_vinecop%pdf_grid`): for every row `i` of `u` and every grid value `g`, the density at
#' `u[i, ]` with column `tgt` replaced by `g`. Edges that do not involve `tgt` are evaluated
#' once per row instead of once per grid point, which makes predictive densities much cheaper
#' than [dvinecop()] at the same points (the values are the same). Relaxed edges are used
#' (exactly when `tgt` is the variable passed to [vinecop_relax()]; their conditioning values
#' are taken from the first grid point otherwise).
#'
#' @param u an `n x d` matrix of points (`n x (d + k)` with discrete variables, as in
#'   [vinecop()]).
#' @param vinecop a [vinecop_dist()] object.
#' @param tgt the variable (column, `1..d`) that varies.
#' @param grid the values of variable `tgt`: a vector (the same grid for every row) or an
#'   `n x nb` matrix (one grid per row).
#' @param grid_minus for a discrete `tgt`: \eqn{F(x^-)} at the grid points, same shape as
#'   `grid`.
#' @param cores number of OpenMP threads.
#' @return An `n x nb` matrix of densities.
#' @examples
#' vc <- vinecop_dist(list(list(bicop_dist("gaussian", 0, 0.6), bicop_dist("clayton", 0, 2))),
#'                    dvine_structure(1:3))
#' u <- rbind(c(0.3, 0.5, 0.6), c(0.8, 0.2, 0.4))
#' g <- c(0.1, 0.5, 0.9)
#' dvinecop_grid(u, vc, tgt = 2, grid = g)
#' dvinecop(cbind(0.3, g, 0.6), vc)                 # the first row, point by point
#' @export
dvinecop_grid <- function(u, vinecop, tgt, grid, grid_minus = NULL, cores = 1) {
  v <- .vine_c(vinecop)
  s <- .split(u, vinecop$var_types %||% rep("c", vinecop$d))
  d <- vinecop$d
  n <- ncol(s$u)
  if (length(tgt) != 1L || !tgt %in% seq_len(d))
    stop("tgt must be one of 1..", d, call. = FALSE)
  shape <- function(g) {
    g <- if (is.null(dim(g))) matrix(as.double(g), n, length(g), byrow = TRUE) else
      as.matrix(g)
    if (nrow(g) != n || ncol(g) < 1L) stop("grid must be a vector or an n x nb matrix",
                                          call. = FALSE)
    storage.mode(g) <- "double"
    g
  }
  G <- shape(grid)
  nb <- ncol(G)
  if (v$disc[tgt] == 1L) {
    if (is.null(grid_minus)) stop("variable ", tgt, " is discrete: give grid_minus",
                                  call. = FALSE)
    Gm <- shape(grid_minus)
    if (!identical(dim(Gm), dim(G))) stop("grid_minus must have the shape of grid",
                                         call. = FALSE)
  } else {
    Gm <- G
  }
  ug <- aperm(array(s$u, c(d, n, nb)), c(1L, 3L, 2L))     # d x nb x n
  ugm <- aperm(array(s$um, c(d, n, nb)), c(1L, 3L, 2L))
  ug[tgt, , ] <- t(G)
  ugm[tgt, , ] <- t(Gm)
  t(.Call(C_vine_pdf_grid, v, ug, ugm, as.integer(nb), as.integer(tgt), .cores(cores)))
}

#' Non-simplified vines
#'
#' Relaxes the simplifying assumption edge by edge (`fv_vinecop%relax`): every parametric edge
#' of a tree \eqn{t \ge 2} whose conditioned pair contains variable `tgt` (`tgt = 0`: any
#' edge) and whose conditioning set does not gets parameters that vary with \eqn{w}, the first
#' principal component of the normal scores of its conditioning variables. The family is
#' refitted with kernel weights at five anchors of \eqn{w} and the parameters are linear in
#' \eqn{w} between them. `mode = 1` relaxes only the edges whose test of constant conditional
#' correlation rejects at 5 %; `mode = 2` relaxes every eligible edge. Later trees are
#' refitted. Not available for vines with discrete variables.
#'
#' The relaxed edges are stored in `$relaxed` (`erel`: edge and tree, `wa`: loadings, `wk`:
#' anchors, `pk`: parameters at the anchors). [dvinecop()] and [dvinecop_grid()] use them;
#' the Rosenblatt transforms and simulation use the middle-anchor parameters, which are also
#' those shown by `summary()`.
#'
#' @param vinecop a fitted [vinecop_dist()] object (its edge log-likelihoods are needed).
#' @param data the data the vine was fitted to (`n x d`).
#' @param tgt the target variable (`0` for any).
#' @param mode `1` (test each edge) or `2` (relax every eligible edge).
#' @param weights optional observation weights.
#' @param cores number of OpenMP threads.
#' `vinecop_rel_params()` returns the parameters of the relaxed edges at given points
#' (`fv_vinecop%rel_params`). `vinecop_unrelax()` deletes the relaxed edges and keeps the
#' parameters `vinecop_relax()` wrote (as `fv_vinecop%unrelax`); its log-likelihood is
#' recomputed from `data`, or `NA` without it.
#'
#' @return `vinecop_relax()` and `vinecop_unrelax()`: a [vinecop_dist()] object, with
#'   `relaxed` and the updated `loglik` and `npars` (each relaxed edge counts its parameters
#'   at every anchor but one) or without them. `vinecop_rel_params()`: an
#'   `n x nrel x 3` array, `[i, j, ]` the parameters of relaxed edge `j` at `u[i, ]`.
#' @examples
#' set.seed(6)
#' u <- matrix(runif(1200), 400, 3)
#' r <- 0.9 * (2 * u[, 2] - 1)                       # the 1-3 dependence varies with u2
#' u[, 3] <- pnorm(r * qnorm(u[, 1]) + sqrt(1 - r^2) * qnorm(u[, 3]))
#' fit <- vinecop(u, family_set = "gaussian", structure = dvine_structure(c(1, 2, 3)))
#' rel <- vinecop_relax(fit, u, tgt = 0, mode = 2)
#' c(simplified = fit$loglik, relaxed = rel$loglik, edges = ncol(rel$relaxed$erel))
#' range(vinecop_rel_params(rel, u)[, 1, 1])        # the (1, 3 | 2) correlation along u2
#' vinecop_unrelax(rel, u)$loglik
#' @export
vinecop_relax <- function(vinecop, data, tgt = 0, mode = 1, weights = numeric(), cores = 1) {
  if (!inherits(vinecop, "vinecop_dist")) stop("vinecop must be a vinecop_dist object",
                                              call. = FALSE)
  if (any(vinecop$var_types == "d"))
    stop("relax is not available for vines with discrete variables", call. = FALSE)
  d <- vinecop$d
  if (length(tgt) != 1L || !tgt %in% 0:d) stop("tgt must be 0 or one of 1..", d, call. = FALSE)
  if (length(mode) != 1L || !mode %in% 1:2) stop("mode must be 1 or 2", call. = FALSE)
  s <- .split(data, vinecop$var_types %||% rep("c", d), check = TRUE)
  ctl <- vinecop$controls %||% .controls()
  ctl$ci[10] <- .cores(cores)
  ell <- vinecop$edge_ll %||% matrix(0, d, d)
  storage.mode(ell) <- "double"
  fs <- unique(c(0L, as.integer(vinecop$fam)))
  r <- .Call(C_vine_relax, .vine_c(vinecop, NULL), ell, s$u, .weights(weights, ncol(s$u)),
             as.integer(tgt), as.integer(mode), fs, ctl$ci, ctl$cr)
  m <- r$nrel
  rel <- if (m > 0L) list(erel = r$erel[, seq_len(m), drop = FALSE],
                          wa = r$wa[, seq_len(m), drop = FALSE],
                          wk = r$wk[, seq_len(m), drop = FALSE],
                          pk = r$pk[, , seq_len(m), drop = FALSE]) else NULL
  out <- .new_vine(vinecop$order, vinecop$trunc, vinecop$struct_array, r$fam, r$rot, r$par,
                   r$edge_ll, if (is.null(vinecop$aux)) NULL else r$aux, vinecop$var_types,
                   relaxed = rel, npars = r$npars)
  for (f in c("nobs", "selcrit", "threshold", "names", "controls")) out[[f]] <- vinecop[[f]]
  out$loglik <- r$loglik
  out
}

#' @rdname vinecop_relax
#' @param u points (`n x d`, or `n x (d + k)` with discrete variables) at which the parameters
#'   of the relaxed edges are wanted.
#' @export
vinecop_rel_params <- function(vinecop, u) {
  rel <- vinecop$relaxed
  if (is.null(rel)) stop("the vine has no relaxed edges", call. = FALSE)
  s <- .split(u, vinecop$var_types %||% rep("c", vinecop$d), check = TRUE)
  for (f in c("wa", "wk", "pk")) storage.mode(rel[[f]]) <- "double"
  aperm(.Call(C_vine_rel_params, rel$wa, rel$wk, rel$pk, s$u), c(3, 2, 1))
}

#' @rdname vinecop_relax
#' @export
vinecop_unrelax <- function(vinecop, data = NULL) {
  rel <- vinecop$relaxed
  if (is.null(rel)) return(vinecop)
  nk <- nrow(rel$wk)
  k <- vapply(seq_len(ncol(rel$erel)), function(j)
    .fam_npars(vinecop$fam[rel$erel[1, j], rel$erel[2, j]]), 0)
  vinecop$relaxed <- NULL
  vinecop$npars <- vinecop$npars - (nk - 1) * sum(k)
  vinecop$loglik <- if (is.null(data)) NA_real_ else sum(log(dvinecop(data, vinecop)))
  vinecop
}

.edge_table <- function(x) {
  lab <- if (is.null(x$names)) as.character(seq_len(x$d)) else x$names
  rel <- x$relaxed
  rows <- list()
  for (t in seq_len(x$trunc)) for (e in seq_len(x$d - t)) {
    b <- x$pair_copulas[[t]][[e]]
    st <- x$struct_array
    row <- data.frame(
      tree = t, edge = e,
      conditioned = paste(lab[x$order[e]], lab[x$order[st[e, t]]], sep = ", "),
      conditioning = paste(lab[x$order[st[e, seq_len(t - 1)]]], collapse = ", "),
      family = b$family, rotation = b$rotation, parameters = .fmt_par(b$parameters),
      npars = b$npars, loglik = if (is.null(b$loglik)) NA_real_ else b$loglik,
      stringsAsFactors = FALSE)
    if (!is.null(rel)) row$relaxed <- any(rel$erel[1, ] == e & rel$erel[2, ] == t)
    rows[[length(rows) + 1L]] <- row
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

#' @export
print.vinecop_dist <- function(x, ...) {
  cat("Vine copula ('vinecop_dist'): d = ", x$d, ", trunc_lvl = ", x$trunc, ", order = ",
      paste(x$order, collapse = " "), "\n", sep = "")
  cat("  npars = ", format(x$npars, digits = 4), if (!is.null(x$loglik))
        sprintf(", loglik = %.2f (n = %d)", x$loglik, x$nobs), "\n", sep = "")
  if (any(x$var_types == "d"))
    cat("  discrete variables: ", paste(which(x$var_types == "d"), collapse = " "), "\n",
        sep = "")
  if (!is.null(x$relaxed))
    cat("  non-simplified (relaxed) edges: ", ncol(x$relaxed$erel), "\n", sep = "")
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
  xm <- as.matrix(if (vec) matrix(x) else x)
  storage.mode(xm) <- "double"
  u <- t(.Call(C_pseudo_obs, t(xm), .cores(cores)))
  if (vec) drop(u) else {
    dimnames(u) <- dimnames(as.matrix(x))
    u
  }
}

#' Fit cache
#'
#' The library can store the result of every pair fit under a hash of its data and options,
#' so that repeated selections on the same data (for instance `vinecop(..., threshold = NA)`,
#' or refits with other vine options) reuse them. `fitcache_enable()` switches the cache on
#' with room for `capacity` fits, `fitcache_clear()` empties it, `fitcache_disable()` frees
#' it, and `fitcache_stats()` returns its counters. The cache is off by default.
#'
#' @param capacity number of fits the cache can hold.
#' @return A named numeric vector: `hits`, `misses` and `used` (slots in use), invisibly
#'   except for `fitcache_stats()`.
#' @examples
#' set.seed(7)
#' u <- rvinecop(300, vinecop_dist(list(list(bicop_dist("frank", 0, 4),
#'                                           bicop_dist("joe", 0, 2))), dvine_structure(1:3)))
#' fitcache_enable(1000)
#' f1 <- vinecop(u, family_set = "parametric")
#' f2 <- vinecop(u, family_set = "parametric")       # the same pair fits: cache hits
#' fitcache_stats()
#' fitcache_disable()
#' @export
fitcache_enable <- function(capacity = 20000L) {
  if (!is.numeric(capacity) || length(capacity) != 1L || !(capacity >= 1))
    stop("capacity must be a positive number", call. = FALSE)
  invisible(.Call(C_fitcache, 1L, as.integer(capacity)))
}

#' @rdname fitcache_enable
#' @export
fitcache_clear <- function() invisible(.Call(C_fitcache, 2L, 0L))

#' @rdname fitcache_enable
#' @export
fitcache_disable <- function() invisible(.Call(C_fitcache, 0L, 0L))

#' @rdname fitcache_enable
#' @export
fitcache_stats <- function() .Call(C_fitcache, 3L, 0L)
