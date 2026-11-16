# Family table and the conversion of family names / codes / sets to the integer codes of the
# C interface. The table itself comes from the Fortran library (fortvine_r.f90).

.fv <- new.env(parent = emptyenv())

#' List the pair-copula families of fortvine
#'
#' Returns the 111 distinct families of the library: the 26 hand-coded families (codes 0-25)
#' and the 85 generic families (codes 10000 and above). The names and parameter counts are
#' read from the Fortran library.
#'
#' @return A data frame with columns `code`, `name` and `npars` (`NA` for the nonparametric
#'   families `tll` and `bernstein`, which the C interface cannot use).
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
  code <- tab$code[i]
  if (any(code %in% 24:25))
    stop("the nonparametric families 'tll' (24) and 'bernstein' (25) need their estimate ",
         "(aux vector), which the fortvine C interface does not pass; use the Fortran API",
         call. = FALSE)
  code
}

# family_set: names, codes and the groups "all" (0-22, fortvine's default set), "parametric"
# (0-11, vinecopulib's parametric set) and "generic" (the 85 generic families)
.family_set <- function(family_set) {
  tab <- fortvine_families()
  groups <- list(all = 0:22, parametric = 0:11, generic = tab$code[tab$code >= 10000])
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
