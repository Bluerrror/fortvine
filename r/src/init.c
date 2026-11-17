/* init.c -- .Call entry points of the R package. Each one forwards to a bind(C) routine of the
 * fortvine library: the fvc_* routines of src/fortvine_capi.f90 (the C interface to the typed
 * API: weights, discrete margins, the nonparametric families, every control, relaxed vines,
 * pdf_grid, the fit cache and the family metadata) and fv_pseudo_obs of src/fortvine.f90.
 * The R side passes arrays already in the C-interface layout (observations in the last
 * Fortran dimension) and with the right storage mode; nothing is computed here.
 *
 * A vine is passed as an R list with the elements order, trunc, st, fam, rot, par, aux, disc,
 * nrel, erel, wa, wk, pk (built by .vine_c() in R/vinecop.R).
 */
#include <string.h>
#include <stdint.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

int fvc_naux(void);
int fvc_nk(void);
void fvc_fam_name(int fam, char *buf);
int fvc_families(int nmax, int *codes, int *npars);
void fvc_par_bounds(int fam, double *lb, double *ub);
void fvc_gen_flip(int *fam, double *par);
void fvc_fitcache(int op, int capacity, int64_t *hits, int64_t *misses, int *used);
void fvc_bicop_eval(int fam, int rot, const double *par, int naux, const double *aux, int n,
                    const double *u, int what, double *out);
void fvc_bicop_fit(int n, const double *u, int usew, const double *w, int sel, int nf,
                   const int *famset, const int *ci, const double *cr, int *fam, int *rot,
                   double *par, double *ll, int naux, double *aux, double *npe);
void fvc_vine_select(int n, int d, const double *u, const double *um, const int *disc,
                     int usew, const double *w, int fixed, int nf, const int *famset,
                     const int *ci, const double *cr, int *trunc, int *order, int *st,
                     int *fam, int *rot, double *par, double *edge_ll, int naux, double *aux,
                     double *ll, double *npe, double *thr, int *info);
void fvc_vine_eval(int d, int trunc, const int *order, const int *st, const int *fam,
                   const int *rot, const double *par, int naux, const double *aux,
                   const int *disc, int nrel, const int *erel, const double *wa,
                   const double *wk, const double *pk, int what, int n, const double *u,
                   const double *um, int seed, int nthreads, double *out);
void fvc_vine_pdf_grid(int d, int trunc, const int *order, const int *st, const int *fam,
                       const int *rot, const double *par, int naux, const double *aux,
                       const int *disc, int nrel, const int *erel, const double *wa,
                       const double *wk, const double *pk, int nb, int n, const double *ug,
                       const double *ugm, int tgt, int nthreads, double *out);
void fvc_vine_relax(int n, int d, const double *u, int usew, const double *w, int tgt,
                    int mode, int nf, const int *famset, const int *ci, const double *cr,
                    const int *trunc, const int *order, const int *st, int *fam, int *rot,
                    double *par, double *edge_ll, int naux, double *aux, const int *disc,
                    int nmax, int *nrel, int *erel, double *wa, double *wk, double *pk,
                    double *ll, double *npe);
void fvc_vine_structure(int kind, int d, const int *order, int trunc, int *st, int *trunc_out);
void fvc_vine_matrix(int d, const int *order, int trunc, const int *st, int *m);
void fvc_vine_rel_params(int d, int nrel, const double *wa, const double *wk, const double *pk,
                         int n, const double *uc, double *prel);
void fv_pseudo_obs(int n, int d, const double *x, double *out, int nthreads);

static SEXP named_list(int k, const char **nm, SEXP *val) {
  SEXP out = PROTECT(allocVector(VECSXP, k)), names = PROTECT(allocVector(STRSXP, k));
  for (int i = 0; i < k; i++) {
    SET_VECTOR_ELT(out, i, val[i]);
    SET_STRING_ELT(names, i, mkChar(nm[i]));
  }
  setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}

static void need(SEXP x, SEXPTYPE type, const char *what) {
  if ((SEXPTYPE) TYPEOF(x) != type)
    error("fortvine: argument '%s' has the wrong storage mode", what);
}

/* element `name` of the list x, checked for its storage mode */
static SEXP el(SEXP x, const char *name, SEXPTYPE type) {
  SEXP nm = getAttrib(x, R_NamesSymbol);
  if (TYPEOF(x) != VECSXP || TYPEOF(nm) != STRSXP) error("fortvine: the vine must be a list");
  for (int i = 0; i < LENGTH(x); i++)
    if (strcmp(CHAR(STRING_ELT(nm, i)), name) == 0) {
      SEXP v = VECTOR_ELT(x, i);
      need(v, type, name);
      return v;
    }
  error("fortvine: vine element '%s' is missing", name);
  return R_NilValue;
}

static SEXP mat_int(int d) { return allocMatrix(INTSXP, d, d); }

static void check_ctl(SEXP famset, SEXP ci, SEXP cr) {
  need(famset, INTSXP, "famset");
  need(ci, INTSXP, "ci");
  need(cr, REALSXP, "cr");
  if (LENGTH(ci) != 12 || LENGTH(cr) != 2 || LENGTH(famset) < 1)
    error("fortvine: controls have the wrong length");
}

/* ------------------------------------------------------------ metadata */
static SEXP C_families(void) {
  int dummy[1], n = fvc_families(0, dummy, dummy);
  SEXP code = PROTECT(allocVector(INTSXP, n)), np = PROTECT(allocVector(INTSXP, n));
  SEXP nm = PROTECT(allocVector(STRSXP, n));
  fvc_families(n, INTEGER(code), INTEGER(np));
  char buf[49];
  for (int i = 0; i < n; i++) {
    int k = 48;
    fvc_fam_name(INTEGER(code)[i], buf);
    while (k > 0 && buf[k - 1] == ' ') k--;
    buf[k] = '\0';
    if (INTEGER(np)[i] < 0) INTEGER(np)[i] = NA_INTEGER;
    SET_STRING_ELT(nm, i, mkChar(buf));
  }
  SEXP naux = PROTECT(ScalarInteger(fvc_naux())), nk = PROTECT(ScalarInteger(fvc_nk()));
  const char *lab[] = {"code", "name", "npars", "naux", "nk"};
  SEXP val[] = {code, nm, np, naux, nk};
  SEXP out = named_list(5, lab, val);
  UNPROTECT(5);
  return out;
}

static SEXP C_par_bounds(SEXP fam) {
  SEXP lb = PROTECT(allocVector(REALSXP, 3)), ub = PROTECT(allocVector(REALSXP, 3));
  fvc_par_bounds(asInteger(fam), REAL(lb), REAL(ub));
  const char *nm[] = {"lower", "upper"};
  SEXP val[] = {lb, ub};
  SEXP out = named_list(2, nm, val);
  UNPROTECT(2);
  return out;
}

static SEXP C_gen_flip(SEXP fam, SEXP par) {
  need(par, REALSXP, "par");
  if (LENGTH(par) != 3) error("fortvine: par must have length 3");
  SEXP f = PROTECT(ScalarInteger(asInteger(fam))), p = PROTECT(duplicate(par));
  fvc_gen_flip(INTEGER(f), REAL(p));
  const char *nm[] = {"family", "par"};
  SEXP val[] = {f, p};
  SEXP out = named_list(2, nm, val);
  UNPROTECT(2);
  return out;
}

/* op: 0 disable, 1 enable(capacity), 2 clear, 3 statistics only */
static SEXP C_fitcache(SEXP op, SEXP capacity) {
  int64_t h = 0, m = 0;
  int used = 0;
  fvc_fitcache(asInteger(op), asInteger(capacity), &h, &m, &used);
  SEXP out = PROTECT(allocVector(REALSXP, 3)), nm = PROTECT(allocVector(STRSXP, 3));
  REAL(out)[0] = (double) h;
  REAL(out)[1] = (double) m;
  REAL(out)[2] = (double) used;
  SET_STRING_ELT(nm, 0, mkChar("hits"));
  SET_STRING_ELT(nm, 1, mkChar("misses"));
  SET_STRING_ELT(nm, 2, mkChar("used"));
  setAttrib(out, R_NamesSymbol, nm);
  UNPROTECT(2);
  return out;
}

/* ------------------------------------------------------------ pair copula */
/* u: 2 x n double matrix; what as in fvc_bicop_eval (0 log density, ..., 6 distribution) */
static SEXP C_bicop_eval(SEXP fam, SEXP rot, SEXP par, SEXP aux, SEXP u, SEXP what) {
  need(par, REALSXP, "par");
  need(aux, REALSXP, "aux");
  need(u, REALSXP, "u");
  if (LENGTH(par) != 3 || LENGTH(aux) < 1) error("fortvine: par / aux have the wrong length");
  int n = LENGTH(u) / 2;
  SEXP out = PROTECT(allocVector(REALSXP, n));
  fvc_bicop_eval(asInteger(fam), asInteger(rot), REAL(par), LENGTH(aux), REAL(aux), n, REAL(u),
                 asInteger(what), REAL(out));
  UNPROTECT(1);
  return out;
}

/* sel = 0: fit fam / rot; sel = 1: select among famset. w: weights (length 0: none) */
static SEXP C_bicop_fit(SEXP u, SEXP w, SEXP sel, SEXP famset, SEXP ci, SEXP cr, SEXP fam,
                        SEXP rot, SEXP naux) {
  need(u, REALSXP, "u");
  need(w, REALSXP, "weights");
  check_ctl(famset, ci, cr);
  int n = LENGTH(u) / 2, usew = LENGTH(w) > 0, na = asInteger(naux);
  if (usew && LENGTH(w) != n) error("fortvine: weights must have one entry per observation");
  if (na < 1) error("fortvine: naux must be positive");
  SEXP f = PROTECT(ScalarInteger(asInteger(fam))), r = PROTECT(ScalarInteger(asInteger(rot)));
  SEXP par = PROTECT(allocVector(REALSXP, 3)), ll = PROTECT(allocVector(REALSXP, 1));
  SEXP ax = PROTECT(allocVector(REALSXP, na)), npe = PROTECT(allocVector(REALSXP, 1));
  fvc_bicop_fit(n, REAL(u), usew, usew ? REAL(w) : REAL(u), asInteger(sel), LENGTH(famset),
                INTEGER(famset), INTEGER(ci), REAL(cr), INTEGER(f), INTEGER(r), REAL(par),
                REAL(ll), na, REAL(ax), REAL(npe));
  const char *nm[] = {"family", "rotation", "par", "loglik", "aux", "npars"};
  SEXP val[] = {f, r, par, ll, ax, npe};
  SEXP out = named_list(6, nm, val);
  UNPROTECT(6);
  return out;
}

/* ------------------------------------------------------------ vine copula */
/* u, um: d x n; disc: d; order: d and st: d x d (copied; read when fixed = 1) */
static SEXP C_vine_select(SEXP u, SEXP um, SEXP disc, SEXP w, SEXP fixed, SEXP famset,
                          SEXP ci, SEXP cr, SEXP trunc, SEXP order, SEXP st_in, SEXP naux) {
  need(u, REALSXP, "u");
  need(um, REALSXP, "u_minus");
  need(disc, INTSXP, "disc");
  need(w, REALSXP, "weights");
  need(order, INTSXP, "order");
  need(st_in, INTSXP, "st");
  check_ctl(famset, ci, cr);
  int d = LENGTH(order), n = LENGTH(u) / d, usew = LENGTH(w) > 0, na = asInteger(naux);
  if (LENGTH(um) != LENGTH(u) || LENGTH(disc) != d || LENGTH(st_in) != d * d || na < 1)
    error("fortvine: inconsistent vine dimensions");
  if (usew && LENGTH(w) != n) error("fortvine: weights must have one entry per observation");
  SEXP tr = PROTECT(ScalarInteger(asInteger(trunc))), ord = PROTECT(duplicate(order));
  SEXP st = PROTECT(duplicate(st_in)), fam = PROTECT(mat_int(d)), rot = PROTECT(mat_int(d));
  SEXP par = PROTECT(alloc3DArray(REALSXP, 3, d, d));
  SEXP ell = PROTECT(allocMatrix(REALSXP, d, d));
  SEXP ax = PROTECT(alloc3DArray(REALSXP, na, d, d)), ll = PROTECT(allocVector(REALSXP, 1));
  SEXP npe = PROTECT(allocVector(REALSXP, 1)), thr = PROTECT(allocVector(REALSXP, 1));
  SEXP info = PROTECT(allocVector(INTSXP, 1));
  fvc_vine_select(n, d, REAL(u), REAL(um), INTEGER(disc), usew, usew ? REAL(w) : REAL(u),
                  asInteger(fixed), LENGTH(famset), INTEGER(famset), INTEGER(ci), REAL(cr),
                  INTEGER(tr), INTEGER(ord), INTEGER(st), INTEGER(fam), INTEGER(rot),
                  REAL(par), REAL(ell), na, REAL(ax), REAL(ll), REAL(npe), REAL(thr),
                  INTEGER(info));
  const char *nm[] = {"trunc", "order", "st", "fam", "rot", "par", "edge_ll", "aux", "loglik",
                      "npars", "threshold", "info"};
  SEXP val[] = {tr, ord, st, fam, rot, par, ell, ax, ll, npe, thr, info};
  SEXP out = named_list(12, nm, val);
  UNPROTECT(12);
  return out;
}

/* the vine arrays of the list v, with their sizes checked */
typedef struct {
  int d, trunc, naux, nrel;
  int *order, *st, *fam, *rot, *disc, *erel;
  double *par, *aux, *wa, *wk, *pk;
} vine_t;

static vine_t vine_of(SEXP v) {
  vine_t x;
  SEXP order = el(v, "order", INTSXP), aux = el(v, "aux", REALSXP);
  SEXP st = el(v, "st", INTSXP), fam = el(v, "fam", INTSXP), rot = el(v, "rot", INTSXP);
  SEXP par = el(v, "par", REALSXP), disc = el(v, "disc", INTSXP), erel = el(v, "erel", INTSXP);
  SEXP wa = el(v, "wa", REALSXP), wk = el(v, "wk", REALSXP), pk = el(v, "pk", REALSXP);
  int d = LENGTH(order), nk = fvc_nk();
  if (d < 1) error("fortvine: empty vine");
  x.d = d;
  x.trunc = asInteger(el(v, "trunc", INTSXP));
  x.nrel = asInteger(el(v, "nrel", INTSXP));
  x.naux = LENGTH(aux) / (d * d);
  int m = x.nrel > 1 ? x.nrel : 1;
  if (LENGTH(st) != d * d || LENGTH(fam) != d * d || LENGTH(rot) != d * d ||
      LENGTH(par) != 3 * d * d || x.naux < 1 || LENGTH(aux) != x.naux * d * d ||
      LENGTH(disc) != d || x.trunc < 0 || x.trunc > d - 1 || x.nrel < 0 ||
      LENGTH(erel) != 2 * m || LENGTH(wa) != d * m || LENGTH(wk) != nk * m ||
      LENGTH(pk) != 3 * nk * m)
    error("fortvine: inconsistent vine arrays");
  x.order = INTEGER(order); x.st = INTEGER(st); x.fam = INTEGER(fam); x.rot = INTEGER(rot);
  x.par = REAL(par); x.aux = REAL(aux); x.disc = INTEGER(disc); x.erel = INTEGER(erel);
  x.wa = REAL(wa); x.wk = REAL(wk); x.pk = REAL(pk);
  return x;
}

/* what: 0 density (row 1 of the result), 1 Rosenblatt, 2 inverse Rosenblatt; u, um: d x n */
static SEXP C_vine_eval(SEXP v, SEXP u, SEXP um, SEXP what, SEXP seed, SEXP nth) {
  vine_t x = vine_of(v);
  need(u, REALSXP, "u");
  need(um, REALSXP, "u_minus");
  if (LENGTH(um) != LENGTH(u) || LENGTH(u) % x.d != 0)
    error("fortvine: u and u_minus do not match the vine");
  int n = LENGTH(u) / x.d;
  SEXP out = PROTECT(allocMatrix(REALSXP, x.d, n));
  fvc_vine_eval(x.d, x.trunc, x.order, x.st, x.fam, x.rot, x.par, x.naux, x.aux, x.disc,
                x.nrel, x.erel, x.wa, x.wk, x.pk, asInteger(what), n, REAL(u), REAL(um),
                asInteger(seed), asInteger(nth), REAL(out));
  UNPROTECT(1);
  return out;
}

/* ug, ugm: d x nb x n; out: nb x n */
static SEXP C_vine_pdf_grid(SEXP v, SEXP ug, SEXP ugm, SEXP nb, SEXP tgt, SEXP nth) {
  vine_t x = vine_of(v);
  need(ug, REALSXP, "ug");
  need(ugm, REALSXP, "ugm");
  int b = asInteger(nb);
  if (b < 1 || LENGTH(ugm) != LENGTH(ug) || LENGTH(ug) % (x.d * b) != 0)
    error("fortvine: inconsistent grid dimensions");
  int n = LENGTH(ug) / (x.d * b);
  SEXP out = PROTECT(allocMatrix(REALSXP, b, n));
  fvc_vine_pdf_grid(x.d, x.trunc, x.order, x.st, x.fam, x.rot, x.par, x.naux, x.aux, x.disc,
                    x.nrel, x.erel, x.wa, x.wk, x.pk, b, n, REAL(ug), REAL(ugm),
                    asInteger(tgt), asInteger(nth), REAL(out));
  UNPROTECT(1);
  return out;
}

/* relax the vine v on u (d x n); its fam, rot, par, edge_ll, aux are copied and updated */
static SEXP C_vine_relax(SEXP v, SEXP edge_ll, SEXP u, SEXP w, SEXP tgt, SEXP mode,
                         SEXP famset, SEXP ci, SEXP cr) {
  vine_t x = vine_of(v);
  need(edge_ll, REALSXP, "edge_ll");
  need(u, REALSXP, "u");
  need(w, REALSXP, "weights");
  check_ctl(famset, ci, cr);
  int d = x.d, n = LENGTH(u) / d, usew = LENGTH(w) > 0, nk = fvc_nk(), nmax = d * d;
  if (LENGTH(edge_ll) != d * d || LENGTH(u) != d * n)
    error("fortvine: inconsistent vine arrays");
  if (usew && LENGTH(w) != n) error("fortvine: weights must have one entry per observation");
  SEXP fam = PROTECT(duplicate(el(v, "fam", INTSXP)));
  SEXP rot = PROTECT(duplicate(el(v, "rot", INTSXP)));
  SEXP par = PROTECT(duplicate(el(v, "par", REALSXP))), ell = PROTECT(duplicate(edge_ll));
  SEXP ax = PROTECT(duplicate(el(v, "aux", REALSXP))), nrel = PROTECT(allocVector(INTSXP, 1));
  SEXP erel = PROTECT(allocMatrix(INTSXP, 2, nmax));
  SEXP wa = PROTECT(allocMatrix(REALSXP, d, nmax)), wk = PROTECT(allocMatrix(REALSXP, nk, nmax));
  SEXP pk = PROTECT(alloc3DArray(REALSXP, 3, nk, nmax));
  SEXP ll = PROTECT(allocVector(REALSXP, 1)), npe = PROTECT(allocVector(REALSXP, 1));
  memset(INTEGER(erel), 0, sizeof(int) * 2 * nmax);
  memset(REAL(wa), 0, sizeof(double) * d * nmax);
  memset(REAL(wk), 0, sizeof(double) * nk * nmax);
  memset(REAL(pk), 0, sizeof(double) * 3 * nk * nmax);
  fvc_vine_relax(n, d, REAL(u), usew, usew ? REAL(w) : REAL(u), asInteger(tgt),
                 asInteger(mode), LENGTH(famset), INTEGER(famset), INTEGER(ci), REAL(cr),
                 &x.trunc, x.order, x.st, INTEGER(fam), INTEGER(rot), REAL(par), REAL(ell),
                 x.naux, REAL(ax), x.disc, nmax, INTEGER(nrel), INTEGER(erel), REAL(wa),
                 REAL(wk), REAL(pk), REAL(ll), REAL(npe));
  const char *nm[] = {"fam", "rot", "par", "edge_ll", "aux", "nrel", "erel", "wa", "wk", "pk",
                      "loglik", "npars"};
  SEXP val[] = {fam, rot, par, ell, ax, nrel, erel, wa, wk, pk, ll, npe};
  SEXP out = named_list(12, nm, val);
  UNPROTECT(12);
  return out;
}

/* kind 0: D-vine, 1: C-vine */
static SEXP C_vine_structure(SEXP kind, SEXP order, SEXP trunc) {
  need(order, INTSXP, "order");
  int d = LENGTH(order);
  SEXP st = PROTECT(mat_int(d)), tr = PROTECT(allocVector(INTSXP, 1));
  fvc_vine_structure(asInteger(kind), d, INTEGER(order), asInteger(trunc), INTEGER(st),
                     INTEGER(tr));
  const char *nm[] = {"st", "trunc"};
  SEXP val[] = {st, tr};
  SEXP out = named_list(2, nm, val);
  UNPROTECT(2);
  return out;
}

static SEXP C_vine_matrix(SEXP order, SEXP trunc, SEXP st) {
  need(order, INTSXP, "order");
  need(st, INTSXP, "st");
  int d = LENGTH(order);
  if (LENGTH(st) != d * d) error("fortvine: inconsistent structure array");
  SEXP m = PROTECT(mat_int(d));
  fvc_vine_matrix(d, INTEGER(order), asInteger(trunc), INTEGER(st), INTEGER(m));
  UNPROTECT(1);
  return m;
}

/* wa: d x nrel, wk: nk x nrel, pk: 3 x nk x nrel, u: d x n -> parameters 3 x nrel x n */
static SEXP C_vine_rel_params(SEXP wa, SEXP wk, SEXP pk, SEXP u) {
  need(wa, REALSXP, "wa");
  need(wk, REALSXP, "wk");
  need(pk, REALSXP, "pk");
  need(u, REALSXP, "u");
  int d = nrows(u), n = ncols(u), nrel = ncols(wa);
  if (nrows(wa) != d || ncols(wk) != nrel || LENGTH(pk) != 3 * nrows(wk) * nrel)
    error("fortvine: inconsistent relaxed-edge arrays");
  SEXP out = PROTECT(alloc3DArray(REALSXP, 3, nrel, n));
  fvc_vine_rel_params(d, nrel, REAL(wa), REAL(wk), REAL(pk), n, REAL(u), REAL(out));
  UNPROTECT(1);
  return out;
}

/* x: d x n */
static SEXP C_pseudo_obs(SEXP x, SEXP nth) {
  need(x, REALSXP, "x");
  int d = nrows(x), n = ncols(x);
  SEXP out = PROTECT(allocMatrix(REALSXP, d, n));
  fv_pseudo_obs(n, d, REAL(x), REAL(out), asInteger(nth));
  UNPROTECT(1);
  return out;
}

static const R_CallMethodDef call_methods[] = {
  {"C_families", (DL_FUNC) &C_families, 0},
  {"C_par_bounds", (DL_FUNC) &C_par_bounds, 1},
  {"C_gen_flip", (DL_FUNC) &C_gen_flip, 2},
  {"C_fitcache", (DL_FUNC) &C_fitcache, 2},
  {"C_bicop_eval", (DL_FUNC) &C_bicop_eval, 6},
  {"C_bicop_fit", (DL_FUNC) &C_bicop_fit, 9},
  {"C_vine_select", (DL_FUNC) &C_vine_select, 12},
  {"C_vine_eval", (DL_FUNC) &C_vine_eval, 6},
  {"C_vine_pdf_grid", (DL_FUNC) &C_vine_pdf_grid, 6},
  {"C_vine_relax", (DL_FUNC) &C_vine_relax, 9},
  {"C_vine_structure", (DL_FUNC) &C_vine_structure, 3},
  {"C_vine_matrix", (DL_FUNC) &C_vine_matrix, 3},
  {"C_vine_rel_params", (DL_FUNC) &C_vine_rel_params, 4},
  {"C_pseudo_obs", (DL_FUNC) &C_pseudo_obs, 2},
  {NULL, NULL, 0}
};

void R_init_fortvine(DllInfo *dll) {
  R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
