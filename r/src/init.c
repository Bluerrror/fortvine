/* init.c -- .Call entry points of the R package. Each one forwards to a bind(C) routine of
 * the fortvine library (src/fortvine.f90, documented in docs/manual/c-api.md) or of the
 * metadata shim fortvine_r.f90. The R side passes arrays already in the C-interface layout
 * (observations in the last Fortran dimension) and with the right storage mode.
 */
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void fv_bicop_eval(int fam, int rot, const double *par, int n, const double *u, int what,
                   double *out, int nthreads);
void fv_bicop_fit(int n, const double *u, int fam, int rot, const int *ci, double *par,
                  double *ll, int nthreads);
void fv_bicop_select(int n, const double *u, int nf, const int *famset, const int *ci,
                     double psi0, int *fam, int *rot, double *par, double *ll, int nthreads);
void fv_vine_select(int n, int d, const double *u, int nf, const int *famset, const int *ci,
                    const double *cr, int *order_io, const int *struct_in, int *struct_out,
                    int *fam_out, int *rot_out, double *par_out, double *ll_out, int nthreads,
                    int *info);
void fv_vine_pdf(int nsamp, int bins, int d, int tgt, const double *u, const int *order,
                 int trunc, const int *st, const int *fam, const int *rot, const double *par,
                 double *out, int nthreads);
void fv_vine_rosenblatt(int n, int d, const double *u, const int *order, int trunc,
                        const int *st, const int *fam, const int *rot, const double *par,
                        double *out, int nthreads);
void fv_vine_inverse_rosenblatt(int n, int d, const double *u, const int *order, int trunc,
                                const int *st, const int *fam, const int *rot,
                                const double *par, double *out, int nthreads);
void fv_pseudo_obs(int n, int d, const double *x, double *out, int nthreads);
void fvr_families(int nmax, int *n, int *codes, int *npars, char *names);
void fvr_set_dvine(int d, const int *order, int trunc, int *st);

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

static SEXP mat_int(int d) { return allocMatrix(INTSXP, d, d); }

/* u: 2 x n double matrix */
static SEXP C_bicop_eval(SEXP fam, SEXP rot, SEXP par, SEXP u, SEXP what, SEXP nth) {
  int n = LENGTH(u) / 2;
  SEXP out = PROTECT(allocVector(REALSXP, n));
  fv_bicop_eval(asInteger(fam), asInteger(rot), REAL(par), n, REAL(u), asInteger(what),
                REAL(out), asInteger(nth));
  UNPROTECT(1);
  return out;
}

static SEXP C_bicop_fit(SEXP u, SEXP fam, SEXP rot, SEXP ci, SEXP nth) {
  SEXP par = PROTECT(allocVector(REALSXP, 3)), ll = PROTECT(allocVector(REALSXP, 1));
  fv_bicop_fit(LENGTH(u) / 2, REAL(u), asInteger(fam), asInteger(rot), INTEGER(ci), REAL(par),
               REAL(ll), asInteger(nth));
  const char *nm[] = {"par", "loglik"};
  SEXP val[] = {par, ll};
  SEXP out = named_list(2, nm, val);
  UNPROTECT(2);
  return out;
}

static SEXP C_bicop_select(SEXP u, SEXP famset, SEXP ci, SEXP psi0, SEXP nth) {
  SEXP fam = PROTECT(allocVector(INTSXP, 1)), rot = PROTECT(allocVector(INTSXP, 1));
  SEXP par = PROTECT(allocVector(REALSXP, 3)), ll = PROTECT(allocVector(REALSXP, 1));
  fv_bicop_select(LENGTH(u) / 2, REAL(u), LENGTH(famset), INTEGER(famset), INTEGER(ci),
                  asReal(psi0), INTEGER(fam), INTEGER(rot), REAL(par), REAL(ll), asInteger(nth));
  const char *nm[] = {"family", "rotation", "par", "loglik"};
  SEXP val[] = {fam, rot, par, ll};
  SEXP out = named_list(4, nm, val);
  UNPROTECT(4);
  return out;
}

/* u: d x n; order: d (in/out, copied); st_in: d x d */
static SEXP C_vine_select(SEXP u, SEXP famset, SEXP ci, SEXP cr, SEXP order, SEXP st_in,
                          SEXP nth) {
  int d = LENGTH(order), n = LENGTH(u) / d;
  SEXP ord = PROTECT(duplicate(order)), st = PROTECT(mat_int(d)), fam = PROTECT(mat_int(d));
  SEXP rot = PROTECT(mat_int(d)), ll = PROTECT(allocMatrix(REALSXP, d, d));
  SEXP par = PROTECT(alloc3DArray(REALSXP, 3, d, d)), info = PROTECT(allocVector(INTSXP, 1));
  fv_vine_select(n, d, REAL(u), LENGTH(famset), INTEGER(famset), INTEGER(ci), REAL(cr),
                 INTEGER(ord), INTEGER(st_in), INTEGER(st), INTEGER(fam), INTEGER(rot),
                 REAL(par), REAL(ll), asInteger(nth), INTEGER(info));
  const char *nm[] = {"order", "st", "fam", "rot", "par", "edge_ll", "info"};
  SEXP val[] = {ord, st, fam, rot, par, ll, info};
  SEXP out = named_list(7, nm, val);
  UNPROTECT(7);
  return out;
}

/* what: 0 density, 1 Rosenblatt, 2 inverse Rosenblatt; u: d x n */
static SEXP C_vine_eval(SEXP u, SEXP order, SEXP trunc, SEXP st, SEXP fam, SEXP rot, SEXP par,
                        SEXP what, SEXP nth) {
  int d = LENGTH(order), n = LENGTH(u) / d, w = asInteger(what);
  SEXP out = PROTECT(w == 0 ? allocVector(REALSXP, n) : allocMatrix(REALSXP, d, n));
  if (w == 0)
    fv_vine_pdf(n, 1, d, -1, REAL(u), INTEGER(order), asInteger(trunc), INTEGER(st),
                INTEGER(fam), INTEGER(rot), REAL(par), REAL(out), asInteger(nth));
  else if (w == 1)
    fv_vine_rosenblatt(n, d, REAL(u), INTEGER(order), asInteger(trunc), INTEGER(st),
                       INTEGER(fam), INTEGER(rot), REAL(par), REAL(out), asInteger(nth));
  else
    fv_vine_inverse_rosenblatt(n, d, REAL(u), INTEGER(order), asInteger(trunc), INTEGER(st),
                               INTEGER(fam), INTEGER(rot), REAL(par), REAL(out),
                               asInteger(nth));
  UNPROTECT(1);
  return out;
}

/* x: d x n */
static SEXP C_pseudo_obs(SEXP x, SEXP nth) {
  int d = nrows(x), n = ncols(x);
  SEXP out = PROTECT(allocMatrix(REALSXP, d, n));
  fv_pseudo_obs(n, d, REAL(x), REAL(out), asInteger(nth));
  UNPROTECT(1);
  return out;
}

static SEXP C_families(void) {
  enum { NMAX = 512, LEN = 48 };
  int n, codes[NMAX], npars[NMAX];
  char *names = R_alloc(LEN * NMAX, 1), buf[LEN + 1];
  fvr_families(NMAX, &n, codes, npars, names);
  if (n < 0) error("fortvine: more than %d families", NMAX);
  SEXP code = PROTECT(allocVector(INTSXP, n)), np = PROTECT(allocVector(INTSXP, n));
  SEXP nm = PROTECT(allocVector(STRSXP, n));
  for (int i = 0; i < n; i++) {
    int k = LEN;
    memcpy(buf, names + LEN * i, LEN);
    while (k > 0 && buf[k - 1] == ' ') k--;
    buf[k] = '\0';
    INTEGER(code)[i] = codes[i];
    INTEGER(np)[i] = npars[i] < 0 ? NA_INTEGER : npars[i];
    SET_STRING_ELT(nm, i, mkChar(buf));
  }
  const char *lab[] = {"code", "name", "npars"};
  SEXP val[] = {code, nm, np};
  SEXP out = named_list(3, lab, val);
  UNPROTECT(3);
  return out;
}

static SEXP C_set_dvine(SEXP order, SEXP trunc) {
  int d = LENGTH(order);
  SEXP st = PROTECT(mat_int(d));
  fvr_set_dvine(d, INTEGER(order), asInteger(trunc), INTEGER(st));
  UNPROTECT(1);
  return st;
}

static const R_CallMethodDef call_methods[] = {
  {"C_bicop_eval", (DL_FUNC) &C_bicop_eval, 6},
  {"C_bicop_fit", (DL_FUNC) &C_bicop_fit, 5},
  {"C_bicop_select", (DL_FUNC) &C_bicop_select, 5},
  {"C_vine_select", (DL_FUNC) &C_vine_select, 7},
  {"C_vine_eval", (DL_FUNC) &C_vine_eval, 9},
  {"C_pseudo_obs", (DL_FUNC) &C_pseudo_obs, 2},
  {"C_families", (DL_FUNC) &C_families, 0},
  {"C_set_dvine", (DL_FUNC) &C_set_dvine, 2},
  {NULL, NULL, 0}
};

void R_init_fortvine(DllInfo *dll) {
  R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
