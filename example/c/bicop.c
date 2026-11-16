/* bicop.c -- one pair copula through the fortvine C API: evaluate, simulate, fit, select.
 *
 * Build the shared library (bash build.sh), then from the repository root:
 *   Linux:   gcc -O2 example/c/bicop.c -Lbuild -lfortvine -Wl,-rpath,"$PWD/build" -o build/bicop_c
 *   macOS:   gcc -O2 example/c/bicop.c -Lbuild -lfortvine -Wl,-rpath,"$PWD/build" -o build/bicop_c
 *   Windows (MinGW-w64): gcc -O2 example/c/bicop.c build/libfortvine.dll -o build/bicop_c.exe
 *
 * Arrays follow the kernels' convention: observations in the last Fortran dimension, so the
 * pair (u1, u2) of observation i sits at u[2*i], u[2*i+1] (a row-major n x 2 array).
 */
#include <stdio.h>
#include <stdlib.h>

/* bind(C) routines of src/fortvine.f90 (integer = int, real(c_double) = double) */
int fv_version(void);
void fv_bicop_eval(int fam, int rot, const double *par, int n, const double *u, int what,
                   double *out, int nthreads);
void fv_bicop_fit(int n, const double *u, int fam, int rot, const int *ci, double *par,
                  double *ll, int nthreads);
void fv_bicop_select(int n, const double *u, int nf, const int *famset, const int *ci,
                     double psi0, int *fam, int *rot, double *par, double *ll, int nthreads);
double fv_kendall_tau(int n, const double *x, const double *y);

enum { CLAYTON = 3, N = 2000, NTH = 4 };
enum { LOGPDF = 0, PDF = 1, HFUNC1 = 2, HFUNC2 = 3, HINV1 = 4, HINV2 = 5 };

/* uniform in (0, 1) from a 64-bit xorshift generator (deterministic) */
static double unif(unsigned long long *s) {
  *s ^= *s << 13; *s ^= *s >> 7; *s ^= *s << 17;
  return ((*s >> 11) + 0.5) / 9007199254740992.0;
}

int main(void) {
  static const char *name[] = {"indep", "gaussian", "student", "clayton", "gumbel", "frank",
    "joe", "bb1", "bb6", "bb7", "bb8", "tawn", "plackett", "fgm", "amh", "galambos",
    "husler_reiss", "bb2", "bb3", "bb5", "bb9", "bb10", "bb4"};
  /* ci: criterion (2 BIC), method (0 MLE), preselect, rotations, fast fits, exact tau */
  const int ci[6] = {2, 0, 1, 1, 0, 1};
  double par[3] = {2.0, 0.0, 0.0}, est[3], ll, diag[6] = {0.1, 0.1, 0.5, 0.5, 0.9, 0.9}, dens[3];
  double *uw = malloc(2 * N * sizeof *uw), *u = malloc(2 * N * sizeof *u), *v = malloc(N * sizeof *v);
  double *x = malloc(N * sizeof *x), *y = malloc(N * sizeof *y);
  int famset[23], fam, rot, i;
  unsigned long long seed = 88172645463325252ULL;

  printf("fortvine C API version %d\n", fv_version());
  fv_bicop_eval(CLAYTON, 0, par, 3, diag, PDF, dens, 1);
  printf("Clayton(2) density on the diagonal: %.4f %.4f %.4f\n", dens[0], dens[1], dens[2]);

  /* simulate: u2 = h1^-1(w | u1) for independent uniforms (u1, w) */
  for (i = 0; i < N; i++) { uw[2 * i] = unif(&seed); uw[2 * i + 1] = unif(&seed); }
  fv_bicop_eval(CLAYTON, 0, par, N, uw, HINV1, v, NTH);
  for (i = 0; i < N; i++) {
    u[2 * i] = x[i] = uw[2 * i];
    u[2 * i + 1] = y[i] = v[i];
  }
  printf("empirical Kendall tau: %.3f\n", fv_kendall_tau(N, x, y));

  fv_bicop_fit(N, u, CLAYTON, 0, ci, est, &ll, NTH);
  printf("fitted theta: %.4f   log-likelihood: %.2f\n", est[0], ll);

  for (i = 0; i < 23; i++) famset[i] = i;          /* families 0..22 */
  fv_bicop_select(N, u, 23, famset, ci, 0.9, &fam, &rot, est, &ll, NTH);
  printf("selected: %s rotation %d par %.3f %.3f %.3f\n", name[fam], rot, est[0], est[1], est[2]);

  free(uw); free(u); free(v); free(x); free(y);
  return 0;
}
