/* example_interop.c -- the C side of the interop figure (tools/example_figures.py interop).
 * Reads interop_specs.txt and interop_points.txt written by tools/example_interop.f90 from the
 * current directory, evaluates pdf, h1 and inverse h1 of every copula through the C API
 * (fv_bicop_eval) and writes interop_c.txt with the layout of interop_fortran.txt.
 *   gcc -O2 tools/example_interop.c build/libfortvine.dll -o build/example_interop_c.exe   (MinGW)
 *   gcc -O2 tools/example_interop.c -Lbuild -lfortvine -Wl,-rpath,"$PWD/build" -o build/example_interop_c
 */
#include <stdio.h>
#include <stdlib.h>

void fv_bicop_eval(int fam, int rot, const double *par, int n, const double *u, int what,
                   double *out, int nthreads);

enum { MAXC = 64, MAXP = 100000 };

int main(void) {
  int fam[MAXC], rot[MAXC], nc = 0, m = 0, k, i, j;
  double par[MAXC][3], *u = malloc(2 * MAXP * sizeof *u), *out;
  FILE *f = fopen("interop_specs.txt", "r");
  if (!f) return 1;
  while (nc < MAXC && fscanf(f, "%d %d %lf %lf %lf", &fam[nc], &rot[nc], &par[nc][0], &par[nc][1], &par[nc][2]) == 5) nc++;
  fclose(f);
  f = fopen("interop_points.txt", "r");
  if (!f) return 1;
  while (m < MAXP && fscanf(f, "%lf %lf", &u[2 * m], &u[2 * m + 1]) == 2) m++;
  fclose(f);
  out = malloc((size_t)m * 3 * nc * sizeof *out);
  for (k = 0; k < nc; k++)
    for (j = 0; j < 3; j++) {                      /* what = 1 pdf, 2 h1, 4 h1 inverse */
      static const int what[3] = {1, 2, 4};
      double *col = malloc(m * sizeof *col);
      fv_bicop_eval(fam[k], rot[k], par[k], m, u, what[j], col, 1);
      for (i = 0; i < m; i++) out[(size_t)i * 3 * nc + 3 * k + j] = col[i];
      free(col);
    }
  f = fopen("interop_c.txt", "w");
  for (i = 0; i < m; i++) {
    for (k = 0; k < 3 * nc; k++) fprintf(f, " %.17e", out[(size_t)i * 3 * nc + k]);
    fprintf(f, "\n");
  }
  fclose(f);
  printf("C: %d copulas x %d points\n", nc, m);
  free(u); free(out);
  return 0;
}
