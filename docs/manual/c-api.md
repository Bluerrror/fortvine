# C interface

The kernel module `fortvine` exports 13 routines with `bind(C)` names, so the library can be
called from C, C++, Python (ctypes, cffi), Julia (`ccall`) or R (`.C` / `.Call`) without a
Fortran compiler on the caller's side. They cover pair-copula evaluation, fitting and selection,
vine selection, vine density and Rosenblatt transforms, and the rank utilities.

The C interface is lower-level than the [typed Fortran API](../api/index.md): no derived
types and no optional arguments. A second module, `fortvine_capi` (`src/fortvine_capi.f90`),
adds C routines over the typed API for everything the 13 kernel routines leave out:
observation weights, discrete margins, TLL and Bernstein, every control, non-simplified vines,
`pdf_grid`, the fit cache and the abort flag (see [The extended interface](#the-extended-interface)).
The Python and R packages are built on both.

## Conventions

**Linking.** `bash build.sh` produces `build/libfortvine.a` and a shared library
(`libfortvine.so`, `.dylib` or `.dll`). With the static archive, link the Fortran and OpenMP
runtimes as well, e.g. `gcc prog.c build/libfortvine.a -lgfortran -lquadmath -fopenmp`. On
Windows the DLL depends on the MinGW runtime DLLs (`libgomp-1.dll`, `libquadmath-0.dll`,
`libwinpthread-1.dll`), which must be on the DLL search path.

**Arguments.** Scalars marked *value* are passed by value; everything else (arrays and scalar
outputs) by pointer. Integers are C `int`, reals C `double`.

**Array layout: observations are the last Fortran dimension.** Fortran arrays are
column-major, so an array declared `u(2, n)` in Fortran is, in C or NumPy (row-major), a
contiguous `n x 2` array in which the two values of one observation are adjacent. In general:

| Fortran declaration | C / NumPy row-major view | meaning |
|---|---|---|
| `u(2, n)` | `double u[n][2]`, NumPy `(n, 2)` | pair data, one row per observation |
| `u(d, n)` | `double u[n][d]`, NumPy `(n, d)` | data matrix, one row per observation |
| `u(d, bins, nsamp)` | NumPy `(nsamp, bins, d)` | grid for `fv_vine_pdf` |
| `st(d, d)`, `fam(d, d)`, `rot(d, d)` = `(edge, tree)` | NumPy `(d, d)` indexed `[tree, edge]` | vine edge arrays |
| `par(3, d, d)` = `(k, edge, tree)` | NumPy `(d, d, 3)` indexed `[tree, edge, k]` | vine parameters |
| `x0(3, njob)` | NumPy `(njob, 3)` | per-job parameter vectors |

So the usual NumPy/PyTorch data matrix of shape `(n, d)` can be passed as is. Column-major
callers (Julia, R) pass the transpose: a `d x n` matrix. Indices stored *in* the arrays (`order`,
`st`) are 1-based variable labels in natural order, exactly as in
[`fv_vinecop`](vinecop.md#how-the-structure-is-stored); `tgt` of `fv_vine_pdf` and the dataset
index `ds` of `fv_fit_jobs` are 0-based.

**Control vectors.** Options that the Fortran API keeps in `fv_controls` are passed as small
integer and real arrays:

| entry | `ci(6)` (pair routines) and `ci(1:6)` of vine selection |
|---|---|
| `ci(1)` | selection criterion: 0 loglik, 1 AIC, 2 BIC, 3 mBIC, 4 mBICv |
| `ci(2)` | method: 0 MLE, 1 tau inversion |
| `ci(3)` | preselect: 0 / 1 |
| `ci(4)` | allow rotations: 0 / 1 |
| `ci(5)` | fast MLE: 0 / 1 |
| `ci(6)` | tau mode: 0 auto, 1 exact, 2 approximate |

| entry | additional entries of `ci(10)` for `fv_vine_select` |
|---|---|
| `ci(7)` | tree criterion: 0 tau, 1 rho, 2 joe, 3 Hoeffding, 4 xi, 5 MI, 6 Ferreira |
| `ci(8)` | tree algorithm: 0 Prim, 1 Kruskal |
| `ci(9)` | truncation level, $-1$ for $d - 1$ |
| `ci(10)` | 1 if a fixed structure is given in `order_io` / `struct_in`, else 0 |

`cr(2)` of `fv_vine_select` holds `[threshold, psi0]`. The defaults of `fv_controls`
correspond to `ci = [2, 0, 1, 1, 0, 1, 0, 0, -1, 0]`, `cr = [0.0, 0.9]`.

## Routines

### `fv_version`

```c
int fv_version(void);
```

Returns the interface version, currently 1.

### `fv_bicop_eval`

```c
void fv_bicop_eval(int fam, int rot, const double par[3], int n, const double *u,
                   int what, double *out, int nthreads);
```

| argument | | description |
|---|---|---|
| `fam`, `rot` | value | family code and rotation |
| `par[3]` | in | parameters |
| `n` | value | number of points |
| `u` | in, `u(2, n)` | points |
| `what` | value | 0 logpdf, 1 pdf, 2 hfunc1, 3 hfunc2, 4 hinv1, 5 hinv2 |
| `out` | out, `out(n)` | result |
| `nthreads` | value | threads (used for $n > 4096$) |

`hinv1` returns $v$ with $h_1(u_1, v) = u_2$, `hinv2` returns $v$ with $h_2(v, u_2) = u_1$,
as in [`fv_bicop`](bicop.md#h-functions-and-their-inverses). Not usable for TLL/Bernstein (no
`aux` argument).

### `fv_bicop_pdf_hfuncs`

```c
void fv_bicop_pdf_hfuncs(int fam, int rot, const double par[3], int n, const double *u,
                         double *pdf, double *h1, double *h2, int nthreads);
```

Density and both h-functions in one pass, `pdf(n)`, `h1(n)`, `h2(n)` (cheaper than three
`fv_bicop_eval` calls). Parametric families only.

### `fv_bicop_fit`

```c
void fv_bicop_fit(int n, const double *u, int fam, int rot, const int ci[6],
                  double par[3], double *ll, int nthreads);
```

Fits family `fam` with rotation `rot` to `u(2, n)` as
[`fv_bicop%fit`](bicop.md#fitting-a-given-family) does; returns `par` and the log-likelihood
`ll`. Of `ci` only the method, fast-MLE and tau-mode entries matter.

### `fv_bicop_select`

```c
void fv_bicop_select(int n, const double *u, int nf, const int *famset, const int ci[6],
                     double psi0, int *fam, int *rot, double par[3], double *ll, int nthreads);
```

Selects among the `nf` family codes in `famset` (all rotations allowed by `ci`), as
[`fv_bicop%select`](bicop.md#selecting-the-family). `psi0` (value) is the mBIC prior. Outputs:
`fam`, `rot`, `par`, `ll`.

### `fv_fit_jobs`

```c
void fv_fit_jobs(int njob, int n, int nds, const double *u, const int *ds, const int *fam,
                 const int *rot, const int *method, const double *x0, const double *lb,
                 const double *ub, int usew, const double *w, double *par, double *ll,
                 int nthreads);
```

Runs `njob` independent fits in parallel, each on one of `nds` datasets of equal length `n`.

| argument | | description |
|---|---|---|
| `njob`, `n`, `nds` | value | number of jobs, observations per dataset, number of datasets |
| `u` | in, `u(2, n, nds)` | the datasets (NumPy `(nds, n, 2)`) |
| `ds(njob)` | in | 0-based dataset index of each job |
| `fam(njob)`, `rot(njob)` | in | family and rotation of each job |
| `method(njob)` | in | 0 evaluate the log-likelihood at `x0`; 1 MLE within `[lb, ub]` starting at `x0`; 2 Student with $\rho$ fixed at `x0(1)`, $\nu$ by MLE |
| `x0`, `lb`, `ub` | in, `(3, njob)` | start (or evaluation point) and box of each job |
| `usew` | value | 1 if `w` holds observation weights (shared by all datasets) |
| `w(n)` | in | weights (ignored if `usew = 0`, but must be a valid pointer) |
| `par`, `ll` | out, `(3, njob)`, `(njob)` | estimates and log-likelihoods |

This is the routine to use when you want to drive the search yourself (custom starts, boxes,
profile likelihoods). Parametric families only.

### `fv_vine_select`

```c
void fv_vine_select(int n, int d, const double *u, int nf, const int *famset,
                    const int ci[10], const double cr[2], int *order_io,
                    const int *struct_in, int *struct_out, int *fam_out, int *rot_out,
                    double *par_out, double *ll_out, int nthreads, int *info);
```

| argument | | description |
|---|---|---|
| `n`, `d` | value | observations, dimension |
| `u` | in, `u(d, n)` | pseudo-observations (NumPy `(n, d)`) |
| `nf`, `famset(nf)` | value, in | candidate families |
| `ci(10)`, `cr(2)` | in | controls, see [above](#conventions) |
| `order_io(d)` | in/out | with `ci(10) = 1`: the fixed order; on return the order of the selected vine |
| `struct_in(d, d)` | in | with `ci(10) = 1`: the fixed natural-order structure array; otherwise ignored (still pass a `d x d` array) |
| `struct_out`, `fam_out`, `rot_out` | out, `(d, d)` | structure, families, rotations, `(edge, tree)` |
| `par_out` | out, `(3, d, d)` | parameters |
| `ll_out` | out, `(d, d)` | edge log-likelihoods |
| `nthreads` | value | threads |
| `info` | out | 0 on success; 2, 3, 4 indicate an internal failure |

The number of trees fitted is `ci(9)` (or $d - 1$ for $-1$); pass the same value as `trunc`
to the evaluation routines below. Rows of `struct_out` beyond that level are 0.

### `fv_vine_pdf`

```c
void fv_vine_pdf(int nsamp, int bins, int d, int tgt, const double *u, const int *order,
                 int trunc, const int *st, const int *fam, const int *rot, const double *par,
                 double *out, int nthreads);
```

Vine density on `u(d, bins, nsamp)` in original variable order; `out(bins, nsamp)`. Variable
`tgt` (0-based; $-1$ for none) is the only one that may vary across the `bins` of a sample,
which lets the routine evaluate the other edges once per sample (as
[`pdf_grid`](vinecop.md#density-and-log-likelihood)). For a plain density call use
`bins = 1`, `tgt = -1`.

### `fv_vine_rosenblatt`, `fv_vine_inverse_rosenblatt`

```c
void fv_vine_rosenblatt(int n, int d, const double *u, const int *order, int trunc,
                        const int *st, const int *fam, const int *rot, const double *par,
                        double *out, int nthreads);
void fv_vine_inverse_rosenblatt(int n, int d, const double *u, const int *order, int trunc,
                                const int *st, const int *fam, const int *rot,
                                const double *par, double *out, int nthreads);
```

Rosenblatt transform and its inverse, `u(d, n)` to `out(d, n)`, see
[Vine copulas](vinecop.md#rosenblatt-transform-and-simulation). To simulate, apply the inverse
to independent uniforms.

### `fv_pseudo_obs`

```c
void fv_pseudo_obs(int n, int d, const double *x, double *out, int nthreads);
```

Column-wise pseudo-observations of `x(d, n)`: average ranks (ties share the mean rank) divided
by $n + 1$.

### `fv_kendall_tau`, `fv_tau_matrix`

```c
double fv_kendall_tau(int n, const double *x, const double *y);
void   fv_tau_matrix(int n, int d, const double *x, double *tm, int nthreads);
```

Kendall's $\tau_b$ (tie-corrected, $O(n \log n)$) of two vectors, and the symmetric `d x d`
matrix of all pairwise values for `x(d, n)` (diagonal 1).

## The extended interface

`fortvine_capi` passes a pair copula or a vine as its arrays and rebuilds the
`fv_bicop` / `fv_vinecop` object inside each call, so the caller owns all state (the Python
and R objects are plain arrays and can be saved). Each routine only converts its arguments and
calls the typed API.

**Controls.** One integer vector `ci(12)` and `cr(2) = [threshold, psi0]` carry every field of
`fv_controls`, together with the family set `famset(nf)`:

| entry | meaning |
|---|---|
| `ci(1:6)` | as above: criterion, method, preselect, allow rotations, fast MLE, tau mode |
| `ci(7:9)` | as above: tree criterion, tree algorithm, truncation level ($-1$: $d - 1$) |
| `ci(10)` | threads (0: `FORTVINE_THREADS` or half the CPUs) |
| `ci(11)` | `select_trunc`: truncation level chosen by mBICv |
| `ci(12)` | `select_threshold`: threshold chosen by mBICv |

**A vine as arrays.** `d`, `trunc`, `order(d)`, `st(d, d)`, `fam(d, d)`, `rot(d, d)`,
`par(3, d, d)`, `aux(naux, d, d)`, `disc(d)` (1 for a discrete variable), and the relaxed
edges `nrel`, `erel(2, *)`, `wa(d, *)`, `wk(nk, *)`, `pk(3, nk, *)` (each sized
$\max(\mathrm{nrel}, 1)$). `naux` is `fvc_naux()` when an edge is TLL or Bernstein, otherwise 1;
`nk = fvc_nk()` is the number of anchors of a relaxed edge.

```c
int  fvc_naux(void);                         /* length of a TLL / Bernstein estimate */
int  fvc_nk(void);                           /* anchors of a relaxed edge */
void fvc_fam_name(int fam, char buf[48]);    /* blank padded */
int  fvc_families(int nmax, int *codes, int *npars);   /* all 111; npars -1 for tll, bernstein */
void fvc_par_bounds(int fam, double lb[3], double ub[3]);
void fvc_gen_flip(int *fam, double par[3]);  /* generic family of C(u2, u1) */
void fvc_set_abort(int flag);                /* fv_abort: stop a running selection */
void fvc_fitcache(int op, int capacity, int64_t *hits, int64_t *misses, int *used);
                                             /* op 0 disable, 1 enable, 2 clear, 3 stats */

void fvc_bicop_eval(int fam, int rot, const double par[3], int naux, const double *aux,
                    int n, const double *u, int what, double *out);
                    /* what 0 logpdf, 1 pdf, 2 h1, 3 h2, 4 hinv1, 5 hinv2, 6 cdf */
void fvc_bicop_fit(int n, const double *u, int usew, const double *w, int sel, int nf,
                   const int *famset, const int ci[12], const double cr[2], int *fam, int *rot,
                   double par[3], double *ll, int naux, double *aux, double *npe);

void fvc_vine_select(int n, int d, const double *u, const double *um, const int *disc,
                     int usew, const double *w, int fixed, int nf, const int *famset,
                     const int ci[12], const double cr[2], int *trunc, int *order, int *st,
                     int *fam, int *rot, double *par, double *edge_ll, int naux, double *aux,
                     double *ll, double *npe, double *thr, int *info);
void fvc_vine_eval(int d, int trunc, const int *order, const int *st, const int *fam,
                   const int *rot, const double *par, int naux, const double *aux,
                   const int *disc, int nrel, const int *erel, const double *wa,
                   const double *wk, const double *pk, int what, int n, const double *u,
                   const double *um, int seed, int nthreads, double *out);
                   /* what 0 density (out(1, :)), 1 Rosenblatt, 2 inverse, 3 simulate */
void fvc_vine_pdf_grid(/* vine arrays as in fvc_vine_eval */ ..., int nb, int n,
                       const double *ug, const double *ugm, int tgt, int nthreads, double *out);
void fvc_vine_relax(int n, int d, const double *u, int usew, const double *w, int tgt, int mode,
                    int nf, const int *famset, const int ci[12], const double cr[2],
                    int trunc, const int *order, const int *st, int *fam, int *rot, double *par,
                    double *edge_ll, int naux, double *aux, const int *disc, int nmax, int *nrel,
                    int *erel, double *wa, double *wk, double *pk, double *ll, double *npe);
void fvc_vine_rel_params(int d, int nrel, const double *wa, const double *wk, const double *pk,
                         int n, const double *uc, double *prel);   /* prel(3, nrel, n) */
void fvc_vine_structure(int kind, int d, const int *order, int trunc, int *st, int *trunc_out);
void fvc_vine_matrix(int d, const int *order, int trunc, const int *st, int *m);
```

`fvc_bicop_fit` fits the given `fam`/`rot` (`sel = 0`) or selects among `famset` (`sel = 1`),
with weights `w(n)` if `usew = 1`; `npe` is the number of parameters (effective degrees of
freedom for TLL and Bernstein). `fvc_vine_select` takes `u(d, n)` and, for discrete variables,
`um(d, n)` $= F(x^-)$; with `fixed = 1` the structure in `order`, `st`, `trunc` is kept; `info`
is 0, or $-2$ if `naux` is too small for the family set. In `fvc_vine_eval` the Rosenblatt
transform of a vine with discrete variables is randomised with `seed` (0: the conditional
distribution function itself), and simulation is reproducible for `seed > 0`.
`fvc_vine_relax` updates the vine in place (later trees are refitted) and writes the relaxed
edges (`nrel = -1` if more than `nmax`). `tgt`, `mode`, `nb` and the grid layout are those of
[`relax` and `pdf_grid`](vinecop.md).

## Limitations

None: every feature of the typed API can be reached from C through `fortvine_capi`. The 13
kernel routines above do not take weights, discrete margins or TLL/Bernstein estimates (do
not include codes 24 or 25 in their family sets); use the `fvc_*` routines for those. The
generic families (codes $\ge 10000$) work through every routine.

## Example in C

```c
#include <stdio.h>
#include <stdlib.h>

/* prototypes of the fortvine C interface (all arrays column-major, observations last) */
int    fv_version(void);
void   fv_bicop_eval(int fam, int rot, const double par[3], int n, const double *u,
                     int what, double *out, int nthreads);
void   fv_bicop_fit(int n, const double *u, int fam, int rot, const int ci[6],
                    double par[3], double *ll, int nthreads);
void   fv_bicop_select(int n, const double *u, int nf, const int *famset, const int ci[6],
                       double psi0, int *fam, int *rot, double par[3], double *ll, int nthreads);
double fv_kendall_tau(int n, const double *x, const double *y);

int main(void) {
  enum { N = 2000 };
  double par[3] = {2.0, 0.0, 0.0}, *u = malloc(2 * N * sizeof(double));
  double *w = malloc(N * sizeof(double)), *x = malloc(N * sizeof(double)), *y = malloc(N * sizeof(double));
  srand(1);
  for (int i = 0; i < N; ++i) {          /* u[2*i] = first variable, u[2*i+1] = level w */
    u[2 * i] = (rand() + 0.5) / (RAND_MAX + 1.0);
    u[2 * i + 1] = (rand() + 0.5) / (RAND_MAX + 1.0);
  }
  /* simulate a Clayton(2) pair: second column = hinv1(u1, w) */
  fv_bicop_eval(3, 0, par, N, u, 4, w, 1);
  for (int i = 0; i < N; ++i) { u[2 * i + 1] = w[i]; x[i] = u[2 * i]; y[i] = w[i]; }
  printf("fortvine C API version %d, tau = %.3f\n", fv_version(), fv_kendall_tau(N, x, y));

  int ci[6] = {2, 0, 1, 1, 0, 1};       /* BIC, MLE, preselect, rotations, no fast_mle, exact tau */
  double p[3], ll;
  fv_bicop_fit(N, u, 3, 0, ci, p, &ll, 4);
  printf("Clayton fit: theta = %.4f, loglik = %.2f\n", p[0], ll);

  int famset[6] = {1, 2, 3, 4, 5, 6}, fam, rot;
  fv_bicop_select(N, u, 6, famset, ci, 0.9, &fam, &rot, p, &ll, 4);
  printf("selected family %d, rotation %d, par %.4f %.4f\n", fam, rot, p[0], p[1]);
  free(u); free(w); free(x); free(y);
  return 0;
}
```

```bash
gcc -O2 capi.c build/libfortvine.a -lgfortran -lquadmath -fopenmp -o capi
```

Output (MinGW-w64 gcc 16.1, Windows; `rand()` sequences differ between C libraries, so other
platforms give slightly different numbers):

```text
fortvine C API version 1, tau = 0.522
Clayton fit: theta = 2.0579, loglik = 923.65
selected family 3, rotation 0, par 2.0579 0.0000
```

## Example in Python (ctypes)

A pair copula: simulate a Gumbel sample with `hinv1`, select a family, evaluate the density.

```python
import os, ctypes as ct
import numpy as np

# Windows only: make the MinGW runtime DLLs (libgomp-1.dll, ...) findable
os.add_dll_directory("C:/path/to/mingw64/bin")
lib = ct.CDLL("path/to/fortvine/build/libfortvine.dll")       # .so / .dylib elsewhere

dptr = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
iptr = np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS")
lib.fv_bicop_eval.argtypes = [ct.c_int, ct.c_int, dptr, ct.c_int, dptr, ct.c_int, dptr, ct.c_int]
lib.fv_bicop_eval.restype = None
lib.fv_bicop_select.argtypes = [ct.c_int, dptr, ct.c_int, iptr, iptr, ct.c_double,
                                ct.POINTER(ct.c_int), ct.POINTER(ct.c_int), dptr,
                                ct.POINTER(ct.c_double), ct.c_int]
lib.fv_bicop_select.restype = None

rng = np.random.default_rng(1)
n = 2000
u = rng.uniform(size=(n, 2))           # row-major (n, 2) == Fortran u(2, n)
par = np.array([2.0, 0.0, 0.0])
v = np.empty(n)
lib.fv_bicop_eval(4, 0, par, n, u, 4, v, 1)   # Gumbel(2): v = hinv1(u1, w)
u[:, 1] = v

ci = np.array([2, 0, 1, 1, 0, 1], dtype=np.int32)       # BIC, MLE, preselect, rotations, -, exact tau
famset = np.arange(0, 12, dtype=np.int32)
fam, rot, ll = ct.c_int(), ct.c_int(), ct.c_double()
p = np.zeros(3)
lib.fv_bicop_select(n, u, famset.size, famset, ci, 0.9, ct.byref(fam), ct.byref(rot), p, ct.byref(ll), 4)
print("family", fam.value, "rotation", rot.value, "par", p.round(4), "loglik", round(ll.value, 2))

dens = np.empty(n)
lib.fv_bicop_eval(fam.value, rot.value, p, n, u, 1, dens, 1)
print("mean log density", np.log(dens).mean().round(4))
```

Output:

```text
family 4 rotation 0 par [1.9756 0.     0.    ] loglik 728.44
mean log density 0.3642
```

A vine: select on a NumPy `(n, d)` matrix and check that the density routine reproduces the
sum of the edge log-likelihoods.

```python
I = ct.c_int
d64, i32 = dptr, iptr
lib.fv_vine_select.argtypes = [I, I, d64, I, i32, i32, d64, i32, i32, i32, i32, i32, d64, d64, I, ct.POINTER(I)]
lib.fv_vine_pdf.argtypes = [I, I, I, I, d64, i32, I, i32, i32, i32, d64, d64, I]
lib.fv_pseudo_obs.argtypes = [I, I, d64, d64, I]

rng = np.random.default_rng(0)
n, d = 1500, 4
z = rng.standard_normal((n, d)); z[:, 1] += z[:, 0]; z[:, 2] += z[:, 1]; z[:, 3] += 0.5 * z[:, 0]
u = np.empty_like(z); lib.fv_pseudo_obs(n, d, z, u, 4)
famset = np.arange(0, 12, dtype=np.int32)
ci = np.array([2, 0, 1, 1, 0, 1, 0, 0, -1, 0], dtype=np.int32)
cr = np.array([0.0, 0.9])
order = np.zeros(d, np.int32); st_in = np.zeros((d, d), np.int32)
st, fam, rot = (np.zeros((d, d), np.int32) for _ in range(3))
par = np.zeros((d, d, 3)); ll = np.zeros((d, d)); info = I()
lib.fv_vine_select(n, d, u, famset.size, famset, ci, cr, order, st_in, st, fam, rot, par, ll, 4, ct.byref(info))
print("info", info.value, "order", order)
print("fam[tree, edge]:\n", fam)
print("sum of edge logliks", ll.sum().round(3))
dens = np.empty((n, 1))
lib.fv_vine_pdf(n, 1, d, -1, u.reshape(n, 1, d).copy(), order, d - 1, st, fam, rot, par, dens, 4)
print("sum log pdf      ", np.log(dens).sum().round(3))
```

Output:

```text
info 0 order [3 2 1 4]
fam[tree, edge]:
 [[1 1 1 0]
 [0 0 0 0]
 [0 0 0 0]
 [0 0 0 0]]
sum of edge logliks 1465.024
sum log pdf       1465.024
```

The data are a Gaussian Markov tree (1-2-3 and 1-4), so the first tree has three Gaussian
edges (family 1) and the higher trees are independent; row `t` of `fam` lists the edges of tree
`t`.
