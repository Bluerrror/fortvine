# Calling fortvine from Python and C

**Programs:** `example/python/fortvine_ctypes.py`, `example/c/bicop.c` &nbsp;·&nbsp;
**Topics:** the shared library, the `bind(C)` interface of `src/fortvine.f90`, array layout,
`ctypes` + numpy, linking a C program.

Every computational kernel in `src/fortvine.f90` has a C-interoperable entry point
(`bind(C, name="fv_...")`). `bash build.sh` links them into the shared library
`build/libfortvine.so` (Linux), `libfortvine.dylib` (macOS) or `libfortvine.dll`
(Windows), which any language with a C foreign-function interface can load. No wrapper
code has to be compiled.

## The conventions

* **Scalars** (`integer(c_int)`, `real(c_double)`) declared with `value` are passed by value;
  output scalars (`ll`, `fam`, `rot`, `info`) are passed by reference.
* **Arrays are passed by pointer** and must be contiguous `double` / `int32` buffers.
* **Observations are the last Fortran dimension.** `u(2, n)` for a pair, `u(d, n)` for a
  vine: one observation's values are adjacent in memory. In C or numpy that is a row-major
  array of shape `(n, 2)` or `(n, d)`, so ordinary C-ordered numpy arrays are passed as they
  are.
* **Vine structure arrays** `st`, `fam`, `rot` are Fortran `(edge, tree)` arrays of size
  `d x d`, which in row-major languages are `[tree][edge]`; parameters `par(3, edge, tree)`
  are `[tree][edge][k]`, i.e. numpy shape `(d, d, 3)`. Indices of trees and edges start at 0
  on the C side; `order` holds 1-based variable labels.
* **`nthreads`** is an explicit OpenMP thread count (use 1 or more).

The C prototypes of the routines used below (see the [C interface](../manual/c-api.md) for
all of them):

```c
int    fv_version(void);
void   fv_bicop_eval(int fam, int rot, const double *par, int n, const double *u,
                     int what, double *out, int nthreads);   /* what: 0 logpdf 1 pdf 2 h1 3 h2 4 hinv1 5 hinv2 */
void   fv_bicop_fit(int n, const double *u, int fam, int rot, const int *ci,
                    double *par, double *ll, int nthreads);
void   fv_bicop_select(int n, const double *u, int nf, const int *famset, const int *ci,
                       double psi0, int *fam, int *rot, double *par, double *ll, int nthreads);
double fv_kendall_tau(int n, const double *x, const double *y);
void   fv_tau_matrix(int n, int d, const double *x, double *tau, int nthreads);
void   fv_pseudo_obs(int n, int d, const double *x, double *u, int nthreads);
void   fv_vine_select(int n, int d, const double *u, int nf, const int *famset,
                      const int *ci, const double *cr, int *order, const int *st_in,
                      int *st, int *fam, int *rot, double *par, double *edge_ll,
                      int nthreads, int *info);
void   fv_vine_pdf(int nsamp, int bins, int d, int tgt, const double *u, const int *order,
                   int trunc, const int *st, const int *fam, const int *rot,
                   const double *par, double *out, int nthreads);
void   fv_vine_rosenblatt(int n, int d, const double *u, const int *order, int trunc,
                          const int *st, const int *fam, const int *rot,
                          const double *par, double *out, int nthreads);
void   fv_vine_inverse_rosenblatt(int n, int d, const double *u, const int *order, int trunc,
                                  const int *st, const int *fam, const int *rot,
                                  const double *par, double *out, int nthreads);
```

The integer control vector `ci` mirrors `fv_controls`:

| index (0-based) | meaning | default |
|---:|---|---:|
| 0 | selection criterion: 0 loglik, 1 AIC, 2 BIC, 3 mBIC, 4 mBICv | 2 |
| 1 | method: 0 MLE, 1 inversion of Kendall's tau | 0 |
| 2 | preselect families (1 = yes) | 1 |
| 3 | allow rotations | 1 |
| 4 | fast one-parameter fits (tau inversion only) | 0 |
| 5 | tau estimator: 0 auto, 1 exact, 2 approximate | 1 |
| 6 | (vine only) tree criterion: 0 tau, 1 rho, 2 Joe, 3 Hoeffding, ... | 0 |
| 7 | (vine only) spanning-tree algorithm: 0 Prim, 1 Kruskal | 0 |
| 8 | (vine only) truncation level, -1 = none | -1 |
| 9 | (vine only) structure given in `order`/`st_in` (1) or to be selected (0) | 0 |

and `cr = {threshold, psi0}` (default `{0, 0.9}`).

The C interface covers unweighted, continuous data with the parametric families. Weights,
discrete variables, the nonparametric families and the mBICv-based truncation search are
reached through the Fortran API (`fortvine_api`).

## Python with ctypes

```python
--8<-- "example/python/fortvine_ctypes.py"
```

Run it from the repository root after `bash build.sh` (it needs numpy):

```bash
python example/python/fortvine_ctypes.py
```

Output:

```text
fortvine C API version 1
Clayton(2) density on the diagonal: [5.3702 1.481  2.1578]
empirical Kendall tau: 0.485
fitted theta: 1.9290   log-likelihood: 805.75
selected: clayton rotation 0 par [1.929 0.    0.   ]

vine order: [1 2 3 4]
  tree 1: gaussian(0.582), gaussian(0.595), gaussian(0.574)
  tree 2: indep(0.000), indep(0.000)
  tree 3: indep(0.000)
log-likelihood: 1865.83  (sum of edge log-likelihoods 1865.83)
Kendall tau, data      : [1.    0.395 0.221 0.125]
Kendall tau, simulated : [1.    0.404 0.232 0.152]
```

### What the script does

* **Loading.** `ctypes.CDLL` opens the library. On Windows the DLL imports the MinGW
  runtime DLLs (`libgomp-1.dll`, `libquadmath-0.dll`, `libwinpthread-1.dll`), so the script
  adds the directory of `gfortran` on `PATH` with `os.add_dll_directory` first.
* **Signatures.** Declaring `argtypes` with `numpy.ctypeslib.ndpointer(..., flags='C_CONTIGUOUS')`
  makes ctypes check every array's type and layout, so a wrongly shaped or Fortran-ordered
  array raises an error instead of producing garbage.
* **One pair.** `fv_bicop_eval` with `what = 1` gives the same Clayton densities as the
  Fortran [pair copula example](pair-copula.md); `what = 4` (inverse $h_1$) simulates;
  `fv_bicop_fit` and `fv_bicop_select` fit and select with the default controls.
* **A vine.** Four-dimensional data with the correlation matrix $0.6^{|i-j|}$ (a Gaussian
  AR(1) process) are converted to pseudo-observations by `fv_pseudo_obs`, and
  `fv_vine_select` chooses the structure. The answer is the one theory predicts: the
  Gaussian AR(1) process is Markov, so the D-vine 1-2-3-4 with Gaussian pairs of
  correlation 0.6 in tree 1 and independence in all higher trees is exact. The estimates
  0.57-0.60 match, and BIC sets the three higher edges to independence.
* **Density and simulation.** `fv_vine_pdf` with one bin per sample and `tgt = -1` returns
  the ordinary vine density (the `bins`/`tgt` arguments serve the fast prediction-grid
  evaluation of [`pdf_grid`](../manual/vinecop.md)); its log-sum equals the sum of the edge
  log-likelihoods returned by the selection. `fv_vine_inverse_rosenblatt` maps uniform noise
  to a new sample whose Kendall taus reproduce those of the data.
* **Truncation.** `fv_vine_select` does not return a truncation level; with no truncation
  requested (`ci[8] = -1`) it is `d - 1`, the value passed to the other routines.

## C

```c
--8<-- "example/c/bicop.c"
```

Compile against the shared library and run (from the repository root):

=== "Linux / macOS"

    ```bash
    gcc -O2 example/c/bicop.c -Lbuild -lfortvine -Wl,-rpath,"$PWD/build" -o build/bicop_c
    ./build/bicop_c
    ```

=== "Windows (MinGW-w64)"

    ```bash
    gcc -O2 example/c/bicop.c build/libfortvine.dll -o build/bicop_c.exe
    PATH="$PWD/build:$PATH" ./build/bicop_c.exe
    ```

    MinGW's linker can link directly against a DLL. At run time the DLL must be found on
    `PATH` (or next to the executable), together with the MinGW runtime DLLs it imports.

Output:

```text
fortvine C API version 1
Clayton(2) density on the diagonal: 5.3702 1.4810 2.1578
empirical Kendall tau: 0.504
fitted theta: 2.0611   log-likelihood: 879.27
selected: clayton rotation 0 par 2.061 0.000 0.000
```

The program declares the prototypes itself, so no header file is needed. The pairs are
stored as a row-major `n x 2` array (`u[2*i]`, `u[2*i+1]`), which is the Fortran `u(2, n)`
the kernel expects. A small xorshift generator makes the run reproducible.

!!! warning "Linking statically from C"
    `build/libfortvine.a` can also be linked into a C program, but then the Fortran and
    OpenMP runtimes must be linked too. The simplest way is to let `gfortran` drive the
    link: `gcc -c prog.c && gfortran -fopenmp prog.o build/libfortvine.a -o prog` (add
    `-static` on Windows for a self-contained executable). With `gcc` as the driver, use
    `gcc -fopenmp prog.o build/libfortvine.a -lgfortran -lquadmath -lm -o prog`.
