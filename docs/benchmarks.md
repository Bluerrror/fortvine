# Benchmarks

This page compares fortvine with two other vine-copula libraries on identical data:

* **pyvinecopulib 1.0.1** — the Python bindings of [vinecopulib](https://github.com/vinecopulib/vinecopulib)
  (C++), the reference implementation that fortvine's conventions follow;
* **torchvine 0.2.2** — a pure-PyTorch reimplementation of the pyvinecopulib API (CPU build of
  PyTorch 2.15).

All numbers on this page come from `bench/run_bench.py`; the raw results are in
`bench/results/results.json`. You can reproduce them on your own machine (see the end of the page).

## Set-up

| | |
|---|---|
| machine | Intel Core (Raptor Lake, 32 logical CPUs), Windows 11 |
| fortvine | gfortran 16.1 (WinLibs), `-O3 -march=native -funroll-loops -fopenmp` |
| Python side | Python 3.14.0, NumPy 2.5.5, pyvinecopulib 1.0.1 (PyPI wheel), torch 2.15.0+cpu, torchvine 0.2.2 |
| family set | the 12 families all libraries share: independence, Gaussian, Student t, Clayton, Gumbel, Frank, Joe, BB1, BB6, BB7, BB8, Tawn (torchvine has no Tawn) — with rotations |
| criterion | BIC everywhere |
| threads | "one thread" runs pin every library to one thread (`FORTVINE_THREADS=1`, `num_threads=1`, `torch.set_num_threads(1)`) |

**Data.** Pair-copula data are simulated with pyvinecopulib from each family at Kendall's
$\tau \approx 0.5$ ($n = 5000$). Vine data come from a heavy-tailed, skewed two-factor model
turned into pseudo-observations ($n = 2000$, $d = 5, 10, 20$), so that no library's own model
class generated them.

Timings are wall-clock means after a warm-up call. Python timings include the (small) call
overhead of the bindings.

## Results

![Vine selection time against dimension, one and eight threads](assets/bench_vine.png)

![Pair-copula fit and selection times per family](assets/bench_pair.png)

Each dot or line is one library; both axes of time are logarithmic. The full tables follow.

--8<-- "bench/results/tables.md"

## What the numbers say

**Same answers.** Where the two libraries fit the same model, fortvine and pyvinecopulib reach the
same maximum: the log-likelihoods of every pair fit, every family selection and every selected
vine agree to the digits shown, and both select the same family for every simulated pair. This
is by design — fortvine reproduces vinecopulib's rotations, structure arrays, candidate
pre-selection and Dissmann algorithm — and it makes the speed comparison a like-for-like one.
torchvine's maxima are slightly lower for several one-parameter families (it starts from Kendall's
$\tau$ inversion), and its selection picks a different family for the Clayton and Joe samples (it has no Tawn family, so the Tawn sample cannot be matched).

**Pair-copula evaluation.** The density is 1.05–3.7× faster than pyvinecopulib's for every family (Student t is the closest).
The largest gaps are in the inverse h-function: fortvine inverts h with closed forms where they
exist and with a safeguarded Newton step that uses the density as the exact derivative elsewhere,
so BB1, BB6, BB7, BB8 and Tawn inversions run in about 0.6–1.2 µs per point, against 10–17 µs.

**Fitting.** A maximum-likelihood fit of a known family is 3–14× faster than pyvinecopulib's
(1–12 ms against 6–243 ms at $n = 5000$). Exact parameter scores from dual numbers let
Newton/BFGS converge in a few steps, and large samples are fitted in two stages (a subsample
first, then a short polish on all data).

**Selection.** Choosing among the 12 families and their rotations takes 50–110 ms in fortvine
and 1–2 s in pyvinecopulib — 9–41× faster. Candidate pre-selection, racing (a cheap bound
decides whether a candidate can still win before it is polished) and the shared closed-form
kernels account for this.

**Vines.** Selecting a regular vine is 11–17× faster on one thread and 15–23× faster on eight,
with identical log-likelihoods. Density evaluation is 3–9× faster and simulation 14–40× faster.
On eight threads a 20-dimensional vine is selected in 1.2 s.

| d | threads | fortvine select | pyvinecopulib select | speed-up |
|---:|---:|---:|---:|---:|
| 5 | 1 | 0.49 s | 6.94 s | 14× |
| 10 | 1 | 2.04 s | 34.3 s | 17× |
| 20 | 1 | 14.1 s | 158.9 s | 11× |
| 5 | 8 | 0.08 s | 1.26 s | 15× |
| 10 | 8 | 0.33 s | 6.70 s | 20× |
| 20 | 8 | 1.19 s | 27.6 s | 23× |

!!! note "Read these numbers with care"
    Benchmarks measure one machine, one compiler and one set of data. pyvinecopulib was run from
    its PyPI wheel on Windows; a native Linux build of vinecopulib may be faster. torchvine is
    designed for GPUs and autograd, and its CPU timings do not show what it is built for. The
    script is in the repository so you can check the comparison on your own hardware.

## Beyond speed

| | fortvine | pyvinecopulib | torchvine |
|---|---|---|---|
| language | Fortran 2008 + OpenMP | C++17 (Python bindings) | Python / PyTorch |
| dependencies | compiler only | Eigen, Boost (header-only), pybind/nanobind | PyTorch |
| parametric families | 108 | 11 | 10 |
| asymmetric families | 32 (Tawn, asymmetric EV, Archimax, Khoudraji) | 1 (Tawn) | 0 |
| nonparametric families | TLL, Bernstein | TLL | TLL |
| family rotations | yes | yes | yes |
| exact parameter gradients | dual numbers, every family | analytic (1.0) | autograd |
| discrete / mixed margins | yes | yes | yes |
| observation weights | yes | yes | yes |
| truncation / threshold by mBICv | yes | yes | yes |
| non-simplified (covariate-dependent) edges | yes (`relax`) | — | — |
| GPU | no | no | yes |
| C ABI | yes | C++ API | — |

## Reproducing

```bash
bash build.sh bench                       # builds build/fv_bench
python -m venv .venv && . .venv/bin/activate
pip install numpy pyvinecopulib torch torchvine pytorch-minimize
python bench/run_bench.py                 # about 10 minutes; --quick for a short run
```

The script writes `bench/results/results.json` and `bench/results/tables.md` (included above).
