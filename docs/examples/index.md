# Examples

Every example is a complete, short program in the repository's `example/` directory. The
first four also come as Python and R scripts (`example/python/`, `example/r/`) that use the
[Python](../python.md) and [R](../r.md) packages; each page shows the three versions in
language tabs. The Fortran programs are compiled by

```bash
bash build.sh examples        # -> build/ex01_pair_copula, build/ex02_vine, ...
```

and run from anywhere (`./build/ex02_vine`); the scripts with `python example/python/ex02_vine.py`
and `Rscript example/r/ex02_vine.R`. The outputs shown on these pages are copied from real runs
on Windows 11 with gfortran 16.1 (WinLibs), `OMP_NUM_THREADS=4`, Python 3.14 and R 4.3.2.

| Page | Program | What it shows |
|---|---|---|
| [Pair copula basics](pair-copula.md) | `example/ex01_pair_copula.f90`, `python/ex01_pair_copula.py`, `r/ex01_pair_copula.R` | density, simulation by inverse h-function, MLE of a given family, family selection by BIC |
| [Vine copula workflow](vine.md) | `example/ex02_vine.f90`, `python/ex02_vine.py`, `r/ex02_vine.R` | simulate from a known D-vine, select a regular vine, density, Rosenblatt transform, new samples |
| [Generic families](generic.md) | `example/ex03_generic_families.f90`, `python/ex03_generic_families.py`, `r/ex03_generic_families.R` | the 85 generic family codes, an Archimax copula, an asymmetric Khoudraji copula, selection among all 109 parametric families |
| [Tail dependence and asymmetry](tails.md) | `example/ex04_tails.f90`, `python/ex04_tails.py`, `r/ex04_tails.R` | tail coefficients from family formulas, from the copula CDF at finite levels and from samples; what BIC selects |
| [Weighted and discrete data](weights-discrete.md) | `example/ex05_weights_discrete.f90` | observation weights (time decay, duplicated rows) and a vine with a Poisson count variable |
| [Calling fortvine from Python and C](interop.md) | `example/python/fortvine_ctypes.py`, `example/c/bicop.c` | the `bind(C)` API through `ctypes` + numpy, and from a C program |

!!! info "Random numbers and reproducibility"
    Every program fixes its seeds (`random_seed(put=...)` or the `seed=` argument of
    `fv_vinecop%simulate` in Fortran, `np.random.default_rng(seed)` and `seed=` in Python,
    `set.seed()` in R), so it prints the same numbers every time with the same compiler and
    package. The three languages draw their uniforms from different generators, though, so
    the numbers that depend on a simulated sample differ between the Fortran, Python and R
    versions of an example; deterministic quantities (densities, h-functions, exact tail
    coefficients) agree to the printed digits.

The remaining programs in `example/` (`bench.f90`, `famspeed.f90`, `fitspeed.f90`,
`selspeed.f90`) are timing drivers used for the [benchmarks](../benchmarks.md).
