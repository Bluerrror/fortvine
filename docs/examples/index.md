# Examples

Every example is a complete, short program in the repository's `example/` directory. The
Fortran ones are compiled by

```bash
bash build.sh examples        # -> build/ex01_pair_copula, build/ex02_vine, ...
```

and run from anywhere (`./build/ex02_vine`). The outputs shown on these pages are copied from
real runs on Windows 11 with gfortran 16.1 (WinLibs), `OMP_NUM_THREADS=4`.

| Page | Program | What it shows |
|---|---|---|
| [Pair copula basics](pair-copula.md) | `example/ex01_pair_copula.f90` | density, simulation by inverse h-function, MLE of a given family, family selection by BIC |
| [Vine copula workflow](vine.md) | `example/ex02_vine.f90` | simulate from a known D-vine, select a regular vine, density, Rosenblatt transform, new samples |
| [Generic families](generic.md) | `example/ex03_generic_families.f90` | the 85 generic family codes, an Archimax copula, an asymmetric Khoudraji copula, selection among all 112 families |
| [Tail dependence and asymmetry](tails.md) | `example/ex04_tails.f90` | tail coefficients from family formulas, from the copula CDF at finite levels and from samples; what BIC selects |
| [Weighted and discrete data](weights-discrete.md) | `example/ex05_weights_discrete.f90` | observation weights (time decay, duplicated rows) and a vine with a Poisson count variable |
| [Calling fortvine from Python and C](interop.md) | `example/python/fortvine_ctypes.py`, `example/c/bicop.c` | the `bind(C)` API through `ctypes` + numpy, and from a C program |

!!! info "Random numbers and reproducibility"
    `ex01` and `ex03` seed the Fortran generator with `call random_seed()`, which in
    gfortran draws a fresh seed from the operating system: their numbers change slightly
    from run to run. `ex02`, `ex04`, `ex05` and the Python and C programs use fixed seeds
    and print the same numbers every time with the same compiler. `fv_vinecop%simulate`
    takes a `seed=` argument for that purpose.

The remaining programs in `example/` (`bench.f90`, `famspeed.f90`, `fitspeed.f90`,
`selspeed.f90`) are timing drivers used for the [benchmarks](../benchmarks.md).
