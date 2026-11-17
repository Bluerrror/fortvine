## Submission

First submission of fortvine 1.0.0: R bindings to the fortvine Fortran vine-copula library.
All computations are done by the bundled Fortran code (src/fortvine/, a copy of the
library sources) through its C interface (fortvine_capi).

## Test environments

* local: Windows 11 x64, R 4.3.2, Rtools43 (gcc / gfortran 12.3.0)
* local, CRAN-like flags: `-Wall -pedantic` (C) and `-Wall -pedantic -std=f2008` (Fortran):
  no warnings from the package's own C and Fortran glue
* GitHub Actions: ubuntu-latest, R release, `R CMD check --as-cran` (.github/workflows/r.yml)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Notes for the reviewers

* The bundled library writes nothing to the console and never stops the process when built
  for R (`-DFORTVINE_NO_IO` in Makevars); errors are returned to R and raised there.
* OpenMP is used through `$(SHLIB_OPENMP_FFLAGS)` / `$(SHLIB_OPENMP_CFLAGS)`; examples and
  tests use at most 2 threads (the default is 1). The Fortran sources live in src/fortvine/
  and are compiled by explicit rules in Makevars, whose flags are in the variable
  `FV_FCFLAGS`; R applies `PKG_FCFLAGS` only to free-form sources at the top of src/.
  `$(FLIBS)` is in `PKG_LIBS` for the same reason (the C compiler links the package).
