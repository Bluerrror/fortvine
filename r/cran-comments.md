## Submission

First submission of fortvine 1.0.0: R bindings to the fortvine Fortran vine-copula library.
All computations are done by the bundled Fortran code (src/fortvine/, a copy of the
library sources) through its C interface.

## Test environments

* local: Windows 11 x64, R 4.3.2, Rtools43 (gcc / gfortran 12.3.0)
* local, CRAN-like flags: the same with `-Wall -pedantic` (C) and `-Wall -pedantic -std=f2008`
  (Fortran): no warnings from the package's own C and Fortran glue
* GitHub Actions: ubuntu-latest, R release (.github/workflows/r.yml)

To add before submitting: win-builder (R-devel, `devtools::check_win_devel()`) and macOS
(`devtools::check_mac_release()`), or R-hub (`rhub::rhub_check()`).

## R CMD check results

`R CMD check --as-cran --no-manual` on the built tarball (Windows, R 4.3.2):

0 errors | 1 warning | 2 notes

* WARNING: "Non-portable flags in variable 'PKG_FCFLAGS': -ffree-line-length-none".
  Three statements of the library sources exceed 132 columns. They are being wrapped
  upstream; the flag is then removed and this warning disappears. **Do not submit before
  that.**
* NOTE: "New submission".
* NOTE: "unable to verify current time" (no time server reachable from the check machine;
  not a package issue).
