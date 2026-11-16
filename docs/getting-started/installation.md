# Installation

fortvine is plain modern Fortran with OpenMP directives. It has no dependencies beyond the
compiler: no BLAS, no LAPACK, no CMake, no Python. You can compile it with the bundled shell
script, with the Fortran Package Manager, or by hand with three compiler calls.

## Requirements

| What | Version | Notes |
|---|---|---|
| gfortran | 10 or newer | Tested up to GCC 16. Other Fortran 2008 compilers with a C preprocessor should work but are not tested. |
| OpenMP | the one that ships with gfortran (libgomp) | Optional. Without it the library runs on one thread and gives the same results. |
| A POSIX shell | bash | Only for `build.sh`. On Windows, Git Bash or the MSYS2 shell. |

Supported platforms:

=== "Linux"

    Any distribution's gfortran works:

    ```bash
    sudo apt-get install gfortran        # Debian / Ubuntu
    sudo dnf install gcc-gfortran        # Fedora / RHEL
    ```

=== "macOS"

    Homebrew's GCC contains gfortran and libgomp:

    ```bash
    brew install gcc
    export FC=gfortran-14                # the versioned name Homebrew installs
    ```

    Apple's `clang` has no Fortran compiler; the shared library is built as
    `libfortvine.dylib`.

=== "Windows"

    Use a MinGW-w64 gfortran, for example one of these:

    * **WinLibs** (standalone GCC + MinGW-w64, unpack and add `mingw64\bin` to `PATH`),
    * **MSYS2** (`pacman -S mingw-w64-ucrt-x86_64-gcc-fortran`),
    * any other MinGW-w64 distribution of GCC 10+.

    Run `build.sh` from Git Bash or the MSYS2 shell. The script detects MinGW / MSYS /
    Cygwin and links the example and test programs statically (see
    [Troubleshooting](#troubleshooting)).

## Building with `build.sh`

From the repository root:

```bash
bash build.sh
```

This compiles the three source files (`src/fortvine_ad.f90`, `src/fortvine.f90`,
`src/fortvine_api.f90`; the `.inc` files are pulled in by the preprocessor) and writes, in
`build/`:

| File | Purpose |
|---|---|
| `libfortvine.a` | static library for Fortran programs |
| `libfortvine.so` / `libfortvine.dylib` / `libfortvine.dll` | shared library for C, C++, Python, Julia, R, ... |
| `fortvine_api.mod`, `fortvine.mod`, `fortvine_ad*.mod` | module files the compiler needs for `use fortvine_api` |

Extra targets can be appended in any combination:

```bash
bash build.sh test                 # build and run test/test_fortvine and test/test_genfam
bash build.sh examples             # build every example/*.f90 into build/<name>
bash build.sh bench                # build bench/fv_bench (driven by bench/run_bench.py)
bash build.sh test examples bench  # all of them
```

`bash build.sh test` must end with the two lines `ALL TESTS PASSED` and
`ALL GENERIC FAMILY TESTS PASSED`; see [Validation](../theory/validation.md) for what the
suites check. A full `bash build.sh test examples` takes about half a minute on a laptop.

### Compiler and flags

Two environment variables control the build:

| Variable | Default | Meaning |
|---|---|---|
| `FC` | `gfortran` | the Fortran compiler |
| `FFLAGS` | `-O3 -march=native -funroll-loops -fno-math-errno -ffree-line-length-none -fopenmp` | compile flags |

`build.sh` always adds `-cpp -I ../src -fPIC` itself, so `FFLAGS` only has to carry
optimisation and OpenMP choices. `-ffree-line-length-none` is required: a few source lines
are longer than the standard 132 columns.

```bash
FC=gfortran-13 bash build.sh test
FFLAGS="-O2 -ffree-line-length-none -fopenmp" bash build.sh test examples   # what CI uses
FFLAGS="-O2 -ffree-line-length-none" bash build.sh                           # serial build
```

!!! tip "Portable binaries"
    `-march=native` tunes the code for the CPU that compiles it. A library built this way
    can crash with *illegal instruction* on an older machine. If you distribute binaries,
    replace it by a baseline such as `-march=x86-64-v2` (or `-march=x86-64-v3` for AVX2
    machines), or leave it out:

    ```bash
    FFLAGS="-O3 -march=x86-64-v2 -funroll-loops -fno-math-errno -ffree-line-length-none -fopenmp" bash build.sh
    ```

## Building with fpm

The repository is also a [Fortran Package Manager](https://fpm.fortran-lang.org) package.
`fpm.toml` declares the library in `src/`, turns on the C preprocessor
(`[preprocess.cpp]`, needed for the `.inc` family files), disables implicit typing, and
discovers the tests in `test/` and the examples in `example/` automatically.

fpm does not know about the long source lines or about OpenMP, so pass the flags explicitly:

```bash
fpm build --profile release --flag "-ffree-line-length-none -fopenmp" --link-flag "-fopenmp"
fpm test  --profile release --flag "-ffree-line-length-none -fopenmp" --link-flag "-fopenmp"
fpm run --example ex01_pair_copula --profile release --flag "-ffree-line-length-none -fopenmp" --link-flag "-fopenmp"
```

To use fortvine as a dependency of your own fpm project:

```toml
[dependencies]
fortvine = { git = "https://github.com/Bluerrror/fortvine" }
```

!!! note
    The fpm route has not been exercised by the continuous integration (which runs
    `build.sh`); if your fpm version does not preprocess `.f90` files, add `-cpp` to `--flag`.

## Linking your own program

A Fortran program needs the module files at compile time and the static library at link
time:

```fortran title="hello_vine.f90"
program hello_vine
  use fortvine_api
  implicit none
  type(fv_bicop) :: cop
  real(kind(1d0)) :: c(1)
  cop%family = FV_FRANK; cop%par(1) = 5.0d0
  call cop%pdf([0.3d0], [0.4d0], c)
  print *, 'Frank(5) density at (0.3, 0.4):', c(1)
end program hello_vine
```

=== "Linux / macOS"

    ```bash
    gfortran -O3 -fopenmp -I path/to/fortvine/build hello_vine.f90 \
             path/to/fortvine/build/libfortvine.a -o hello_vine
    ```

=== "Windows (MinGW-w64)"

    ```bash
    gfortran -O3 -fopenmp -I path/to/fortvine/build hello_vine.f90 \
             path/to/fortvine/build/libfortvine.a -static -o hello_vine.exe
    ```

Rules of thumb:

* Compile your program with the **same compiler** (and preferably the same version) as the
  library: `.mod` files are not portable between compilers or major gfortran versions.
* If the library was built with OpenMP, **link with `-fopenmp`** too, or the linker reports
  undefined `GOMP_*` / `omp_*` symbols. A serial library links without it.
* The library must come **after** your source or object files on the command line.

## The shared library (C, C++, Python, ...)

`build.sh` also produces `build/libfortvine.so` (Linux), `.dylib` (macOS) or `.dll`
(Windows). Every kernel of `src/fortvine.f90` with a `bind(C, name="fv_...")` interface is
exported under that name, e.g. `fv_bicop_eval`, `fv_bicop_fit`, `fv_bicop_select`,
`fv_vine_select`, `fv_vine_pdf`, `fv_vine_rosenblatt`, `fv_vine_inverse_rosenblatt`,
`fv_kendall_tau`, `fv_tau_matrix`, `fv_pseudo_obs`.

* C and C++: declare the prototypes and link against the library, see
  [Calling fortvine from Python and C](../examples/interop.md) and the
  [C interface](../manual/c-api.md).
* Python: load it with `ctypes` (or `cffi`) and pass `numpy` arrays; a complete script is in
  `example/python/fortvine_ctypes.py`.

On Windows the DLL is linked with `-static-libgfortran -static-libgcc`, so it does not need
`libgfortran`, but it still imports `libgomp-1.dll`, `libquadmath-0.dll` and
`libwinpthread-1.dll`. Those ship with every MinGW-w64 gfortran: keep its `bin` directory on
`PATH` (or call `os.add_dll_directory` in Python), or copy the three DLLs next to
`libfortvine.dll`.

## Troubleshooting

??? question "`Error: Line truncated` or `Unterminated character constant`"
    `-ffree-line-length-none` is missing from `FFLAGS`. If you set `FFLAGS` yourself, keep
    that flag.

??? question "`undefined reference to GOMP_parallel` / `omp_get_num_procs`"
    The library was compiled with `-fopenmp` but your program was linked without it. Add
    `-fopenmp` to the link command, or rebuild the library serially
    (`FFLAGS="-O3 -ffree-line-length-none" bash build.sh`).

??? question "OpenMP is not available on my compiler"
    Build without `-fopenmp` as above. Every routine still works, on one thread. Densities,
    h-functions and log-likelihood sums do not depend on the number of threads (sums are
    accumulated over a fixed set of blocks in a fixed order). Family *selection* races
    candidate fits against each other when several threads are available, so a near-tie
    between two candidates can occasionally be decided differently; see
    [Performance and threads](../manual/performance.md).

??? question "*Illegal instruction* on another machine"
    The library was built with `-march=native` on a newer CPU. Rebuild with a baseline
    architecture (see the tip under [Compiler and flags](#compiler-and-flags)).

??? question "Windows: the program starts and exits silently, or complains about missing DLLs"
    A dynamically linked executable needs `libgfortran-5.dll`, `libgomp-1.dll`,
    `libquadmath-0.dll`, `libwinpthread-1.dll` and `libgcc_s_seh-1.dll` at run time.
    Either link with `-static` (what `build.sh` does for the tests and examples) or put the
    MinGW-w64 `bin` directory on `PATH`.

??? question "Windows: a program with big arrays exits silently (exit code 127 in Git Bash)"
    With `-fopenmp`, gfortran also implies `-frecursive`, which places local arrays, even
    those of the main program, on the stack, and the default Windows stack is 2 MB.
    Make large arrays `allocatable` (as `example/ex04_tails.f90` does) or enlarge the stack
    at link time with `-Wl,--stack,268435456`.

??? question "`Fatal Error: Cannot read module file 'fortvine_api.mod' ... created by a different version`"
    Rebuild fortvine with the compiler you use for your program.
