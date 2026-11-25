#!/usr/bin/env bash
# Build fortvine with gfortran (no fpm / CMake needed).
#   bash build.sh            -> build/libfortvine.a, build/libfortvine.so|.dll, build/*.mod
#   bash build.sh test       -> also builds and runs test/test_fortvine and test/test_genfam
#   bash build.sh examples   -> also builds the programs in example/ into build/
#   bash build.sh bench      -> also builds bench/fv_bench (driven by bench/run_bench.py)
# Environment: FC (default gfortran), FFLAGS (default: -O3 -march=native ... -fopenmp)
set -euo pipefail
cd "$(dirname "$0")"
FC=${FC:-gfortran}
FFLAGS=${FFLAGS:-"-O3 -march=native -funroll-loops -fno-math-errno -ffree-line-length-none -fopenmp"}
mkdir -p build
cd build
OBJS=""
for f in fv_util fortvine_ad fortvine fortvine_api fortvine_capi; do
  $FC $FFLAGS -cpp -I ../src -fPIC -J . -c ../src/$f.f90 -o $f.o
  OBJS="$OBJS $f.o"
done
ar rcs libfortvine.a $OBJS
case "$(uname -s)" in
  MINGW*|MSYS*|CYGWIN*) $FC -shared -static-libgfortran -static-libgcc -fopenmp $OBJS -o libfortvine.dll; STATIC="-static" ;;
  Darwin*) $FC -shared -fopenmp $OBJS -o libfortvine.dylib; STATIC="" ;;
  *) $FC -shared -fopenmp $OBJS -o libfortvine.so; STATIC="" ;;
esac
for target in "$@"; do
  case $target in
    test)
      for t in test_fortvine test_genfam; do
        $FC $FFLAGS -I . ../test/$t.f90 libfortvine.a $STATIC -o $t && ./$t
      done ;;
    bench)
      $FC $FFLAGS -I . ../bench/fv_bench.f90 libfortvine.a $STATIC -o fv_bench ;;
    examples)
      for e in ../example/*.f90; do
        b=$(basename "$e" .f90); $FC $FFLAGS -I . -I ../example "$e" libfortvine.a $STATIC -o "$b"
      done ;;
  esac
done
echo "fortvine built in $(pwd)"
