#!/usr/bin/env bash
# Copy the Fortran sources of the library into the R package, which must be self-contained
# because remotes::install_github(subdir = "r") builds only that folder. The copies go to
# r/src/fortvine/ with the same content; the .f90 files get the extension .F90 so that every
# Fortran compiler runs the C preprocessor on them (the .inc files are #included) without a
# compiler-specific flag. R CMD check accepts any file in a subdirectory of src/.
# CI runs this script and then `git diff --exit-code r/src` to check the copies are current.
set -euo pipefail
cd "$(dirname "$0")/.."
dst=r/src/fortvine
mkdir -p "$dst"
rm -f "$dst"/*.F90 "$dst"/*.inc
for f in src/*.f90; do
  cp "$f" "$dst/$(basename "$f" .f90).F90"
done
cp src/*.inc "$dst/"
echo "copied $(ls src/*.f90 src/*.inc | wc -l) files from src/ to $dst/"
