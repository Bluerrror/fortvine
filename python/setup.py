"""Build step: compile the Fortran sources of the repository (../src) plus the small
bind(C) shim fortvine_py.f90 into fortvine/_lib/libfortvine.{so,dylib,dll}.

An sdist carries a copy of ../src in fsrc/ and of ../LICENSE (made by the sdist command),
so it builds outside the repository.  Environment: FC (default gfortran), FFLAGS (default: portable
-O3, no -march=native), FORTVINE_SRC (overrides the source directory).
"""
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from setuptools import setup
from setuptools.command.bdist_wheel import bdist_wheel
from setuptools.command.build_py import build_py
from setuptools.command.sdist import sdist
from setuptools.dist import Distribution

HERE = Path(__file__).resolve().parent
SOURCES = ["fortvine_ad.f90", "fortvine.f90", "fortvine_api.f90"]
LIBNAME = {"win32": "libfortvine.dll", "darwin": "libfortvine.dylib"}.get(sys.platform, "libfortvine.so")
FFLAGS = "-O3 -funroll-loops -fno-math-errno -ffree-line-length-none -fopenmp"


def fortran_src():
    if "FORTVINE_SRC" in os.environ:
        return Path(os.environ["FORTVINE_SRC"]).resolve()
    return HERE / "fsrc" if (HERE / "fsrc").is_dir() else HERE.parent / "src"


def compile_library(target):
    src = fortran_src()
    if not (src / "fortvine.f90").is_file():
        sys.exit(f"fortvine: Fortran sources not found in {src} (build from a clone of the "
                 "repository or from the sdist, or set FORTVINE_SRC)")
    fc = os.environ.get("FC", "gfortran")
    flags = shlex.split(os.environ.get("FFLAGS", FFLAGS))
    # Windows: link the gfortran, gomp, quadmath and winpthread runtimes into the DLL so it
    # loads in any environment; elsewhere the compiler's shared runtimes are used (wheels
    # get them bundled by auditwheel / delocate)
    link = ["-static"] if sys.platform == "win32" else []
    files = [src / f for f in SOURCES] + [HERE / "fortvine_py.f90"]
    with tempfile.TemporaryDirectory() as tmp:
        objs = []
        for f in files:
            obj = Path(tmp) / (f.stem + ".o")
            run([fc, *flags, "-cpp", "-fPIC", "-I", str(src), "-J", tmp, "-c", str(f), "-o", str(obj)])
            objs.append(str(obj))
        run([fc, "-shared", *flags, *link, *objs, "-o", str(target)])


def run(cmd):
    print(" ".join(cmd), flush=True)
    try:
        subprocess.check_call(cmd)
    except FileNotFoundError:
        sys.exit(f"fortvine: Fortran compiler '{cmd[0]}' not found; install gfortran or set FC")


class BuildPy(build_py):
    def run(self):
        super().run()
        out = Path(self.build_lib) / "fortvine" / "_lib"
        out.mkdir(parents=True, exist_ok=True)
        compile_library(out / LIBNAME)


class BdistWheel(bdist_wheel):
    """Platform wheel, but no CPython extension (ctypes): py3-none-<platform>."""

    def get_tag(self):
        return ("py3", "none", super().get_tag()[2])


class BinaryDistribution(Distribution):
    def has_ext_modules(self):           # platform-specific build and wheel
        return True


class Sdist(sdist):
    """Copy the Fortran sources (../src) into the sdist as fsrc/, and the repository LICENSE."""

    def make_release_tree(self, base_dir, files):
        super().make_release_tree(base_dir, files)
        src = fortran_src()
        dest = Path(base_dir) / "fsrc"
        dest.mkdir(exist_ok=True)
        for f in sorted(src.glob("*.f90")) + sorted(src.glob("*.inc")):
            shutil.copy2(f, dest / f.name)
        if (HERE.parent / "LICENSE").is_file():
            shutil.copy2(HERE.parent / "LICENSE", Path(base_dir) / "LICENSE")


setup(cmdclass={"build_py": BuildPy, "bdist_wheel": BdistWheel, "sdist": Sdist}, distclass=BinaryDistribution)
