"""Load the compiled fortvine library and declare the C signatures (docs/manual/c-api.md).

Observations are the last Fortran dimension, so a C-order NumPy (n, k) array is the
Fortran array (k, n) the routines expect.
"""
import ctypes as ct
import os
import sys
from pathlib import Path

import numpy as np

_NAME = {"win32": "libfortvine.dll", "darwin": "libfortvine.dylib"}.get(sys.platform, "libfortvine.so")


def _load():
    path = Path(os.environ.get("FORTVINE_LIB", Path(__file__).parent / "_lib" / _NAME))
    if not path.is_file():
        raise ImportError(f"fortvine: compiled library not found at {path}; reinstall with "
                          "`pip install ./python` (needs gfortran)")
    if sys.platform == "win32":
        os.add_dll_directory(str(path.parent))
    return ct.CDLL(str(path))


lib = _load()

I, D = ct.c_int, ct.c_double
Iref, Dref = ct.POINTER(I), ct.POINTER(D)
dptr = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
iptr = np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS")

_SIG = {
    # src/fortvine.f90
    "fv_version": ([], I),
    "fv_kendall_tau": ([I, dptr, dptr], D),
    "fv_bicop_eval": ([I, I, dptr, I, dptr, I, dptr, I], None),
    "fv_bicop_fit": ([I, dptr, I, I, iptr, dptr, Dref, I], None),
    "fv_bicop_select": ([I, dptr, I, iptr, iptr, D, Iref, Iref, dptr, Dref, I], None),
    "fv_pseudo_obs": ([I, I, dptr, dptr, I], None),
    "fv_vine_select": ([I, I, dptr, I, iptr, iptr, dptr, iptr, iptr, iptr, iptr, iptr, dptr, dptr, I, Iref], None),
    "fv_vine_pdf": ([I, I, I, I, dptr, iptr, I, iptr, iptr, iptr, dptr, dptr, I], None),
    "fv_vine_rosenblatt": ([I, I, dptr, iptr, I, iptr, iptr, iptr, dptr, dptr, I], None),
    "fv_vine_inverse_rosenblatt": ([I, I, dptr, iptr, I, iptr, iptr, iptr, dptr, dptr, I], None),
    # python/fortvine_py.f90
    "fvpy_fam_name": ([I, ct.c_char_p], None),
    "fvpy_npars": ([I], I),
    "fvpy_families": ([I, iptr], I),
    "fvpy_bicop_cdf": ([I, I, dptr, I, dptr, dptr], None),
    "fvpy_vine_structure": ([I, I, iptr, I, iptr, Iref], None),
    "fvpy_vine_matrix": ([I, iptr, I, iptr, iptr], None),
}
for _n, (_a, _r) in _SIG.items():
    getattr(lib, _n).argtypes = _a
    getattr(lib, _n).restype = _r
