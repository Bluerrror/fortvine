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

I64ref = ct.POINTER(ct.c_int64)

_SIG = {
    # src/fortvine.f90
    "fv_version": ([], I),
    "fv_kendall_tau": ([I, dptr, dptr], D),
    "fv_bicop_eval": ([I, I, dptr, I, dptr, I, dptr, I], None),
    "fv_pseudo_obs": ([I, I, dptr, dptr, I], None),
    # src/fortvine_capi.f90: C interface to the typed API (fortvine_api)
    "fvc_naux": ([], I),
    "fvc_nk": ([], I),
    "fvc_fam_name": ([I, ct.c_char_p], None),
    "fvc_families": ([I, iptr, iptr], I),
    "fvc_par_bounds": ([I, dptr, dptr], None),
    "fvc_gen_flip": ([Iref, dptr], None),
    "fvc_set_abort": ([I], None),
    "fvc_fitcache": ([I, I, I64ref, I64ref, Iref], None),
    "fvc_bicop_eval": ([I, I, dptr, I, dptr, I, dptr, I, dptr], None),
    "fvc_bicop_fit": ([I, dptr, I, dptr, I, I, iptr, iptr, dptr, Iref, Iref, dptr, Dref, I, dptr, Dref], None),
    "fvc_vine_select": ([I, I, dptr, dptr, iptr, I, dptr, I, I, iptr, iptr, dptr, Iref, iptr, iptr, iptr, iptr,
                         dptr, dptr, I, dptr, Dref, Dref, Dref, Iref], None),
    "fvc_vine_eval": ([I, I, iptr, iptr, iptr, iptr, dptr, I, dptr, iptr, I, iptr, dptr, dptr, dptr,
                       I, I, dptr, dptr, I, I, dptr], None),
    "fvc_vine_pdf_grid": ([I, I, iptr, iptr, iptr, iptr, dptr, I, dptr, iptr, I, iptr, dptr, dptr, dptr,
                           I, I, dptr, dptr, I, I, dptr], None),
    # (fvc_vine_relax takes trunc by reference)
    "fvc_vine_relax": ([I, I, dptr, I, dptr, I, I, I, iptr, iptr, dptr, Iref, iptr, iptr, iptr, iptr, dptr, dptr,
                        I, dptr, iptr, I, Iref, iptr, dptr, dptr, dptr, Dref, Dref], None),
    "fvc_vine_rel_params": ([I, I, dptr, dptr, dptr, I, dptr, dptr], None),
    "fvc_vine_structure": ([I, I, iptr, I, iptr, Iref], None),
    "fvc_vine_matrix": ([I, iptr, I, iptr, iptr], None),
}
for _n, (_a, _r) in _SIG.items():
    getattr(lib, _n).argtypes = _a
    getattr(lib, _n).restype = _r
