"""fortvine_ctypes.py -- call the fortvine shared library from Python (ctypes + numpy).

Build the library first (bash build.sh), then run

    python example/python/fortvine_ctypes.py

The library is looked up in build/ next to this repository (or at $FORTVINE_LIB).

Array convention of the C API: observations are the LAST Fortran dimension, i.e. one
observation's values are contiguous.  A numpy array of shape (n, k) in C order is exactly
the Fortran array (k, n) the kernels expect -- no transposes, no copies.
"""
import ctypes as ct
import os
import shutil
import sys
from pathlib import Path

import numpy as np

FAMILY = ['indep', 'gaussian', 'student', 'clayton', 'gumbel', 'frank', 'joe', 'bb1', 'bb6', 'bb7',
          'bb8', 'tawn', 'plackett', 'fgm', 'amh', 'galambos', 'husler_reiss', 'bb2', 'bb3', 'bb5',
          'bb9', 'bb10', 'bb4', 'tev', 'tll', 'bernstein']


def load_library():
    root = Path(__file__).resolve().parents[2]
    name = {'win32': 'libfortvine.dll', 'darwin': 'libfortvine.dylib'}.get(sys.platform, 'libfortvine.so')
    path = Path(os.environ.get('FORTVINE_LIB', root / 'build' / name))
    if sys.platform == 'win32':
        # the DLL needs the MinGW runtime (libgomp-1.dll, libquadmath-0.dll, libwinpthread-1.dll)
        gfortran = shutil.which('gfortran')
        if gfortran:
            os.add_dll_directory(str(Path(gfortran).parent))
    return ct.CDLL(str(path))


lib = load_library()
I, D = ct.c_int, ct.c_double
dptr = np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')
iptr = np.ctypeslib.ndpointer(dtype=np.int32, flags='C_CONTIGUOUS')
Iref, Dref = ct.POINTER(I), ct.POINTER(D)

# signatures of the bind(C) routines in src/fortvine.f90 (scalars by value unless marked out)
lib.fv_version.restype = I
lib.fv_kendall_tau.argtypes = [I, dptr, dptr]
lib.fv_kendall_tau.restype = D
lib.fv_bicop_eval.argtypes = [I, I, dptr, I, dptr, I, dptr, I]           # fam, rot, par(3), n, u(n,2), what, out(n), nthreads
lib.fv_bicop_fit.argtypes = [I, dptr, I, I, iptr, dptr, Dref, I]          # n, u, fam, rot, ci(6), par(3), ll, nthreads
lib.fv_bicop_select.argtypes = [I, dptr, I, iptr, iptr, D, Iref, Iref, dptr, Dref, I]
lib.fv_pseudo_obs.argtypes = [I, I, dptr, dptr, I]                        # n, d, x(n,d), out(n,d), nthreads
lib.fv_tau_matrix.argtypes = [I, I, dptr, dptr, I]                        # n, d, x(n,d), tau(d,d), nthreads
lib.fv_vine_select.argtypes = [I, I, dptr, I, iptr, iptr, dptr, iptr, iptr, iptr, iptr, iptr, dptr, dptr, I, Iref]
lib.fv_vine_pdf.argtypes = [I, I, I, I, dptr, iptr, I, iptr, iptr, iptr, dptr, dptr, I]
lib.fv_vine_inverse_rosenblatt.argtypes = [I, I, dptr, iptr, I, iptr, iptr, iptr, dptr, dptr, I]
for f in (lib.fv_bicop_eval, lib.fv_bicop_fit, lib.fv_bicop_select, lib.fv_pseudo_obs, lib.fv_tau_matrix,
          lib.fv_vine_select, lib.fv_vine_pdf, lib.fv_vine_inverse_rosenblatt):
    f.restype = None

NTH = 4
# ci(1:6): selection criterion (2 = BIC), method (0 = MLE), preselect, allow rotations,
#          fast one-parameter fits, tau estimator (1 = exact) -- the fv_controls defaults
CI6 = np.array([2, 0, 1, 1, 0, 1], dtype=np.int32)
PARAMETRIC = np.arange(23, dtype=np.int32)          # families 0..22 (the default set: no tev, tll, bernstein)


def bicop_eval(fam, rot, par, u, what):
    """what: 0 log-density, 1 density, 2 h1, 3 h2, 4 h1 inverse, 5 h2 inverse; u has shape (n, 2)."""
    u = np.ascontiguousarray(u, dtype=np.float64)
    out = np.empty(len(u))
    lib.fv_bicop_eval(fam, rot, np.asarray(par, dtype=np.float64), len(u), u, what, out, NTH)
    return out


# ------------------------------------------------------------------------------ one pair
print('fortvine C API version', lib.fv_version())
par = np.array([2.0, 0.0, 0.0])                      # Clayton, theta = 2
g = np.array([0.1, 0.5, 0.9])
print('Clayton(2) density on the diagonal:', bicop_eval(3, 0, par, np.column_stack([g, g]), 1).round(4))

rng = np.random.default_rng(1)
n = 2000
u1, w = rng.random(n), rng.random(n)
u2 = bicop_eval(3, 0, par, np.column_stack([u1, w]), 4)      # u2 = h1^-1(w | u1)
u = np.column_stack([u1, u2])
print('empirical Kendall tau: %.3f' % lib.fv_kendall_tau(n, u1, u2))

est, ll = np.zeros(3), D()
lib.fv_bicop_fit(n, u, 3, 0, CI6, est, ct.byref(ll), NTH)
print('fitted theta: %.4f   log-likelihood: %.2f' % (est[0], ll.value))

fam, rot = I(), I()
lib.fv_bicop_select(n, u, len(PARAMETRIC), PARAMETRIC, CI6, 0.9, ct.byref(fam), ct.byref(rot), est, ct.byref(ll), NTH)
print('selected:', FAMILY[fam.value], 'rotation', rot.value, 'par', est.round(3))

# ------------------------------------------------------------------------------ a vine
# 4-dimensional data with a Gaussian dependence, turned into pseudo-observations
d, n = 4, 3000
corr = 0.6 ** np.abs(np.subtract.outer(np.arange(d), np.arange(d)))   # AR(1) correlation
x = rng.multivariate_normal(np.zeros(d), corr, size=n)
u = np.empty((n, d))
lib.fv_pseudo_obs(n, d, x, u, NTH)

# ci(1:10) = CI6, tree criterion (0 tau), tree algorithm (0 Prim), truncation (-1 = none),
# structure given (0 = select it); cr = threshold, psi0
ci10 = np.concatenate([CI6, [0, 0, -1, 0]]).astype(np.int32)
cr = np.array([0.0, 0.9])
order = np.zeros(d, dtype=np.int32)
st_in = np.zeros((d, d), dtype=np.int32)
st, fam, rot = (np.zeros((d, d), dtype=np.int32) for _ in range(3))    # [tree, edge]
vpar, edge_ll = np.zeros((d, d, 3)), np.zeros((d, d))                   # [tree, edge, k], [tree, edge]
info = I()
lib.fv_vine_select(n, d, u, len(PARAMETRIC), PARAMETRIC, ci10, cr, order, st_in, st, fam, rot, vpar, edge_ll,
                   NTH, ct.byref(info))
assert info.value == 0
print('\nvine order:', order)
for t in range(d - 1):
    print('  tree %d:' % (t + 1), ', '.join('%s(%.3f)' % (FAMILY[fam[t, e]], vpar[t, e, 0]) for e in range(d - 1 - t)))

trunc = d - 1
dens = np.empty(n)                                   # density: one bin per sample, tgt = -1
lib.fv_vine_pdf(n, 1, d, -1, u, order, trunc, st, fam, rot, vpar, dens, NTH)
print('log-likelihood: %.2f  (sum of edge log-likelihoods %.2f)' % (np.log(dens).sum(), edge_ll.sum()))

# simulation = inverse Rosenblatt transform of independent uniforms
sim = np.empty((n, d))
lib.fv_vine_inverse_rosenblatt(n, d, rng.random((n, d)), order, trunc, st, fam, rot, vpar, sim, NTH)
t_data, t_sim = np.empty((d, d)), np.empty((d, d))
lib.fv_tau_matrix(n, d, u, t_data, NTH)
lib.fv_tau_matrix(n, d, sim, t_sim, NTH)
print('Kendall tau, data      :', t_data[0].round(3))
print('Kendall tau, simulated :', t_sim[0].round(3))
