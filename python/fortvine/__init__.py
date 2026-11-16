"""fortvine: Python binding to the pure-Fortran vine-copula library fortvine.

Every density, h-function, fit, selection and vine operation runs in the compiled Fortran
library (through its C interface); this module only converts arguments and wraps results.
Names follow pyvinecopulib where the C interface offers the same operation.
"""
import ctypes as ct
import os
from collections import namedtuple
from functools import lru_cache

import numpy as np

from ._capi import lib


def _nt(num_threads):
    """Thread count passed to the library: as given, else FORTVINE_THREADS, else half the CPUs
    (the same rule as fortvine's fv_threads)."""
    if num_threads:
        return int(num_threads)
    try:
        v = int(os.environ.get("FORTVINE_THREADS", "0"))
    except ValueError:
        v = 0
    return v if v > 0 else max(1, (os.cpu_count() or 2) // 2)

__version__ = "1.0.0"
__all__ = ["Bicop", "Vinecop", "families", "to_pseudo_obs", "kendall_tau", "cite", "Family", "__version__"]

Family = namedtuple("Family", "code name npars")

# fv_controls default family set: every parametric hand-coded family except tev (codes 0-22)
DEFAULT_FAMILY_SET = tuple(range(23))
_NOAUX = (24, 25)
_CRIT = {"loglik": 0, "aic": 1, "bic": 2, "mbic": 3, "mbicv": 4}
_METHOD = {"mle": 0, "itau": 1}
_TAU = {"auto": 0, "exact": 1, "approx": 2}
_TREE_CRIT = {"tau": 0, "rho": 1, "joe": 2, "hoeffd": 3, "xi": 4, "mi": 5, "ferreira": 6}
_TREE_ALGO = {"prim": 0, "mst_prim": 0, "kruskal": 1, "mst_kruskal": 1}


# --------------------------------------------------------------------------- families
def _fam_name(code):
    buf = ct.create_string_buffer(48)
    lib.fvpy_fam_name(int(code), buf)
    return buf.raw.decode("ascii").strip()


@lru_cache(maxsize=None)
def _families():
    n = lib.fvpy_families(0, np.zeros(1, np.int32))
    codes = np.zeros(n, np.int32)
    lib.fvpy_families(n, codes)
    return tuple(Family(int(c), _fam_name(c), int(lib.fvpy_npars(int(c)))) for c in codes)


def families():
    """All distinct families of the library: list of Family(code, name, npars)."""
    return list(_families())


@lru_cache(maxsize=None)
def _by_name():
    return {f.name.lower(): f.code for f in _families()}


def _code(family):
    """Family name (case-insensitive, as fv_fam_name) or integer code -> code."""
    if isinstance(family, str):
        try:
            code = _by_name()[family.lower()]
        except KeyError:
            raise ValueError(f"unknown family {family!r}; see fortvine.families()") from None
    else:
        code = int(family)
        if _fam_name(code) == "invalid":
            raise ValueError(f"invalid family code {code}")
    if code in _NOAUX:
        raise ValueError("tll and bernstein (codes 24, 25) need their nonparametric estimate, which the "
                         "C interface does not pass; use the Fortran API for them")
    return code


def _family_set(family_set):
    codes = DEFAULT_FAMILY_SET if family_set is None else [_code(f) for f in family_set]
    return np.ascontiguousarray(codes, dtype=np.int32)


def _ci6(criterion="bic", method="mle", preselect=True, allow_rotations=True, fast_mle=False, tau_mode="exact"):
    return [_CRIT[criterion.lower()], _METHOD[method.lower()], int(bool(preselect)), int(bool(allow_rotations)),
            int(bool(fast_mle)), _TAU[tau_mode.lower()]]


def _data(u, k=None):
    u = np.ascontiguousarray(u, dtype=np.float64)
    if u.ndim != 2 or (k is not None and u.shape[1] != k):
        raise ValueError(f"expected an (n, {k or 'd'}) array, got shape {u.shape}")
    return u


# --------------------------------------------------------------------------- utilities
def to_pseudo_obs(x, num_threads=None):
    """Column-wise pseudo-observations: average ranks / (n + 1) (fv_pseudo_obs)."""
    x = np.asarray(x, dtype=np.float64)
    x2 = _data(x.reshape(-1, 1) if x.ndim == 1 else x)
    out = np.empty_like(x2)
    lib.fv_pseudo_obs(x2.shape[0], x2.shape[1], x2, out, _nt(num_threads))
    return out.reshape(x.shape)


def kendall_tau(x, y):
    """Kendall's tau_b of two vectors (fv_kendall_tau)."""
    x, y = np.ascontiguousarray(x, np.float64), np.ascontiguousarray(y, np.float64)
    return lib.fv_kendall_tau(x.size, x, y)


_CITATION = """Shojaeezadeh, S. A. (2026). fortvine: a complete compendium to work on vine copula modeling. In preparation.

@unpublished{shojaeezadeh2026fortvine,
  author = {Shojaeezadeh, Shahab Aldin},
  title  = {fortvine: a complete compendium to work on vine copula modeling},
  year   = {2026},
  note   = {In preparation. Software: https://github.com/Bluerrror/fortvine}
}"""


def cite():
    """Print and return the preferred citation of fortvine, with BibTeX."""
    print(_CITATION)
    return _CITATION


# --------------------------------------------------------------------------- pair copula
class Bicop:
    """Bivariate copula: Bicop(family="clayton", rotation=0, parameters=[2.0])."""

    def __init__(self, family="indep", rotation=0, parameters=None):
        self._fam = _code(family)
        if rotation not in (0, 90, 180, 270):
            raise ValueError("rotation must be 0, 90, 180 or 270")
        self.rotation = int(rotation)
        self.parameters = parameters
        self._loglik = None

    family = property(lambda self: _fam_name(self._fam), doc="family name (fv_fam_name)")
    family_code = property(lambda self: self._fam, doc="integer family code")
    npars = property(lambda self: int(lib.fvpy_npars(self._fam)), doc="number of parameters")

    @property
    def parameters(self):
        return self._par[:self.npars].copy()

    @parameters.setter
    def parameters(self, p):
        p = np.zeros(self.npars) if p is None else np.asarray(p, dtype=np.float64).ravel()
        if p.size != self.npars:
            raise ValueError(f"{self.family} has {self.npars} parameter(s), got {p.size}")
        self._par = np.zeros(3)
        self._par[:p.size] = p

    def _eval(self, u, what, num_threads):
        u = _data(u, 2)
        out = np.empty(len(u))
        lib.fv_bicop_eval(self._fam, self.rotation, self._par, len(u), u, what, out, _nt(num_threads))
        return out

    def pdf(self, u, num_threads=None):
        return self._eval(u, 1, num_threads)

    def logpdf(self, u, num_threads=None):
        return self._eval(u, 0, num_threads)

    def hfunc1(self, u, num_threads=None):
        """h1(u1, u2) = P(U2 <= u2 | U1 = u1)."""
        return self._eval(u, 2, num_threads)

    def hfunc2(self, u, num_threads=None):
        """h2(u1, u2) = P(U1 <= u1 | U2 = u2)."""
        return self._eval(u, 3, num_threads)

    def hinv1(self, u, num_threads=None):
        """v with h1(u1, v) = u2."""
        return self._eval(u, 4, num_threads)

    def hinv2(self, u, num_threads=None):
        """v with h2(v, u2) = u1."""
        return self._eval(u, 5, num_threads)

    def cdf(self, u):
        u = _data(u, 2)
        out = np.empty(len(u))
        lib.fvpy_bicop_cdf(self._fam, self.rotation, self._par, len(u), u, out)
        return out

    def loglik(self, u=None, num_threads=None):
        """Log-likelihood on u, or (u=None) the one of the last fit / select."""
        if u is None:
            if self._loglik is None:
                raise ValueError("no fit yet: pass data")
            return self._loglik
        return float(self.logpdf(u, num_threads).sum())

    def fit(self, u, method="mle", fast_mle=False, tau_mode="exact", num_threads=None):
        """Estimate the parameters of the current family and rotation (fv_bicop_fit)."""
        u = _data(u, 2)
        par, ll = np.zeros(3), ct.c_double()
        ci = np.array(_ci6(method=method, fast_mle=fast_mle, tau_mode=tau_mode), np.int32)
        lib.fv_bicop_fit(len(u), u, self._fam, self.rotation, ci, par, ct.byref(ll), _nt(num_threads))
        self._par, self._loglik = par, ll.value
        return self

    def select(self, u, family_set=None, criterion="bic", method="mle", preselect=True, allow_rotations=True,
               fast_mle=False, tau_mode="exact", psi0=0.9, num_threads=None):
        """Select family, rotation and parameters (fv_bicop_select)."""
        u = _data(u, 2)
        fs = _family_set(family_set)
        ci = np.array(_ci6(criterion, method, preselect, allow_rotations, fast_mle, tau_mode), np.int32)
        fam, rot, par, ll = ct.c_int(), ct.c_int(), np.zeros(3), ct.c_double()
        lib.fv_bicop_select(len(u), u, fs.size, fs, ci, psi0, ct.byref(fam), ct.byref(rot), par, ct.byref(ll),
                            _nt(num_threads))
        self._fam, self.rotation, self._par, self._loglik = fam.value, rot.value, par, ll.value
        return self

    def simulate(self, n, seed=None, num_threads=None):
        """n draws (n, 2): u2 = hinv1(u1, w) for independent uniforms u1, w."""
        w = np.random.default_rng(seed).random((n, 2))
        return np.column_stack([w[:, 0], self.hinv1(w, num_threads)])

    def __repr__(self):
        return f"Bicop(family={self.family!r}, rotation={self.rotation}, parameters={self.parameters.tolist()})"


# --------------------------------------------------------------------------- vine copula
class Vinecop:
    """Vine copula in natural order (as fv_vinecop).

    order (d,): variable labels (1-based); st, fam, rot (d, d) and par (d, d, 3) are indexed
    [tree, edge] (0-based): edge e of tree t joins order[e] and order[st[t, e] - 1].
    """

    def __init__(self, d=0):
        self.d, self.trunc = int(d), 0
        self.order = np.arange(1, d + 1, dtype=np.int32)
        self.st, self.fam, self.rot = (np.zeros((d, d), np.int32) for _ in range(3))
        self.par, self.edge_ll = np.zeros((d, d, 3)), np.zeros((d, d))
        self._fitted = False

    # structure --------------------------------------------------------------
    def _set_structure(self, kind, order, trunc_lvl):
        order = np.ascontiguousarray(order, np.int32)
        d = order.size
        self.__init__(d)
        st, tr = np.zeros((d, d), np.int32), ct.c_int()
        lib.fvpy_vine_structure(kind, d, order, d - 1 if trunc_lvl is None else trunc_lvl, st, ct.byref(tr))
        self.order, self.st, self.trunc = order, st, tr.value
        return self

    def set_dvine(self, order, trunc_lvl=None):
        """D-vine (path) in the given order (1-based labels); families are set by select."""
        return self._set_structure(0, order, trunc_lvl)

    def set_cvine(self, order, trunc_lvl=None):
        """C-vine; the root of tree t is order[d - t] (the last element is the first root)."""
        return self._set_structure(1, order, trunc_lvl)

    @property
    def trunc_lvl(self):
        return self.trunc

    @property
    def matrix(self):
        """R-vine matrix in original labels (fv_vinecop%matrix)."""
        m = np.zeros((self.d, self.d), np.int32)
        lib.fvpy_vine_matrix(self.d, self._i(self.order), self.trunc, self._i(self.st), m)
        return m.T.copy()

    def _edges(self):
        return [(t, e) for t in range(self.trunc) for e in range(self.d - 1 - t)]

    @property
    def families(self):
        """Family names per tree, [tree][edge]."""
        return [[_fam_name(self.fam[t, e]) for e in range(self.d - 1 - t)] for t in range(self.trunc)]

    @property
    def rotations(self):
        return [[int(self.rot[t, e]) for e in range(self.d - 1 - t)] for t in range(self.trunc)]

    @property
    def parameters(self):
        return [[self.get_pair_copula(t, e).parameters for e in range(self.d - 1 - t)] for t in range(self.trunc)]

    @property
    def npars(self):
        return sum(int(lib.fvpy_npars(int(self.fam[t, e]))) for t, e in self._edges())

    def get_pair_copula(self, tree, edge):
        n = int(lib.fvpy_npars(int(self.fam[tree, edge])))
        return Bicop(int(self.fam[tree, edge]), int(self.rot[tree, edge]), self.par[tree, edge, :n])

    def set_pair_copula(self, tree, edge, bicop):
        self.fam[tree, edge], self.rot[tree, edge] = bicop.family_code, bicop.rotation
        self.par[tree, edge] = bicop._par

    # selection ----------------------------------------------------------------
    def select(self, u, family_set=None, criterion="bic", method="mle", preselect=True, allow_rotations=True,
               fast_mle=False, tau_mode="exact", tree_criterion="tau", tree_algorithm="prim", trunc_lvl=None,
               threshold=0.0, psi0=0.9, fixed=None, num_threads=None):
        """Select structure (unless fixed), families and parameters (fv_vine_select).

        fixed defaults to True when a structure is set (set_dvine / set_cvine).
        """
        u = _data(u)
        n, d = u.shape
        if self.d != d:
            self.__init__(d)
        fx = self.trunc > 0 if fixed is None else bool(fixed)
        T = d - 1 if trunc_lvl is None or trunc_lvl < 0 else min(d - 1, int(trunc_lvl))
        if fx:
            T = min(T, self.trunc)
        fs = _family_set(family_set)
        ci = np.array(_ci6(criterion, method, preselect, allow_rotations, fast_mle, tau_mode)
                      + [_TREE_CRIT[tree_criterion.lower()], _TREE_ALGO[tree_algorithm.lower()], T, int(fx)], np.int32)
        cr = np.array([threshold, psi0], np.float64)
        order = self._i(self.order).copy()
        st, fam, rot = (np.zeros((d, d), np.int32) for _ in range(3))
        par, ll, info = np.zeros((d, d, 3)), np.zeros((d, d)), ct.c_int()
        lib.fv_vine_select(n, d, u, fs.size, fs, ci, cr, order, self._i(self.st), st, fam, rot, par, ll,
                           _nt(num_threads), ct.byref(info))
        if info.value != 0:
            raise RuntimeError(f"fv_vine_select failed (info = {info.value})")
        self.order, self.st, self.fam, self.rot, self.par, self.edge_ll = order, st, fam, rot, par, ll
        self.trunc, self._fitted = T, True
        return self

    @classmethod
    def from_data(cls, u, **kwargs):
        """Select a vine (structure, families, parameters) for pseudo-observations u (n, d)."""
        return cls(np.shape(u)[1]).select(u, fixed=False, **kwargs)

    # evaluation -----------------------------------------------------------------
    @staticmethod
    def _i(a):
        return np.ascontiguousarray(a, np.int32)

    def _args(self):
        fam = self._i(self.fam)
        if np.isin(fam, _NOAUX).any():
            raise ValueError("tll / bernstein edges are not usable through the C interface")
        return (self._i(self.order), self.trunc, self._i(self.st), fam, self._i(self.rot),
                np.ascontiguousarray(self.par, np.float64))

    def pdf(self, u, num_threads=None):
        u = _data(u, self.d)
        out = np.empty(len(u))
        lib.fv_vine_pdf(len(u), 1, self.d, -1, u, *self._args(), out, _nt(num_threads))
        return out

    def loglik(self, u=None, num_threads=None):
        """Log-likelihood on u, or (u=None) the sum of the edge log-likelihoods of the last select."""
        if u is None:
            if not self._fitted:
                raise ValueError("no fit yet: pass data")
            return float(self.edge_ll.sum())
        return float(np.log(self.pdf(u, num_threads)).sum())

    def _transform(self, f, u, num_threads):
        u = _data(u, self.d)
        out = np.empty_like(u)
        f(len(u), self.d, u, *self._args(), out, _nt(num_threads))
        return out

    def rosenblatt(self, u, num_threads=None):
        return self._transform(lib.fv_vine_rosenblatt, u, num_threads)

    def inverse_rosenblatt(self, w, num_threads=None):
        return self._transform(lib.fv_vine_inverse_rosenblatt, w, num_threads)

    def simulate(self, n, seed=None, num_threads=None):
        """n draws (n, d): inverse Rosenblatt transform of independent uniforms."""
        return self.inverse_rosenblatt(np.random.default_rng(seed).random((n, self.d)), num_threads)

    def __str__(self):
        lines = [f"fortvine vine copula: d = {self.d}, trunc = {self.trunc}",
                 "order: " + " ".join(map(str, self.order))]
        for t, e in self._edges():
            b = self.get_pair_copula(t, e)
            lines.append(f"  tree {t + 1} edge {e + 1}: {b.family:<12} rot {b.rotation:3d} par "
                         + " ".join(f"{p:.4f}" for p in b.parameters))
        return "\n".join(lines)

    def __repr__(self):
        return f"<fortvine.Vinecop d={self.d} trunc={self.trunc} npars={self.npars}>"
