"""fortvine: Python binding to the pure-Fortran vine-copula library fortvine.

Every density, h-function, fit, selection and vine operation runs in the compiled Fortran
library (through its C interfaces, src/fortvine.f90 and src/fortvine_capi.f90); this module
only converts arguments and wraps results.  Names follow pyvinecopulib where the library
offers the same operation.
"""
import ctypes as ct
import os
import threading
from concurrent.futures import ThreadPoolExecutor, TimeoutError as _FutureTimeout
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
__all__ = ["Bicop", "Vinecop", "families", "family_group", "to_pseudo_obs", "kendall_tau", "par_bounds",
           "gen_flip", "set_abort", "fitcache_enable", "fitcache_clear", "fitcache_disable", "fitcache_stats",
           "cite", "Family", "FitCacheStats", "NAUX", "NK", "__version__"]

Family = namedtuple("Family", "code name npars")
FitCacheStats = namedtuple("FitCacheStats", "hits misses used")

NAUX = int(lib.fvc_naux())        # length of the estimate of a tll / bernstein pair copula
NK = int(lib.fvc_nk())            # anchors of a non-simplified (relaxed) edge

# fv_controls default family set: every parametric hand-coded family except tev (codes 0-22)
DEFAULT_FAMILY_SET = tuple(range(23))
_NP = (24, 25)                    # tll, bernstein
_GROUPS = {
    "default": DEFAULT_FAMILY_SET,           # fv_controls default
    "parametric": tuple(range(12)),          # FV_PARAMETRIC: vinecopulib's parametric families
    "all_parametric": tuple(range(24)),      # FV_ALL: every hand-coded parametric family (with tev)
    "nonparametric": _NP,                    # FV_NONPAR: tll, bernstein
    "all": tuple(range(26)),                 # every hand-coded family, codes 0-25
    "one_par": (1, 3, 4, 5, 6),              # FV_ONE_PAR
}
_CRIT = {"loglik": 0, "aic": 1, "bic": 2, "mbic": 3, "mbicv": 4}
_METHOD = {"mle": 0, "itau": 1}
_TAU = {"auto": 0, "exact": 1, "approx": 2}
_TREE_CRIT = {"tau": 0, "rho": 1, "joe": 2, "hoeffd": 3, "xi": 4, "mi": 5, "ferreira": 6}
_TREE_ALGO = {"prim": 0, "mst_prim": 0, "kruskal": 1, "mst_kruskal": 1}
_RELAX_MODE = {"none": 0, "test": 1, "all": 2}
# pyvinecopulib (FitControls*) and fortvine (fv_controls) spellings of the options
_ALIASES = {"selection_criterion": "criterion", "selcrit": "criterion", "parametric_method": "method",
            "tree_crit": "tree_criterion", "tree_algo": "tree_algorithm", "preselect_families": "preselect",
            "select_trunc": "select_trunc_lvl", "tau": "tau_mode"}


def _aliases(kw, allowed):
    out = {}
    for k, v in kw.items():
        name = _ALIASES.get(k, k)
        if name not in allowed:
            raise TypeError(f"unexpected keyword argument {k!r}")
        out[name] = v
    return out


def _key(table, value, what):
    if isinstance(value, (int, np.integer)):
        if int(value) not in table.values():
            raise ValueError(f"invalid {what} {value}")
        return int(value)
    try:
        return table[str(value).lower()]
    except KeyError:
        raise ValueError(f"unknown {what} {value!r}; one of {sorted(table)}") from None


# --------------------------------------------------------------------------- families
def _fam_name(code):
    buf = ct.create_string_buffer(48)
    lib.fvc_fam_name(int(code), buf)
    return buf.raw.decode("ascii").strip()


@lru_cache(maxsize=None)
def _families():
    n = lib.fvc_families(0, np.zeros(1, np.int32), np.zeros(1, np.int32))
    codes, npars = np.zeros(n, np.int32), np.zeros(n, np.int32)
    lib.fvc_families(n, codes, npars)
    return tuple(Family(int(c), _fam_name(c), int(k)) for c, k in zip(codes, npars))


def families():
    """All distinct families of the library: list of Family(code, name, npars); npars is -1 for
    tll and bernstein, whose number of parameters is the effective one of each estimate."""
    return list(_families())


def family_group(name):
    """Family codes of a group: 'default' (0-22), 'parametric' (0-11), 'all_parametric' (0-23),
    'nonparametric' (24, 25), 'all' (0-25), 'one_par' (1, 3, 4, 5, 6) or 'generic' (10000+)."""
    name = name.lower()
    if name == "generic":
        return tuple(f.code for f in _families() if f.code >= 10000)
    try:
        return _GROUPS[name]
    except KeyError:
        raise ValueError(f"unknown family group {name!r}; one of {sorted(_GROUPS) + ['generic']}") from None


@lru_cache(maxsize=None)
def _by_name():
    return {f.name.lower(): f.code for f in _families()}


def _code(family):
    """Family name (case-insensitive, as fv_fam_name) or integer code -> code."""
    if isinstance(family, str):
        try:
            return _by_name()[family.lower()]
        except KeyError:
            raise ValueError(f"unknown family {family!r}; see fortvine.families()") from None
    code = int(family)
    if _fam_name(code) == "invalid":
        raise ValueError(f"invalid family code {code}")
    return code


def _bounds(code):
    lb, ub = np.zeros(3), np.zeros(3)
    lib.fvc_par_bounds(int(code), lb, ub)
    return lb, ub


@lru_cache(maxsize=None)
def _npars(code):
    """Number of parameters of a parametric family (0 for indep, tll, bernstein)."""
    for f in _families():
        if f.code == code:
            return max(f.npars, 0)
    lb, ub = _bounds(code)            # a generic code that duplicates another one
    return int(np.sum(ub > lb))


def par_bounds(family):
    """(lower, upper) parameter bounds of a family (the MLE search box, fv_par_bounds)."""
    code = _code(family)
    k = _npars(code)
    lb, ub = _bounds(code)
    return lb[:k], ub[:k]


def gen_flip(family, parameters):
    """Generic family code and parameters of the transposed copula C(u2, u1) (gen_flip);
    other families are returned unchanged."""
    fam = ct.c_int(_code(family))
    p = np.zeros(3)
    q = np.asarray(parameters, np.float64).ravel()
    p[:q.size] = q
    lib.fvc_gen_flip(ct.byref(fam), p)
    return fam.value, p[:_npars(fam.value)]


def _family_set(family_set):
    """None, a family / group name, a code, or a list of them -> int32 codes (no duplicates)."""
    if family_set is None:
        items = DEFAULT_FAMILY_SET
    elif isinstance(family_set, (str, int, np.integer)):
        items = [family_set]
    else:
        items = family_set
    codes = []
    for f in items:
        if isinstance(f, str) and (f.lower() in _GROUPS or f.lower() == "generic"):
            codes.extend(family_group(f))
        else:
            codes.append(_code(f))
    codes = list(dict.fromkeys(codes))
    if not codes or len(codes) > 192:
        raise ValueError("a family set needs 1 to 192 families")
    return np.ascontiguousarray(codes, dtype=np.int32)


def _controls(family_set=None, criterion="bic", method="mle", tree_criterion="tau", tree_algorithm="prim",
              tau_mode="exact", trunc_lvl=None, num_threads=None, preselect=True, allow_rotations=True,
              fast_mle=False, select_trunc_lvl=False, select_threshold=False, threshold=0.0, psi0=0.9):
    """fv_controls as the (famset, ci(12), cr(2)) arrays of fortvine_capi."""
    tl = -1 if trunc_lvl is None or trunc_lvl < 0 else int(trunc_lvl)
    ci = [_key(_CRIT, criterion, "criterion"), _key(_METHOD, method, "method"), int(bool(preselect)),
          int(bool(allow_rotations)), int(bool(fast_mle)), _key(_TAU, tau_mode, "tau mode"),
          _key(_TREE_CRIT, tree_criterion, "tree criterion"), _key(_TREE_ALGO, tree_algorithm, "tree algorithm"),
          tl, _nt(num_threads), int(bool(select_trunc_lvl)), int(bool(select_threshold))]
    return (_family_set(family_set), np.array(ci, np.int32), np.array([threshold, psi0], np.float64))


def _data(u, k=None):
    u = np.ascontiguousarray(u, dtype=np.float64)
    if u.ndim != 2 or (k is not None and u.shape[1] != k):
        raise ValueError(f"expected an (n, {k or 'd'}) array, got shape {u.shape}")
    return u


def _weights(weights, n):
    if weights is None:
        return 0, np.zeros(1)
    w = np.ascontiguousarray(weights, dtype=np.float64).ravel()
    if w.size != n:
        raise ValueError(f"weights: expected {n} values, got {w.size}")
    return 1, w


def _split(u, disc):
    """Data of a vine with discreteness flags disc(d) -> (u, u_minus), both (n, d).

    pyvinecopulib convention: with k discrete variables u is (n, d + k), the last k columns
    being F(x-) of the discrete variables in their order, or (n, 2d) with F(x-) of every
    variable in columns d..2d-1 (only the discrete ones are read)."""
    d, k = disc.size, int(disc.sum())
    u = np.asarray(u, dtype=np.float64)
    if k == 0:
        x = _data(u, d)
        return x, x
    if u.ndim != 2 or u.shape[1] not in (d + k, 2 * d):
        raise ValueError(f"{k} discrete variable(s): expected an (n, {d + k}) array (F(x) for the {d} variables, "
                         f"then F(x-) for the discrete ones) or (n, {2 * d}), got shape {u.shape}")
    x = np.ascontiguousarray(u[:, :d])
    if u.shape[1] == 2 * d:
        return x, np.ascontiguousarray(u[:, d:])
    xm = x.copy()
    xm[:, disc] = u[:, d:]
    return x, xm


def _types(var_types, d=None):
    disc = np.array([str(t).lower() == "d" for t in var_types], dtype=bool)
    if not all(str(t).lower() in ("c", "d") for t in var_types):
        raise ValueError("var_types: 'c' (continuous) or 'd' (discrete) per variable")
    if d is not None and disc.size != d:
        raise ValueError(f"var_types: expected {d} entries, got {disc.size}")
    return disc


# --------------------------------------------------------------------------- run-time switches
def set_abort(flag=True):
    """Set (or clear) the library's abort flag: running fits and selections return early and their
    result is discarded.  Vine selection and relax set it themselves on Ctrl-C (KeyboardInterrupt)."""
    lib.fvc_set_abort(int(bool(flag)))


_worker = []


def _interruptible(fn):
    """Run a long library call so that Ctrl-C aborts it: from the main thread the call runs in a
    persistent worker thread (ctypes releases the GIL; one thread, so its OpenMP pool is reused)
    and a KeyboardInterrupt sets the abort flag, waits for the call, clears the flag and re-raises."""
    if threading.current_thread() is not threading.main_thread():
        return fn()
    if not _worker:
        _worker.append(ThreadPoolExecutor(max_workers=1, thread_name_prefix="fortvine"))
    fut = _worker[0].submit(fn)
    try:
        while True:
            try:
                return fut.result(timeout=0.1)
            except _FutureTimeout:
                continue
    except KeyboardInterrupt:
        set_abort(True)
        try:
            fut.result()
        finally:
            set_abort(False)
        raise


def _fitcache(op, capacity=0):
    h, m, used = ct.c_int64(), ct.c_int64(), ct.c_int()
    lib.fvc_fitcache(op, int(capacity), ct.byref(h), ct.byref(m), ct.byref(used))
    return FitCacheStats(h.value, m.value, used.value)


def fitcache_enable(capacity=20000):
    """Enable the pair-fit memo cache shared by every fit and select (fv_fitcache_enable)."""
    return _fitcache(1, capacity)


def fitcache_clear():
    """Empty the fit cache and reset its counters (fv_fitcache_clear)."""
    return _fitcache(2)


def fitcache_disable():
    """Disable the fit cache and free its table (fv_fitcache_disable)."""
    return _fitcache(0)


def fitcache_stats():
    """FitCacheStats(hits, misses, used) of the fit cache (fv_fitcache_stats)."""
    return _fitcache(3)


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
    """Bivariate copula: Bicop(family="clayton", rotation=0, parameters=[2.0]).

    tll (24) and bernstein (25) carry their estimate in `aux` (NAUX values, set by fit /
    select); their `parameters` are empty and `npars` is the effective number of parameters.
    """

    def __init__(self, family="indep", rotation=0, parameters=None, aux=None):
        self._fam = _code(family)
        if rotation not in (0, 90, 180, 270):
            raise ValueError("rotation must be 0, 90, 180 or 270")
        self.rotation = int(rotation)
        self.parameters = parameters
        self.aux = None
        if aux is not None:
            a = np.ascontiguousarray(aux, np.float64).ravel()
            if a.size != NAUX:
                raise ValueError(f"aux: expected {NAUX} values, got {a.size}")
            self.aux = a.copy()
        self._loglik = None
        self._nobs = None

    family = property(lambda self: _fam_name(self._fam), doc="family name (fv_fam_name)")
    family_code = property(lambda self: self._fam, doc="integer family code")
    is_nonparametric = property(lambda self: self._fam in _NP, doc="tll or bernstein")

    @property
    def npars(self):
        """Number of parameters; for tll / bernstein the effective degrees of freedom of the
        estimate (npars_eff, 0.0 before a fit), the penalty used by AIC / BIC."""
        if self._fam in _NP:
            return 0.0 if self.aux is None else float(self.aux[1])
        return _npars(self._fam)

    @property
    def parameters(self):
        return self._par[:_npars(self._fam)].copy()

    @parameters.setter
    def parameters(self, p):
        k = _npars(self._fam)
        p = np.zeros(k) if p is None else np.asarray(p, dtype=np.float64).ravel()
        if p.size != k:
            raise ValueError(f"{self.family} has {k} parameter(s), got {p.size}")
        self._par = np.zeros(3)
        self._par[:p.size] = p

    @property
    def parameters_lower_bounds(self):
        return par_bounds(self._fam)[0]

    @property
    def parameters_upper_bounds(self):
        return par_bounds(self._fam)[1]

    def _aux(self):
        if self._fam not in _NP:
            return 1, np.zeros(1)
        if self.aux is None:
            raise ValueError(f"{self.family}: no estimate yet; call fit or select first")
        return NAUX, self.aux

    def _eval(self, u, what, num_threads):
        u = _data(u, 2)
        out = np.empty(len(u))
        if self._fam in _NP or what == 6:
            na, aux = self._aux()
            lib.fvc_bicop_eval(self._fam, self.rotation, self._par, na, aux, len(u), u, what, out)
        else:
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
        return self._eval(u, 6, None)

    def loglik(self, u=None, num_threads=None):
        """Log-likelihood on u, or (u=None) the one of the last fit / select (weighted if
        weights were given)."""
        if u is None:
            if self._loglik is None:
                raise ValueError("no fit yet: pass data")
            return self._loglik
        return float(self.logpdf(u, num_threads).sum())

    def aic(self, u=None):
        """-2 loglik + 2 npars (u=None: the last fit)."""
        return -2.0 * self.loglik(u) + 2.0 * self.npars

    def bic(self, u=None):
        """-2 loglik + log(n) npars (u=None: the last fit and its n)."""
        n = self._nobs if u is None else len(u)
        return -2.0 * self.loglik(u) + np.log(n) * self.npars

    def _fit(self, u, sel, controls, weights):
        u = _data(u, 2)
        usew, w = _weights(weights, len(u))
        fs, ci, cr = controls
        if not sel:
            fs = np.array([self._fam], np.int32)
        na = NAUX if np.isin(fs, _NP).any() else 1
        fam, rot = ct.c_int(self._fam), ct.c_int(self.rotation)
        par, ll, aux, npe = np.zeros(3), ct.c_double(), np.zeros(na), ct.c_double()
        lib.fvc_bicop_fit(len(u), u, usew, w, int(sel), fs.size, fs, ci, cr, ct.byref(fam), ct.byref(rot), par,
                          ct.byref(ll), na, aux, ct.byref(npe))
        self._fam, self.rotation, self._par, self._loglik = fam.value, rot.value, par, ll.value
        self.aux = aux if self._fam in _NP else None
        self._nobs = len(u)
        return self

    def fit(self, u, method="mle", fast_mle=False, tau_mode="exact", num_threads=None, *, weights=None):
        """Estimate the parameters (or, for tll / bernstein, the estimate `aux`) of the current
        family and rotation; `weights` (n,): observation weights."""
        return self._fit(u, False, _controls(method=method, fast_mle=fast_mle, tau_mode=tau_mode,
                                             num_threads=num_threads), weights)

    def select(self, u, family_set=None, criterion="bic", method="mle", preselect=True, allow_rotations=True,
               fast_mle=False, tau_mode="exact", psi0=0.9, num_threads=None, *, weights=None, **aliases):
        """Select family, rotation and parameters among family_set (names, codes or group names,
        see family_group; default: codes 0-22); `weights` (n,): observation weights."""
        opts = dict(family_set=family_set, criterion=criterion, method=method, preselect=preselect,
                    allow_rotations=allow_rotations, fast_mle=fast_mle, tau_mode=tau_mode, psi0=psi0,
                    num_threads=num_threads)
        opts.update(_aliases(aliases, opts))
        return self._fit(u, True, _controls(**opts), weights)

    def simulate(self, n, seed=None, num_threads=None):
        """n draws (n, 2): u2 = hinv1(u1, w) for independent uniforms u1, w."""
        w = np.random.default_rng(seed).random((n, 2))
        return np.column_stack([w[:, 0], self.hinv1(w, num_threads)])

    def __repr__(self):
        if self._fam in _NP:
            return f"Bicop(family={self.family!r}, npars={self.npars:.2f})"
        return f"Bicop(family={self.family!r}, rotation={self.rotation}, parameters={self.parameters.tolist()})"


# --------------------------------------------------------------------------- vine copula
class Vinecop:
    """Vine copula in natural order (as fv_vinecop).

    order (d,): variable labels (1-based); st, fam, rot (d, d) and par (d, d, 3) are indexed
    [tree, edge] (0-based): edge e of tree t joins order[e] and order[st[t, e] - 1].
    aux (d, d, naux): tll / bernstein estimates (naux = NAUX when the vine has such edges, else 1).
    disc (d,): discrete variables (original labels - 1).  Non-simplified edges set by relax:
    nrel, erel (nrel, 2) = (edge, tree) 1-based, wa (nrel, d), wk (nrel, NK), pk (nrel, NK, 3).
    """

    def __init__(self, d=0, var_types=None):
        self.d, self.trunc = int(d), 0
        self.order = np.arange(1, d + 1, dtype=np.int32)
        self.st, self.fam, self.rot = (np.zeros((d, d), np.int32) for _ in range(3))
        self.par, self.edge_ll = np.zeros((d, d, 3)), np.zeros((d, d))
        self.aux = np.zeros((d, d, 1))
        self.disc = np.zeros(d, bool) if var_types is None else _types(var_types, d)
        self.threshold = 0.0
        self.unrelax()
        self._fitted = False
        self._nobs = None

    # structure --------------------------------------------------------------
    def _set_structure(self, kind, order, trunc_lvl):
        order = np.ascontiguousarray(order, np.int32)
        d = order.size
        disc = self.disc if self.d == d else None
        self.__init__(d)
        if disc is not None:
            self.disc = disc
        st, tr = np.zeros((d, d), np.int32), ct.c_int()
        lib.fvc_vine_structure(kind, d, order, d - 1 if trunc_lvl is None else trunc_lvl, st, ct.byref(tr))
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
    def var_types(self):
        """'c' / 'd' per variable (pyvinecopulib var_types)."""
        return ["d" if x else "c" for x in self.disc]

    @var_types.setter
    def var_types(self, value):
        self.disc = _types(value, self.d)

    @property
    def matrix(self):
        """R-vine matrix in original labels (fv_vinecop%matrix)."""
        m = np.zeros((self.d, self.d), np.int32)
        lib.fvc_vine_matrix(self.d, self._i(self.order), self.trunc, self._i(self.st), m)
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
        """Number of parameters (fv_vinecop%npars_eff): effective degrees of freedom for tll /
        bernstein edges (then a float) and NK - 1 extra parameter sets per relaxed edge."""
        k, np_edge = 0.0, False
        for t, e in self._edges():
            f = int(self.fam[t, e])
            if f in _NP:
                k += float(self.aux[t, e, 1]) if self.aux.shape[2] > 1 else 0.0
                np_edge = True
            else:
                k += _npars(f)
        for e, t in self.erel:
            k += (NK - 1) * _npars(int(self.fam[t - 1, e - 1]))
        return k if np_edge else int(k)

    @property
    def nrel(self):
        """Number of non-simplified (relaxed) edges."""
        return len(self.erel)

    @property
    def relaxed_edges(self):
        """(tree, edge) 0-based of the relaxed edges."""
        return [(int(t) - 1, int(e) - 1) for e, t in self.erel]

    def get_pair_copula(self, tree, edge):
        f = int(self.fam[tree, edge])
        if f in _NP:
            if self.aux.shape[2] != NAUX:
                return Bicop(f)
            return Bicop(f, aux=self.aux[tree, edge])
        return Bicop(f, int(self.rot[tree, edge]), self.par[tree, edge, :_npars(f)])

    def set_pair_copula(self, tree, edge, bicop):
        if bicop.is_nonparametric:
            if bicop.aux is None:
                raise ValueError("set_pair_copula: the tll / bernstein copula has no estimate (fit it first)")
            if self.aux.shape[2] != NAUX:
                self.aux = np.zeros((self.d, self.d, NAUX))
            self.aux[tree, edge] = bicop.aux
        elif self.aux.shape[2] > 1:
            self.aux[tree, edge] = 0.0
        self.fam[tree, edge], self.rot[tree, edge] = bicop.family_code, bicop.rotation
        self.par[tree, edge] = bicop._par

    def rel_params(self, u):
        """Parameters (n, nrel, 3) of the relaxed edges at the samples u (fv_vinecop%rel_params):
        row i holds, for each edge of relaxed_edges, the pair-copula parameters at u[i]."""
        x, _ = _split(u, self.disc)
        out = np.zeros((len(x), max(self.nrel, 1), 3))
        a = self._args()
        lib.fvc_vine_rel_params(self.d, self.nrel, a[12], a[13], a[14], len(x), x, out)
        return out[:, :self.nrel].copy()

    def unrelax(self):
        """Delete the non-simplified edges (fv_vinecop%unrelax; the parameters relax wrote stay)."""
        self.erel = np.zeros((0, 2), np.int32)
        self.wa, self.wk, self.pk = np.zeros((0, self.d)), np.zeros((0, NK)), np.zeros((0, NK, 3))
        return self

    # selection ----------------------------------------------------------------
    def select(self, u, family_set=None, criterion="bic", method="mle", preselect=True, allow_rotations=True,
               fast_mle=False, tau_mode="exact", tree_criterion="tau", tree_algorithm="prim", trunc_lvl=None,
               threshold=0.0, psi0=0.9, fixed=None, num_threads=None, *, weights=None, var_types=None,
               select_trunc_lvl=False, select_threshold=False, **aliases):
        """Select structure (unless fixed), families and parameters (fv_vinecop%select).

        fixed defaults to True when a structure is set (set_dvine / set_cvine).
        weights (n,): observation weights.  var_types: 'c' / 'd' per variable (default: the
        vine's, all continuous for a new vine); with k discrete variables u is (n, d + k), F(x)
        of all variables followed by F(x-) of the discrete ones, or (n, 2d).
        select_trunc_lvl / select_threshold: truncation level / threshold chosen by mBICv.
        """
        opts = dict(family_set=family_set, criterion=criterion, method=method, preselect=preselect,
                    allow_rotations=allow_rotations, fast_mle=fast_mle, tau_mode=tau_mode,
                    tree_criterion=tree_criterion, tree_algorithm=tree_algorithm, trunc_lvl=trunc_lvl,
                    threshold=threshold, psi0=psi0, num_threads=num_threads, select_trunc_lvl=select_trunc_lvl,
                    select_threshold=select_threshold)
        opts.update(_aliases(aliases, opts))
        if var_types is not None:
            disc = _types(var_types)
        elif self.d > 0 and self.disc.any():
            disc = self.disc
        else:
            disc = np.zeros(np.shape(u)[1], bool)
        x, xm = _split(u, disc)
        n, d = x.shape
        if self.d != d:
            self.__init__(d)
        fx = self.trunc > 0 if fixed is None else bool(fixed)
        fs, ci, cr = _controls(**opts)
        usew, w = _weights(weights, n)
        na = NAUX if np.isin(fs, _NP).any() else 1
        order, st, trunc = self._i(self.order).copy(), self._i(self.st).copy(), ct.c_int(self.trunc)
        fam, rot = np.zeros((d, d), np.int32), np.zeros((d, d), np.int32)
        par, ell, aux = np.zeros((d, d, 3)), np.zeros((d, d)), np.zeros((d, d, na))
        ll, npe, thr, info = ct.c_double(), ct.c_double(), ct.c_double(), ct.c_int()
        di = disc.astype(np.int32)
        _interruptible(lambda: lib.fvc_vine_select(
            n, d, x, xm, di, usew, w, int(fx), fs.size, fs, ci, cr, ct.byref(trunc), order, st, fam, rot, par, ell,
            na, aux, ct.byref(ll), ct.byref(npe), ct.byref(thr), ct.byref(info)))
        if info.value != 0:
            raise RuntimeError(f"fvc_vine_select failed (info = {info.value})")
        self.order, self.st, self.fam, self.rot, self.par, self.edge_ll = order, st, fam, rot, par, ell
        self.aux, self.disc, self.trunc, self.threshold = aux, disc.copy(), trunc.value, thr.value
        self.unrelax()
        self._fitted, self._nobs, self._npe = True, n, npe.value
        return self

    @classmethod
    def from_data(cls, u, **kwargs):
        """Select a vine (structure, families, parameters) for pseudo-observations u (n, d), or
        (n, d + k) with var_types= for k discrete variables."""
        return cls().select(u, fixed=False, **kwargs)

    def relax(self, u, target=0, mode="test", method="mle", fast_mle=False, tau_mode="exact", num_threads=None, *,
              weights=None):
        """Non-simplified edges (fv_vinecop%relax): parameters of the eligible edges of trees >= 2
        vary with the first principal component of the normal scores of their conditioning set.

        target: variable label (1-based) whose edges are relaxed, 0 = any.  mode: 'test' (1, edges
        whose constant-conditional-correlation test rejects at 5 %) or 'all' (2).  Later trees are
        refitted.  pdf, loglik and pdf_grid use the relaxed edges; rosenblatt, inverse_rosenblatt
        and simulate use par (the middle-anchor values).  Not available with discrete variables.
        """
        if self.disc.any():
            raise ValueError("relax is not available for vines with discrete variables")
        x = _data(u, self.d)
        n, d = x.shape
        md = _key(_RELAX_MODE, mode, "relax mode")
        fs, ci, cr = _controls(method=method, fast_mle=fast_mle, tau_mode=tau_mode, num_threads=num_threads)
        usew, w = _weights(weights, n)
        na = self.aux.shape[2]
        fam, rot = self._i(self.fam).copy(), self._i(self.rot).copy()
        par, ell = np.ascontiguousarray(self.par, np.float64).copy(), np.ascontiguousarray(self.edge_ll).copy()
        aux = np.ascontiguousarray(self.aux, np.float64).copy()
        nmax = max(1, d * d)
        nrel, erel = ct.c_int(), np.zeros((nmax, 2), np.int32)
        wa, wk, pk = np.zeros((nmax, d)), np.zeros((nmax, NK)), np.zeros((nmax, NK, 3))
        ll, npe = ct.c_double(), ct.c_double()
        _interruptible(lambda: lib.fvc_vine_relax(
            n, d, x, usew, w, int(target), md, fs.size, fs, ci, cr, ct.byref(ct.c_int(self.trunc)),
            self._i(self.order), self._i(self.st), fam, rot, par, ell, na, aux, self.disc.astype(np.int32), nmax,
            ct.byref(nrel), erel, wa, wk, pk, ct.byref(ll), ct.byref(npe)))
        if nrel.value < 0:
            raise RuntimeError("fvc_vine_relax: too many relaxed edges")
        k = nrel.value
        self.fam, self.rot, self.par, self.edge_ll, self.aux = fam, rot, par, ell, aux
        self.erel, self.wa, self.wk, self.pk = erel[:k].copy(), wa[:k].copy(), wk[:k].copy(), pk[:k].copy()
        self._npe = npe.value
        return self

    # evaluation -----------------------------------------------------------------
    @staticmethod
    def _i(a):
        return np.ascontiguousarray(a, np.int32)

    def _args(self):
        """The vine arrays of fvc_vine_eval / fvc_vine_pdf_grid (d ... pk)."""
        na = self.aux.shape[2]
        if na != NAUX and any(int(self.fam[t, e]) in _NP for t, e in self._edges()):
            raise ValueError("tll / bernstein edge without its estimate (aux)")
        k = self.nrel

        def pad(a, shape):
            out = np.zeros((max(k, 1),) + shape, a.dtype)
            out[:k] = a
            return out

        return (self.d, self.trunc, self._i(self.order), self._i(self.st), self._i(self.fam), self._i(self.rot),
                np.ascontiguousarray(self.par, np.float64), na, np.ascontiguousarray(self.aux, np.float64),
                self.disc.astype(np.int32), k, pad(self._i(self.erel), (2,)), pad(self.wa, (self.d,)),
                pad(self.wk, (NK,)), pad(self.pk, (NK, 3)))

    def _eval(self, what, x, xm, seed, num_threads):
        out = np.empty_like(x)
        lib.fvc_vine_eval(*self._args(), what, len(x), x, xm, int(seed), _nt(num_threads), out)
        return out

    def pdf(self, u, num_threads=None):
        """Vine density at u (n, d); with k discrete variables u is (n, d + k) or (n, 2d)."""
        x, xm = _split(u, self.disc)
        return self._eval(0, x, xm, 0, num_threads)[:, 0].copy()

    def loglik(self, u=None, num_threads=None):
        """Log-likelihood on u, or (u=None) the sum of the edge log-likelihoods of the last select /
        relax (weighted if weights were given)."""
        if u is None:
            if not self._fitted:
                raise ValueError("no fit yet: pass data")
            return float(self.edge_ll[:self.trunc].sum())
        return float(np.log(self.pdf(u, num_threads)).sum())

    def aic(self, u=None):
        """-2 loglik + 2 npars (u=None: the last fit)."""
        return -2.0 * self.loglik(u) + 2.0 * self.npars

    def bic(self, u=None):
        """-2 loglik + log(n) npars (u=None: the last fit and its n)."""
        n = self._nobs if u is None else len(u)
        return -2.0 * self.loglik(u) + np.log(n) * self.npars

    def pdf_grid(self, ug, target, num_threads=None):
        """Density on a prediction grid (fv_vinecop%pdf_grid): ug (n, bins, d) holds, per sample,
        bins points that differ only in variable `target` (1-based label); returns (n, bins).
        With k discrete variables the last axis is d + k (or 2d), as for pdf."""
        ug = np.asarray(ug, dtype=np.float64)
        if ug.ndim != 3:
            raise ValueError(f"expected an (n, bins, d) array, got shape {ug.shape}")
        n, nb = ug.shape[:2]
        x, xm = _split(ug.reshape(n * nb, -1), self.disc)
        x, xm = (np.ascontiguousarray(a.reshape(n, nb, self.d)) for a in (x, xm))
        if not 1 <= int(target) <= self.d:
            raise ValueError("target: a variable label 1..d")
        out = np.empty((n, nb))
        lib.fvc_vine_pdf_grid(*self._args(), nb, n, x, xm, int(target), _nt(num_threads), out)
        return out

    def rosenblatt(self, u, num_threads=None, randomize_discrete=True, seed=1):
        """Rosenblatt transform (n, d).  Discrete variables: Brockwell's randomised transform
        (reproducible from seed != 0), or with randomize_discrete=False F(x | ...)."""
        x, xm = _split(u, self.disc)
        if randomize_discrete and self.disc.any() and int(seed) == 0:
            raise ValueError("seed must be non-zero for the randomised transform")
        return self._eval(1, x, xm, seed if randomize_discrete else 0, num_threads)

    def inverse_rosenblatt(self, w, num_threads=None):
        x = _data(w, self.d)
        return self._eval(2, x, x, 0, num_threads)

    def simulate(self, n, seed=None, num_threads=None):
        """n draws (n, d): inverse Rosenblatt transform of independent uniforms (continuous
        scale; apply the marginal quantile functions for discrete variables)."""
        return self.inverse_rosenblatt(np.random.default_rng(seed).random((n, self.d)), num_threads)

    def __str__(self):
        lines = [f"fortvine vine copula: d = {self.d}, trunc = {self.trunc}",
                 "order: " + " ".join(map(str, self.order))]
        if self.disc.any():
            lines.append("discrete: " + " ".join(str(j + 1) for j in np.flatnonzero(self.disc)))
        for t, e in self._edges():
            b = self.get_pair_copula(t, e)
            if b.is_nonparametric:
                lines.append(f"  tree {t + 1} edge {e + 1}: {b.family:<12} edf {b.npars:.2f}")
            else:
                lines.append(f"  tree {t + 1} edge {e + 1}: {b.family:<12} rot {b.rotation:3d} par "
                             + " ".join(f"{p:.4f}" for p in b.parameters))
        for r, (t, e) in enumerate(self.relaxed_edges):
            lines.append(f"  non-simplified tree {t + 1} edge {e + 1}: first parameter at the w anchors "
                         + " ".join(f"{p:.4f}" for p in self.pk[r, :, 0]))
        return "\n".join(lines)

    def __repr__(self):
        return f"<fortvine.Vinecop d={self.d} trunc={self.trunc} npars={self.npars}>"
