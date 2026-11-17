import os
import pickle

import numpy as np
import pytest

if int(os.environ.get("FORTVINE_THREADS", "0") or 0) not in range(1, 9):
    os.environ["FORTVINE_THREADS"] = "4"          # tests: at most 8 threads

import fortvine as fv  # noqa: E402


@pytest.fixture(scope="module")
def clayton_sample():
    return fv.Bicop("clayton", 0, [2.0]).simulate(2000, seed=1)


def test_families():
    fams = fv.families()
    assert len(fams) == 111 and len({f.code for f in fams}) == 111
    assert len({f.name.lower() for f in fams}) == 111            # names are unambiguous
    assert fams[3] == fv.Family(3, "clayton", 1)
    assert fv.Bicop("Student").npars == 2 and fv.Bicop("TAWN").npars == 3


def test_family_groups_and_npars():
    assert fv.family_group("nonparametric") == (24, 25) and len(fv.family_group("generic")) == 85
    assert list(fv._family_set(["parametric", "tll", "gaussian"])) == list(range(12)) + [24]
    for f in fv.families():                        # parameter counts agree with the bounds table
        if f.npars >= 0:
            lb, ub = fv.par_bounds(f.code)
            assert len(lb) == f.npars and np.all(ub > lb)
    assert [f.npars for f in fv.families()][24:26] == [-1, -1]
    lb, ub = fv.par_bounds("bb1")
    assert len(lb) == 2 and ub[1] > lb[1] >= 1.0


def test_gen_flip():
    u = np.random.default_rng(1).uniform(0.05, 0.95, (50, 2))
    for name, p in (("tawn-asymmixed", [0.5, 0.5]), ("tawn-type1", [3.0, 0.5])):
        code, par = fv.gen_flip(name, p)
        b, f = fv.Bicop(name, 0, p), fv.Bicop(code, 0, par)
        assert code != b.family_code or not np.allclose(par, p)
        np.testing.assert_allclose(f.pdf(u[:, ::-1]), b.pdf(u), rtol=1e-8)   # c'(u2, u1) = c(u1, u2)
        assert not np.allclose(b.pdf(u[:, ::-1]), b.pdf(u))                    # (not exchangeable)


@pytest.fixture(scope="module")
def mixture_sample():
    rng = np.random.default_rng(5)
    a = fv.Bicop("clayton", 0, [4.0]).simulate(1500, seed=1)
    b = fv.Bicop("gumbel", 90, [3.0]).simulate(1500, seed=2)
    return np.where(rng.random((1500, 1)) < 0.5, a, b)


def test_tll_bernstein(mixture_sample):
    u = mixture_sample
    m = 200
    g = (np.arange(m) + 0.5) / m
    grid = np.column_stack([np.repeat(g, m), np.tile(g, m)])
    for name in ("tll", "bernstein"):
        b = fv.Bicop(name).fit(u)
        assert b.aux.shape == (fv.NAUX,) and b.parameters.size == 0
        assert b.npars == pytest.approx(b.aux[1]) and b.npars > 1
        dens = b.pdf(grid)
        assert dens.min() >= 0 and dens.mean() == pytest.approx(1.0, abs=5e-3)
        assert b.loglik() == pytest.approx(b.loglik(u), rel=1e-8)
        h = b.hfunc1(u[:200])
        np.testing.assert_allclose(b.hinv1(np.column_stack([u[:200, 0], h])), u[:200, 1], atol=1e-6)
        assert b.aic() == pytest.approx(-2 * b.loglik() + 2 * b.npars)
        c = pickle.loads(pickle.dumps(b))
        np.testing.assert_array_equal(c.pdf(u[:20]), b.pdf(u[:20]))
    with pytest.raises(ValueError, match="fit or select"):
        fv.Bicop("tll").pdf(u[:3])
    # AIC prefers the nonparametric estimate on this mixture, BIC a parametric family
    assert fv.Bicop().select(u, family_set=["parametric", "nonparametric"], criterion="aic").family == "tll"
    assert not fv.Bicop().select(u, family_set=["parametric", "nonparametric"]).is_nonparametric


def test_weights_equal_duplication(clayton_sample):
    w = np.random.default_rng(2).integers(1, 3, len(clayton_sample))
    dup = np.repeat(clayton_sample, w, axis=0)
    a, b = fv.Bicop("clayton").fit(clayton_sample, weights=w), fv.Bicop("clayton").fit(dup)
    np.testing.assert_allclose(a.parameters, b.parameters, rtol=1e-10)
    assert a.loglik() == pytest.approx(b.loglik(), rel=1e-10)
    a, b = fv.Bicop().select(clayton_sample, weights=w), fv.Bicop().select(dup)
    assert (a.family, a.rotation) == (b.family, b.rotation)
    np.testing.assert_allclose(a.parameters, b.parameters, rtol=1e-8)


def test_hinv_inverts_hfunc():
    rng = np.random.default_rng(0)
    u = rng.uniform(0.02, 0.98, (500, 2))
    for b in (fv.Bicop("clayton", 90, [3.0]), fv.Bicop("bb1", 0, [0.5, 1.5]), fv.Bicop("tawn", 0, [0.4, 0.8, 2.0])):
        h = b.hfunc1(u)
        np.testing.assert_allclose(b.hinv1(np.column_stack([u[:, 0], h])), u[:, 1], atol=1e-6)
        h = b.hfunc2(u)
        np.testing.assert_allclose(b.hinv2(np.column_stack([h, u[:, 1]])), u[:, 0], atol=1e-6)


def test_pdf_integrates_to_one():
    m = 400
    g = (np.arange(m) + 0.5) / m
    grid = np.column_stack([np.repeat(g, m), np.tile(g, m)])
    for b in (fv.Bicop("frank", 0, [5.0]), fv.Bicop("gaussian", 0, [0.5])):
        assert b.pdf(grid).mean() == pytest.approx(1.0, abs=5e-3)


def test_cdf_and_loglik():
    b = fv.Bicop("gumbel", 0, [2.0])
    u = np.array([[0.3, 0.7], [0.5, 0.5]])
    assert np.all(np.diff(b.cdf(np.array([[0.2, 0.5], [0.4, 0.5], [0.8, 0.5]]))) > 0)
    assert b.cdf(np.array([[0.4, 1.0]]))[0] == pytest.approx(0.4, abs=1e-8)
    assert b.loglik(u) == pytest.approx(np.log(b.pdf(u)).sum())


def test_fit_recovers_parameter(clayton_sample):
    b = fv.Bicop("clayton").fit(clayton_sample)
    assert b.parameters[0] == pytest.approx(2.0, abs=0.15)
    assert b.loglik() == pytest.approx(b.loglik(clayton_sample), rel=1e-8)


def test_select_picks_true_family(clayton_sample):
    b = fv.Bicop().select(clayton_sample, family_set=["gaussian", "clayton", "gumbel", "frank", "joe"])
    assert (b.family, b.rotation) == ("clayton", 0)
    u = fv.Bicop("gumbel", 180, [2.5]).simulate(2000, seed=2)
    b = fv.Bicop().select(u)
    assert b.family in ("gumbel", "galambos", "husler_reiss") and b.rotation == 180


def test_simulate_seed():
    b = fv.Bicop("joe", 0, [2.0])
    np.testing.assert_array_equal(b.simulate(50, seed=7), b.simulate(50, seed=7))
    assert fv.kendall_tau(*b.simulate(3000, seed=3).T) == pytest.approx(0.36, abs=0.04)


@pytest.fixture(scope="module")
def vine_data():
    rng = np.random.default_rng(0)
    z = rng.standard_normal((1500, 4))
    z[:, 1] += z[:, 0]; z[:, 2] += z[:, 1]; z[:, 3] += 0.5 * z[:, 0]
    return fv.to_pseudo_obs(z)


def test_pseudo_obs(vine_data):
    assert vine_data.shape == (1500, 4)
    np.testing.assert_allclose(np.sort(vine_data[:, 0]), np.arange(1, 1501) / 1501)


def test_vine_loglik_equals_sum_log_pdf(vine_data):
    vc = fv.Vinecop.from_data(vine_data, family_set=["indep", "gaussian", "clayton", "gumbel", "frank"])
    assert vc.loglik() == pytest.approx(np.log(vc.pdf(vine_data)).sum(), rel=1e-8)
    assert vc.families[0] == ["gaussian"] * 3 and vc.npars == 3
    assert sorted(vc.matrix[::-1].diagonal()) == [1, 2, 3, 4]


def test_rosenblatt_round_trip(vine_data):
    vc = fv.Vinecop.from_data(vine_data, trunc_lvl=2)
    w = vc.rosenblatt(vine_data[:200])
    np.testing.assert_allclose(vc.inverse_rosenblatt(w), vine_data[:200], atol=1e-6)
    np.testing.assert_array_equal(vc.simulate(30, seed=5), vc.simulate(30, seed=5))


def test_dvine_fixed_structure(vine_data):
    vc = fv.Vinecop().set_dvine([1, 2, 3, 4]).select(vine_data)
    assert list(vc.order) in ([1, 2, 3, 4], [4, 3, 2, 1]) and vc.trunc == 3   # the same path
    vc2 = fv.Vinecop(4).set_dvine([1, 2, 3, 4])
    vc2.set_pair_copula(0, 0, fv.Bicop("clayton", 0, [2.0]))
    assert vc2.npars == 1 and vc2.get_pair_copula(0, 0).family == "clayton"


def test_generic_family():
    b = fv.Bicop("archimax(gamma|huslerreiss)")           # generic code 10203, two parameters
    assert b.family_code == 10203 and b.npars == 2
    b.fit(fv.Bicop("bb1", 0, [0.8, 1.6]).simulate(2000, seed=4))
    u = b.simulate(4000, seed=5)
    refit = fv.Bicop(10203).fit(u)
    np.testing.assert_allclose(refit.parameters, b.parameters, rtol=0.2)
    h = b.hfunc1(u[:200])
    np.testing.assert_allclose(b.hinv1(np.column_stack([u[:200, 0], h])), u[:200, 1], atol=1e-6)


def test_vine_weights(vine_data):
    w = np.full(len(vine_data), 2.0)
    a = fv.Vinecop.from_data(vine_data, family_set=["gaussian", "clayton", "frank"])
    b = fv.Vinecop.from_data(vine_data, family_set=["gaussian", "clayton", "frank"], weights=w)
    assert a.families == b.families and list(a.order) == list(b.order)
    assert b.loglik() == pytest.approx(2 * a.loglik(), rel=1e-6)


def test_controls(vine_data):
    a = fv.Vinecop.from_data(vine_data, select_trunc_lvl=True)
    assert 1 <= a.trunc <= 3
    rng = np.random.default_rng(4)                 # Markov chain: tree 1 carries everything
    z = np.zeros((2000, 6))
    z[:, 0] = rng.standard_normal(2000)
    for j in range(1, 6):
        z[:, j] = 0.7 * z[:, j - 1] + np.sqrt(0.51) * rng.standard_normal(2000)
    m = fv.Vinecop().set_dvine(range(1, 7)).select(fv.to_pseudo_obs(z), family_set=["indep", "gaussian"],
                                                    select_trunc_lvl=True, criterion="mbicv")
    assert m.trunc < 5
    t = fv.Vinecop.from_data(vine_data, select_threshold=True, criterion="mbicv")
    assert 0.0 < t.threshold <= 1.0
    # pyvinecopulib / fv_controls spellings and the other options are accepted
    k = fv.Vinecop.from_data(vine_data, tree_crit="rho", tree_algo="kruskal", parametric_method="itau",
                             selection_criterion="aic", preselect_families=False, tau_mode="approx",
                             allow_rotations=False, fast_mle=True, psi0=0.5, threshold=0.05, trunc_lvl=2)
    assert k.trunc == 2 and np.all(k.rot[:2] == 0)
    with pytest.raises(TypeError):
        fv.Vinecop.from_data(vine_data, no_such_option=1)


def test_discrete_vine_matches_pyvinecopulib():
    pv = pytest.importorskip("pyvinecopulib")
    stats = pytest.importorskip("scipy.stats")
    truth = fv.Vinecop().set_dvine([1, 2, 3])
    truth.fam[0, :2] = 1
    truth.par[0, :2, 0] = [0.7, 0.5]
    truth.fam[1, 0], truth.par[1, 0, 0] = 5, 2.0
    v = truth.simulate(2000, seed=1)
    k = stats.poisson.ppf(v[:, 1], 2.0)            # variable 2: a Poisson(2) count
    data = np.column_stack([v[:, 0], stats.poisson.cdf(k, 2.0), v[:, 2], stats.poisson.cdf(k - 1, 2.0)])
    vc = fv.Vinecop.from_data(data, var_types=["c", "d", "c"], family_set=["indep", "gaussian", "frank"])
    assert vc.var_types == ["c", "d", "c"] and vc.fam[0, 0] == 1 and vc.fam[0, 1] == 1
    assert vc.loglik() == pytest.approx(vc.loglik(data), rel=1e-8)
    pcs = [[pv.Bicop(family=getattr(pv.BicopFamily, b.family), parameters=b.parameters.reshape(-1, 1))
            if b.npars else pv.Bicop() for b in (vc.get_pair_copula(t, e) for e in range(2 - t))] for t in range(2)]
    ref = pv.Vinecop.from_structure(matrix=vc.matrix.astype(np.uint64), pair_copulas=pcs, var_types=vc.var_types)
    np.testing.assert_allclose(vc.pdf(data), ref.pdf(data), rtol=1e-6)
    full = np.column_stack([data[:, :3], data[:, :1], data[:, 3:], data[:, 2:3]])        # (n, 2d) form
    np.testing.assert_allclose(vc.pdf(full), vc.pdf(data), rtol=1e-12)
    r = vc.rosenblatt(data, seed=3)
    assert np.all(np.abs(r.mean(axis=0) - 0.5) < 0.03) and np.all(np.abs(r.var(axis=0) - 1 / 12) < 0.01)
    r0 = vc.rosenblatt(data, randomize_discrete=False)
    np.testing.assert_allclose(r0[:, [0, 2]], r[:, [0, 2]])
    assert np.all(r0[:, 1] >= r[:, 1] - 1e-12)        # F(x | .) >= the randomised value
    np.testing.assert_array_equal(pickle.loads(pickle.dumps(vc)).pdf(data[:10]), vc.pdf(data[:10]))
    with pytest.raises(ValueError, match="discrete"):
        vc.pdf(data[:, :3])


def test_nonparametric_vine(vine_data):
    vc = fv.Vinecop.from_data(vine_data[:800], family_set=["indep", "nonparametric"])
    assert vc.aux.shape == (4, 4, fv.NAUX) and isinstance(vc.npars, float)
    assert vc.npars == pytest.approx(vc._npe, rel=1e-12)
    assert vc.loglik() == pytest.approx(vc.loglik(vine_data[:800]), rel=1e-8)
    w = vc.rosenblatt(vine_data[:200])
    np.testing.assert_allclose(vc.inverse_rosenblatt(w), vine_data[:200], atol=1e-6)
    b = vc.get_pair_copula(0, 0)
    assert b.is_nonparametric and b.npars == pytest.approx(vc.aux[0, 0, 1])
    vc2 = fv.Vinecop(4).set_dvine([1, 2, 3, 4])
    vc2.set_pair_copula(0, 0, fv.Bicop("tll").fit(vine_data[:, :2]))
    assert vc2.aux.shape[2] == fv.NAUX and np.all(np.isfinite(vc2.pdf(vine_data[:20])))


@pytest.fixture(scope="module")
def nonsimplified():
    stats = pytest.importorskip("scipy.stats")
    rng = np.random.default_rng(1)
    n = 2000
    u2, w1, w3 = rng.random(n), rng.random(n), rng.random(n)
    x1 = fv.Bicop("clayton", 0, [2.0]).hinv1(np.column_stack([u2, w1]))
    rho = -0.6 + 1.4 * u2                          # the (1, 3 | 2) correlation varies with u2
    v3 = stats.norm.cdf(rho * stats.norm.ppf(w1) + np.sqrt(1 - rho ** 2) * stats.norm.ppf(w3))
    x3 = fv.Bicop("frank", 0, [5.0]).hinv1(np.column_stack([u2, v3]))
    return np.column_stack([x1, u2, x3])


def test_relax(nonsimplified):
    u = nonsimplified
    vc = fv.Vinecop().set_dvine([1, 2, 3]).select(u, family_set=["indep", "gaussian", "clayton", "frank"])
    ll0, k0 = vc.loglik(), vc.npars
    r = pickle.loads(pickle.dumps(vc)).relax(u, mode="test")
    assert r.nrel == 1 and r.relaxed_edges == [(1, 0)] and r.pk.shape == (1, fv.NK, 3)
    assert r.loglik() > ll0 + 50 and r.npars == k0 + (fv.NK - 1) == r._npe
    assert r.loglik() == pytest.approx(r.loglik(u), rel=1e-8)         # pdf uses the relaxed edges
    g = np.repeat(u[:5, None, :], 4, axis=1)
    g[:, :, 0] = [0.2, 0.4, 0.6, 0.8]          # variable 1 does not condition the relaxed edge: exact
    np.testing.assert_allclose(r.pdf_grid(g, 1), r.pdf(g.reshape(-1, 3)).reshape(5, 4), rtol=1e-10)
    q = pickle.loads(pickle.dumps(r))
    np.testing.assert_array_equal(q.pdf(u[:50]), r.pdf(u[:50]))
    assert "non-simplified tree 2 edge 1" in str(r)
    p = r.rel_params(u)                        # piecewise linear between the anchor parameters
    assert p.shape == (len(u), 1, 3)
    assert r.pk[0, :, 0].min() - 1e-12 <= p[:, 0, 0].min() and p[:, 0, 0].max() <= r.pk[0, :, 0].max() + 1e-12
    assert np.corrcoef(p[:, 0, 0], u[:, 1])[0, 1] > 0.9   # follows the conditioning variable u2
    assert r.unrelax().nrel == 0
    assert fv.Vinecop().set_dvine([1, 2, 3]).select(u).relax(u, mode="none").nrel == 0


def test_pdf_grid(vine_data):
    vc = fv.Vinecop.from_data(vine_data, family_set=["gaussian", "clayton", "frank"])
    g = np.repeat(vine_data[:30, None, :], 9, axis=1)
    g[:, :, 2] = np.linspace(0.05, 0.95, 9)
    np.testing.assert_allclose(vc.pdf_grid(g, 3), vc.pdf(g.reshape(-1, 4)).reshape(30, 9), rtol=1e-12)


def test_fitcache_and_abort(vine_data):
    fv.fitcache_enable(50000)
    try:
        fv.fitcache_clear()
        a = fv.Vinecop.from_data(vine_data)
        s1 = fv.fitcache_stats()
        b = fv.Vinecop.from_data(vine_data)
        s2 = fv.fitcache_stats()
    finally:
        fv.fitcache_disable()
    assert s1.misses > 0 and s1.hits == 0 and s2.hits >= s1.misses and s2.misses == s1.misses
    assert a.loglik() == b.loglik() and fv.fitcache_stats().used == 0
    fv.set_abort(True)
    try:
        aborted = fv.Vinecop.from_data(vine_data)
    finally:
        fv.set_abort(False)
    assert aborted.trunc == 0
    assert fv.Vinecop.from_data(vine_data).trunc == 3


def test_keyboard_interrupt_aborts():
    import _thread
    import threading
    import time
    u = fv.to_pseudo_obs(np.random.default_rng(0).standard_normal((20000, 10)).cumsum(axis=1))
    timer = threading.Timer(0.3, _thread.interrupt_main)
    t0 = time.perf_counter()
    timer.start()
    with pytest.raises(KeyboardInterrupt):
        fv.Vinecop.from_data(u, family_set=["all", "generic"])
    assert time.perf_counter() - t0 < 30
    assert fv.Vinecop.from_data(u[:500, :3]).trunc == 2           # the abort flag was cleared


def test_against_pyvinecopulib():
    pv = pytest.importorskip("pyvinecopulib")
    rng = np.random.default_rng(3)
    u = rng.uniform(0.05, 0.95, (200, 2))
    for name, par in (("clayton", [2.0]), ("gumbel", [1.7]), ("frank", [4.0]), ("joe", [2.2]), ("bb1", [0.6, 1.4])):
        ref = pv.Bicop(family=getattr(pv.BicopFamily, name), parameters=np.array(par).reshape(-1, 1))
        np.testing.assert_allclose(fv.Bicop(name, 0, par).pdf(u), ref.pdf(u), rtol=1e-6)
        np.testing.assert_allclose(fv.Bicop(name, 0, par).hfunc1(u), ref.hfunc1(u), rtol=1e-6, atol=1e-10)


def test_cite(capsys):
    text = fv.cite()
    assert "@unpublished{shojaeezadeh2026fortvine" in text and text in capsys.readouterr().out
