import numpy as np
import pytest

import fortvine as fv


@pytest.fixture(scope="module")
def clayton_sample():
    return fv.Bicop("clayton", 0, [2.0]).simulate(2000, seed=1)


def test_families():
    fams = fv.families()
    assert len(fams) == 111 and len({f.code for f in fams}) == 111
    assert len({f.name.lower() for f in fams}) == 111            # names are unambiguous
    assert fams[3] == fv.Family(3, "clayton", 1)
    assert fv.Bicop("Student").npars == 2 and fv.Bicop("TAWN").npars == 3


def test_tll_bernstein_rejected():
    for f in ("tll", "bernstein", 24, 25):
        with pytest.raises(ValueError, match="C interface"):
            fv.Bicop(f)


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
