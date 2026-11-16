"""Cross-library benchmark: fortvine vs pyvinecopulib (vinecopulib C++) vs torchvine (PyTorch).

Usage:  python bench/run_bench.py [--fv build/fv_bench] [--quick]

Writes bench/results/results.json and bench/results/tables.md (pasted into docs/benchmarks.md).
Every library gets the same data, the same 12 common families (indep, gaussian, student,
clayton, gumbel, frank, joe, bb1, bb6, bb7, bb8, tawn) and BIC.  Single-thread runs pin
every library to one thread.
"""
import argparse, json, os, platform, re, subprocess, sys, tempfile, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")

# family, fortvine code, parameters (Kendall tau ~ 0.5)
FAMS = [("gaussian", 1, [0.7071]), ("student", 2, [0.7071, 4.0]), ("clayton", 3, [2.0]),
        ("gumbel", 4, [2.0]), ("frank", 5, [5.74]), ("joe", 6, [2.86]), ("bb1", 7, [0.5, 1.5]),
        ("bb6", 8, [1.5, 1.6]), ("bb7", 9, [1.8, 1.2]), ("bb8", 10, [4.0, 0.8]),
        ("tawn", 11, [0.6, 0.9, 3.0])]
COMMON = ["indep", "gaussian", "student", "clayton", "gumbel", "frank", "joe", "bb1", "bb6", "bb7",
          "bb8", "tawn"]
COMMON_TV = [f for f in COMMON if f != "tawn"]   # torchvine has no Tawn family


def timed(f, reps=1):
    f()                                       # warm-up
    t0 = time.perf_counter()
    for _ in range(reps):
        r = f()
    return (time.perf_counter() - t0) / reps, r


def run_fv(fv, args, threads=1):
    env = dict(os.environ, FORTVINE_THREADS=str(threads), OMP_NUM_THREADS=str(threads))
    out = subprocess.run([fv] + [str(a) for a in args], capture_output=True, text=True, env=env, check=True).stdout
    kv = {}
    for k, v in re.findall(r"(\w+)=\s*(\S+)", out):
        try:
            kv[k] = float(v)
        except ValueError:
            kv[k] = v
    return kv


def factor_data(n, d, seed):
    """t-factor model with skewed loadings -> pseudo-observations (a neutral, non-Gaussian vine)."""
    rng = np.random.default_rng(seed)
    f = rng.standard_t(5, size=(n, 2))
    load = rng.uniform(0.3, 0.9, size=(2, d)) * np.array([[1.0], [0.6]])
    x = f @ load + 0.6 * rng.standard_t(4, size=(n, d)) + 0.3 * np.abs(f[:, :1]) ** 1.5
    return (np.argsort(np.argsort(x, axis=0), axis=0) + 1.0) / (n + 1.0)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser()
    ap.add_argument("--fv", default=os.path.join(HERE, "..", "build", "fv_bench"))
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    fv = a.fv + (".exe" if sys.platform == "win32" and not a.fv.endswith(".exe") else "")
    os.makedirs(OUT, exist_ok=True)
    tmp = tempfile.mkdtemp()

    import pyvinecopulib as pv
    try:
        import torch, torchvine as tv
        torch.set_num_threads(1)
    except Exception as e:                    # torchvine is optional
        print("torchvine not available:", e); tv = None

    res = {"machine": {"platform": platform.platform(), "processor": platform.processor(),
                       "cpus": os.cpu_count(), "python": sys.version.split()[0]}, "eval": [], "fit": [],
           "select": [], "vine": []}
    n_eval = 200_000 if a.quick else 1_000_000
    n_fit = 5000
    pvfam = lambda name: getattr(pv.BicopFamily, name)

    # ---------------------------------------------------------------- 1. evaluation
    for name, code, par in FAMS:
        rng = np.random.default_rng(1)
        u = rng.uniform(1e-6, 1 - 1e-6, size=(n_eval, 2))
        row = {"family": name}
        p3 = (par + [0.0, 0.0])[:3]
        r = run_fv(fv, ["eval", code, *p3, n_eval, 3])
        row["fortvine"] = {k: r[k + "_ns"] for k in ("logpdf", "hfunc1", "hinv1")}
        c = pv.Bicop(family=pvfam(name), parameters=np.array(par, float).reshape(-1, 1))
        row["pyvinecopulib"] = {"logpdf": 1e9 * timed(lambda: c.pdf(u), 3)[0] / n_eval,
                                "hfunc1": 1e9 * timed(lambda: c.hfunc1(u), 3)[0] / n_eval,
                                "hinv1": 1e9 * timed(lambda: c.hinv1(u), 3)[0] / n_eval}
        if tv is not None:
            try:
                ct = tv.Bicop(getattr(tv.BicopFamily, name), parameters=torch.tensor(par, dtype=torch.float64))
                ut = torch.tensor(u)
                row["torchvine"] = {"logpdf": 1e9 * timed(lambda: ct.pdf(ut), 3)[0] / n_eval,
                                    "hfunc1": 1e9 * timed(lambda: ct.hfunc1(ut), 3)[0] / n_eval,
                                    "hinv1": 1e9 * timed(lambda: ct.hinv1(ut), 3)[0] / n_eval}
            except Exception as e:
                row["torchvine"] = {"error": str(e)[:120]}
        res["eval"].append(row); print("eval", row)

    # ---------------------------------------------------------- 2./3. fit and select
    for name, code, par in FAMS:
        c = pv.Bicop(family=pvfam(name), parameters=np.array(par, float).reshape(-1, 1))
        u = c.simulate(n_fit, seeds=[11, 12])
        fn = os.path.join(tmp, f"pair_{name}.csv"); np.savetxt(fn, u, delimiter=",", fmt="%.15f")
        reps = 3 if a.quick else 10
        r = run_fv(fv, ["fit", fn, code, reps])
        rowf = {"family": name, "fortvine": {"ms": r["ms"], "loglik": r["loglik"]}}
        ctl = pv.FitControlsBicop(family_set=[pvfam(name)], num_threads=1)
        def pvfit():
            b = pv.Bicop(family=pvfam(name)); b.fit(u, controls=ctl); return b
        t, b = timed(pvfit, reps)
        rowf["pyvinecopulib"] = {"ms": 1e3 * t, "loglik": float(b.loglik(u))}
        if tv is not None:
            try:
                ut = torch.tensor(u)
                ctlt = tv.FitControlsBicop(family_set=[getattr(tv.BicopFamily, name)])
                def tvfit():
                    bt = tv.Bicop(getattr(tv.BicopFamily, name)); bt.fit(ut, controls=ctlt); return bt
                t, bt = timed(tvfit, 1)
                rowf["torchvine"] = {"ms": 1e3 * t, "loglik": float(bt.loglik(ut))}
            except Exception as e:
                rowf["torchvine"] = {"error": str(e)[:120]}
        res["fit"].append(rowf); print("fit", rowf)

        r = run_fv(fv, ["select", fn, 1 if a.quick else 3])
        rows = {"truth": name, "fortvine": {"ms": r["ms"], "family": r["family"], "loglik": r["loglik"]}}
        ctls = pv.FitControlsBicop(family_set=[pvfam(f) for f in COMMON], selection_criterion="bic", num_threads=1)
        def pvsel():
            b = pv.Bicop(); b.select(u, controls=ctls); return b
        t, b = timed(pvsel, 1)
        rows["pyvinecopulib"] = {"ms": 1e3 * t, "family": str(b.family).split(".")[-1], "loglik": float(b.loglik(u))}
        if tv is not None:
            try:
                ctlt = tv.FitControlsBicop(family_set=[getattr(tv.BicopFamily, f) for f in COMMON_TV],
                                           selection_criterion="bic")
                def tvsel():
                    bt = tv.Bicop(); bt.select(ut, controls=ctlt); return bt
                t, bt = timed(tvsel, 1)
                rows["torchvine"] = {"ms": 1e3 * t, "family": str(bt.family).split(".")[-1],
                                     "loglik": float(bt.loglik(ut))}
            except Exception as e:
                rows["torchvine"] = {"error": str(e)[:120]}
        res["select"].append(rows); print("select", rows)

    # ------------------------------------------------------------------- 4./5. vines
    n_vine = 2000
    for d in ([5, 10] if a.quick else [5, 10, 20]):
        u = factor_data(n_vine, d, seed=d)
        fn = os.path.join(tmp, f"vine_{d}.csv"); np.savetxt(fn, u, delimiter=",", fmt="%.15f")
        for thr in (1, 8):
            row = {"d": d, "n": n_vine, "threads": thr}
            r = run_fv(fv, ["vine", fn, thr, 1], threads=thr)
            row["fortvine"] = {k: r[k] for k in ("select_s", "loglik", "npars", "pdf_s", "simulate_s")}
            ctl = pv.FitControlsVinecop(family_set=[pvfam(f) for f in COMMON], selection_criterion="bic",
                                        num_threads=thr)
            t, v = timed(lambda: pv.Vinecop.from_data(u, controls=ctl), 1)
            tp, _ = timed(lambda: v.pdf(u, num_threads=thr), 1)
            tm, _ = timed(lambda: v.simulate(n_vine, num_threads=thr), 1)
            row["pyvinecopulib"] = {"select_s": t, "loglik": float(v.loglik(u)), "npars": float(v.npars),
                                    "pdf_s": tp, "simulate_s": tm}
            if tv is not None and thr == 1 and d <= 10:
                try:
                    ut = torch.tensor(u)
                    ctlt = tv.FitControlsVinecop(family_set=[getattr(tv.BicopFamily, f) for f in COMMON_TV],
                                                 selection_criterion="bic")
                    t, vt = timed(lambda: tv.Vinecop.from_data(ut, controls=ctlt), 1)
                    tp, _ = timed(lambda: vt.pdf(ut), 1)
                    tm, _ = timed(lambda: vt.simulate(n_vine), 1)
                    row["torchvine"] = {"select_s": t, "loglik": float(vt.loglik(ut)), "npars": float(vt.npars),
                                        "pdf_s": tp, "simulate_s": tm}
                except Exception as e:
                    row["torchvine"] = {"error": str(e)[:160]}
            res["vine"].append(row); print("vine", row)

    with open(os.path.join(OUT, "results.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    write_tables(res)


def write_tables(res):
    libs = ["fortvine", "pyvinecopulib", "torchvine"]
    def cell(row, lib, key, fmt):
        v = row.get(lib, {}).get(key)
        return "n/a" if v is None else (fmt % v if isinstance(v, (int, float)) else str(v))
    L = []
    L.append("### Pair-copula evaluation (ns per point, one thread, n = 10^6)\n")
    L.append("| family | " + " | ".join(f"{l} pdf | {l} h | {l} h⁻¹" for l in libs) + " |")
    L.append("|---" * (1 + 3 * len(libs)) + "|")
    for r in res["eval"]:
        L.append(f"| {r['family']} | " + " | ".join(" | ".join(cell(r, l, k, "%.0f") for k in ("logpdf", "hfunc1", "hinv1"))
                                                for l in libs) + " |")
    L.append("\n### Pair-copula MLE fit of a known family (ms, one thread, n = 5000)\n")
    L.append("| family | " + " | ".join(f"{l} ms | {l} loglik" for l in libs) + " |")
    L.append("|---" * (1 + 2 * len(libs)) + "|")
    for r in res["fit"]:
        L.append(f"| {r['family']} | " + " | ".join(f"{cell(r, l, 'ms', '%.2f')} | {cell(r, l, 'loglik', '%.1f')}"
                                                for l in libs) + " |")
    L.append("\n### Family selection among the 12 common families + rotations (ms, one thread, n = 5000)\n")
    L.append("| truth | " + " | ".join(f"{l} ms | {l} chosen" for l in libs) + " |")
    L.append("|---" * (1 + 2 * len(libs)) + "|")
    for r in res["select"]:
        L.append(f"| {r['truth']} | " + " | ".join(f"{cell(r, l, 'ms', '%.1f')} | {cell(r, l, 'family', '%s')}"
                                               for l in libs) + " |")
    L.append("\n### Vine copula selection, density and simulation (n = 2000)\n")
    L.append("| d | threads | " + " | ".join(f"{l} select s | {l} loglik | {l} pdf s | {l} simulate s" for l in libs) + " |")
    L.append("|---" * (2 + 4 * len(libs)) + "|")
    for r in res["vine"]:
        L.append(f"| {r['d']} | {r['threads']} | " + " | ".join(
            " | ".join(cell(r, l, k, f) for k, f in (("select_s", "%.3f"), ("loglik", "%.1f"), ("pdf_s", "%.4f"),
                                                      ("simulate_s", "%.4f"))) for l in libs) + " |")
    with open(os.path.join(OUT, "tables.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
