"""Cross-library benchmark: fortvine (Fortran, and its Python and R packages) vs pyvinecopulib and
rvinecopulib (vinecopulib C++) vs torchvine (PyTorch).

Usage:  python bench/run_bench.py [--fv build/fv_bench] [--rscript Rscript] [--quick]

The R part (fortvine R package, rvinecopulib) runs in bench/bench_r.R on the same data files;
it is skipped if Rscript is not found.

Writes bench/results/results.json and bench/results/tables.md (pasted into docs/benchmarks.md).
Every library gets the same data, the same 12 common families (indep, gaussian, student,
clayton, gumbel, frank, joe, bb1, bb6, bb7, bb8, tawn) and BIC.  Single-thread runs pin
every library to one thread.
"""
import argparse, json, os, platform, re, shutil, subprocess, sys, tempfile, time
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
    t0 = time.perf_counter()
    r = f()                                   # warm-up; a call over 2 s is timed by itself
    if time.perf_counter() - t0 > 2.0:
        return time.perf_counter() - t0, r
    t0 = time.perf_counter()
    for _ in range(reps):
        r = f()
    return (time.perf_counter() - t0) / reps, r


def unthrottle(handle=None):
    """Windows: opt a process (default: this one) out of power throttling. Windows 11 moves a
    windowless process to the efficiency cores after ~2 s of work (EcoQoS), which slows long
    timings by ~1.6x and would favour libraries timed in short subprocesses."""
    if sys.platform != "win32":
        return
    import ctypes as ct

    class State(ct.Structure):                # PROCESS_POWER_THROTTLING_STATE
        _fields_ = [("version", ct.c_ulong), ("control", ct.c_ulong), ("state", ct.c_ulong)]
    k = ct.WinDLL("kernel32")
    k.GetCurrentProcess.restype = ct.c_void_p
    k.SetProcessInformation.argtypes = [ct.c_void_p, ct.c_int, ct.c_void_p, ct.c_ulong]
    s = State(1, 1, 0)                         # control EXECUTION_SPEED, state off
    h = ct.c_void_p(int(handle)) if handle is not None else k.GetCurrentProcess()
    if not k.SetProcessInformation(h, 4, ct.byref(s), ct.sizeof(s)):   # 4 = ProcessPowerThrottling
        print("warning: could not disable power throttling")


def run(cmd, env=None):
    """Run cmd unthrottled; returns its stdout."""
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, env=env)
    unthrottle(getattr(p, "_handle", None))
    out, _ = p.communicate()
    if p.returncode:
        raise subprocess.CalledProcessError(p.returncode, cmd)
    return out


def run_fv(fv, args, threads=1):
    env = dict(os.environ, FORTVINE_THREADS=str(threads), OMP_NUM_THREADS=str(threads))
    out = run([fv] + [str(a) for a in args], env=env)
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
    ap.add_argument("--rscript", default=shutil.which("Rscript") or "Rscript")
    a = ap.parse_args()
    unthrottle()
    fv = a.fv + (".exe" if sys.platform == "win32" and not a.fv.endswith(".exe") else "")
    os.makedirs(OUT, exist_ok=True)
    tmp = tempfile.mkdtemp()

    import pyvinecopulib as pv
    import fortvine as fvp
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
    dims = [5, 10] if a.quick else [5, 10, 20]
    n_vine = 2000

    # data files (shared with bench_r.R), then the R part
    for name, code, par in FAMS:
        c = pv.Bicop(family=pvfam(name), parameters=np.array(par, float).reshape(-1, 1))
        np.savetxt(os.path.join(tmp, f"pair_{name}.csv"), c.simulate(n_fit, seeds=[11, 12]), delimiter=",",
                   fmt="%.15f")
    for d in dims:
        np.savetxt(os.path.join(tmp, f"vine_{d}.csv"), factor_data(n_vine, d, seed=d), delimiter=",", fmt="%.15f")
    rres = run_r(a.rscript, tmp, n_eval, a.quick)

    # ---------------------------------------------------------------- 1. evaluation
    for name, code, par in FAMS:
        rng = np.random.default_rng(1)
        u = rng.uniform(1e-6, 1 - 1e-6, size=(n_eval, 2))
        row = {"family": name}
        p3 = (par + [0.0, 0.0])[:3]
        r = run_fv(fv, ["eval", code, *p3, n_eval, 3])
        row["fortvine"] = {k: r[k + "_ns"] for k in ("logpdf", "hfunc1", "hinv1")}
        cf = fvp.Bicop(name, 0, par)
        row["fortvine_py"] = {"logpdf": 1e9 * timed(lambda: cf.logpdf(u, num_threads=1), 3)[0] / n_eval,
                              "hfunc1": 1e9 * timed(lambda: cf.hfunc1(u, num_threads=1), 3)[0] / n_eval,
                              "hinv1": 1e9 * timed(lambda: cf.hinv1(u, num_threads=1), 3)[0] / n_eval}
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
        fn = os.path.join(tmp, f"pair_{name}.csv")
        u = np.loadtxt(fn, delimiter=",")
        reps = 3 if a.quick else 10
        r = run_fv(fv, ["fit", fn, code, reps])
        rowf = {"family": name, "fortvine": {"ms": r["ms"], "loglik": r["loglik"]}}
        t, b = timed(lambda: fvp.Bicop(name).fit(u, num_threads=1), reps)
        rowf["fortvine_py"] = {"ms": 1e3 * t, "loglik": b.loglik()}
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
        t, b = timed(lambda: fvp.Bicop().select(u, family_set=COMMON, criterion="bic", num_threads=1), 1)
        rows["fortvine_py"] = {"ms": 1e3 * t, "family": b.family, "loglik": b.loglik()}
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
    for d in dims:
        fn = os.path.join(tmp, f"vine_{d}.csv")
        u = np.loadtxt(fn, delimiter=",")
        for thr in (1, 8):
            row = {"d": d, "n": n_vine, "threads": thr}
            r = run_fv(fv, ["vine", fn, thr, 1], threads=thr)
            row["fortvine"] = {k: r[k] for k in ("select_s", "loglik", "npars", "pdf_s", "simulate_s")}
            t, v = timed(lambda: fvp.Vinecop.from_data(u, family_set=COMMON, criterion="bic", num_threads=thr), 1)
            row["fortvine_py"] = {"select_s": t, "loglik": v.loglik(u, num_threads=thr), "npars": float(v.npars),
                                  "pdf_s": timed(lambda: v.pdf(u, num_threads=thr), 3)[0],
                                  "simulate_s": timed(lambda: v.simulate(n_vine, num_threads=thr), 3)[0]}
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

    if rres:                                  # merge the R results row by row (same order)
        res["machine"]["R"] = rres["versions"]
        for sec in ("eval", "fit", "select", "vine"):
            for row, rr in zip(res[sec], rres[sec]):
                row.update({k: v for k, v in rr.items() if k in ("fortvine_r", "rvinecopulib")})
    with open(os.path.join(OUT, "results.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    write_tables(res)


def run_r(rscript, tmp, n_eval, quick):
    """bench_r.R on the files in tmp; {} if R or the packages are missing."""
    json.dump({"fams": [{"name": n, "par": p} for n, _, p in FAMS], "common": COMMON,
               "dims": [5, 10] if quick else [5, 10, 20]}, open(os.path.join(tmp, "spec.json"), "w"))
    try:
        args = [os.path.join(HERE, "bench_r.R"), tmp, str(n_eval), "1" if quick else "0"]
        # Windows: Rscript.exe starts Rterm.exe as a child, which would stay throttled; run Rterm itself
        rterm = os.path.join(os.path.dirname(rscript), "x64", "Rterm.exe")
        if sys.platform == "win32" and os.path.isfile(rterm):
            run([rterm, "--no-save", "--no-restore", "--slave", "-f", args[0], "--args", *args[1:]])
        else:
            run([rscript, *args])
        return json.load(open(os.path.join(tmp, "r.json")))
    except (OSError, subprocess.CalledProcessError) as e:
        print("R part skipped:", e)
        return {}


LIBS = ["fortvine", "fortvine_py", "fortvine_r", "pyvinecopulib", "rvinecopulib", "torchvine"]
HEAD = {"fortvine": "fortvine (Fortran)", "fortvine_py": "fortvine (Python)", "fortvine_r": "fortvine (R)",
        "pyvinecopulib": "pyvinecopulib", "rvinecopulib": "rvinecopulib", "torchvine": "torchvine"}
NAMES = {"gaussian": "Gaussian", "student": "Student $t$", "clayton": "Clayton", "gumbel": "Gumbel",
         "frank": "Frank", "joe": "Joe", "bb1": "BB1", "bb6": "BB6", "bb7": "BB7", "bb8": "BB8", "tawn": "Tawn"}


def _v(row, lib, key):
    v = row.get(lib)
    return v.get(key) if isinstance(v, dict) else None


def _f(v, fmt):
    if v is None or isinstance(v, str):
        return "—"
    s = fmt % v
    return s.lstrip("+-") if float(s) == 0 else s          # no signed zeros


def write_tables(res):
    """Paper-style tables: one measure per table, one column per library, math in the headers."""
    L = []
    hdr = "| | " + " | ".join(HEAD[l] for l in LIBS) + " |\n|---|" + "---:|" * len(LIBS)

    def table(title, rows, label, fmt):
        L.append(f"**{title}**\n")
        L.append(hdr)
        for r in rows:
            L.append(f"| {label(r)} | " + " | ".join(_f(fmt[0](r, l), fmt[1]) for l in LIBS) + " |")
        L.append("")

    fam = lambda r: NAMES.get(r.get("family", r.get("truth")), r.get("family", r.get("truth")))
    for k, (key, what) in enumerate([("logpdf", r"density $c(u,v)$"), ("hfunc1", r"$h$-function $h_1(v\mid u)$"),
                                     ("hinv1", r"inverse $h_1^{-1}(w\mid u)$")], start=1):
        table(f"Table {k}. Evaluation of the {what}: time per point [ns], one thread, $n = 10^6$.",
              res["eval"], fam, (lambda r, l, key=key: _v(r, l, key), "%.0f"))
    table("Table 4. Maximum-likelihood fit of the true family: time [ms], one thread, $n = 5000$.",
          res["fit"], fam, (lambda r, l: _v(r, l, "ms"), "%.2f"))

    def dll(r, l):          # log-likelihood relative to fortvine's maximum
        a, b = _v(r, l, "loglik"), _v(r, "fortvine", "loglik")
        return None if a is None or b is None else a - b
    L.append(r"**Table 5. Maximised log-likelihood $\hat\ell$ of the fit in Table 4, and the difference "
             r"$\Delta\ell = \hat\ell_{\mathrm{lib}} - \hat\ell_{\mathrm{fortvine}}$.**" + "\n")
    L.append(r"| | $\hat\ell$ fortvine | " + " | ".join(rf"$\Delta\ell$ {HEAD[l]}" for l in LIBS[1:]) + " |\n|---|"
             + "---:|" * len(LIBS))
    for r in res["fit"]:
        L.append(f"| {fam(r)} | {_f(_v(r, 'fortvine', 'loglik'), '%.2f')} | "
                 + " | ".join(_f(dll(r, l), '%+.2f') for l in LIBS[1:]) + " |")
    L.append("")
    table("Table 6. Selection among the 12 common families and their rotations (BIC): time [ms], one thread, "
          "$n = 5000$.", res["select"], fam, (lambda r, l: _v(r, l, "ms"), "%.1f"))
    L.append("**Table 7. Family chosen by BIC in Table 6 (true family in the first column).**\n")
    L.append(hdr)
    for r in res["select"]:
        L.append(f"| {fam(r)} | " + " | ".join(NAMES.get(_v(r, l, "family"), _v(r, l, "family") or "—") for l in LIBS) + " |")
    L.append("")
    vl = lambda r: f"$d = {r['d']}$, {r['threads']} thread" + ("s" if r["threads"] > 1 else "")
    table("Table 8. Regular-vine selection (Dissmann, BIC): time [s], $n = 2000$.", res["vine"], vl,
          (lambda r, l: _v(r, l, "select_s"), "%.3f"))
    table(r"Table 9. Log-likelihood $\hat\ell$ of the selected vine, $n = 2000$.", res["vine"], vl,
          (lambda r, l: _v(r, l, "loglik"), "%.1f"))
    table("Table 10. Vine density on the $n = 2000$ sample: time [ms].", res["vine"], vl,
          (lambda r, l: None if _v(r, l, "pdf_s") is None else 1e3 * _v(r, l, "pdf_s"), "%.2f"))
    table("Table 11. Simulation of $n = 2000$ draws from the selected vine: time [ms].", res["vine"], vl,
          (lambda r, l: None if _v(r, l, "simulate_s") is None else 1e3 * _v(r, l, "simulate_s"), "%.2f"))
    with open(os.path.join(OUT, "tables.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    if sys.argv[1:] == ["--tables"]:          # rebuild tables.md from results.json without rerunning
        write_tables(json.load(open(os.path.join(OUT, "results.json"))))
    else:
        main()
