"""Every figure of the documentation site, from the library itself.

    bash build.sh            # library
    gfortran ... tools/density_grid.f90 -> build/density_grid ; tools/vine_sample.f90 -> build/vine_sample
    python tools/figures.py  [families | hero | bench | vine]   (default: all)

Densities are drawn in normal scores (z = Phi^-1(u)), where the shape of a copula is easiest to read.
"""
import json, os, subprocess, sys, tempfile
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXE = ".exe" if sys.platform == "win32" else ""
GRID = os.path.join(ROOT, "build", "density_grid" + EXE)
VSAMPLE = os.path.join(ROOT, "build", "vine_sample" + EXE)
TMP = tempfile.mkdtemp()
ENV = dict(os.environ, OMP_NUM_THREADS="8")

plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"],
                     "mathtext.fontset": "stix", "font.size": 9, "axes.edgecolor": "#8A8FA3",
                     "axes.linewidth": 0.6, "xtick.color": "#5A5F73", "ytick.color": "#5A5F73",
                     "axes.labelcolor": "#1E2230", "text.color": "#1E2230"})
INK, MUTED, GRIDC = "#1E2230", "#5A5F73", "#E3E4EC"
BLUE, GREEN = "#1F5FA6", "#2E8B57"
LIBC = {"fortvine": "#1F5FA6", "pyvinecopulib": "#2E8B57", "torchvine": "#C7922B"}   # validated, light surface
CMAP = LinearSegmentedColormap.from_list("fv", ["#FFFFFF", "#E3EDF8", "#AFCBEA", "#5F95D0", "#1F5FA6", "#0D2C52"])
TAU, NU, RHO = r"$\tau$", r"$\nu$", r"$\rho$"

z = np.linspace(-3, 3, 121)
Z1, Z2 = np.meshgrid(z, z)
PHI = np.exp(-0.5 * (Z1**2 + Z2**2)) / (2 * np.pi)


def grids(specs):
    """Run density_grid on a list of specs -> [(par, tau, density grid)]."""
    lines = []
    for i, s in enumerate(specs):
        p = (s["p"] + [0, 0, 0])[:3]
        lines.append(f"{s['fam']} {s.get('rot', 0)} {p[0]} {p[1]} {p[2]} {s.get('k', 0)} {s.get('lo', 0)} "
                     f"{s.get('hi', 0)} {s.get('tau', 0)} g_{i}.txt")
    r = subprocess.run([GRID], input="\n".join(lines) + "\n", capture_output=True, text=True, cwd=TMP, env=ENV,
                       check=True)
    out = []
    for i, ln in enumerate(x for x in r.stdout.split("\n") if x.strip()):
        f = ln.split()
        out.append((np.array([float(x) for x in f[2:5]]), float(f[5]), np.loadtxt(os.path.join(TMP, f"g_{i}.txt"))))
    return out


def save(f, path, **kw):
    """savefig with a few retries: on Windows a scanner or preview pane can briefly lock a fresh PNG."""
    import time
    for k in range(5):
        try:
            return f.savefig(path, **kw)
        except OSError:
            if k == 4:
                raise
            time.sleep(1.0)


def panel(a, g, ticks=True):
    d = g * PHI
    lev = np.linspace(0, max(d.max(), 1e-9) * 1.001, 11)
    a.contourf(Z1, Z2, d, levels=lev[1:], cmap=CMAP, vmin=0)          # zero-density band left white
    a.contour(Z1, Z2, d, levels=lev[1:], colors="#0D2C52", linewidths=0.35, alpha=0.55)
    a.set_xlim(-3, 3); a.set_ylim(-3, 3); a.set_aspect("equal")
    for sp in a.spines.values():
        sp.set_visible(False)
    if ticks:
        a.set_xticks([-2, 0, 2]); a.set_yticks([-2, 0, 2]); a.tick_params(labelsize=8, length=0, pad=2)
    else:
        a.set_xticks([]); a.set_yticks([])


# ------------------------------------------------------------------------------------------ family pages
FIGS = {
    "elliptical.png": ([
        dict(name="Gaussian (1)", fam=1, p=[0.5], k=1, lo=0, hi=0.99, tau=0.5),
        dict(name=f"Student t (2), {NU} = 4", fam=2, p=[0.5, 4], k=1, lo=0, hi=0.99, tau=0.5),
        dict(name=f"Student t (2), {RHO} = 0", fam=2, p=[0.0, 3]),
    ], 3),
    "archimedean.png": ([
        dict(name="Clayton (3)", fam=3, p=[2], k=1, lo=0.01, hi=20, tau=0.5),
        dict(name="Gumbel (4)", fam=4, p=[2], k=1, lo=1, hi=20, tau=0.5),
        dict(name="Frank (5)", fam=5, p=[5], k=1, lo=0.1, hi=30, tau=0.5),
        dict(name="Joe (6)", fam=6, p=[2], k=1, lo=1, hi=20, tau=0.5),
        dict(name="BB1 (7)", fam=7, p=[1, 1.5], k=1, lo=0.01, hi=6, tau=0.5),
        dict(name="BB6 (8)", fam=8, p=[1.5, 1], k=2, lo=1, hi=6, tau=0.5),
        dict(name="BB7 (9)", fam=9, p=[1.5, 1], k=2, lo=0.01, hi=20, tau=0.5),
        dict(name="BB8 (10)", fam=10, p=[4, 0.5], k=2, lo=0.01, hi=1, tau=0.5),
    ], 4),
    "archimedean-new.png": ([
        dict(name="PS.Sibuya (11300)", fam=11300, p=[1, 1.5], k=1, lo=1, hi=10, tau=0.5),
        dict(name="gamma.Sibuya (11700)", fam=11700, p=[1, 1.5], k=1, lo=0.001, hi=15, tau=0.5),
        dict(name="gamma.logseries (11800)", fam=11800, p=[1, 5], k=1, lo=0.001, hi=15, tau=0.5),
        dict(name="logseries.PS (12300)", fam=12300, p=[1, 1.5], k=1, lo=0.001, hi=30, tau=0.5),
        dict(name="logseries.gamma (12400)", fam=12400, p=[1, 1], k=1, lo=0.001, hi=30, tau=0.5),
        dict(name="logseries.Sibuya (12500)", fam=12500, p=[1, 1.5], k=1, lo=0.001, hi=30, tau=0.5),
    ], 3),
    "extreme-value.png": ([
        dict(name="Gumbel (4)", fam=4, p=[2], k=1, lo=1, hi=20, tau=0.5),
        dict(name="Galambos (15)", fam=15, p=[1], k=1, lo=0.05, hi=20, tau=0.5),
        dict(name="Husler-Reiss (16)", fam=16, p=[1], k=1, lo=0.05, hi=20, tau=0.5),
        dict(name=f"t-EV (23), {NU} = 4", fam=23, p=[0.5, 4], k=1, lo=-0.9, hi=0.99, tau=0.5),
        dict(name="Tawn (11)", fam=11, p=[0.5, 1.0, 2], k=3, lo=1, hi=40, tau=0.4),
        dict(name="Tawn-type2 (10007)", fam=10007, p=[2, 0.5], k=1, lo=1, hi=20, tau=0.4),
        dict(name="asym. neg. logistic (10008)", fam=10008, p=[1.0, 0.5, 1], k=3, lo=0.05, hi=20, tau=0.4),
        dict(name="Tawn-asymmixed (10011)", fam=10011, p=[0.9, 1.0]),
    ], 4),
    "archimax.png": ([
        dict(name="BB4 (22)", fam=22, p=[1, 1], k=2, lo=0.05, hi=10, tau=0.5),
        dict(name="gamma | HR (10203)", fam=10203, p=[1, 1], k=2, lo=0.05, hi=10, tau=0.5),
        dict(name="gamma | asymlog. 1 (10206)", fam=10206, p=[1, 2, 0.4], k=2, lo=1, hi=20, tau=0.5),
        dict(name="Sibuya | neglog. (10302)", fam=10302, p=[1.5, 1], k=2, lo=0.05, hi=10, tau=0.5),
        dict(name="PS | asymmixed (10111)", fam=10111, p=[1.5, 1, 0]),
        dict(name="negbin | HR (10703)", fam=10703, p=[2, 0.5, 1], k=3, lo=0.05, hi=10, tau=0.5),
    ], 3),
    "khoudraji.png": ([
        dict(name="Khoudraji-Clayton (20200)", fam=20200, p=[3, 0.95, 0.5], k=1, lo=0.01, hi=20, tau=0.35),
        dict(name="Khoudraji-Joe (20300)", fam=20300, p=[3, 0.5, 0.95], k=1, lo=1, hi=20, tau=0.35),
        dict(name="Khoudraji-Frank (20400)", fam=20400, p=[5, 0.6, 1.0], k=1, lo=0.01, hi=35, tau=0.35),
    ], 3),
    "other.png": ([
        dict(name="Plackett (12)", fam=12, p=[5], k=1, lo=1, hi=500, tau=0.5),
        dict(name="FGM (13)", fam=13, p=[1.0]),
        dict(name="AMH (14)", fam=14, p=[0.95]),
    ], 3),
}


def families():
    out = os.path.join(ROOT, "docs", "families", "img")
    os.makedirs(out, exist_ok=True)
    for name, (specs, ncol) in FIGS.items():
        res = grids(specs)
        nrow = (len(specs) + ncol - 1) // ncol
        f, ax = plt.subplots(nrow, ncol, figsize=(1.9 * ncol, 2.2 * nrow), squeeze=False)
        for i, (s, (par, tau, g)) in enumerate(zip(specs, res)):
            a = ax[i // ncol][i % ncol]
            panel(a, g)
            ps = ", ".join(f"{v:.3g}" for v in par[:len(s["p"])])
            a.set_title(f"{s['name']}\n({ps}),  {TAU} = {tau:.2f}", fontsize=9.5, color=INK)
        for j in range(len(specs), nrow * ncol):
            ax[j // ncol][j % ncol].axis("off")
        f.tight_layout(pad=0.3, w_pad=0.6, h_pad=1.0)
        save(f, os.path.join(out, name), dpi=130, bbox_inches="tight", pad_inches=0.05, facecolor="white")
        plt.close(f)
        print(name, os.path.getsize(os.path.join(out, name)) // 1024, "KB")


# ------------------------------------------------------------------------------------------ home page strip
def hero():
    specs = [dict(name="Gaussian", fam=1, p=[0.5], k=1, lo=0, hi=0.99, tau=0.5),
             dict(name="Student t", fam=2, p=[0.5, 3], k=1, lo=0, hi=0.99, tau=0.5),
             dict(name="Clayton", fam=3, p=[2], k=1, lo=0.01, hi=20, tau=0.5),
             dict(name="Gumbel", fam=4, p=[2], k=1, lo=1, hi=20, tau=0.5),
             dict(name="Frank", fam=5, p=[5], k=1, lo=0.1, hi=30, tau=0.5),
             dict(name="Tawn", fam=11, p=[0.6, 1.0, 4.0]),
             dict(name="Archimax", fam=10206, p=[1, 2, 0.4], k=2, lo=1, hi=20, tau=0.5),
             dict(name="Khoudraji", fam=20300, p=[3, 0.5, 0.95], k=1, lo=1, hi=20, tau=0.35)]
    res = grids(specs)
    f, ax = plt.subplots(1, len(specs), figsize=(1.45 * len(specs), 1.75))
    for a, s, (_, tau, g) in zip(ax, specs, res):
        panel(a, g, ticks=False)
        a.set_xlabel(s["name"], fontsize=10.5, color=INK, labelpad=4)
    f.tight_layout(pad=0.2, w_pad=0.5)
    p = os.path.join(ROOT, "docs", "assets", "hero.png")
    save(f, p, dpi=170, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    plt.close(f)
    print("hero.png", os.path.getsize(p) // 1024, "KB")


# ------------------------------------------------------------------------------------------ benchmarks
def _libs(rows, key):
    return [l for l in LIBC if any(isinstance(r.get(l), dict) and key in r[l] for r in rows)]


def bench():
    res = json.load(open(os.path.join(ROOT, "bench", "results", "results.json")))
    out = os.path.join(ROOT, "docs", "assets")
    # 1. vine selection time against dimension, one panel per thread count (one axis each, same scale)
    f, ax = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
    for a, thr in zip(ax, (1, 8)):
        rows = [r for r in res["vine"] if r["threads"] == thr]
        for lib in _libs(rows, "select_s"):
            pts = [(r["d"], r[lib]["select_s"]) for r in rows if isinstance(r.get(lib), dict) and "select_s" in r[lib]]
            d, t = zip(*pts)
            a.plot(d, t, color=LIBC[lib], lw=2, marker="o", ms=5, mec="white", mew=1.2, zorder=3)
            a.annotate(f"{lib}", (d[-1], t[-1]), xytext=(6, 0), textcoords="offset points", va="center",
                       fontsize=8.5, color=INK)
        a.set_yscale("log"); a.set_xticks([5, 10, 20]); a.set_xlim(3.5, 26)
        a.grid(True, which="major", color=GRIDC, lw=0.6); a.set_axisbelow(True)
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
        a.set_title(f"{thr} thread" + ("s" if thr > 1 else ""), fontsize=10, color=INK)
        a.set_xlabel("dimension d  (n = 2000)")
    ax[0].set_ylabel("vine selection time [s]")
    f.tight_layout()
    save(f, os.path.join(out, "bench_vine.png"), dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(f)
    # 2. pair-copula fit and selection, one dot row per family (Cleveland dot plot, log scale)
    f, ax = plt.subplots(1, 2, figsize=(7.2, 3.6), sharey=True)
    fams = [r["family"] for r in res["fit"]]
    yy = np.arange(len(fams))[::-1]
    for a, (sec, key, lab) in zip(ax, (("fit", "ms", "MLE fit of the true family [ms]"),
                                       ("select", "ms", "selection among 12 families [ms]"))):
        rows = res[sec]
        for y in yy:
            a.axhline(y, color=GRIDC, lw=0.6, zorder=0)
        for lib in _libs(rows, key):
            x = [r[lib][key] if isinstance(r.get(lib), dict) and key in r[lib] else np.nan for r in rows]
            a.scatter(x, yy, s=34, color=LIBC[lib], edgecolor="white", linewidth=1.0, zorder=3, label=lib)
        a.set_xscale("log"); a.set_xlabel(lab)
        for sp in ("top", "right", "left"):
            a.spines[sp].set_visible(False)
        a.tick_params(axis="y", length=0)
    pretty = {"student": "Student $t$", "bb1": "BB1", "bb6": "BB6", "bb7": "BB7", "bb8": "BB8"}
    ax[0].set_yticks(yy); ax[0].set_yticklabels([pretty.get(f, f.capitalize()) for f in fams])
    h, l = ax[0].get_legend_handles_labels()
    f.legend(h, l, loc="upper center", ncol=3, frameon=False, fontsize=9, bbox_to_anchor=(0.55, 1.03))
    f.tight_layout(rect=(0, 0, 1, 0.94))
    save(f, os.path.join(out, "bench_pair.png"), dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(f)
    print("bench_vine.png, bench_pair.png")


# ------------------------------------------------------------------------------------------ vine pairs plot
def vine():
    from scipy.stats import norm, kendalltau        # noqa: E402  (only needed here)
    fn = os.path.join(TMP, "vine.csv")
    subprocess.run([VSAMPLE, fn], check=True, env=ENV)
    zz = norm.ppf(np.loadtxt(fn, delimiter=","))
    d = zz.shape[1]
    f, ax = plt.subplots(d, d, figsize=(6.4, 6.4))
    for i in range(d):
        for j in range(d):
            a = ax[i][j]
            a.set_xticks([]); a.set_yticks([])
            for sp in a.spines.values():
                sp.set_color(GRIDC)
            if i == j:
                a.hist(zz[:, i], bins=30, color=BLUE, alpha=0.85, edgecolor="white", linewidth=0.4)
                a.text(0.05, 0.85, f"$u_{i + 1}$", transform=a.transAxes, fontsize=11)
            elif i > j:
                a.scatter(zz[:, j], zz[:, i], s=2.5, color=BLUE, alpha=0.45, linewidths=0)
                a.set_xlim(-3.5, 3.5); a.set_ylim(-3.5, 3.5)
            else:
                t = kendalltau(zz[:, j], zz[:, i])[0]
                a.text(0.5, 0.5, f"{TAU} = {t:.2f}", transform=a.transAxes, ha="center", va="center",
                       fontsize=11 + 8 * abs(t), color=INK)
    f.tight_layout(pad=0.3, w_pad=0.25, h_pad=0.25)
    p = os.path.join(ROOT, "docs", "assets", "vine_pairs.png")
    save(f, p, dpi=140, bbox_inches="tight", facecolor="white")
    plt.close(f)
    print("vine_pairs.png", os.path.getsize(p) // 1024, "KB")


if __name__ == "__main__":
    os.makedirs(os.path.join(ROOT, "docs", "assets"), exist_ok=True)
    for job in (sys.argv[1:] or ["families", "hero", "bench", "vine"]):
        globals()[job]()
