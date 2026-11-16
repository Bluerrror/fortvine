"""Figures of the example pages (docs/examples/*.md), from the library itself.

    bash build.sh                                   # library (build/libfortvine.a, .dll/.so)
    python tools/example_figures.py [pair | vine | generic | tails | weights | interop]   (default: all)

Every number comes from fortvine: the helper programs tools/example_*.f90 re-run the example
programs with their fixed seeds and write the data; this script only plots.  A helper is compiled
on first use (or when its source is newer than the executable) with

    gfortran -O2 -ffree-line-length-none -fopenmp -I build tools/X.f90 build/libfortvine.a -static -o build/X.exe

(gfortran and, for `interop`, gcc must be on PATH; `interop` also loads build/libfortvine.dll|.so
through ctypes and runs example/python/fortvine_ctypes.py).  PNGs go to docs/examples/img/.
Densities are drawn in normal scores (z = Phi^-1(u)), as in tools/figures.py, whose style and
helpers (panel, save, CMAP, colours) are reused.
"""
import contextlib, ctypes as ct, io, os, runpy, shutil, subprocess, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figures as F                                       # noqa: E402  (rcParams, panel, save, CMAP, grid)
import matplotlib.pyplot as plt                           # noqa: E402
from matplotlib.colors import LinearSegmentedColormap      # noqa: E402
from scipy.stats import norm                              # noqa: E402

ROOT, EXE, TMP, ENV = F.ROOT, F.EXE, F.TMP, F.ENV
BUILD = os.path.join(ROOT, "build")
OUT = os.path.join(ROOT, "docs", "examples", "img")
INK, MUTED, GRIDC = "#1E2230", "#5A5F73", "#E3E4EC"
BLUE, GREEN, GOLD = "#1F5FA6", "#2E8B57", "#C7922B"          # primary, second, third
RED, GREY = "#A4493D", "#8A8FA3"                            # fifth series, reference lines
CMAP = LinearSegmentedColormap.from_list("fv", ["#FFFFFF", "#E3EDF8", "#AFCBEA", "#5F95D0", "#1F5FA6", "#0D2C52"])
F.CMAP = CMAP                                               # F.panel draws with this density ramp
Z1, Z2, PHI = F.Z1, F.Z2, F.PHI
WCMAP = LinearSegmentedColormap.from_list("fvw", ["#AFCBEA", "#5F95D0", "#1F5FA6", "#0D2C52"])
plt.rcParams.update({"lines.linewidth": 1.2, "legend.frameon": False, "legend.fontsize": 8.5,
                     "axes.titlesize": 9.5, "axes.labelsize": 9.5})


# ------------------------------------------------------------------------------------------ helpers
def helper(name):
    """Path of build/<name>.exe, compiled from tools/<name>.f90 if missing or out of date."""
    src, exe = os.path.join(ROOT, "tools", name + ".f90"), os.path.join(BUILD, name + EXE)
    if not os.path.exists(exe) or os.path.getmtime(exe) < os.path.getmtime(src):
        flags = ["-static"] if sys.platform == "win32" else []
        subprocess.run(["gfortran", "-O2", "-ffree-line-length-none", "-fopenmp", "-I", BUILD, src,
                        os.path.join(BUILD, "libfortvine.a")] + flags + ["-o", exe], check=True)
    return exe


def run(name):
    """Run a helper in the scratch directory; returns its stdout."""
    r = subprocess.run([helper(name)], cwd=TMP, env=ENV, capture_output=True, text=True, check=True)
    return r.stdout


def load(fn, **kw):
    return np.loadtxt(os.path.join(TMP, fn), **kw)


def blocks(fn, nb, ng=121):
    g = load(fn)
    return [g[k * ng:(k + 1) * ng] for k in range(nb)]


def style(a, grid=True):
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)
    a.tick_params(labelsize=8.5, length=2.5, width=0.6)
    if grid:
        a.grid(True, color=GRIDC, lw=0.6)
        a.set_axisbelow(True)


def tag(a, s):
    a.text(-0.02, 1.02, s, transform=a.transAxes, ha="right", va="bottom", fontsize=10, fontweight="bold",
           color=INK)


def finish(f, name, dpi=140):
    os.makedirs(OUT, exist_ok=True)
    p = os.path.join(OUT, name)
    F.save(f, p, dpi=dpi, bbox_inches="tight", pad_inches=0.05, facecolor="white")
    plt.close(f)
    print(f"{name:28s} {os.path.getsize(p) // 1024:4d} KB")


def lines(a, g, color, lw=0.8, nlev=8, ls="-"):
    """Density contour lines (normal-score grid) on axis a."""
    d = g * PHI
    lev = np.linspace(0, d.max(), nlev + 2)[1:-1]
    return a.contour(Z1, Z2, d, levels=lev, colors=color, linewidths=lw, linestyles=ls)


# ------------------------------------------------------------------------------------------ pair copula
def pair():
    out = run("example_pair")
    th = float(out.split()[1])
    u = load("pair_sample.txt")
    gt, gf = blocks("pair_grid.txt", 2)
    zz = norm.ppf(u)
    # 1. the sample in copula scale and in normal scores, fitted density contours on top
    f, ax = plt.subplots(1, 2, figsize=(7.0, 3.35))
    a = ax[0]
    a.scatter(u[:, 0], u[:, 1], s=2.2, color=BLUE, alpha=0.45, linewidths=0)
    a.set_xlim(0, 1); a.set_ylim(0, 1); a.set_aspect("equal")
    a.set_xlabel("$u_1$"); a.set_ylabel("$u_2$")
    a.set_title("sample, copula scale", color=INK)
    style(a)
    a = ax[1]
    a.scatter(zz[:, 0], zz[:, 1], s=2.2, color=BLUE, alpha=0.35, linewidths=0, label="sample, $n = 2000$")
    lines(a, gf, INK, lw=0.8)
    lines(a, gt, GOLD, lw=0.9, ls="--")
    a.plot([], [], color=INK, lw=0.8, label=rf"fitted Clayton, $\hat\theta = {th:.3f}$")
    a.plot([], [], color=GOLD, lw=0.9, ls="--", label=r"true Clayton, $\theta = 2$")
    a.set_xlim(-3.2, 3.2); a.set_ylim(-3.2, 3.2); a.set_aspect("equal")
    a.set_xlabel(r"$z_1 = \Phi^{-1}(u_1)$"); a.set_ylabel(r"$z_2 = \Phi^{-1}(u_2)$")
    a.set_title("normal scores with density contours", color=INK)
    a.legend(loc="upper left", fontsize=7.5, handlelength=1.6, borderaxespad=0.2, frameon=True, framealpha=0.9,
             facecolor="white", edgecolor="none", markerscale=3)
    style(a)
    tag(ax[0], "a"); tag(ax[1], "b")
    f.tight_layout(w_pad=2.0)
    finish(f, "pair_sample.png")

    # 2. h-function, its inverse and the conditional density for three values of u
    h, hi = load("pair_h.txt"), load("pair_hinv.txt")
    f, ax = plt.subplots(1, 2, figsize=(7.0, 2.9))
    for uc, col in zip((0.1, 0.5, 0.9), (BLUE, GREEN, GOLD)):
        m, mi = h[:, 0] == uc, hi[:, 0] == uc
        ax[0].plot(h[m, 1], h[m, 2], color=col, lw=1.4, label=f"$u = {uc}$")
        ax[0].plot(hi[mi, 2], hi[mi, 1], "o", ms=3.6, mfc="white", mec=col, mew=1.0)
        ax[1].plot(h[m, 1], h[m, 3], color=col, lw=1.4)
    ax[0].plot([], [], "o", ms=3.6, mfc="white", mec=INK, mew=1.0,
               label=r"$h_1^{-1}(w \mid u)$ at $w = 0.05, 0.10, \ldots, 0.95$")
    ax[0].set_xlabel("$v$"); ax[0].set_ylabel(r"$w = h_1(v \mid u)$")
    ax[0].set_title(r"$h$-function and its inverse", color=INK)
    ax[0].set_xlim(0, 1); ax[0].set_ylim(0, 1.0)
    ax[1].set_xlabel("$v$"); ax[1].set_ylabel(r"$c(u, v)$")
    ax[1].set_title(r"conditional density of $V$ given $U = u$", color=INK)
    ax[1].set_xlim(0, 1); ax[1].set_ylim(0, None)
    for a in ax:
        style(a)
    tag(ax[0], "a"); tag(ax[1], "b")
    hh, ll = ax[0].get_legend_handles_labels()
    f.legend(hh, ll, loc="lower center", ncol=4, fontsize=8.5, bbox_to_anchor=(0.5, -0.01), handlelength=1.6)
    f.tight_layout(w_pad=2.0, rect=(0, 0.07, 1, 1))
    finish(f, "pair_hfunc.png")


# ------------------------------------------------------------------------------------------ vine
FAMN = {0: "indep.", 1: "Gaussian", 2: "Student t", 3: "Clayton", 4: "Gumbel", 5: "Frank", 6: "Joe", 7: "BB1",
        8: "BB6", 9: "BB7", 10: "BB8", 11: "Tawn", 12: "Plackett", 15: "Galambos", 16: "Hüsler-Reiss"}


def read_vine(fn):
    edges = []
    for ln in open(os.path.join(TMP, fn)):
        f = ln.split()
        t, e, a, b, fam, rot = map(int, f[:6])
        nc = int(f[10])
        edges.append(dict(t=t, a=a, b=b, fam=fam, rot=rot, par=[float(x) for x in f[6:9]], tau=float(f[9]),
                          cond=tuple(sorted(int(x) for x in f[11:11 + nc]))))
    return edges


def node_label(e):
    s = f"{e['a']},{e['b']}"
    return s + (" | " + ",".join(map(str, e["cond"])) if e["cond"] else "")


def tree_layout(nodes, links):
    """x-positions for a small tree: the longest path on a line, other nodes hung below it."""
    adj = {n: set() for n in nodes}
    for p, q in links:
        adj[p].add(q); adj[q].add(p)

    def far(s):
        dist, prev, todo = {s: 0}, {s: None}, [s]
        for x in todo:
            for y in adj[x]:
                if y not in dist:
                    dist[y], prev[y] = dist[x] + 1, x
                    todo.append(y)
        end = max(dist, key=dist.get)
        path = [end]
        while prev[path[-1]] is not None:
            path.append(prev[path[-1]])
        return end, path
    a, _ = far(nodes[0])
    _, path = far(a)
    if len(path) > 1 and str(path[0]) > str(path[-1]):
        path = path[::-1]
    pos = {n: (float(i), 0.0) for i, n in enumerate(path)}
    rest = [n for n in nodes if n not in pos]
    while rest:
        for n in list(rest):
            for y in adj[n]:
                if y in pos:
                    pos[n] = (pos[y][0] + 0.5, -0.9); rest.remove(n); break
    return pos


def draw_tree(a, edges, t, ref=None, xspan=4.0):
    """Tree t of a vine: nodes = variables (t = 1) or edges of tree t-1; edge labels family and tau."""
    if t == 1:
        nodes = sorted({e["a"] for e in edges if e["t"] == 1} | {e["b"] for e in edges if e["t"] == 1})
        sets = {n: frozenset([n]) for n in nodes}
        labels = {n: str(n) for n in nodes}
    else:
        prev = [e for e in edges if e["t"] == t - 1]
        nodes = list(range(len(prev)))
        sets = {i: frozenset([e["a"], e["b"], *e["cond"]]) for i, e in enumerate(prev)}
        labels = {i: node_label(e) for i, e in enumerate(prev)}
    links = []
    for e in (e for e in edges if e["t"] == t):
        sa, sb = frozenset([e["a"], *e["cond"]]), frozenset([e["b"], *e["cond"]])
        p = next(n for n in nodes if sets[n] == sa); q = next(n for n in nodes if sets[n] == sb)
        links.append((p, q, e))
    pos = tree_layout(nodes, [(p, q) for p, q, _ in links])
    sc = xspan / max(1, max(x for x, _ in pos.values()))
    pos = {n: (x * sc, y) for n, (x, y) in pos.items()}
    for p, q, e in links:
        (x1, y1), (x2, y2) = pos[p], pos[q]
        indep = e["fam"] == 0
        new = ref is not None and t == 1 and frozenset([e["a"], e["b"]]) not in ref
        col = GREY if indep else (GOLD if new else INK)
        a.plot([x1, x2], [y1, y2], color=col, lw=0.9 if indep else 1.4, ls=":" if indep else "-", zorder=1)
        xm, ym = 0.5 * (x1 + x2), 0.5 * (y1 + y2)
        name = FAMN.get(e["fam"], str(e["fam"])) + (f" {e['rot']}°" if e["rot"] else "")
        a.text(xm, ym + 0.17, name, ha="center", va="bottom", fontsize=7.6, color=INK if not indep else MUTED)
        if not indep:
            a.text(xm, ym - 0.17, rf"$\tau = {e['tau']:.2f}$", ha="center", va="top", fontsize=7.4, color=MUTED)
    for n, (x, y) in pos.items():
        a.text(x, y, labels[n], ha="center", va="center", fontsize=8 if t == 1 else 7.2, color=INK, zorder=3,
               bbox=dict(boxstyle="round,pad=0.25" if t > 1 else "circle,pad=0.25", fc="white", ec=BLUE, lw=0.9))
    a.set_xlim(-0.55, xspan + 0.55); a.set_ylim(-0.75, 0.75)
    a.axis("off")


def vine():
    out = run("example_vine")
    ll = {ln.split()[0]: (float(ln.split()[1]), int(ln.split()[2])) for ln in out.strip().split("\n")}
    truth, sel = read_vine("vine_truth.txt"), read_vine("vine_select.txt")
    ref = {frozenset([e["a"], e["b"]]) for e in truth if e["t"] == 1}
    # 1. the true and the selected vine, tree by tree
    f, ax = plt.subplots(4, 2, figsize=(7.2, 4.6), gridspec_kw=dict(hspace=0.25, wspace=0.08))
    for c, (edges, title) in enumerate(((truth, f"true D-vine: {ll['truth'][1]} parameters, "
                                                 rf"$\ell = {ll['truth'][0]:.1f}$"),
                                        (sel, f"selected (Dissmann): {ll['select'][1]} parameters, "
                                              rf"$\ell = {ll['select'][0]:.1f}$"))):
        for t in range(1, 5):
            draw_tree(ax[t - 1][c], edges, t, ref=ref if c == 1 else None, xspan=4.0 - 0.0 * t)
            ax[t - 1][c].text(-0.55, 0.62, f"$T_{t}$", fontsize=9, color=MUTED, ha="left", va="center")
        ax[0][c].set_title(title, fontsize=9.5, color=INK, pad=6)
    h = [plt.Line2D([], [], color=INK, lw=1.4), plt.Line2D([], [], color=GOLD, lw=1.4),
         plt.Line2D([], [], color=GREY, lw=0.9, ls=":")]
    f.legend(h, ["pair copula", "tree-1 edge not in the true vine", "independence"], loc="lower center", ncol=3,
             fontsize=8.5, bbox_to_anchor=(0.5, -0.02))
    f.subplots_adjust(left=0.01, right=0.99, top=0.93, bottom=0.07)
    finish(f, "vine_trees.png")

    # 2. why: the empirical Kendall tau of every pair, and which pairs each first tree uses
    tm = load("vine_taumat.txt")
    d = tm.shape[0]
    selt1 = {frozenset([e["a"], e["b"]]) for e in sel if e["t"] == 1}
    f, ax = plt.subplots(1, 2, figsize=(7.0, 2.9), gridspec_kw=dict(width_ratios=[1, 1.35]))
    a = ax[0]
    m = np.where(np.eye(d) == 1, np.nan, tm)
    a.imshow(m, cmap=CMAP, vmin=0.25, vmax=0.6)
    for i in range(d):
        for j in range(d):
            if i != j:
                a.text(j, i, f"{tm[i, j]:.3f}", ha="center", va="center", fontsize=7.4,
                       color="white" if tm[i, j] > 0.47 else INK)
    for i in range(d):
        for j in range(i + 1, d):
            pr = frozenset([i + 1, j + 1])
            if pr in ref:
                for (r, c) in ((i, j), (j, i)):
                    a.add_patch(plt.Rectangle((c - 0.46, r - 0.46), 0.92, 0.92, fill=False, ec=GREEN, lw=1.6))
            if pr in selt1:
                for (r, c) in ((i, j), (j, i)):
                    a.add_patch(plt.Rectangle((c - 0.36, r - 0.36), 0.72, 0.72, fill=False, ec=GOLD, lw=1.4,
                                              ls="--"))
    a.set_xticks(range(d)); a.set_yticks(range(d))
    a.set_xticklabels([f"$u_{k + 1}$" for k in range(d)]); a.set_yticklabels([f"$u_{k + 1}$" for k in range(d)])
    a.tick_params(length=0, labelsize=8.5)
    for sp in a.spines.values():
        sp.set_visible(False)
    a.set_title(r"empirical Kendall $\hat\tau_{ij}$", color=INK)
    a = ax[1]
    pairs = sorted(((tm[i, j], i + 1, j + 1) for i in range(d) for j in range(i + 1, d)), reverse=True)
    yy = np.arange(len(pairs))[::-1]
    for y, (tv, i, j) in zip(yy, pairs):
        pr = frozenset([i, j])
        a.axhline(y, color=GRIDC, lw=0.6, zorder=0)
        a.scatter(tv, y, s=30, color=GREEN if pr in ref else "white", edgecolor=GREEN if pr in ref else GREY,
                  lw=1.0, zorder=3)
        if pr in selt1:
            a.scatter(tv, y, s=95, facecolor="none", edgecolor=GOLD, lw=1.3, ls="--", zorder=4)
    a.set_yticks(yy); a.set_yticklabels([f"({i},{j})" for _, i, j in pairs])
    a.set_xlabel(r"$\hat\tau_{ij}$")
    a.set_title("pairs ranked by $\\hat\\tau$: the spanning tree takes the top four", color=INK)
    for sp in ("top", "right", "left"):
        a.spines[sp].set_visible(False)
    a.tick_params(axis="y", length=0, labelsize=8.5); a.tick_params(axis="x", labelsize=8.5)
    a.set_xlim(0.34, 0.53)
    h = [plt.Line2D([], [], marker="o", color="none", mfc=GREEN, mec=GREEN, ms=6),
         plt.Line2D([], [], marker="o", color="none", mfc="none", mec=GOLD, ms=10, mew=1.3)]
    a.legend(h, ["in the true tree 1", "in the selected tree 1"], loc="lower right", fontsize=8)
    tag(ax[0], "a"); tag(ax[1], "b")
    f.tight_layout(w_pad=1.5)
    finish(f, "vine_tau.png")


# ------------------------------------------------------------------------------------------ generic families
def generic():
    out = run("example_generic")
    vals = {ln.split()[0]: [float(x) for x in ln.split()[1:]] for ln in out.strip().split("\n") if
            ln.startswith(("archimax", "khoudraji"))}
    pa, pk = vals["archimax_hat"], vals["khoudraji_hat"]
    g = blocks("gen_grid.txt", 4)
    titles = ["Archimax 10203, true\n" + r"$(\theta, \lambda) = (1, 1.5)$",
              "Archimax 10203, fitted\n" + rf"$(\hat\theta, \hat\lambda) = ({pa[0]:.3f}, {pa[1]:.3f})$",
              "Khoudraji-Clayton 20200, true\n" + r"$(\theta, a, b) = (2, 0.6, 0.95)$",
              "Khoudraji-Clayton 20200, fitted\n" + rf"$({pk[0]:.2f}, {pk[1]:.2f}, {pk[2]:.2f})$"]
    f, ax = plt.subplots(1, 4, figsize=(7.2, 2.25))
    for a, gg, ti in zip(ax, g, titles):
        F.panel(a, gg)
        a.plot([-3, 3], [-3, 3], color=GREY, lw=0.5, ls=":")
        a.set_title(ti, fontsize=8.6, color=INK)
    f.tight_layout(pad=0.3, w_pad=0.8)
    finish(f, "generic_density.png")

    p = load("gen_pickands.txt")
    w = p[:, 0]
    f, ax = plt.subplots(1, 2, figsize=(7.0, 2.9), sharey=True)
    for a in ax:
        a.fill_between(w, np.maximum(w, 1 - w), 1, color="#F1F3F8", lw=0, zorder=0)
        a.plot(w, np.maximum(w, 1 - w), color=GREY, lw=0.6)
        a.plot([0, 1], [1, 1], color=GREY, lw=0.6)
        a.set_xlim(0, 1); a.set_ylim(0.48, 1.02); a.set_xlabel("$w$")
        style(a)
    a = ax[0]
    a.plot(w, p[:, 1], color=BLUE, lw=1.6, label=r"from $C_{10203}$, $(\theta, \lambda) = (1, 1.5)$")
    a.plot(w, p[:, 2], color=GOLD, lw=1.1, ls="--", label=rf"from $C_{{10203}}$, fitted $\hat\lambda = {pa[1]:.3f}$")
    a.plot(w[::10], p[::10, 3], "o", ms=3.6, mfc="white", mec=INK, mew=0.9, label=r"Hüsler-Reiss (16), $\lambda = 1.5$")
    a.set_ylabel("$A(w)$")
    a.set_title("inner Pickands function of the Archimax copula", color=INK)
    a.legend(loc="upper center", fontsize=7.6, bbox_to_anchor=(0.5, 0.97), handlelength=1.6, labelspacing=0.3)
    a = ax[1]
    names = [r"Tawn-mixed 10004, $\theta = 0.8$", r"Tawn type 1 10006, $(3, 0.5)$", r"Tawn type 2 10007, $(3, 0.5)$",
             r"asym. neg. logistic 10008, $(0.6, 0.9, 2)$", r"asym. neg. log. type 1 10009, $(2, 0.5)$"]
    for k, (nm, col, ls) in enumerate(zip(names, (BLUE, GREEN, GOLD, INK, RED), ("-", "-", "--", "-", "-"))):
        a.plot(w, p[:, 4 + k], color=col, lw=1.3, ls=ls, label=nm)
    a.set_title("the five extreme-value codes listed by the program", color=INK)
    a.legend(loc="lower center", fontsize=7.2, handlelength=1.6, labelspacing=0.25)
    a.text(0.5, 0.995, "independence $A = 1$", ha="center", va="top", fontsize=7.4, color=MUTED)
    ax[0].text(0.5, 0.505, r"comonotone $\max(w, 1 - w)$", ha="center", va="bottom", fontsize=7.4, color=MUTED,
               bbox=dict(fc="white", ec="none", pad=0.5))
    tag(ax[0], "a"); tag(ax[1], "b")
    f.tight_layout(w_pad=1.2)
    finish(f, "generic_pickands.png")


# ------------------------------------------------------------------------------------------ tails
def tails():
    out = run("example_tails")
    c, e = load("tails_curve.txt"), load("tails_emp.txt")
    q = c[:, 0]
    names = ["Gumbel, $\\theta = 2$", "Clayton, $\\theta = 2$", "survival Clayton, $\\theta = 2$",
             "Tawn, $(0.4, 0.9, 3)$", "Khoudraji-Joe, $(3, 0.5, 0.9)$"]
    cols = [BLUE, GREEN, GOLD, INK, RED]
    lim = [(0, 2 - 2 ** 0.5), (2 ** -0.5, 0), (0, 2 ** -0.5),
           (0, 0.4 + 0.9 - (0.4 ** 3 + 0.9 ** 3) ** (1 / 3)), (0, 0.5 + 0.9 - (0.5 ** 3 + 0.9 ** 3) ** (1 / 3))]
    f, ax = plt.subplots(1, 2, figsize=(7.2, 3.3), sharey=True)
    for s, a in enumerate(ax):
        for k in range(5):
            a.plot(q, c[:, 1 + 2 * k + s], color=cols[k], lw=1.4, label=names[k], zorder=3)
            a.plot(e[:, 0], e[:, 1 + 2 * k + s], "o", ms=3.4, mfc="white", mec=cols[k], mew=0.9, zorder=4)
            if lim[k][s] > 0:
                a.axhline(lim[k][s], color=cols[k], lw=0.9, ls=(0, (4, 3)), zorder=2)
        a.plot(q, c[:, 11 + s], color=GREY, lw=1.1, ls=":", label=r"Gaussian, $\rho = 0.707$ ($\tau = 0.5$)", zorder=2)
        a.axhline(0, color=GREY, lw=0.9, ls=(0, (4, 3)), zorder=2)
        a.set_xscale("log"); a.set_xlim(1e-6, 0.5); a.set_ylim(-0.02, 0.8)
        a.set_xlabel("$q$ (distance from the corner)")
        style(a)
    ax[0].set_ylabel(r"$\lambda(q)$")
    ax[0].set_title(r"lower tail $\lambda_L(q) = C(q, q)/q$", color=INK)
    ax[1].set_title(r"upper tail $\lambda_U(q) = \{2q - 1 + C(1 - q, 1 - q)\}/q$", color=INK)
    h, l = ax[0].get_legend_handles_labels()
    h += [plt.Line2D([], [], color=INK, lw=0.9, ls=(0, (4, 3))),
          plt.Line2D([], [], marker="o", color="none", mfc="white", mec=INK, ms=4)]
    l += [r"limit $\lambda_L$ / $\lambda_U$", "empirical, 200 000 draws"]
    f.legend(h, l, loc="lower center", ncol=4, fontsize=8, bbox_to_anchor=(0.5, -0.12), columnspacing=1.2,
             handlelength=1.8)
    tag(ax[0], "a"); tag(ax[1], "b")
    f.tight_layout(w_pad=1.0)
    finish(f, "tails_lambda.png")


# ------------------------------------------------------------------------------------------ weights / discrete
def weights():
    out = run("example_weights")
    th = {ln.split()[0]: float(ln.split()[1]) for ln in out.strip().split("\n") if ln.startswith("theta")}
    s, pr = load("w_sample.txt"), load("w_profile.txt")
    zz = norm.ppf(s[:, :2])
    # 1. weighted vs unweighted fit
    f, ax = plt.subplots(1, 2, figsize=(7.0, 3.1), gridspec_kw=dict(width_ratios=[1, 1.25]))
    a = ax[0]
    sc = a.scatter(zz[:, 0], zz[:, 1], c=s[:, 2], cmap=WCMAP, vmin=0, vmax=1, s=2.6, linewidths=0, alpha=0.8)
    cb = f.colorbar(sc, ax=a, fraction=0.046, pad=0.03)
    cb.set_label(r"weight $w_i = 0.998^{\,n - i}$", fontsize=8.5); cb.ax.tick_params(labelsize=8, length=2)
    cb.outline.set_linewidth(0.4)
    a.set_xlim(-3.5, 3.5); a.set_ylim(-3.5, 3.5); a.set_aspect("equal")
    a.set_xlabel("$z_1$"); a.set_ylabel("$z_2$")
    a.set_title("two-regime sample, coloured by weight", color=INK)
    style(a)
    a = ax[1]
    for col, k, lab, key in ((GREY, 1, "equal weights", "theta_equal"), (BLUE, 2, "time-decay weights", "theta_decay")):
        y = pr[:, k] - pr[:, k].max()
        a.plot(pr[:, 0], y, color=col, lw=1.5, label=lab + rf", $\hat\theta = {th[key]:.2f}$")
        a.plot(th[key], 0, "v", color=col, ms=6, mec="white", mew=0.6, zorder=4, clip_on=False)
    for tv in (1, 4):
        a.axvline(tv, color=GREEN, lw=0.8, ls=(0, (4, 3)))
    a.text(1.06, 0.075, r"old regime, $\theta = 1$", color=GREEN, fontsize=8, ha="left", va="center")
    a.text(4.06, 0.075, r"recent regime, $\theta = 4$", color=GREEN, fontsize=8, ha="left", va="center")
    a.set_xlim(0.5, 6.5); a.set_ylim(-0.6, 0.11)
    a.set_xlabel(r"Clayton $\theta$")
    a.set_ylabel(r"$\ell_w(\theta)/\sum w_i$, relative to its maximum")
    a.set_title("weighted log-likelihood per unit weight", color=INK)
    a.legend(loc="lower center", fontsize=8, bbox_to_anchor=(0.47, 0.0))
    style(a)
    tag(ax[0], "a"); tag(ax[1], "b")
    f.tight_layout(w_pad=1.5)
    finish(f, "weights_fit.png")

    # 2. the mixed data: a Poisson count against two continuous variables
    dd = load("d_data.txt")
    z1, x2, z3, up, um = dd.T
    rng = np.random.default_rng(3)                      # display jitter only
    xj = x2 + rng.uniform(-0.3, 0.3, len(x2))
    f, ax = plt.subplots(1, 3, figsize=(7.2, 2.6))
    for a, zc, lab, ti in ((ax[0], z1, "$z_1$", r"$(x_2, z_1)$: Gaussian, $\rho = 0.7$"),
                           (ax[1], z3, "$z_3$", r"$(x_2, z_3)$: Clayton, $\theta = 2$")):
        a.scatter(xj, zc, s=2.0, color=BLUE, alpha=0.4, linewidths=0)
        a.set_xlabel("$x_2$ (Poisson count, jittered)"); a.set_ylabel(lab)
        a.set_xticks(range(0, 10, 2)); a.set_ylim(-3.6, 3.6)
        a.set_title(ti, color=INK)
        style(a)
    a = ax[2]
    m = 150
    u3 = norm.cdf(z3[:m])
    a.hlines(u3, um[:m], up[:m], color=BLUE, lw=0.6, alpha=0.7)
    a.plot(0.5 * (um[:m] + up[:m]), u3, "o", ms=1.6, color=GOLD, alpha=0.9, mew=0)
    a.plot([], [], color=BLUE, lw=1.0, label=r"$[F(x_2^-), F(x_2)]$")
    a.plot([], [], "o", ms=3, color=GOLD, label="naive mid-point")
    a.set_xlim(0, 1); a.set_ylim(0, 1)
    a.set_xlabel("$u_2$"); a.set_ylabel("$u_3$")
    a.set_title(f"copula scale, first {m} rows", color=INK)
    a.legend(loc="upper left", fontsize=7.4, handlelength=1.2, borderaxespad=0.1, frameon=True, facecolor="white",
             edgecolor="none", framealpha=0.95)
    style(a)
    for a, t in zip(ax, "abc"):
        tag(a, t)
    f.tight_layout(w_pad=1.0)
    finish(f, "discrete_data.png")

    # 3. randomised Rosenblatt residuals
    r = load("d_ros.txt")
    n = len(r)
    nb = 20
    edges = np.linspace(0, 1, nb + 1)
    se = np.sqrt((1 / nb) * (1 - 1 / nb) / n) * nb        # sd of a bin's density estimate under uniformity
    ytop = max(2.0, np.histogram(r[:, 4], bins=edges, density=True)[0].max() + 0.35)
    f, ax = plt.subplots(1, 3, figsize=(7.2, 2.45), sharey=True)
    for k, a in enumerate(ax):
        a.fill_between([0, 1], 1 - 2 * se, 1 + 2 * se, color="#E3EDF8", lw=0, zorder=0)
        a.hist(r[:, k], bins=edges, density=True, color=BLUE, alpha=0.85, edgecolor="white", lw=0.4,
               label="randomised" if k == 1 else None)
        if k == 1:
            a.hist(r[:, 3 + k], bins=edges, density=True, histtype="step", color=GOLD, lw=1.3,
                   label=r"not randomised, $F(x_2 \mid x_1)$")
        a.axhline(1, color=INK, lw=0.8)
        a.set_xlim(0, 1); a.set_ylim(0, ytop)
        a.set_xlabel(f"$r_{k + 1}$" + (" (discrete variable)" if k == 1 else ""))
        style(a, grid=False)
        a.text(0.5, ytop - 0.08, rf"mean {r[:, k].mean():.3f}, var {r[:, k].var():.4f}", ha="center", va="top", fontsize=7.8,
               color=MUTED)
    ax[0].set_ylabel("density")
    ax[1].legend(loc="upper left", fontsize=7.6, bbox_to_anchor=(0.0, 0.88), handlelength=1.4)
    f.tight_layout(w_pad=0.8)
    finish(f, "rosenblatt_residuals.png")


# ------------------------------------------------------------------------------------------ interop
def load_dll():
    name = {"win32": "libfortvine.dll", "darwin": "libfortvine.dylib"}.get(sys.platform, "libfortvine.so")
    if sys.platform == "win32" and shutil.which("gfortran"):
        os.add_dll_directory(os.path.dirname(shutil.which("gfortran")))
    return ct.CDLL(os.path.join(BUILD, name))


def interop():
    run("example_interop")
    src, cexe = os.path.join(ROOT, "tools", "example_interop.c"), os.path.join(BUILD, "example_interop_c" + EXE)
    if not os.path.exists(cexe) or os.path.getmtime(cexe) < os.path.getmtime(src):
        lnk = [os.path.join(BUILD, "libfortvine.dll")] if sys.platform == "win32" else \
            ["-L" + BUILD, "-lfortvine", "-Wl,-rpath," + BUILD]
        subprocess.run(["gcc", "-O2", src] + lnk + ["-o", cexe], check=True)
    env = dict(ENV, PATH=BUILD + os.pathsep + ENV.get("PATH", ""))
    subprocess.run([cexe], cwd=TMP, env=env, check=True, capture_output=True)
    specs, pts = load("interop_specs.txt", ndmin=2), load("interop_points.txt")
    vf, vc = load("interop_fortran.txt"), load("interop_c.txt")
    lib = load_dll()
    dptr = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
    lib.fv_bicop_eval.argtypes = [ct.c_int, ct.c_int, dptr, ct.c_int, dptr, ct.c_int, dptr, ct.c_int]
    lib.fv_bicop_eval.restype = None
    vp = np.empty_like(vf)
    u = np.ascontiguousarray(pts)
    for k, sp in enumerate(specs):
        for j, what in enumerate((1, 2, 4)):
            o = np.empty(len(u))
            lib.fv_bicop_eval(int(sp[0]), int(sp[1]), np.ascontiguousarray(sp[2:5]), len(u), u, what, o, 1)
            vp[:, 3 * k + j] = o
    dc, dpy = np.abs(vc - vf).max(), np.abs(vp - vf).max()
    print(f"max |C - Fortran| = {dc:.3g}, max |Python - Fortran| = {dpy:.3g} over {vf.size} values")
    # the Python demo itself (its numbers are those printed on the page)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        g = runpy.run_path(os.path.join(ROOT, "example", "python", "fortvine_ctypes.py"), run_name="__main__")
    td, ts, d = g["t_data"], g["t_sim"], g["d"]

    f, ax = plt.subplots(1, 2, figsize=(7.0, 3.1))
    a = ax[0]
    pf = vf[:, 0::3].ravel()
    a.plot([1e-3, 30], [1e-3, 30], color=GREY, lw=0.8, zorder=1)
    a.scatter(pf, vp[:, 0::3].ravel(), s=12, facecolor="none", edgecolor=BLUE, lw=0.6, label="Python (ctypes)", zorder=2)
    a.scatter(pf, vc[:, 0::3].ravel(), s=7, marker="x", color=GREEN, lw=0.6, label="C", zorder=3)
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(1e-2, 30); a.set_ylim(1e-2, 30); a.set_aspect("equal")
    a.set_xlabel(r"$c(u, v)$ from Fortran (fv_bicop%pdf)"); a.set_ylabel(r"$c(u, v)$ from Python / C")
    a.set_title(f"{len(specs)} copulas $\\times$ {len(pts)} points", color=INK)
    a.legend(loc="upper left", fontsize=8)
    a.text(0.97, 0.04, f"largest difference over all\n{vf.size} values ($c$, $h_1$, $h_1^{{-1}}$): {max(dc, dpy):.0f}",
           transform=a.transAxes, ha="right", va="bottom", fontsize=7.8, color=MUTED)
    style(a)
    a = ax[1]
    iu = np.triu_indices(d, 1)
    a.plot([0, 0.5], [0, 0.5], color=GREY, lw=0.8, zorder=1)
    a.scatter(td[iu], ts[iu], s=30, color=BLUE, edgecolor="white", lw=0.8, zorder=3)
    # label offsets (points) for the seeded demo: three pairs cluster near 0.4, two near 0.23
    off = {(0, 1): (-6, 3), (2, 3): (-6, -10), (1, 2): (6, -3), (0, 2): (-6, -6)}
    for i, j in zip(*iu):
        dx, dy = off.get((i, j), (6, -3))
        a.annotate(f"({i + 1},{j + 1})", (td[i, j], ts[i, j]), xytext=(dx, dy), textcoords="offset points",
                   ha="right" if dx < 0 else "left", fontsize=7.8, color=INK)
    a.set_xlim(0, 0.5); a.set_ylim(0, 0.5); a.set_aspect("equal")
    a.set_xlabel(r"$\hat\tau_{ij}$ of the data"); a.set_ylabel(r"$\hat\tau_{ij}$ of the simulated sample")
    a.set_title(r"Python demo: $\hat\tau$ of data and simulation", color=INK)
    style(a)
    tag(ax[0], "a"); tag(ax[1], "b")
    f.tight_layout(w_pad=2.0)
    finish(f, "interop_agree.png")


JOBS = {"pair": pair, "vine": vine, "generic": generic, "tails": tails, "weights": weights, "interop": interop}

if __name__ == "__main__":
    for job in (sys.argv[1:] or list(JOBS)):
        JOBS[job]()
