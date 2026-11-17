"""Figures of the documentation pages outside the manual and the examples, from the library itself.

    pip install ./python scipy matplotlib pyvinecopulib     # fortvine binding (+ pyvinecopulib for `validation`)
    bash build.sh                                            # build/density_grid (parameters at a given tau)
    python tools/page_figures.py [concepts | quickstart | numerics | scores | validation |
                                  families | nonparametric | python]                 (default: all)

Every copula number comes from fortvine through its Python binding (densities, h-functions and
their inverses, distribution functions, fits, selections, simulation, vines); this script only
transforms margins with scipy, counts, and plots.  `validation` also runs pyvinecopulib for the
cross-library comparison; `scores` differentiates fortvine's own log-likelihood numerically.
Seeds are fixed.  PNGs go to the img/ folder next to each page (docs/assets/ for python.md).
Style and helpers (fonts, colours, density panels in normal scores, tree drawings) are those of
tools/figures.py and tools/example_figures.py.
"""
import os, sys, time
os.environ.setdefault("FORTVINE_THREADS", "8")             # at most 8 threads, as the other figure scripts
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figures as F                                       # noqa: E402  (rcParams, panel, save, grids)
import example_figures as E                               # noqa: E402  (style, tag, lines, tree drawing)
import matplotlib.pyplot as plt                           # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, LogNorm   # noqa: E402
from scipy.stats import norm                              # noqa: E402
import fortvine as fv                                     # noqa: E402

ROOT = F.ROOT
DOCS = os.path.join(ROOT, "docs")
INK, MUTED, GRIDC = E.INK, E.MUTED, E.GRIDC
BLUE, GREEN, GOLD, GREY, RED = E.BLUE, E.GREEN, E.GOLD, E.GREY, E.RED
CMAP = E.CMAP
Z1, Z2, PHI = F.Z1, F.Z2, F.PHI
UG = norm.cdf(np.column_stack([Z1.ravel(), Z2.ravel()]))  # the normal-score grid on the copula scale
RAMP = ["#AFCBEA", "#7FAEDD", "#4F86C6", "#1F5FA6", "#0D2C52"]  # five shades of the blue ramp
EPS = np.finfo(float).eps
style, tag, lines = E.style, E.tag, E.lines


def finish(f, folder, name, dpi=150):
    out = os.path.join(DOCS, folder)
    os.makedirs(out, exist_ok=True)
    p = os.path.join(out, name)
    F.save(f, p, dpi=dpi, bbox_inches="tight", pad_inches=0.05, facecolor="white")
    plt.close(f)
    print(f"{folder}/{name:26s} {os.path.getsize(p) // 1024:4d} KB")


def zgrid(cop):
    """Copula density of `cop` on the 121 x 121 normal-score grid of tools/figures.py."""
    return cop.pdf(UG).reshape(Z1.shape)


def at_tau(specs):
    """Parameters that give the requested Kendall tau (fortvine's density_grid helper inverts tau)."""
    return [par[:len(s["p"])] for s, (par, _, _) in zip(specs, F.grids(specs))]


def tau_of(cop, n=40000, seed=11):
    """Kendall's tau of a pair copula, from a large fortvine sample."""
    s = cop.simulate(n, seed=seed)
    return fv.kendall_tau(s[:, 0], s[:, 1])


NAMES = {**E.FAMN, **{13: "FGM", 14: "AMH", 17: "BB2", 18: "BB3", 19: "BB5", 20: "BB9", 21: "BB10", 22: "BB4",
                       23: "t-EV", 24: "TLL", 25: "Bernstein"}}


def fam_label(b, rot=True):
    """Display name of a fitted Bicop (with its rotation)."""
    s = NAMES.get(b.family_code, b.family)
    return s + (f" {b.rotation}°" if rot and b.rotation else "")


def zlim(a, lim=3.2):
    a.set_xlim(-lim, lim); a.set_ylim(-lim, lim); a.set_aspect("equal")


def vine_edges(vc, taus=True):
    """Edges of a fortvine vine as dicts for E.draw_tree (conditioned pair, conditioning set, family, tau)."""
    out = []
    for t in range(vc.trunc):
        for e in range(vc.d - 1 - t):
            cop = vc.get_pair_copula(t, e)
            tau = tau_of(cop, 20000) if (taus and cop.family_code != 0) else 0.0
            out.append(dict(t=t + 1, a=int(vc.order[e]), b=int(vc.order[vc.st[t, e] - 1]), fam=int(vc.fam[t, e]),
                            rot=int(vc.rot[t, e]), tau=tau, cond=tuple(sorted(int(vc.order[x - 1]) for x in vc.st[:t, e]))))
    return out


# ================================================================================ getting-started/concepts
def concepts():
    from scipy.stats import lognorm, gamma
    folder = os.path.join("getting-started", "img")
    # 1. raw data -> pseudo-observations -> normal scores
    n = 800
    v = fv.Bicop("gumbel", parameters=[2.0]).simulate(n, seed=3)
    x = np.column_stack([lognorm.ppf(v[:, 0], 0.8, scale=20.0), gamma.ppf(v[:, 1], 2.0, scale=3.0)])
    u = fv.to_pseudo_obs(x)
    z = norm.ppf(u)
    fit = fv.Bicop("gumbel").fit(u)
    f, ax = plt.subplots(1, 3, figsize=(7.2, 2.65))
    a = ax[0]
    a.scatter(x[:, 0], x[:, 1], s=3, color=BLUE, alpha=0.5, linewidths=0)
    a.set_xlabel("$x_1$ (log-normal margin)"); a.set_ylabel("$x_2$ (gamma margin)")
    a.set_title("raw data", color=INK)
    a = ax[1]
    a.scatter(u[:, 0], u[:, 1], s=3, color=BLUE, alpha=0.5, linewidths=0)
    a.set_xlim(0, 1); a.set_ylim(0, 1); a.set_aspect("equal")
    a.set_xlabel(r"$\hat u_1 = R_1 / (n + 1)$"); a.set_ylabel(r"$\hat u_2 = R_2 / (n + 1)$")
    a.set_title("pseudo-observations", color=INK)
    a = ax[2]
    a.scatter(z[:, 0], z[:, 1], s=3, color=BLUE, alpha=0.4, linewidths=0, label="data")
    lines(a, zgrid(fit), INK, lw=0.75)
    a.plot([], [], color=INK, lw=0.75, label=rf"fitted Gumbel, $\hat\theta = {fit.parameters[0]:.2f}$")
    zlim(a)
    a.set_xlabel(r"$z_1 = \Phi^{-1}(\hat u_1)$"); a.set_ylabel(r"$z_2 = \Phi^{-1}(\hat u_2)$")
    a.set_title("normal scores", color=INK)
    a.legend(loc="upper left", fontsize=7.2, handlelength=1.4, borderaxespad=0.1, frameon=True, framealpha=0.9,
             facecolor="white", edgecolor="none", markerscale=3)
    for k, a in enumerate(ax):
        style(a); tag(a, "abc"[k])
    f.tight_layout(w_pad=1.2)
    finish(f, folder, "pseudo_obs.png")

    # 2. h-functions, the density as their derivative, and the inverse
    cop = fv.Bicop("clayton", parameters=[2.0])
    vv = np.linspace(0.0005, 0.9995, 400)
    us = [0.05, 0.25, 0.5, 0.75, 0.95]
    f, ax = plt.subplots(1, 3, figsize=(7.2, 2.6))
    for col, u1 in zip(RAMP, us):
        uu = np.column_stack([np.full_like(vv, u1), vv])
        ax[0].plot(vv, cop.hfunc1(uu), color=col, lw=1.3, label=f"$u = {u1}$")
        ax[1].plot(vv, cop.pdf(uu), color=col, lw=1.3)
        ax[2].plot(vv, cop.hinv1(uu), color=col, lw=1.3)
    ax[0].set_xlabel("$v$"); ax[0].set_ylabel(r"$h_1(v \mid u) = P(V \leq v \mid U = u)$")
    ax[0].set_title("h-function", color=INK)
    ax[1].set_xlabel("$v$"); ax[1].set_ylabel(r"$c(u, v)$")
    ax[1].set_yscale("log"); ax[1].set_ylim(1e-2, 50)
    ax[1].set_title(r"density $= \partial h_1 / \partial v$", color=INK)
    ax[2].set_xlabel("$w$"); ax[2].set_ylabel(r"$h_1^{-1}(w \mid u)$")
    ax[2].set_title("inverse h-function", color=INK)
    for k, a in enumerate(ax):
        style(a); tag(a, "abc"[k])
        if k != 1:
            a.set_xlim(0, 1); a.set_ylim(0, 1)
    hh, ll = ax[0].get_legend_handles_labels()
    f.legend(hh, ll, loc="lower center", ncol=5, fontsize=8.5, bbox_to_anchor=(0.5, -0.02), handlelength=1.8,
             title=r"Clayton, $\theta = 2$ ($\tau = 0.5$), conditioning value", title_fontsize=8.5)
    f.tight_layout(w_pad=1.0, rect=(0, 0.13, 1, 1))
    finish(f, folder, "hfunctions.png")

    # 3. the four rotations of a Gumbel copula
    f, ax = plt.subplots(1, 4, figsize=(7.2, 2.1))
    for a, r in zip(ax, (0, 90, 180, 270)):
        F.panel(a, zgrid(fv.Bicop("gumbel", r, [2.0])))
        a.set_title(f"rotation = {r}" + {0: "", 90: "", 180: " (survival)", 270: ""}[r], fontsize=9.5, color=INK)
        a.set_xlabel("$z_1$", fontsize=8.5, labelpad=1)
    ax[0].set_ylabel("$z_2$", fontsize=8.5, labelpad=1)
    f.tight_layout(pad=0.3, w_pad=0.8)
    finish(f, folder, "rotations.png")

    # 4. C-vine and D-vine structures on five variables (from set_cvine / set_dvine)
    d = 5
    cv = fv.Vinecop().set_cvine([5, 4, 3, 2, 1])
    dv = fv.Vinecop().set_dvine([1, 2, 3, 4, 5])
    f, ax = plt.subplots(d - 1, 2, figsize=(7.2, 4.1), gridspec_kw=dict(hspace=0.15, wspace=0.08))
    for c, (vc, title) in enumerate(((cv, "C-vine: every tree is a star"), (dv, "D-vine: every tree is a path"))):
        edges = vine_edges(vc, taus=False)
        for t in range(1, d):
            structure_tree(ax[t - 1][c], edges, t)
            ax[t - 1][c].text(-0.75, 0.55, f"$T_{t}$", fontsize=9, color=MUTED, ha="left", va="center")
        ax[0][c].set_title(title, fontsize=10, color=INK)
    f.tight_layout(pad=0.3)
    finish(f, folder, "vine_shapes.png")


def draw_tree(a, edges, t, ref=None):
    """Tree t of a fitted vine: family and Kendall tau per edge (as E.draw_tree, but stars drawn as stars)."""
    nodes, labels, links, pos, star = tree_geometry(edges, t)
    for p, q, e in links:
        (x1, y1), (x2, y2) = pos[p], pos[q]
        indep = e["fam"] == 0
        new = ref is not None and t == 1 and frozenset([e["a"], e["b"]]) not in ref
        col = GREY if indep else (GOLD if new else INK)
        a.plot([x1, x2], [y1, y2], color=col, lw=0.9 if indep else 1.4, ls=":" if indep else "-", zorder=1)
        name = E.FAMN.get(e["fam"], str(e["fam"])) + (f" {e['rot']}°" if e["rot"] else "")
        if star:                                    # label under the leaf
            txt = name if indep else name + rf", $\tau = {e['tau']:.2f}$"
            a.text(x2, y2 - 0.2, txt, ha="center", va="top", fontsize=7.2, color=INK)
        else:
            xm, ym = 0.5 * (x1 + x2), 0.5 * (y1 + y2)
            a.text(xm, ym + 0.17, name, ha="center", va="bottom", fontsize=7.6, color=INK if not indep else MUTED)
            if not indep:
                a.text(xm, ym - 0.17, rf"$\tau = {e['tau']:.2f}$", ha="center", va="top", fontsize=7.4, color=MUTED)
    for n, (x, y) in pos.items():
        a.text(x, y, labels[n], ha="center", va="center", fontsize=8 if t == 1 else 7.2, color=INK, zorder=3,
               bbox=dict(boxstyle="round,pad=0.25" if t > 1 else "circle,pad=0.25", fc="white", ec=BLUE, lw=0.9))
    a.set_xlim(-0.75, 4.75); a.set_ylim(-0.75, 0.75)
    a.axis("off")


def tree_geometry(edges, t):
    """Nodes, node labels, links and positions of tree t; stars (one hub, >= 3 links) get a hub-over-leaves layout."""
    if t == 1:
        nodes = sorted({e["a"] for e in edges if e["t"] == 1} | {e["b"] for e in edges if e["t"] == 1})
        sets = {n: frozenset([n]) for n in nodes}
        labels = {n: str(n) for n in nodes}
    else:
        prev = [e for e in edges if e["t"] == t - 1]
        nodes = list(range(len(prev)))
        sets = {i: frozenset([e["a"], e["b"], *e["cond"]]) for i, e in enumerate(prev)}
        labels = {i: E.node_label(e) for i, e in enumerate(prev)}
    links = []
    for e in (e for e in edges if e["t"] == t):
        sa, sb = frozenset([e["a"], *e["cond"]]), frozenset([e["b"], *e["cond"]])
        p = next(n for n in nodes if sets[n] == sa); q = next(n for n in nodes if sets[n] == sb)
        links.append((p, q, e))
    deg = {n: sum(n in (p, q) for p, q, _ in links) for n in nodes}
    star = len(links) > 2 and max(deg.values()) == len(links)
    if star:
        hub = max(deg, key=deg.get)
        leaves = [n for n in nodes if n != hub]
        pos = {hub: (2.0, 0.5)}
        for i, n in enumerate(leaves):
            pos[n] = (4.0 * i / (len(leaves) - 1), -0.25)
        links = [(hub, q if p == hub else p, e) for p, q, e in links]
    else:
        pos = E.tree_layout(nodes, [(p, q) for p, q, _ in links])
        xmax = max(1.0, max(x for x, _ in pos.values()))
        pos = {n: (4.0 * x / xmax, y) for n, (x, y) in pos.items()}
    return nodes, labels, links, pos, star


def structure_tree(a, edges, t):
    """Tree t of a vine as a plain structure diagram: nodes and the conditioned pair | conditioning set."""
    nodes, labels, links, pos, star = tree_geometry(edges, t)
    for p, q, e in links:
        (x1, y1), (x2, y2) = pos[p], pos[q]
        a.plot([x1, x2], [y1, y2], color=INK, lw=1.1, zorder=1)
        if not star:
            a.text(0.5 * (x1 + x2), 0.16, E.node_label(e), ha="center", va="bottom", fontsize=7.0, color=MUTED)
    for n, (x, y) in pos.items():
        a.text(x, y, labels[n], ha="center", va="center", fontsize=8 if t == 1 else 7.2, color=INK, zorder=3,
               bbox=dict(boxstyle="round,pad=0.25" if t > 1 else "circle,pad=0.25", fc="white", ec=BLUE, lw=0.9))
    a.set_xlim(-0.75, 4.75); a.set_ylim(-0.75, 0.75)
    a.axis("off")


# ================================================================================ getting-started/quickstart
def quickstart():
    """The Python run of the quick-start program (same seeds as docs/getting-started/quickstart.md)."""
    folder = os.path.join("getting-started", "img")
    n = 1000
    cop = fv.Bicop("gumbel", parameters=[1.5])
    u12 = cop.simulate(n, seed=42)
    est = fv.Bicop("gumbel").fit(u12)
    rows = []                                       # every family of the default set, best rotation
    for code in fv.family_group("default"):
        b = fv.Bicop().select(u12, family_set=[code])
        rows.append((fam_label(b), b.bic()))
    rows.sort(key=lambda r: r[1])
    rows = rows[:12]
    f, ax = plt.subplots(1, 2, figsize=(7.2, 3.3), gridspec_kw=dict(width_ratios=[1, 1.05]))
    a = ax[0]
    z = norm.ppf(u12)
    a.scatter(z[:, 0], z[:, 1], s=2.5, color=BLUE, alpha=0.35, linewidths=0, label=f"sample, $n = {n}$")
    lines(a, zgrid(est), INK, lw=0.8)
    lines(a, zgrid(cop), GOLD, lw=0.9, ls="--")
    a.plot([], [], color=INK, lw=0.8, label=rf"fitted Gumbel, $\hat\theta = {est.parameters[0]:.3f}$")
    a.plot([], [], color=GOLD, lw=0.9, ls="--", label=r"true Gumbel, $\theta = 1.5$")
    zlim(a)
    a.set_xlabel(r"$z_1 = \Phi^{-1}(u_1)$"); a.set_ylabel(r"$z_2 = \Phi^{-1}(u_2)$")
    a.set_title("the simulated pair, normal scores", color=INK)
    a.legend(loc="upper left", fontsize=7.3, handlelength=1.6, borderaxespad=0.15, frameon=True, framealpha=0.9,
             facecolor="white", edgecolor="none", markerscale=3)
    a = ax[1]
    yy = np.arange(len(rows))[::-1]
    d0 = rows[0][1]
    for y, (name, b) in zip(yy, rows):
        a.plot([0, b - d0], [y, y], color=GRIDC, lw=1.0, zorder=1)
        a.scatter([b - d0], [y], s=30, color=BLUE if y != yy[0] else GREEN, edgecolor="white", lw=0.8, zorder=3)
    a.set_yticks(yy); a.set_yticklabels([r[0] for r in rows], fontsize=8.5)
    a.set_xlabel(r"$\mathrm{BIC} - \mathrm{BIC}_{\min}$")
    a.set_title("BIC of each family (best rotation)", color=INK)
    a.text(1.0, yy[0], "  selected", va="center", ha="left", fontsize=8, color=INK)
    style(a, grid=False); a.grid(True, axis="x", color=GRIDC, lw=0.6); a.set_axisbelow(True)
    a.spines["left"].set_visible(False); a.tick_params(axis="y", length=0)
    style(ax[0]); tag(ax[0], "a"); tag(ax[1], "b")
    f.tight_layout(w_pad=1.5)
    finish(f, folder, "quickstart_pair.png")

    # the four-dimensional vine: truth and selection, tree by tree
    truth = fv.Vinecop().set_cvine([1, 2, 3, 4])
    truth.fam[0, :3] = 3; truth.par[0, :3, 0] = 3.0
    truth.fam[1, :2] = 3; truth.par[1, :2, 0] = 1.0
    u = truth.simulate(n, seed=1)
    vc = fv.Vinecop.from_data(u)
    te, se = vine_edges(truth), vine_edges(vc)
    ref = {frozenset([e["a"], e["b"]]) for e in te if e["t"] == 1}
    f, ax = plt.subplots(3, 2, figsize=(7.2, 3.6), gridspec_kw=dict(hspace=0.25, wspace=0.08))
    for c, (edges, title) in enumerate(((te, f"true C-vine: {truth.npars} parameters, "
                                             rf"$\ell = {truth.loglik(u):.1f}$"),
                                        (se, f"selected from the sample: {vc.npars} parameters, "
                                             rf"$\ell = {vc.loglik(u):.1f}$"))):
        for t in range(1, 4):
            draw_tree(ax[t - 1][c], edges, t, ref=ref if c == 1 else None)
            ax[t - 1][c].text(-0.55, 0.62, f"$T_{t}$", fontsize=9, color=MUTED, ha="left", va="center")
        ax[0][c].set_title(title, fontsize=10, color=INK)
    h = [plt.Line2D([], [], color=INK, lw=1.4), plt.Line2D([], [], color=GOLD, lw=1.4),
         plt.Line2D([], [], color=GREY, lw=0.9, ls=":")]
    f.legend(h, ["pair copula", "tree-1 edge not in the true vine", "independence"], loc="lower center", ncol=3,
             fontsize=8.5, bbox_to_anchor=(0.5, -0.03))
    f.subplots_adjust(left=0.01, right=0.99, top=0.9, bottom=0.1)
    finish(f, folder, "quickstart_vine.png")


# ================================================================================ theory/numerics
def numerics():
    folder = os.path.join("theory", "img")
    # 1. dynamic range of a strongly dependent density
    zz = np.linspace(-5, 5, 241)
    A, B = np.meshgrid(zz, zz)
    U = norm.cdf(np.column_stack([A.ravel(), B.ravel()]))
    cop = fv.Bicop("gumbel", parameters=[5.0])
    lc = cop.logpdf(U).reshape(A.shape) / np.log(10)
    f, ax = plt.subplots(1, 2, figsize=(7.2, 3.1), gridspec_kw=dict(width_ratios=[1, 1.25]))
    a = ax[0]
    lev = np.arange(np.floor(lc.min() / 10) * 10, lc.max() + 10, 5)
    cf = a.contourf(A, B, lc, levels=lev, cmap=CMAP)
    cs = a.contour(A, B, lc, levels=[-40, -20, -10, 0], colors=INK, linewidths=0.5)
    a.clabel(cs, fmt=lambda v: f"{v:.0f}", fontsize=7, inline_spacing=2)
    zlim(a, 5)
    a.set_xlabel(r"$z_1 = \Phi^{-1}(u_1)$"); a.set_ylabel(r"$z_2 = \Phi^{-1}(u_2)$")
    a.set_title(r"$\log_{10} c(u_1, u_2)$, Gumbel $\theta = 5$", color=INK)
    cb = f.colorbar(cf, ax=a, fraction=0.046, pad=0.03)
    cb.ax.tick_params(labelsize=7.5); cb.outline.set_visible(False)
    a = ax[1]
    s = np.linspace(-6.3, 6.3, 401)                 # u from 1.5e-10 to 1 - 1.5e-10 (inside the clamp)
    for (lab, col, sgn, ls) in (("diagonal $z_2 = z_1$", BLUE, 1, "-"), ("anti-diagonal $z_2 = -z_1$", GREEN, -1, "-")):
        uu = norm.cdf(np.column_stack([s, sgn * s]))
        a.plot(s, cop.logpdf(uu) / np.log(10), color=col, lw=1.4, ls=ls, label=lab)
    for th, ls in ((2.0, "--"),):
        c2 = fv.Bicop("gumbel", parameters=[th])
        uu = norm.cdf(np.column_stack([s, -s]))
        a.plot(s, c2.logpdf(uu) / np.log(10), color=GOLD, lw=1.2, ls=ls, label=rf"anti-diagonal, $\theta = {th:g}$")
    a.set_xlabel("$z_1$"); a.set_ylabel(r"$\log_{10} c$")
    a.set_title("density along the two diagonals", color=INK)
    a.legend(loc="lower center", fontsize=7.8)
    a.set_xlim(-6.3, 6.3)
    style(a); style(ax[0], grid=False)
    tag(ax[0], "a"); tag(ax[1], "b")
    f.tight_layout(w_pad=1.6)
    finish(f, folder, "density_range.png")

    # 2. residual of the inverse h-functions over the square
    fams = [("Gumbel (4), $\\theta = 5$", fv.Bicop("gumbel", parameters=[5.0]), "monotone Newton"),
            ("BB7 (9), $(2, 3)$", fv.Bicop("bb7", parameters=[2.0, 3.0]), "safeguarded Newton"),
            ("Archimax 10206, $(1, 2, 0.4)$", fv.Bicop(10206, 0, [1.0, 2.0, 0.4]), "safeguarded Newton")]
    zz = np.linspace(-4.5, 4.5, 181)
    A, B = np.meshgrid(zz, zz)
    U1, W = norm.cdf(A.ravel()), norm.cdf(B.ravel())
    f, ax = plt.subplots(1, 3, figsize=(7.2, 2.75))
    norm_ = LogNorm(1e-18, 1e-10)
    rcm = LinearSegmentedColormap.from_list("res", ["#FFFFFF", "#E3EDF8", "#AFCBEA", "#5F95D0", "#1F5FA6",
                                                    "#0D2C52"])
    for a, (name, c, how) in zip(ax, fams):
        v = c.hinv1(np.column_stack([U1, W]))
        r = np.abs(c.hfunc1(np.column_stack([U1, v])) - W).reshape(A.shape)
        im = a.pcolormesh(A, B, np.maximum(r, 1e-18), norm=norm_, cmap=rcm, shading="auto", rasterized=True)
        zlim(a, 4.5)
        a.set_title(f"{name}\nmax {r.max():.1e}, {how}", fontsize=8.8, color=INK)
        a.set_xlabel(r"$\Phi^{-1}(u)$", fontsize=8.5)
        a.set_xticks([-4, -2, 0, 2, 4]); a.set_yticks([-4, -2, 0, 2, 4]); a.tick_params(labelsize=7.5, length=2)
    ax[0].set_ylabel(r"$\Phi^{-1}(w)$", fontsize=8.5)
    f.tight_layout(w_pad=0.6, rect=(0, 0, 0.92, 1))
    cax = f.add_axes([0.93, 0.2, 0.015, 0.6])
    cb = f.colorbar(im, cax=cax, extend="both")
    cb.set_label(r"$|h_1(h_1^{-1}(w \mid u) \mid u) - w|$", fontsize=8.5)
    cb.ax.tick_params(labelsize=7.5); cb.outline.set_visible(False)
    finish(f, folder, "hinv_residual.png")

    # 3. cost per point of density, h-function and inverse h-function (one thread)
    cases = [("Gaussian", fv.Bicop("gaussian", parameters=[0.7])),
             ("Student t", fv.Bicop("student", parameters=[0.7, 4.0])),
             ("Clayton", fv.Bicop("clayton", parameters=[2.0])),
             ("Gumbel", fv.Bicop("gumbel", parameters=[2.0])),
             ("Frank", fv.Bicop("frank", parameters=[5.7])),
             ("Joe", fv.Bicop("joe", parameters=[2.9])),
             ("BB1", fv.Bicop("bb1", parameters=[0.5, 1.5])),
             ("BB7", fv.Bicop("bb7", parameters=[1.5, 1.5])),
             ("Tawn", fv.Bicop("tawn", parameters=[0.5, 0.9, 3.0])),
             ("Hüsler-Reiss", fv.Bicop("husler_reiss", parameters=[2.0])),
             ("Archimax 10206", fv.Bicop(10206, 0, [1.0, 2.0, 0.4]))]
    n = 200000
    uu = np.random.default_rng(0).random((n, 2))
    res = []
    for name, c in cases:
        row = []
        for fn in (c.pdf, c.hfunc1, c.hinv1):
            best = np.inf
            for _ in range(3):
                t0 = time.perf_counter(); fn(uu, num_threads=1); best = min(best, time.perf_counter() - t0)
            row.append(best / n * 1e9)
        res.append(row)
    res = np.array(res)
    f, a = plt.subplots(figsize=(7.2, 3.2))
    yy = np.arange(len(cases))[::-1]
    for k, (lab, col) in enumerate((("density", BLUE), ("h-function $h_1$", GREEN),
                                    ("inverse $h_1^{-1}$", GOLD))):
        a.barh(yy + (1 - k) * 0.27, res[:, k], height=0.25, color=col, label=lab)
    for y, r in zip(yy, res):
        a.text(r.max() * 1.08, y - 0.27, f"{r[2] / r[1]:.1f}×", va="center", fontsize=7.5, color=MUTED)
    a.set_xscale("log")
    a.set_yticks(yy); a.set_yticklabels([c[0] for c in cases])
    a.set_xlabel("time per point [ns], one thread")
    a.legend(loc="lower center", fontsize=8.5, ncol=3, bbox_to_anchor=(0.45, 1.0))
    a.set_xlim(right=res.max() * 1.6)
    style(a, grid=False); a.grid(True, axis="x", color=GRIDC, lw=0.6, which="both"); a.set_axisbelow(True)
    a.tick_params(axis="y", length=0); a.spines["left"].set_visible(False)
    f.tight_layout()
    finish(f, folder, "eval_cost.png")
    print("  ns/point (pdf, h1, hinv1):", {c[0]: np.round(r, 1).tolist() for c, r in zip(cases, res)})


# ================================================================================ theory/scores
def ridders(fun, x, h):
    """Ridders' extrapolated central difference (Numerical Recipes dfridr): derivative and error estimate."""
    con, con2, ntab = 1.4, 1.96, 12
    a = np.zeros((ntab, ntab))
    hh = h
    a[0, 0] = (fun(x + hh) - fun(x - hh)) / (2 * hh)
    err, ans = np.inf, a[0, 0]
    for i in range(1, ntab):
        hh /= con
        a[0, i] = (fun(x + hh) - fun(x - hh)) / (2 * hh)
        fac = con2
        for j in range(1, i + 1):
            a[j, i] = (a[j - 1, i] * fac - a[j - 1, i - 1]) / (fac - 1)
            fac *= con2
            errt = max(abs(a[j, i] - a[j - 1, i]), abs(a[j, i] - a[j - 1, i - 1]))
            if errt <= err:
                err, ans = errt, a[j, i]
        if abs(a[i, i] - a[i - 1, i - 1]) >= 2 * err:
            break
    return ans, err


def scores():
    folder = os.path.join("theory", "img")
    cases = [("Frank, $\\theta$", "frank", [5.7], 0),
             ("Hüsler-Reiss, $\\lambda$", "husler_reiss", [2.0], 0),
             ("BB1, $\\delta$ (second parameter)", "bb1", [0.6, 1.4], 1)]
    n = 2000
    hs = np.logspace(-13, -1, 61)
    f, ax = plt.subplots(1, 3, figsize=(7.2, 2.8), sharey=True)
    for a, (title, fam, p, k) in zip(ax, cases):
        u = fv.Bicop(fam, 0, p).simulate(n, seed=5)
        x0 = 0.85 * p[k]                             # away from the MLE: a sizeable score

        def ll(x):
            q = list(p); q[k] = x
            return float(np.sum(fv.Bicop(fam, 0, q).logpdf(u)))
        ref, err = ridders(ll, x0, 0.05 * x0)
        l0 = ll(x0)
        fwd = np.array([(ll(x0 + h * x0) - l0) / (h * x0) for h in hs])
        cen = np.array([(ll(x0 + h * x0) - ll(x0 - h * x0)) / (2 * h * x0) for h in hs])
        a.loglog(hs, np.maximum(np.abs(fwd - ref) / abs(ref), 1e-16), color=GOLD, lw=1.3, label="forward difference")
        a.loglog(hs, np.maximum(np.abs(cen - ref) / abs(ref), 1e-16), color=BLUE, lw=1.3, label="central difference")
        a.axhline(max(err / abs(ref), 1e-16), color=GREY, lw=0.8, ls=":")
        for hv, lab in ((np.sqrt(EPS), r"$\sqrt{\varepsilon}$"), (EPS ** (1 / 3), r"$\varepsilon^{1/3}$")):
            a.axvline(hv, color=GREY, lw=0.7, ls="--")
            a.text(hv * 1.25, 3e-1, lab, fontsize=8.5, color=MUTED, va="center")
        a.set_title(title, fontsize=9.5, color=INK)
        a.set_xlabel("relative step $h / \\theta$")
        a.set_ylim(1e-15, 1e1); a.set_xlim(1e-13, 1e-1)
        style(a)
        print(f"  {fam}: score {ref:.6g} (Ridders error {err:.1e}), best central "
              f"{np.min(np.abs(cen - ref)) / abs(ref):.1e}, best forward {np.min(np.abs(fwd - ref)) / abs(ref):.1e}")
    ax[0].set_ylabel("relative error of the score")
    for k, a in enumerate(ax):
        tag(a, "abc"[k])
    hh, lab = ax[0].get_legend_handles_labels()
    hh.append(plt.Line2D([], [], color=GREY, lw=0.8, ls=":")); lab.append("error estimate of the reference")
    f.legend(hh, lab, loc="lower center", ncol=3, fontsize=8.5, bbox_to_anchor=(0.5, -0.02))
    f.tight_layout(w_pad=0.8, rect=(0, 0.08, 1, 1))
    finish(f, folder, "fd_error.png")


# ================================================================================ theory/validation
PVC_FAMS = {  # fortvine name: (parameter sampler, rotations)
    "gaussian": (lambda r: [r.uniform(-0.9, 0.9)], (0,)),
    "student": (lambda r: [r.uniform(-0.9, 0.9), r.uniform(2.5, 20)], (0,)),
    "clayton": (lambda r: [r.uniform(0.1, 8)], (0, 90, 180, 270)),
    "gumbel": (lambda r: [r.uniform(1.05, 6)], (0, 90, 180, 270)),
    "frank": (lambda r: [r.choice([-1, 1]) * r.uniform(0.5, 20)], (0,)),
    "joe": (lambda r: [r.uniform(1.05, 6)], (0, 90, 180, 270)),
    "bb1": (lambda r: [r.uniform(0.1, 3), r.uniform(1, 3)], (0, 90, 180, 270)),
    "bb6": (lambda r: [r.uniform(1, 3), r.uniform(1, 3)], (0, 90, 180, 270)),
    "bb7": (lambda r: [r.uniform(1, 3), r.uniform(0.1, 3)], (0, 90, 180, 270)),
    "bb8": (lambda r: [r.uniform(1, 5), r.uniform(0.2, 1)], (0, 90, 180, 270)),
    "tawn": (lambda r: [r.uniform(0.1, 1), r.uniform(0.1, 1), r.uniform(1.2, 8)], (0, 90, 180, 270)),
}
PRETTY = {"gaussian": "Gaussian", "student": "Student t", "clayton": "Clayton", "gumbel": "Gumbel",
          "frank": "Frank", "joe": "Joe", "bb1": "BB1", "bb6": "BB6", "bb7": "BB7", "bb8": "BB8", "tawn": "Tawn"}


def validation():
    import pyvinecopulib as pv
    folder = os.path.join("theory", "img")
    rng = np.random.default_rng(2026)
    names = list(PVC_FAMS)
    dl, dh, dfit = {}, {}, {}
    for nm in names:
        samp, rots = PVC_FAMS[nm]
        dl[nm], dh[nm], dfit[nm] = [], [], []
        for k in range(20):
            p, r = samp(rng), int(rng.choice(rots))
            u = rng.uniform(0.001, 0.999, (500, 2))
            a = fv.Bicop(nm, r, p)
            b = pv.Bicop(family=getattr(pv.BicopFamily, nm), rotation=r, parameters=np.array(p).reshape(-1, 1))
            dl[nm].append(np.max(np.abs(a.logpdf(u) - np.log(b.pdf(u)))))
            dh[nm].append(max(np.max(np.abs(a.hfunc1(u) - b.hfunc1(u))), np.max(np.abs(a.hfunc2(u) - b.hfunc2(u)))))
            if k < 10:                                   # maximum likelihood on a fortvine sample
                x = a.simulate(1000, seed=100 + k)
                ea = fv.Bicop(nm, r).fit(x)
                eb = pv.Bicop(family=getattr(pv.BicopFamily, nm), rotation=r)
                eb.fit(x, controls=pv.FitControlsBicop(parametric_method="mle"))
                dfit[nm].append(ea.loglik() - eb.loglik(x))
    f, ax = plt.subplots(1, 3, figsize=(7.2, 3.4), sharey=True)
    yy = np.arange(len(names))[::-1]
    jit = np.random.default_rng(1)
    for a, data, lab, col in ((ax[0], dl, r"max $|\Delta \log c|$ on 500 points", BLUE),
                              (ax[1], dh, r"max $|\Delta h_1|, |\Delta h_2|$", GREEN)):
        for y, nm in zip(yy, names):
            v = np.maximum(np.array(data[nm]), 1e-17)
            a.scatter(v, y + jit.uniform(-0.18, 0.18, v.size), s=11, color=col, alpha=0.7, linewidths=0)
        a.set_xscale("log"); a.set_xlabel(lab); a.set_xlim(1e-17, 1e-3)
    a = ax[2]
    for y, nm in zip(yy, names):
        v = np.array(dfit[nm])
        a.scatter(v, y + jit.uniform(-0.18, 0.18, v.size), s=11, color=GOLD, alpha=0.8, linewidths=0)
    a.set_xscale("symlog", linthresh=1e-6); a.axvline(0, color=GREY, lw=0.8)
    a.set_xlabel(r"$\hat\ell_{\mathrm{fortvine}} - \hat\ell_{\mathrm{pyvinecopulib}}$ (MLE)")
    a.set_xlim(-1.5, 1e-3)
    a.set_xticks([-1, -1e-2, -1e-4, -1e-6, 0, 1e-6])
    notes = {"student": r"Student: $\nu \leq 30$ vs 50"}      # the parameter boxes differ
    for nm in names:
        if min(dfit[nm]) > -1e-3:
            continue
        txt, i = notes.get(nm, PRETTY[nm]), names.index(nm)
        a.annotate(txt, (min(dfit[nm]), yy[i]), xytext=(7, 0), textcoords="offset points", fontsize=7.2,
                   color=INK, va="center")
    a.tick_params(axis="x", labelsize=7, rotation=45)
    ax[0].set_yticks(yy); ax[0].set_yticklabels([PRETTY[n] for n in names])
    for k, a in enumerate(ax):
        style(a, grid=False); a.grid(True, axis="x", color=GRIDC, lw=0.6); a.set_axisbelow(True)
        a.tick_params(axis="y", length=0)
        for y in yy:
            a.axhline(y, color=GRIDC, lw=0.5, zorder=0)
        tag(a, "abc"[k])
    f.tight_layout(w_pad=0.8)
    finish(f, folder, "pyvinecopulib_agreement.png")
    allfit = np.concatenate([dfit[n] for n in names])
    print("  max |dlogc|", max(max(v) for v in dl.values()), " max |dh|", max(max(v) for v in dh.values()),
          " fit diff min/median/max", allfit.min(), np.median(allfit), allfit.max())

    # 2. generic engine against the hand-coded formula, over the square
    pairs = [("Gumbel", 4, 10100, [3.0]), ("Galambos", 15, 10002, [1.3]),
             ("BB5", 19, 10102, [1.4, 1.2]), ("BB4", 22, 10202, [1.1, 1.3])]
    zz = np.linspace(norm.ppf(0.001), norm.ppf(0.999), 161)
    A, B = np.meshgrid(zz, zz)
    U = norm.cdf(np.column_stack([A.ravel(), B.ravel()]))
    f, ax = plt.subplots(1, 4, figsize=(7.2, 2.35))
    nrm = LogNorm(1e-16, 1e-9)
    for a, (nm, hc, gc, p) in zip(ax, pairs):
        dd = np.abs(fv.Bicop(hc, 0, p).logpdf(U) - fv.Bicop(gc, 0, p).logpdf(U)).reshape(A.shape)
        im = a.pcolormesh(A, B, np.maximum(dd, 1e-17), norm=nrm, cmap=CMAP, shading="auto", rasterized=True)
        a.set_title(f"{nm}: {hc} vs {gc}\nmax {dd.max():.1e}", fontsize=8.6, color=INK)
        zlim(a, zz[-1])
        a.set_xticks([-2, 0, 2]); a.set_yticks([-2, 0, 2]); a.tick_params(labelsize=7.5, length=2)
        a.set_xlabel("$z_1$", fontsize=8.5, labelpad=1)
    ax[0].set_ylabel("$z_2$", fontsize=8.5)
    f.tight_layout(w_pad=0.5, rect=(0, 0, 0.91, 1))
    cax = f.add_axes([0.925, 0.2, 0.013, 0.6])
    cb = f.colorbar(im, cax=cax, extend="both")
    cb.set_label(r"$|\log c_{\mathrm{hand}} - \log c_{\mathrm{generic}}|$", fontsize=8)
    cb.ax.tick_params(labelsize=7); cb.outline.set_visible(False)
    finish(f, folder, "generic_vs_hand.png")


# ================================================================================ families/index
TAIL_SPECS = [
    dict(name="Gaussian", fam=1, p=[0.5], k=1, lo=0, hi=0.99, tau=0.5),
    dict(name=r"Student t, $\nu = 4$", fam=2, p=[0.5, 4], k=1, lo=0, hi=0.99, tau=0.5),
    dict(name="Clayton", fam=3, p=[2], k=1, lo=0.01, hi=20, tau=0.5),
    dict(name="Gumbel", fam=4, p=[2], k=1, lo=1, hi=20, tau=0.5),
    dict(name="Frank", fam=5, p=[5], k=1, lo=0.1, hi=30, tau=0.5),
    dict(name="Joe", fam=6, p=[2], k=1, lo=1, hi=20, tau=0.5),
    dict(name=r"BB1, $\delta = 1.3$", fam=7, p=[1, 1.3], k=1, lo=0.01, hi=6, tau=0.5),
    dict(name=r"BB7, $\theta = 1.4$", fam=9, p=[1.4, 1], k=2, lo=0.01, hi=20, tau=0.5),
    dict(name="Galambos", fam=15, p=[1], k=1, lo=0.05, hi=20, tau=0.5),
    dict(name="Hüsler-Reiss", fam=16, p=[1], k=1, lo=0.05, hi=20, tau=0.5),
    dict(name=r"BB4, $\theta = 0.6$", fam=22, p=[0.6, 1], k=2, lo=0.05, hi=10, tau=0.5),
    dict(name=r"Tawn, $(0.5, 1)$", fam=11, p=[0.5, 1.0, 2], k=3, lo=1, hi=40, tau=0.4),
]


def families():
    folder = os.path.join("families", "img")
    # 1. finite-level tail dependence (lambda_L(q), lambda_U(q)) as q -> 0, at Kendall's tau = 0.5
    pars = at_tau(TAIL_SPECS)
    qs = np.logspace(-1, -8, 15)
    f, a = plt.subplots(figsize=(7.2, 4.6))
    a.plot([0, 1], [0, 1], color=GRIDC, lw=0.8, zorder=0)
    pts = {}
    for s, p in zip(TAIL_SPECS, pars):
        c = fv.Bicop(s["fam"], 0, list(p))
        lo = c.cdf(np.column_stack([qs, qs])) / qs
        hi = (2 * qs - 1 + c.cdf(np.column_stack([1 - qs, 1 - qs]))) / qs
        a.plot(lo, hi, color=GREY, lw=0.9, zorder=1)
        a.scatter(lo[0], hi[0], s=14, facecolor="white", edgecolor=GREY, lw=0.8, zorder=2)
        a.scatter(lo[-1], hi[-1], s=30, color=BLUE, edgecolor="white", lw=0.8, zorder=3)
        pts[s["name"]] = (lo[-1], hi[-1])
    off = {"Gaussian": (6, 2), r"Student t, $\nu = 4$": (-64, 9), "Clayton": (-16, 9), "Gumbel": (8, -7),
           "Frank": (6, -8), "Joe": (8, -9), r"BB1, $\delta = 1.3$": (6, 0), r"BB7, $\theta = 1.4$": (6, 0),
           "Galambos": (8, 4), "Hüsler-Reiss": (8, -2), r"BB4, $\theta = 0.6$": (6, 0),
           r"Tawn, $(0.5, 1)$": (8, -9)}
    for nm, (x, y) in pts.items():
        if nm in ("Galambos", "Hüsler-Reiss"):
            continue                                # drawn with Gumbel: the three coincide at q = 1e-8
        lab = "Gumbel, Galambos, Hüsler-Reiss" if nm == "Gumbel" else nm
        a.annotate(lab, (x, y), xytext=off.get(nm, (6, 0)), textcoords="offset points", fontsize=8, color=INK,
                   va="center")
    a.set_xlim(-0.03, 0.85); a.set_ylim(-0.03, 0.75)
    a.set_xlabel(r"lower tail $\lambda_L(q) = C(q, q) / q$")
    a.set_ylabel(r"upper tail $\lambda_U(q) = \{2q - 1 + C(1 - q, 1 - q)\} / q$")
    h = [plt.Line2D([], [], color=GREY, marker="o", mfc="white", mec=GREY, ms=4, lw=0.9),
         plt.Line2D([], [], color=BLUE, marker="o", ls="", ms=5.5)]
    a.legend(h, [r"path from $q = 10^{-1}$", r"$q = 10^{-8}$"], loc="upper right", fontsize=8.5)
    style(a)
    f.tight_layout()
    finish(f, folder, "tail_map.png")
    for nm, (x, y) in pts.items():
        print(f"  {nm:28s} lambda_L {x:.3f} lambda_U {y:.3f}")

    # 2. look-alikes: which family BIC picks for samples of a Gumbel copula
    rep, n = 100, 1000
    truth = fv.Bicop("gumbel", parameters=[2.0])
    picks, dbic = {}, {nm: [] for nm in ("galambos", "husler_reiss", "tev", "bb1", "joe")}
    for r in range(rep):
        u = truth.simulate(n, seed=500 + r)
        b = fv.Bicop().select(u, family_set="all_parametric")
        key = E.FAMN.get(b.family_code, {23: "t-EV", 19: "BB5", 22: "BB4", 17: "BB2", 18: "BB3",
                                         20: "BB9", 21: "BB10"}.get(b.family_code, b.family))
        key += f" {b.rotation}°" if b.rotation else ""
        picks[key] = picks.get(key, 0) + 1
        g = fv.Bicop("gumbel").fit(u).bic()
        for nm in dbic:
            dbic[nm].append(fv.Bicop(nm).fit(u).bic() - g)
    f, ax = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw=dict(width_ratios=[1, 1.2]))
    a = ax[0]
    items = sorted(picks.items(), key=lambda kv: kv[1])
    a.barh(range(len(items)), [v for _, v in items], color=BLUE, height=0.6)
    for i, (_, v) in enumerate(items):
        a.text(v + 0.8, i, str(v), va="center", fontsize=8, color=INK)
    a.set_yticks(range(len(items))); a.set_yticklabels([k for k, _ in items])
    a.set_xlabel(f"times selected (of {rep} samples)")
    a.set_title(r"BIC choice, truth Gumbel $\theta = 2$", color=INK)
    style(a, grid=False); a.tick_params(axis="y", length=0); a.spines["left"].set_visible(False)
    a = ax[1]
    lab = {"galambos": "Galambos", "husler_reiss": "Hüsler-Reiss", "tev": "t-EV", "bb1": "BB1", "joe": "Joe"}
    yy = np.arange(len(dbic))[::-1]
    jit = np.random.default_rng(3)
    for y, nm in zip(yy, dbic):
        v = np.array(dbic[nm])
        a.scatter(v, y + jit.uniform(-0.18, 0.18, v.size), s=9, color=GREEN, alpha=0.6, linewidths=0)
        a.plot([np.median(v)] * 2, [y - 0.3, y + 0.3], color=INK, lw=1.2)
    a.axvline(0, color=GREY, lw=0.8)
    a.set_xscale("symlog", linthresh=10)
    a.set_yticks(yy); a.set_yticklabels([lab[k] for k in dbic])
    a.set_xlabel(r"$\mathrm{BIC} - \mathrm{BIC}_{\mathrm{Gumbel}}$ (median: black bar)")
    a.set_title("the same samples, one family at a time", color=INK)
    style(a, grid=False); a.grid(True, axis="x", color=GRIDC, lw=0.6); a.set_axisbelow(True)
    a.tick_params(axis="y", length=0); a.spines["left"].set_visible(False)
    tag(ax[0], "a"); tag(ax[1], "b")
    f.tight_layout(w_pad=1.5)
    finish(f, folder, "gumbel_lookalikes.png")
    print("  picks", picks)


# ================================================================================ families/nonparametric
def mixture(n, seed):
    """Equal mixture of Clayton(4) and Gumbel(3) rotated 90 degrees (the python.md example)."""
    rng = np.random.default_rng(seed)
    a = fv.Bicop("clayton", parameters=[4.0]).simulate(n, seed=seed + 1)
    b = fv.Bicop("gumbel", 90, [3.0]).simulate(n, seed=seed + 2)
    return np.where(rng.random((n, 1)) < 0.5, a, b)


def mixture_pdf(u):
    return 0.5 * fv.Bicop("clayton", parameters=[4.0]).pdf(u) + 0.5 * fv.Bicop("gumbel", 90, [3.0]).pdf(u)


def nonparametric():
    folder = os.path.join("families", "img")
    n = 2000
    u = mixture(n, 5)
    tll = fv.Bicop("tll").fit(u)
    bern = fv.Bicop("bernstein").fit(u)
    par = fv.Bicop().select(u, family_set="all_parametric")
    gt = mixture_pdf(UG).reshape(Z1.shape)
    f, ax = plt.subplots(1, 4, figsize=(7.2, 2.25))
    for a, (g, ttl) in zip(ax, ((gt, "truth: mixture"), (zgrid(tll), f"TLL, edf {tll.npars:.0f}"),
                                (zgrid(bern), f"Bernstein, {bern.npars:.0f} weights"),
                                (zgrid(par), f"best parametric (BIC):\n{E.FAMN.get(par.family_code, par.family)}"))):
        F.panel(a, g)
        a.set_title(ttl, fontsize=9, color=INK)
    f.tight_layout(pad=0.3, w_pad=0.8)
    finish(f, folder, "np_estimates.png")

    # 2. held-out log-likelihood per observation against n
    ns = [250, 500, 1000, 2000, 4000, 8000]
    test = mixture(20000, 99)
    true_ll = np.mean(np.log(mixture_pdf(test)))
    res = {k: [] for k in ("TLL", "Bernstein", "parametric (BIC)")}
    edf = {k: [] for k in ("TLL", "Bernstein")}
    for m in ns:
        acc = {k: [] for k in res}
        acc_e = {k: [] for k in edf}
        for r in range(4):
            x = mixture(m, 1000 + 17 * r + m)
            t, b = fv.Bicop("tll").fit(x), fv.Bicop("bernstein").fit(x)
            p = fv.Bicop().select(x, family_set="all_parametric")
            acc["TLL"].append(np.mean(t.logpdf(test))); acc["Bernstein"].append(np.mean(b.logpdf(test)))
            acc["parametric (BIC)"].append(np.mean(p.logpdf(test)))
            acc_e["TLL"].append(t.npars); acc_e["Bernstein"].append(b.npars)
        for k in res:
            res[k].append(np.mean(acc[k]))
        for k in edf:
            edf[k].append(np.mean(acc_e[k]))
    f, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
    a = ax[0]
    for (k, col, mk) in (("TLL", BLUE, "o"), ("Bernstein", GREEN, "s"), ("parametric (BIC)", GOLD, "^")):
        a.plot(ns, res[k], color=col, marker=mk, ms=4.5, mec="white", mew=0.8, lw=1.3, label=k)
    a.axhline(true_ll, color=INK, lw=0.9, ls="--", label="true density")
    a.set_xscale("log"); a.set_xlabel("training sample size $n$")
    a.set_ylabel("held-out log-likelihood per observation")
    a = ax[1]
    for (k, col, mk) in (("TLL", BLUE, "o"), ("Bernstein", GREEN, "s")):
        a.plot(ns, edf[k], color=col, marker=mk, ms=4.5, mec="white", mew=0.8, lw=1.3, label=k)
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel("training sample size $n$")
    a.set_ylabel("effective number of parameters (npars)")
    for k, a in enumerate(ax):
        style(a); tag(a, "ab"[k])
        a.set_xticks(ns); a.set_xticklabels([str(m) for m in ns]); a.minorticks_off()
    ax[1].set_yticks([20, 50, 100, 200, 400]); ax[1].set_yticklabels(["20", "50", "100", "200", "400"])
    hh, ll = ax[0].get_legend_handles_labels()
    f.legend(hh, ll, loc="lower center", ncol=4, fontsize=8.5, bbox_to_anchor=(0.5, -0.02))
    f.tight_layout(w_pad=1.5, rect=(0, 0.08, 1, 1))
    finish(f, folder, "np_convergence.png")
    print("  heldout", {k: np.round(v, 3).tolist() for k, v in res.items()}, "truth", round(true_ll, 3))


# ================================================================================ python.md
def python():
    folder = "assets"
    # 1. observation weights: a dependence that strengthens over time, time-decay weights
    n = 2000
    rng = np.random.default_rng(8)
    th_t = np.linspace(0.5, 6.0, n)                 # Clayton parameter drifting upwards
    w1, w2 = rng.random(n), rng.random(n)
    u = np.empty((n, 2)); u[:, 0] = w1
    for i in range(n):
        u[i, 1] = fv.Bicop("clayton", parameters=[th_t[i]]).hinv1(np.array([[w1[i], w2[i]]]))[0]
    halfl = [50, 100, 200, 400, 800, 1600, 3200]
    est = []
    for hl in halfl:
        w = 0.5 ** (np.arange(n)[::-1] / hl)
        est.append(fv.Bicop("clayton").fit(u, weights=w).parameters[0])
    plain = fv.Bicop("clayton").fit(u).parameters[0]
    f, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
    a = ax[0]
    z = norm.ppf(u)
    wcol = 0.5 ** (np.arange(n)[::-1] / 200)
    sc = a.scatter(z[:, 0], z[:, 1], c=wcol, cmap=E.WCMAP, s=3, linewidths=0, vmin=0, vmax=1)
    zlim(a, 3.4)
    a.set_xlabel(r"$\Phi^{-1}(u_1)$"); a.set_ylabel(r"$\Phi^{-1}(u_2)$")
    a.set_title("colour: weight, half-life 200", color=INK, fontsize=9.3)
    cb = f.colorbar(sc, ax=a, fraction=0.046, pad=0.03); cb.ax.tick_params(labelsize=7.5); cb.outline.set_visible(False)
    a = ax[1]
    a.plot(halfl, est, color=BLUE, marker="o", ms=4.5, mec="white", lw=1.3, label="weighted fit")
    a.axhline(plain, color=GREEN, lw=1.2, ls="--", label="unweighted fit")
    a.axhline(th_t[-1], color=INK, lw=0.9, ls=":", label=rf"latest true $\theta = {th_t[-1]:g}$")
    a.axhline(th_t.mean(), color=GREY, lw=0.9, ls="-.", label=rf"average true $\theta = {th_t.mean():.2f}$")
    a.set_xscale("log"); a.set_xlabel("half-life of the weights [observations]")
    a.set_ylabel(r"Clayton $\hat\theta$")
    a.set_title(r"$w_i = 2^{-(n - i)/\mathrm{half\,life}}$", color=INK, fontsize=9.3)
    a.legend(loc="upper right", fontsize=7.8)
    style(ax[0], grid=False); style(a)
    tag(ax[0], "a"); tag(a, "b")
    f.tight_layout(w_pad=1.2)
    finish(f, folder, "python_weights.png")

    # 2. discrete margins: the python.md example (Gaussian D-vine, variable 2 a Poisson count)
    from scipy.stats import poisson
    truth = fv.Vinecop().set_dvine([1, 2, 3])
    truth.fam[0, :2] = 1
    truth.par[0, :2, 0] = [0.7, 0.5]
    v = truth.simulate(2000, seed=1)
    k = poisson.ppf(v[:, 1], 2.0)
    ud = np.column_stack([v[:, 0], poisson.cdf(k, 2.0), v[:, 2], poisson.cdf(k - 1, 2.0)])
    vc = fv.Vinecop.from_data(ud, var_types=["c", "d", "c"], family_set=["indep", "gaussian", "frank"])
    r_rand = vc.rosenblatt(ud, seed=3)
    r_plain = vc.rosenblatt(ud, randomize_discrete=False)
    col = list(vc.order).index(2)                  # the column of variable 2 in the transform
    f, ax = plt.subplots(1, 3, figsize=(7.2, 2.7))
    a = ax[0]
    jx = np.random.default_rng(4).uniform(-0.25, 0.25, len(k))
    a.scatter(norm.ppf(v[:, 0]), k + jx, s=2.5, color=BLUE, alpha=0.4, linewidths=0)
    a.set_xlabel(r"variable 1, $\Phi^{-1}(u_1)$"); a.set_ylabel("variable 2: count (jittered)")
    a.set_title("the data", color=INK)
    bins = np.linspace(0, 1, 21)
    for a, rr, ttl, colr in ((ax[1], r_plain, r"$F(x_2 \mid x_1)$, not randomised", GOLD),
                             (ax[2], r_rand, "randomised (Brockwell 2007)", GREEN)):
        a.hist(rr[:, col], bins=bins, density=True, color=colr, edgecolor="white", lw=0.5)
        a.axhline(1, color=INK, lw=0.9, ls="--")
        a.set_xlim(0, 1); a.set_xlabel("Rosenblatt transform, variable 2"); a.set_title(ttl, color=INK, fontsize=9.3)
    ax[1].set_ylabel("density")
    ax[2].set_ylim(ax[1].get_ylim())
    for i, a in enumerate(ax):
        style(a); tag(a, "abc"[i])
    f.tight_layout(w_pad=1.0)
    finish(f, folder, "python_discrete.png")

    # 3. tll on the mixture of python.md
    rng = np.random.default_rng(5)
    a_ = fv.Bicop("clayton", parameters=[4.0]).simulate(3000, seed=1)
    b_ = fv.Bicop("gumbel", 90, [3.0]).simulate(3000, seed=2)
    um = np.where(rng.random((3000, 1)) < 0.5, a_, b_)
    tll = fv.Bicop("tll").fit(um)
    st = fv.Bicop().select(um, family_set=["parametric", "nonparametric"], criterion="bic")
    f, ax = plt.subplots(1, 2, figsize=(7.2, 3.3))
    zm = norm.ppf(um)
    for a, (c, ttl, colr) in zip(ax, ((tll, f"TLL: loglik {tll.loglik():.0f}, edf {tll.npars:.0f} (AIC choice)", INK),
                                      (st, f"{E.FAMN.get(st.family_code, st.family)}: loglik {st.loglik():.0f}, "
                                           f"{st.npars:.0f} parameters (BIC choice)", INK))):
        a.scatter(zm[:, 0], zm[:, 1], s=2, color=BLUE, alpha=0.3, linewidths=0)
        lines(a, zgrid(c), colr, lw=0.8)
        zlim(a, 3.3)
        a.set_title(ttl, fontsize=9, color=INK)
        a.set_xlabel(r"$\Phi^{-1}(u_1)$")
        style(a)
    ax[0].set_ylabel(r"$\Phi^{-1}(u_2)$")
    tag(ax[0], "a"); tag(ax[1], "b")
    f.tight_layout(w_pad=1.5)
    finish(f, folder, "python_tll.png")

    # 4. truncation: log-likelihood and criteria against the truncation level
    rng = np.random.default_rng(4)
    d, n = 8, 2000
    zz = np.zeros((n, d))
    zz[:, 0] = rng.standard_normal(n); zz[:, 1] = 0.6 * zz[:, 0] + 0.8 * rng.standard_normal(n)
    for j in range(2, d):
        zz[:, j] = 0.5 * zz[:, j - 1] + 0.35 * zz[:, j - 2] + 0.6 * rng.standard_normal(n)
    uu = fv.to_pseudo_obs(zz)
    lvls = list(range(1, d))
    ll, bic, npar = [], [], []
    for t in lvls:
        vt = fv.Vinecop.from_data(uu, trunc_lvl=t)
        ll.append(vt.loglik()); bic.append(vt.bic()); npar.append(vt.npars)
    auto = fv.Vinecop.from_data(uu, select_trunc_lvl=True)
    f, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
    a = ax[0]
    a.plot(lvls, ll, color=BLUE, marker="o", ms=4.5, mec="white", lw=1.3, label="log-likelihood")
    a.set_xlabel("truncation level"); a.set_ylabel("log-likelihood")
    a.set_title("fit improves until tree 2, then stalls", color=INK, fontsize=9.3)
    a = ax[1]
    a.plot(lvls, bic, color=GREEN, marker="s", ms=4.5, mec="white", lw=1.3, label="BIC")
    a.axvline(auto.trunc, color=INK, lw=0.9, ls="--", label=f"select_trunc_lvl (mBICv): {auto.trunc}")
    a.set_xlabel("truncation level"); a.set_ylabel("BIC")
    a.set_title("BIC against the truncation level", color=INK, fontsize=9.3)
    for (x, y, k_) in zip(lvls, bic, npar):
        a.annotate(str(k_), (x, y), xytext=(5, 6), textcoords="offset points", ha="left", fontsize=7.2,
                   color=MUTED)
    a.legend(loc="upper right", fontsize=8)
    for i, a in enumerate(ax):
        style(a); tag(a, "ab"[i]); a.set_xticks(lvls); a.set_xlim(0.7, 7.4)
    f.tight_layout(w_pad=1.5)
    finish(f, folder, "python_trunc.png")
    print("  trunc ll", np.round(ll, 1), "bic", np.round(bic, 1), "npars", npar, "auto", auto.trunc)

    # 5. non-simplified vine: the relaxed edge's parameter as a function of the conditioning value
    rng = np.random.default_rng(1)
    n = 2000
    u2, w1, w3 = rng.random(n), rng.random(n), rng.random(n)
    x1 = fv.Bicop("clayton", parameters=[2.0]).hinv1(np.column_stack([u2, w1]))
    rho = -0.6 + 1.4 * u2
    v3 = norm.cdf(rho * norm.ppf(w1) + np.sqrt(1 - rho**2) * norm.ppf(w3))
    x3 = fv.Bicop("frank", parameters=[5.0]).hinv1(np.column_stack([u2, v3]))
    ur = np.column_stack([x1, u2, x3])
    vc = fv.Vinecop().set_dvine([1, 2, 3]).select(ur, family_set=["indep", "gaussian", "clayton", "frank"])
    simp = vc.get_pair_copula(1, 0)
    ll_s = vc.loglik()
    vc.relax(ur, target=0, mode="test")
    t_, e_ = vc.relaxed_edges[0]
    fam = int(vc.fam[t_, e_])
    wk, pk = vc.wk[0], vc.pk[0, :, 0]
    uk = norm.cdf(wk)                               # loadings (0, 1, 0): w = Phi^-1(u2)
    tau_true = lambda uu_: 2 / np.pi * np.arcsin(-0.6 + 1.4 * uu_)
    tau_anchor = [tau_of(fv.Bicop(fam, 0, [p_]), 20000) for p_ in pk]
    tau_simp = tau_of(simp, 20000)
    f, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
    a = ax[0]
    gg = np.linspace(0.001, 0.999, 200)
    a.plot(gg, tau_true(gg), color=INK, lw=1.0, ls="--", label=r"truth: Gaussian, $\rho = -0.6 + 1.4\,u_2$")
    a.axhline(tau_simp, color=GREEN, lw=1.3, label=f"simplified ({E.FAMN.get(simp.family_code, simp.family)})")
    a.plot(np.r_[0, uk, 1], np.r_[tau_anchor[0], tau_anchor, tau_anchor[-1]], color=BLUE, lw=1.3,
           marker="o", ms=4.5, mec="white", markevery=list(range(1, len(uk) + 1)),
           label=f"relaxed {E.FAMN.get(fam, fam)} at the anchors")
    a.set_xlabel(r"conditioning value $u_2$"); a.set_ylabel(r"conditional Kendall $\tau$ of $(U_1, U_3) \mid U_2$")
    a.set_xlim(0, 1)
    a.legend(loc="upper left", fontsize=7.6)
    # predictive density of variable 3 at three values of u2 (pdf_grid), simplified vs relaxed
    a = ax[1]
    vs = fv.Vinecop().set_dvine([1, 2, 3]).select(ur, family_set=["indep", "gaussian", "clayton", "frank"])
    grid3 = np.linspace(0.005, 0.995, 199)
    for u2v, colr in ((0.1, RAMP[1]), (0.5, RAMP[3]), (0.9, RAMP[4])):
        g = np.zeros((1, grid3.size, 3)); g[0, :, 0] = 0.8; g[0, :, 1] = u2v; g[0, :, 2] = grid3
        a.plot(grid3, vc.pdf_grid(g, target=3)[0], color=colr, lw=1.3, label=f"$u_2 = {u2v}$, relaxed")
        a.plot(grid3, vs.pdf_grid(g, target=3)[0], color=colr, lw=1.0, ls=":")
    a.plot([], [], color=GREY, lw=1.0, ls=":", label="simplified (same colours)")
    a.set_xlabel("$u_3$"); a.set_ylabel(r"vine density at $(0.8, u_2, u_3)$")
    a.legend(loc="upper left", fontsize=7.6)
    for i, a in enumerate(ax):
        style(a); tag(a, "ab"[i])
    f.tight_layout(w_pad=1.5)
    finish(f, folder, "python_relax.png")
    print("  relax: loglik", round(ll_s, 2), "->", round(vc.loglik(), 2), "anchors", np.round(pk, 3),
          "tau", np.round(tau_anchor, 3), "simplified tau", round(tau_simp, 3))


JOBS = {"concepts": concepts, "quickstart": quickstart, "numerics": numerics, "scores": scores,
        "validation": validation, "families": families, "nonparametric": nonparametric, "python": python}

if __name__ == "__main__":
    for job in (sys.argv[1:] or list(JOBS)):
        JOBS[job]()
