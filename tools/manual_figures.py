"""Figures of the user-manual pages (docs/manual/*.md), computed by fortvine itself.

    pip install ./python                     # the fortvine Python package (or FORTVINE_LIB=...)
    python tools/manual_figures.py [bicop | rotations | fitting | controls | vinecop | discrete |
                                    nonparametric | generic | performance]      (default: all)

Every density, h-function, fit, selection, Rosenblatt transform, simulation and timing comes from
the compiled library through the `fortvine` Python package; this script only arranges the calls
and plots.  Normal scores z = Phi^-1(u) (scipy) are used for display, as in tools/figures.py,
whose style and helpers (rcParams, panel, save, CMAP, colours) are reused, together with the
small plotting helpers of tools/example_figures.py.  Seeds are fixed; timings depend on the
machine.  No job uses more than 16 threads.  PNGs go to docs/manual/img/.
"""
import os, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figures as F                                       # noqa: E402  (rcParams, panel, save, CMAP)
import example_figures as E                               # noqa: E402  (style, tag, draw_tree)
import matplotlib.pyplot as plt                           # noqa: E402
from scipy.stats import norm, poisson                     # noqa: E402
import fortvine as fv                                     # noqa: E402

ROOT = F.ROOT
OUT = os.path.join(ROOT, "docs", "manual", "img")
INK, MUTED, GRIDC = E.INK, E.MUTED, E.GRIDC
BLUE, GREEN, GOLD, GREY, RED = E.BLUE, E.GREEN, E.GOLD, E.GREY, E.RED
LBLUE = "#7FA8D6"                                         # light blue, a fourth series
Z1, Z2, PHI = F.Z1, F.Z2, F.PHI
UG = np.column_stack([norm.cdf(Z1).ravel(), norm.cdf(Z2).ravel()])     # the 121 x 121 normal-score grid in u
MAXT = 16                                                 # thread cap of every call
NT = 8                                                    # threads of the non-timing jobs
os.environ.setdefault("FORTVINE_THREADS", str(NT))
style, tag = E.style, E.tag


def finish(f, name, dpi=150):
    os.makedirs(OUT, exist_ok=True)
    p = os.path.join(OUT, name)
    F.save(f, p, dpi=dpi, bbox_inches="tight", pad_inches=0.05, facecolor="white")
    plt.close(f)
    print(f"{name:30s} {os.path.getsize(p) // 1024:4d} KB")


def zdens(cop):
    """Copula density of a Bicop on the normal-score grid (F.panel multiplies by phi(z1) phi(z2))."""
    return cop.pdf(UG).reshape(Z1.shape)


def zlines(a, g, color, lw=0.8, nlev=7, ls="-", top=None):
    d = g * PHI
    lev = np.linspace(0, d.max() if top is None else top, nlev + 2)[1:-1]
    a.contour(Z1, Z2, d, levels=lev, colors=color, linewidths=lw, linestyles=ls)


def zaxes(a, lim=3.0, labels=True):
    a.set_xlim(-lim, lim); a.set_ylim(-lim, lim); a.set_aspect("equal")
    if labels:
        a.set_xlabel(r"$z_1 = \Phi^{-1}(u_1)$"); a.set_ylabel(r"$z_2 = \Phi^{-1}(u_2)$")


def timeit(fn, reps=1):
    """Best of `reps` wall-clock times of fn() in seconds, and the last result."""
    best, r = np.inf, None
    for _ in range(reps):
        t0 = time.perf_counter()
        r = fn()
        best = min(best, time.perf_counter() - t0)
    return best, r


def cvine_example():
    """The hand-written C-vine of the vinecop page (roots 1, 2, 3)."""
    cv = fv.Vinecop().set_cvine([4, 3, 2, 1])
    for e in range(3):
        cv.set_pair_copula(0, e, fv.Bicop("clayton", 0, [2.0]))
    for e in range(2):
        cv.set_pair_copula(1, e, fv.Bicop("frank", 0, [-3.0]))
    cv.set_pair_copula(2, 0, fv.Bicop("gumbel", 90, [1.4]))
    return cv


# ------------------------------------------------------------------------------------------ bicop
def bicop():
    cop = fv.Bicop("gumbel", 180, [1.8])                       # the copula of the page's example
    u = cop.simulate(3000, seed=3)
    # 1. density and distribution function
    f, ax = plt.subplots(1, 2, figsize=(7.0, 3.3))
    a = ax[0]
    F.panel(a, zdens(cop))
    zz = norm.ppf(u[:400])
    a.scatter(zz[:, 0], zz[:, 1], s=1.6, color=INK, alpha=0.45, linewidths=0, zorder=4)
    zaxes(a)
    a.set_title(r"density $c(u_1, u_2)$ in normal scores, 400 draws", color=INK)
    a = ax[1]
    g = np.linspace(0.0005, 0.9995, 161)
    G1, G2 = np.meshgrid(g, g)
    C = cop.cdf(np.column_stack([G1.ravel(), G2.ravel()])).reshape(G1.shape)
    lev = np.arange(0.1, 1.0, 0.1)
    a.contour(G1, G2, C, levels=lev, colors=BLUE, linewidths=1.1)
    a.contour(G1, G2, G1 * G2, levels=lev, colors=GREY, linewidths=0.7, linestyles="--")
    pts = np.array([[0.2, 0.3], [0.5, 0.5], [0.8, 0.9]])
    a.plot(pts[:, 0], pts[:, 1], "o", ms=4.5, mfc=GOLD, mec="white", mew=0.8, zorder=5)
    for (x, y), c in zip(pts, cop.cdf(pts)):
        a.annotate(f"$C = {c:.4f}$", (x, y), xytext=(5, -10), textcoords="offset points", fontsize=7.6, color=INK,
                   ha="left" if x < 0.7 else "right", bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.9))
    a.plot([], [], color=BLUE, lw=1.1, label=r"$C(u_1, u_2)$, rotated Gumbel")
    a.plot([], [], color=GREY, lw=0.7, ls=(0, (4, 3)), label=r"independence $u_1 u_2$")
    a.plot([], [], "o", ms=4.5, mfc=GOLD, mec="white", label="points of the example")
    a.set_xlim(0, 1); a.set_ylim(0, 1); a.set_aspect("equal")
    a.set_xlabel("$u_1$"); a.set_ylabel("$u_2$")
    a.set_title("distribution function: contours at 0.1, ..., 0.9", color=INK)
    a.legend(loc="upper left", fontsize=7.3, frameon=True, facecolor="white", edgecolor="none", framealpha=0.95,
             borderaxespad=0.2, handlelength=1.6)
    style(a, grid=False)
    tag(ax[0], "a"); tag(ax[1], "b")
    f.tight_layout(w_pad=2.0)
    finish(f, "bicop_density_cdf.png")

    # 2. h-functions, their inverses, and the round-trip error over many families
    v = np.linspace(0.0005, 0.9995, 400)
    ws = np.arange(0.1, 0.91, 0.1)
    f = plt.figure(figsize=(7.2, 3.15))
    gs = f.add_gridspec(1, 3, width_ratios=[1, 1, 1.15], wspace=0.42)
    ax = [f.add_subplot(gs[0, k]) for k in range(3)]
    cols = (BLUE, GREEN, GOLD)
    for c0, col in zip((0.1, 0.5, 0.9), cols):
        h1 = cop.hfunc1(np.column_stack([np.full_like(v, c0), v]))
        x1 = cop.hinv1(np.column_stack([np.full_like(ws, c0), ws]))
        ax[0].plot(v, h1, color=col, lw=1.4, label=f"$u_1 = {c0}$")
        ax[0].plot(x1, ws, "o", ms=3.4, mfc="white", mec=col, mew=0.9)
        h2 = cop.hfunc2(np.column_stack([v, np.full_like(v, c0)]))
        x2 = cop.hinv2(np.column_stack([ws, np.full_like(ws, c0)]))
        ax[1].plot(v, h2, color=col, lw=1.4, label=f"$u_2 = {c0}$")
        ax[1].plot(x2, ws, "o", ms=3.4, mfc="white", mec=col, mew=0.9)
    ax[0].set_xlabel("$u_2$"); ax[0].set_ylabel(r"$h_1(u_1, u_2) = P(U_2 \leq u_2 \mid U_1 = u_1)$", fontsize=8.6)
    ax[1].set_xlabel("$u_1$"); ax[1].set_ylabel(r"$h_2(u_1, u_2) = P(U_1 \leq u_1 \mid U_2 = u_2)$", fontsize=8.6)
    ax[0].set_title(r"$h_1$ (lines) and hinv1 (dots)", color=INK, fontsize=8.8)
    ax[1].set_title(r"$h_2$ (lines) and hinv2 (dots)", color=INK, fontsize=8.8)
    for a in ax[:2]:
        a.set_xlim(0, 1); a.set_ylim(0, 1)
        a.legend(loc="lower right", fontsize=7.4, handlelength=1.4, labelspacing=0.25, borderaxespad=0.2)
        style(a)
    specs = [("gaussian", 0, [0.7]), ("student", 0, [0.7, 4]), ("clayton", 0, [3.0]), ("gumbel", 180, [1.8]),
             ("frank", 0, [8.0]), ("joe", 0, [3.0]), ("bb1", 0, [0.8, 1.6]), ("bb6", 0, [1.5, 1.5]),
             ("bb7", 0, [1.8, 1.2]), ("bb8", 0, [4.0, 0.7]), ("tawn", 0, [0.6, 0.9, 4.0]),
             ("plackett", 0, [10.0]), ("galambos", 0, [1.5]), ("husler_reiss", 0, [2.0]),
             (10203, 0, [1.0, 1.5]), (20200, 0, [3.0, 0.95, 0.5])]
    g = np.linspace(0.001, 0.999, 60)
    A, W = np.meshgrid(g, g)
    pts = np.column_stack([A.ravel(), W.ravel()])
    names, e1, e2 = [], [], []
    for fam, rot, p in specs:
        c = fv.Bicop(fam, rot, p)
        x = c.hinv1(pts)
        e1.append(np.abs(c.hfunc1(np.column_stack([pts[:, 0], x])) - pts[:, 1]).max())
        y = c.hinv2(pts[:, ::-1].copy())
        e2.append(np.abs(c.hfunc2(np.column_stack([y, pts[:, 0]])) - pts[:, 1]).max())
        nm = c.family if not isinstance(fam, int) else str(fam)
        names.append(nm.replace("husler_reiss", "Hüsler-Reiss").replace("_", " ") + (f" {rot}°" if rot else ""))
    a = ax[2]
    yy = np.arange(len(specs))[::-1]
    for y in yy:
        a.axhline(y, color=GRIDC, lw=0.6, zorder=0)
    a.scatter(np.maximum(e1, 1e-17), yy, s=20, color=BLUE, zorder=3, label=r"$|h_1(u_1, h_1^{-1}(u_1, w)) - w|$")
    a.scatter(np.maximum(e2, 1e-17), yy, s=20, marker="D", facecolor="none", edgecolor=GREEN, lw=0.9, zorder=4,
              label=r"$|h_2(h_2^{-1}(w, u_2), u_2) - w|$")
    a.set_xscale("log"); a.set_xlim(1e-17, 1e-8)
    a.set_yticks(yy); a.set_yticklabels(names, fontsize=7.4)
    a.tick_params(axis="y", length=0); a.tick_params(axis="x", labelsize=8)
    for sp in ("top", "right", "left"):
        a.spines[sp].set_visible(False)
    a.set_xlabel(r"largest round-trip error, $60 \times 60$ grid", fontsize=8.6)
    a.set_title("inverses are exact to rounding", color=INK, fontsize=8.8)
    a.legend(loc="upper center", bbox_to_anchor=(0.45, -0.17), fontsize=7.1, handlelength=1.0, labelspacing=0.2)
    for a, t in zip(ax, "abc"):
        tag(a, t)
    finish(f, "bicop_hfunc.png")

    # 3. selection: the criterion of every (family, rotation) candidate on the example's sample
    rng = np.random.default_rng(3)                         # the Python tab of the page's example
    u1, w = rng.random(3000), rng.random(3000)
    u = np.column_stack([u1, cop.hinv1(np.column_stack([u1, w]))])
    rotless = {0, 1, 2, 5, 12, 13, 14}
    cand = []
    for code in fv.family_group("default"):
        for r in ((0,) if code in rotless else (0, 180)):
            c = fv.Bicop(code, r).fit(u, num_threads=NT)
            cand.append((c.family, r, c.loglik(), c.npars))
    n = len(u)
    bic = np.array([-2 * ll + np.log(n) * k for _, _, ll, k in cand])
    aic = np.array([-2 * ll + 2 * k for _, _, ll, k in cand])
    sb = fv.Bicop().select(u, num_threads=NT)
    sa = fv.Bicop().select(u, criterion="aic", num_threads=NT)
    order = np.argsort(bic)[:16]
    f, ax = plt.subplots(1, 1, figsize=(7.0, 3.3))
    a = ax
    yy = np.arange(len(order))[::-1]
    for y in yy:
        a.axhline(y, color=GRIDC, lw=0.6, zorder=0)
    a.scatter(bic[order] - bic.min(), yy, s=26, color=BLUE, zorder=3, label=r"BIC $-$ smallest BIC")
    a.scatter(aic[order] - aic.min(), yy, s=26, marker="D", facecolor="none", edgecolor=GREEN, lw=1.0, zorder=4,
              label=r"AIC $-$ smallest AIC")
    lab = []
    for i in order:
        nm, r, ll, k = cand[i]
        nm = nm.replace("husler_reiss", "Hüsler-Reiss")
        s = f"{nm} {r}°" if r else nm
        mark = []
        if nm == sb.family and r == sb.rotation:
            mark.append("BIC")
        if nm == sa.family and r == sa.rotation:
            mark.append("AIC")
        lab.append(s + (f"  ({' and '.join(mark)} choice)" if mark else "") + f"   k = {k}")
    a.set_yticks(yy); a.set_yticklabels(lab, fontsize=7.8)
    a.set_xscale("symlog", linthresh=1.0); a.set_xlim(-0.2, 3000)
    a.set_xlabel("criterion difference to the best candidate (symmetric log scale)")
    a.tick_params(axis="y", length=0)
    for sp in ("top", "right", "left"):
        a.spines[sp].set_visible(False)
    a.legend(loc="lower right", fontsize=8)
    a.set_title(f"the 16 best of {len(cand)} candidates, $n = {n}$ draws of the 180° Gumbel, "
                r"$\theta = 1.8$", color=INK)
    f.tight_layout()
    finish(f, "bicop_select.png")


# ------------------------------------------------------------------------------------------ rotations
def rotations():
    rng = np.random.default_rng(1)
    x1, w = rng.random(2000), rng.random(2000)
    f, ax = plt.subplots(2, 4, figsize=(7.2, 3.9))
    for i, (fam, th) in enumerate((("clayton", 3.0), ("gumbel", 2.5))):
        for j, r in enumerate((0, 90, 180, 270)):
            c = fv.Bicop(fam, r, [th])
            a = ax[i][j]
            F.panel(a, zdens(c), ticks=False)
            x2 = c.hinv1(np.column_stack([x1, w]))
            t = fv.kendall_tau(x1, x2)
            pos = r in (0, 180)                                # positive tau: the upper-left corner is empty
            a.text(0.04, 0.96 if pos else 0.04, rf"$\hat\tau = {t:+.2f}$", transform=a.transAxes, fontsize=7.8,
                   color=INK, va="top" if pos else "bottom")
            if i == 0:
                a.set_title(f"rotation {r}°", color=INK, fontsize=9.5)
            if j == 0:
                a.set_ylabel(f"{fam.capitalize()}, $\\theta = {th:g}$", fontsize=9.5, color=INK)
    f.tight_layout(pad=0.3, w_pad=0.5, h_pad=0.6)
    finish(f, "rotations_grid.png")

    # transposition of asymmetric generic families with gen_flip
    cases = [(10006, 0, [3.0, 0.3]), (10011, 90, [0.9, 0.9]), ("tawn", 0, [0.3, 0.9, 4.0])]
    f, ax = plt.subplots(2, 3, figsize=(6.4, 4.5))
    uv = np.random.default_rng(0).random((1000, 2))
    for j, (fam, r, p) in enumerate(cases):
        c = fv.Bicop(fam, r, p)
        if c.family_code >= 10000:
            code, q = fv.gen_flip(c.family_code, c.parameters)
        else:                                              # Tawn: swap psi1 and psi2 (the page's table)
            code, q = c.family_code, c.parameters[[1, 0, 2]]
        rt = {90: 270, 270: 90}.get(r, r)
        t = fv.Bicop(int(code), rt, q)
        err = np.abs(c.pdf(uv) - t.pdf(uv[:, ::-1].copy())).max()
        for i, (b, ttl) in enumerate(((c, "original"), (t, "transpose"))):
            a = ax[i][j]
            F.panel(a, zdens(b), ticks=False)
            a.plot([-3, 3], [-3, 3], color=GREY, lw=0.6, ls=":")
            ps = ", ".join(f"{x:.3g}" for x in b.parameters)
            a.set_title(f"{ttl}: {b.family} ({b.family_code})\n" + (f"rotation {b.rotation}°, " if b.rotation else "")
                        + f"par ({ps})", fontsize=8.4, color=INK)
        ax[1][j].set_xlabel(f"max |c(u,v) - c_T(v,u)| = {err:.1e}", fontsize=7.8, color=MUTED)
    f.tight_layout(pad=0.3, w_pad=0.8, h_pad=0.8)
    finish(f, "rotations_transpose.png")


# ------------------------------------------------------------------------------------------ fitting
def fitting():
    # 1. the log-likelihood around the estimate: a one-parameter profile and a two-parameter surface
    u = fv.Bicop("gumbel", 0, [1.6]).simulate(500, seed=4)
    mle = fv.Bicop("gumbel").fit(u, num_threads=1)
    itau = fv.Bicop("gumbel").fit(u, method="itau", num_threads=1)
    tau = fv.kendall_tau(u[:, 0], u[:, 1])
    th = np.linspace(1.15, 2.2, 300)
    ll = np.array([fv.Bicop("gumbel", 0, [x]).loglik(u) for x in th])
    f, ax = plt.subplots(1, 2, figsize=(7.2, 3.2), gridspec_kw=dict(width_ratios=[1.1, 1]))
    a = ax[0]
    lo, hi = 1 / (1 - (tau - 0.1)), 1 / (1 - (tau + 0.1))      # Gumbel tau inversion (table of this page)
    a.axvspan(lo, hi, color="#E3EDF8", lw=0, zorder=0, label=r"search box: $\theta(\hat\tau \pm 0.1)$")
    a.plot(th, ll - mle.loglik(), color=BLUE, lw=1.5, label=r"$\ell(\theta) - \ell(\hat\theta)$")
    a.axvline(1.6, color=GREY, lw=0.8, ls=":", label=r"true $\theta = 1.6$")
    a.plot(mle.parameters[0], 0, "o", ms=6, color=BLUE, mec="white", zorder=4,
           label=rf"MLE $\hat\theta = {mle.parameters[0]:.4f}$")
    a.plot(itau.parameters[0], itau.loglik() - mle.loglik(), "D", ms=5.5, mfc="white", mec=GREEN, mew=1.2,
           zorder=5, label=rf"$\tau$-inversion $= {itau.parameters[0]:.4f}$")
    a.set_xlim(th[0], th[-1]); a.set_ylim(-60, 6)
    a.set_xlabel(r"Gumbel $\theta$"); a.set_ylabel("log-likelihood relative to the MLE")
    a.set_title(rf"one-parameter family, $n = 500$, $\hat\tau = {tau:.3f}$", color=INK)
    a.legend(loc="lower center", fontsize=7.4, handlelength=1.5, labelspacing=0.25, frameon=True, facecolor="white",
             edgecolor="none", framealpha=0.95)
    style(a)
    a = ax[1]
    ub = fv.Bicop("bb1", 0, [0.6, 1.4]).simulate(1000, seed=6)
    b = fv.Bicop("bb1").fit(ub, num_threads=1)
    t1, t2 = np.linspace(0.15, 1.25, 70), np.linspace(1.1, 1.8, 70)
    T1, T2 = np.meshgrid(t1, t2)
    L = np.array([fv.Bicop("bb1", 0, [x, y]).loglik(ub) for x, y in zip(T1.ravel(), T2.ravel())]).reshape(T1.shape)
    L -= b.loglik()
    lev = [-40, -20, -10, -5, -2, -0.5]
    a.contourf(T1, T2, L, levels=[-1e9] + lev + [1], colors=[F.CMAP(x) for x in np.linspace(0.0, 0.8, len(lev) + 1)])
    cs = a.contour(T1, T2, L, levels=lev, colors="#0D2C52", linewidths=0.5)
    a.clabel(cs, fmt="%g", fontsize=6.8)
    a.plot(*b.parameters, "o", ms=6, color=GOLD, mec="white", mew=0.9, zorder=4,
           label=rf"MLE $({b.parameters[0]:.3f}, {b.parameters[1]:.3f})$")
    a.plot(0.6, 1.4, "+", ms=9, color=INK, mew=1.3, zorder=4, label=r"true $(0.6, 1.4)$")
    a.set_xlabel(r"BB1 $\theta$"); a.set_ylabel(r"BB1 $\delta$")
    a.set_title(r"two parameters: $\ell(\theta, \delta) - \ell_{\max}$, $n = 1000$", color=INK)
    a.legend(loc="upper right", fontsize=7.4, frameon=True, facecolor="white", edgecolor="none", framealpha=0.9)
    style(a, grid=False)
    tag(ax[0], "a"); tag(ax[1], "b")
    f.tight_layout(w_pad=1.8)
    finish(f, "fitting_loglik.png")

    # 2. accuracy of MLE and tau inversion, and the cost of one fit against n
    ns = [100, 200, 500, 1000, 2000, 5000]
    R = 200
    rmse = {}
    for fam, th0 in (("clayton", 2.0), ("gumbel", 2.0)):
        truth = fv.Bicop(fam, 0, [th0])
        for meth in ("mle", "itau"):
            r = []
            for n in ns:
                est = [fv.Bicop(fam).fit(truth.simulate(n, seed=1000 * n + k), method=meth, num_threads=1)
                       .parameters[0] for k in range(R)]
                r.append(np.sqrt(np.mean((np.array(est) - th0) ** 2)))
            rmse[fam, meth] = np.array(r)
    nt = [250, 500, 1000, 1500, 1800, 1801, 2500, 5000, 10000, 20000, 50000]
    big = fv.Bicop("bb1", 0, [0.6, 1.4]).simulate(50000, seed=8)
    tcases = [("clayton", "mle", BLUE, "-", "Clayton, MLE (Newton, closed form)"),
              ("clayton", "itau", BLUE, "--", r"Clayton, $\tau$-inversion"),
              ("student", "mle", GREEN, "-", "Student, MLE (profile in $\\nu$)"),
              ("bb1", "mle", GOLD, "-", "BB1, MLE (Fisher scoring)"),
              ("tawn", "mle", INK, "-", "Tawn, MLE")]
    times = {}
    for fam, meth, *_ in tcases:
        times[fam, meth] = [timeit(lambda: fv.Bicop(fam).fit(big[:n], method=meth, num_threads=1), reps=5)[0]
                            for n in nt]
    f, ax = plt.subplots(1, 2, figsize=(7.2, 3.2))
    a = ax[0]
    for fam, col in (("clayton", BLUE), ("gumbel", GREEN)):
        a.plot(ns, rmse[fam, "mle"], "-o", color=col, ms=4, mec="white", label=f"{fam.capitalize()}, MLE")
        a.plot(ns, rmse[fam, "itau"], "--s", color=col, ms=3.6, mfc="white", label=rf"{fam.capitalize()}, $\tau$-inv.")
    a.plot(ns, 2.2 * np.array(ns, float) ** -0.5, color=GREY, lw=0.7, ls=":", label=r"slope $n^{-1/2}$")
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel("$n$"); a.set_ylabel(r"RMSE of $\hat\theta$")
    a.set_title(rf"true $\theta = 2$, {R} samples per $n$", color=INK)
    a.legend(fontsize=7.4, loc="lower left", handlelength=2.0, labelspacing=0.25)
    style(a)
    a = ax[1]
    a.axvspan(1800, 60000, color="#F1F3F8", lw=0, zorder=0)
    a.text(2100, 0.04, "two-stage fits\n(600-row search, full polish)", fontsize=7.4, color=MUTED, va="bottom",
           transform=a.get_xaxis_transform())
    for fam, meth, col, ls, lab in tcases:
        a.plot(nt, 1e3 * np.array(times[fam, meth]), color=col, ls=ls, marker="o", ms=3, mec="white", label=lab)
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(200, 60000)
    a.set_xlabel("$n$"); a.set_ylabel("time of one fit [ms], 1 thread")
    a.set_title("cost of a fit (BB1 sample)", color=INK)
    a.legend(fontsize=7.0, loc="upper left", handlelength=2.0, labelspacing=0.25)
    style(a)
    tag(ax[0], "a"); tag(ax[1], "b")
    f.tight_layout(w_pad=1.8)
    finish(f, "fitting_accuracy_cost.png")

    # 3. the page's example: what the options cost in a vine selection (4 threads)
    d = 5
    truth = fv.Vinecop().set_dvine([1, 2, 3, 4, 5])
    for e in range(d - 1):
        truth.set_pair_copula(0, e, fv.Bicop("bb1", 0, [0.6, 1.4]))
    for e in range(d - 2):
        truth.set_pair_copula(1, e, fv.Bicop("student", 0, [0.3, 5.0]))
    u = truth.simulate(5000, seed=5)
    runs = [("MLE, BIC (default)", {}), ("fast_mle", dict(fast_mle=True)), ("itau", dict(method="itau")),
            ("MLE, mBICv", dict(criterion="mbicv")), ("fit cache, first pass", None),
            ("fit cache, second pass", None)]
    res = []
    for lab, opt in runs:
        if lab.startswith("fit cache, first"):
            fv.fitcache_enable(20000)
        t, vc = timeit(lambda: fv.Vinecop.from_data(u, num_threads=4, **(opt or {})), reps=1 if opt is None else 3)
        res.append((lab, t, vc.loglik(u), vc.npars))
    st = fv.fitcache_stats()
    fv.fitcache_disable()
    f, ax = plt.subplots(1, 2, figsize=(7.2, 2.6), sharey=True)
    yy = np.arange(len(res))[::-1]
    cols = [BLUE, BLUE, GREEN, BLUE, GOLD, GOLD]
    ax[0].barh(yy, [r[1] for r in res], color=cols, height=0.62)
    for y, r in zip(yy, res):
        ax[0].text(r[1], y, f" {r[1]:.3f} s", va="center", fontsize=7.6, color=INK)
    ax[0].set_xlabel("selection time [s], 4 threads"); ax[0].set_xlim(0, max(r[1] for r in res) * 1.3)
    best = max(r[2] for r in res)
    ax[1].barh(yy, [r[2] - best for r in res], color=cols, height=0.62)
    for y, r in zip(yy, res):
        ax[1].text(r[2] - best, y, f" {r[2]:.1f} ({r[3]} par.)" if r[2] < best - 1 else f"  {r[2]:.1f} ({r[3]} par.)",
                   va="center", ha="right" if r[2] < best - 1 else "left", fontsize=7.6, color=INK)
    ax[1].set_xlabel("log-likelihood $-$ the largest"); ax[1].set_xlim((min(r[2] for r in res) - best) * 1.15, 60)
    ax[0].set_yticks(yy); ax[0].set_yticklabels([r[0] for r in res])
    for a in ax:
        style(a); a.tick_params(axis="y", length=0); a.grid(False, axis="y")
    ax[0].set_title(f"$d = {d}$, $n = 5000$; cache: {st.hits} hits, {st.misses} misses", color=INK, fontsize=9)
    ax[1].set_title("fit of the selected vine", color=INK, fontsize=9)
    tag(ax[0], "a"); tag(ax[1], "b")
    f.tight_layout(w_pad=1.0)
    finish(f, "fitting_options.png")


# ------------------------------------------------------------------------------------------ controls
def decaying_dvine(d=10, n=1000, seed=7):
    """A D-vine whose dependence fades with the tree: tau 0.5, 0.25, 0.12, 0.06, then 0.03."""
    truth = fv.Vinecop().set_dvine(list(range(1, d + 1)))
    fams = ("clayton", "gumbel", "gaussian")
    taus = [0.5, 0.25, 0.12, 0.06] + [0.03] * (d - 5)
    for t in range(d - 1):
        for e in range(d - 1 - t):
            fam, tau = fams[(t + e) % len(fams)], taus[t]
            p = {"clayton": 2 * tau / (1 - tau), "gumbel": 1 / (1 - tau),      # the inversions of the fitting page
                 "gaussian": np.sin(np.pi * tau / 2)}[fam]
            truth.set_pair_copula(t, e, fv.Bicop(fam, 0, [p]))
    return truth, truth.simulate(n, seed=seed)


def controls():
    d, n = 10, 1000
    truth, u = decaying_dvine(d, n)
    # 1. what each criterion keeps
    crits = [("loglik", "log-likelihood", GREY), ("aic", "AIC", GOLD), ("bic", "BIC", BLUE),
             ("mbic", r"mBIC ($\psi_0 = 0.9$)", GREEN), ("mbicv", "mBICv", INK)]
    fits = {c: fv.Vinecop.from_data(u, criterion=c, num_threads=NT) for c, _, _ in crits}
    f, ax = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw=dict(width_ratios=[1.25, 1]))
    a = ax[0]
    tt = np.arange(1, d)
    a.plot(tt, d - tt, color=GRIDC, lw=6, solid_capstyle="round", zorder=0, label="edges in the tree")
    for k, (c, lab, col) in enumerate(crits):
        vc = fits[c]
        q = [(vc.fam[t, :d - 1 - t] != 0).sum() for t in range(d - 1)]
        a.plot(tt + (k - 2) * 0.09, q, "-o", color=col, ms=3.6, mec="white", mew=0.6, lw=1.1, label=lab)
    a.set_xticks(tt); a.set_xlabel("tree $t$"); a.set_ylabel("non-independence edges")
    a.set_title(r"edges kept; true $\tau$ per tree: 0.5, 0.25, 0.12, 0.06, 0.03, ...", color=INK, fontsize=9)
    a.legend(fontsize=7.4, loc="upper right", labelspacing=0.25, handlelength=1.8)
    style(a)
    a = ax[1]
    pts = {}
    for c, lab, col in crits:
        vc = fits[c]
        key = (vc.npars, round(vc.loglik(u), 6))
        a.plot(vc.npars, vc.loglik(u), "o", ms=6.5 + 2.5 * (key in pts), color=col, mec="white", mew=0.8,
               zorder=3 - 0.1 * (key in pts))
        pts.setdefault(key, []).append(lab)
    for (k, ll), labs in pts.items():
        a.annotate(", ".join(labs), (k, ll), xytext=(7, -3), textcoords="offset points", fontsize=7.6, color=INK)
    a.plot(truth.npars, truth.loglik(u), "*", ms=10, color=RED, mec="white", mew=0.5, zorder=3)
    a.annotate(f"true vine ({truth.npars} par.)", (truth.npars, truth.loglik(u)), xytext=(6, -3),
               textcoords="offset points", fontsize=7.6, color=INK)
    a.set_xlabel("number of parameters"); a.set_ylabel("log-likelihood on the sample")
    a.set_xlim(None, a.get_xlim()[1] + 14)
    a.set_title(f"$d = {d}$, $n = {n}$", color=INK)
    style(a)
    tag(ax[0], "a"); tag(ax[1], "b")
    f.tight_layout(w_pad=1.8)
    finish(f, "controls_selcrit.png")

    # 2. truncation level and threshold: log-likelihood against its BIC penalty
    ks = list(range(1, d))
    tr = [fv.Vinecop.from_data(u, trunc_lvl=k, num_threads=NT) for k in ks]
    st = fv.Vinecop.from_data(u, select_trunc_lvl=True, num_threads=NT)
    cs = np.round(np.arange(0.0, 0.301, 0.02), 3)
    th = [fv.Vinecop.from_data(u, threshold=c, num_threads=NT) for c in cs]
    sth = fv.Vinecop.from_data(u, select_threshold=True, criterion="mbicv", num_threads=NT)
    f, ax = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
    for a, xs, vs, chosen, xl, cl in ((ax[0], ks, tr, st.trunc, "truncation level trunc_lvl",
                                       f"select_trunc_lvl: {st.trunc}"),
                                      (ax[1], cs, th, sth.threshold, "threshold",
                                       f"select_threshold: {sth.threshold:.3f}")):
        ll = np.array([v.loglik(u) for v in vs])
        pen = -0.5 * np.array([v.bic(u) for v in vs])
        a.plot(xs, ll, "-o", color=BLUE, ms=3.6, mec="white", label=r"log-likelihood $\ell$")
        a.plot(xs, pen, "-s", color=GREEN, ms=3.4, mec="white", label=r"$-\mathrm{BIC}/2 = \ell - \frac{1}{2} k \log n$")
        a.axvline(chosen, color=GOLD, lw=1.2, ls="--", label=cl + " (mBICv)")
        for x, v, y in zip(xs, vs, ll):
            if x in (xs[0], xs[len(xs) // 2], xs[-1]) or abs(x - chosen) < 1e-9:
                a.annotate(f"{v.npars}", (x, y), xytext=(0, 5), textcoords="offset points", ha="center",
                           fontsize=7, color=MUTED)
        a.set_xlabel(xl)
        a.legend(fontsize=7.4, loc="lower right" if a is ax[0] else "lower left", labelspacing=0.25)
        style(a)
    ax[0].set_xticks(ks)
    ax[0].set_ylabel("log-likelihood units")
    ax[0].set_title("trees fitted (grey numbers: parameters)", color=INK, fontsize=9)
    ax[1].set_title(r"edges with $|\hat\tau| <$ threshold set to independence", color=INK, fontsize=9)
    tag(ax[0], "a"); tag(ax[1], "b")
    f.tight_layout(w_pad=1.2)
    finish(f, "controls_trunc_threshold.png")

    # 3. tree criteria on data with a non-monotone pair
    rng = np.random.default_rng(11)
    m = 2000
    u1, w = rng.random(m), rng.random(m)
    a_ = fv.Bicop("clayton", 0, [4.0]).hinv1(np.column_stack([u1, w]))
    b_ = fv.Bicop("clayton", 90, [4.0]).hinv1(np.column_stack([u1, w]))
    u2 = np.where(rng.random(m) < 0.5, a_, b_)               # an X-shaped pair: tau near 0, strong dependence
    g = fv.Bicop("gaussian", 0, [0.6])
    u3 = g.hinv1(np.column_stack([u2, rng.random(m)]))
    u4 = g.hinv1(np.column_stack([u3, rng.random(m)]))
    w5 = norm.cdf(0.45 * norm.ppf(u4) + np.sqrt(1 - 0.45 ** 2) * rng.standard_normal(m))   # a weak monotone link
    u5 = fv.Bicop("gumbel", 0, [1.5]).hinv1(np.column_stack([u1, w5]))
    x = np.column_stack([u1, u2, u3, u4, u5])
    tcs = [("tau", r"Kendall $\tau$"), ("rho", r"Spearman $\rho$"), ("joe", "Joe (Gaussian MI)"),
           ("hoeffd", r"Hoeffding $D$"), ("xi", r"Chatterjee $\xi$"), ("mi", "mutual inf. (KSG)"),
           ("ferreira", "Ferreira tail")]
    out = []
    for c, lab in tcs:
        t, vc = timeit(lambda: fv.Vinecop.from_data(x, tree_criterion=c, num_threads=NT), reps=3)
        t1 = {frozenset((int(vc.order[e]), int(vc.order[vc.st[0, e] - 1]))) for e in range(4)}
        out.append((lab, vc.loglik(x), t, frozenset((1, 2)) in t1))
    f, ax = plt.subplots(1, 3, figsize=(7.2, 2.7), gridspec_kw=dict(width_ratios=[0.85, 1.2, 1.0]))
    a = ax[0]
    zz = norm.ppf(x[:800, :2])
    a.scatter(zz[:, 0], zz[:, 1], s=1.8, color=BLUE, alpha=0.5, linewidths=0)
    zaxes(a, 3.2, labels=False)
    a.set_xlabel("$z_1$"); a.set_ylabel("$z_2$")
    a.set_title(rf"pair (1, 2): $\hat\tau = {fv.kendall_tau(x[:, 0], x[:, 1]):.3f}$", color=INK, fontsize=9)
    style(a, grid=False)
    yy = np.arange(len(out))[::-1]
    a = ax[1]
    best = max(o[1] for o in out)
    for y, o in zip(yy, out):
        a.axhline(y, color=GRIDC, lw=0.6, zorder=0)
        a.plot(o[1] - best, y, "o", ms=6, color=GREEN if o[3] else GREY, mec="white", zorder=3)
    a.set_yticks(yy); a.set_yticklabels([o[0] for o in out], fontsize=8)
    a.set_xlabel("log-likelihood $-$ the largest")
    a.set_title("selected vine", color=INK, fontsize=9)
    a.plot([], [], "o", color=GREEN, label="tree 1 links 1 and 2")
    a.plot([], [], "o", color=GREY, label="it does not")
    a.legend(fontsize=7.4, loc="center", handlelength=1.0, frameon=True, facecolor="white", edgecolor="none")
    lo = min(o[1] for o in out) - best
    a.set_xlim(1.12 * lo - 2, -0.08 * lo + 2)
    a = ax[2]
    for y, o in zip(yy, out):
        a.axhline(y, color=GRIDC, lw=0.6, zorder=0)
        a.plot(o[2], y, "o", ms=6, color=BLUE, mec="white", zorder=3)
    a.set_yticks(yy); a.set_yticklabels([])
    a.set_xlabel(f"selection time [s], {NT} threads")
    a.set_title(f"$d = 5$, $n = {m}$", color=INK, fontsize=9)
    a.set_xlim(0, 1.15 * max(o[2] for o in out))
    for a in ax[1:]:
        for sp in ("top", "right", "left"):
            a.spines[sp].set_visible(False)
        a.tick_params(axis="y", length=0)
    for a, t in zip(ax, "abc"):
        tag(a, t)
    f.tight_layout(w_pad=0.8)
    finish(f, "controls_tree_crit.png")


# ------------------------------------------------------------------------------------------ vinecop
def vine_edges(vc):
    """The edges of a Vinecop as dicts for example_figures.draw_tree; tau of each pair copula from
    20 000 of its own draws."""
    out = []
    for t in range(vc.trunc):
        for e in range(vc.d - 1 - t):
            b = vc.get_pair_copula(t, e)
            E.FAMN.setdefault(b.family_code, b.family.capitalize())
            tau = fv.kendall_tau(*b.simulate(20000, seed=1).T) if b.family_code else 0.0
            out.append(dict(t=t + 1, a=int(vc.order[e]), b=int(vc.order[vc.st[t, e] - 1]), fam=b.family_code,
                            rot=b.rotation, par=list(b.parameters), tau=tau,
                            cond=tuple(sorted(int(vc.order[vc.st[k, e] - 1]) for k in range(t)))))
    return out


def draw_tree(a, edges, t, ref=None, xspan=3.0):
    """example_figures.draw_tree with side nodes hung straight below their neighbour and the labels of
    such vertical edges beside them, so that star-shaped trees (C-vines) stay legible."""
    if t == 1:
        nodes = sorted({e["a"] for e in edges if e["t"] == 1} | {e["b"] for e in edges if e["t"] == 1})
        sets = {n_: frozenset([n_]) for n_ in nodes}
        labels = {n_: str(n_) for n_ in nodes}
    else:
        prev = [e for e in edges if e["t"] == t - 1]
        nodes = list(range(len(prev)))
        sets = {i: frozenset([e["a"], e["b"], *e["cond"]]) for i, e in enumerate(prev)}
        labels = {i: E.node_label(e) for i, e in enumerate(prev)}
    links = []
    for e in (e for e in edges if e["t"] == t):
        sa, sb = frozenset([e["a"], *e["cond"]]), frozenset([e["b"], *e["cond"]])
        p = next(n_ for n_ in nodes if sets[n_] == sa); q = next(n_ for n_ in nodes if sets[n_] == sb)
        links.append((p, q, e))
    pos = E.tree_layout(nodes, [(p, q) for p, q, _ in links])
    line = {n_: x for n_, (x, y) in pos.items() if y == 0}
    sc = xspan / max(1, max(line.values()))
    for n_, (x, y) in list(pos.items()):
        if y != 0:                                        # hang below its neighbour on the path
            nb = next(q if p == n_ else p for p, q, _ in links if n_ in (p, q) and (q if p == n_ else p) in line)
            pos[n_] = (line[nb], -1.0)
    pos = {n_: (x * sc if y == 0 else x * sc, y) for n_, (x, y) in pos.items()}
    for p, q, e in links:
        (x1, y1), (x2, y2) = pos[p], pos[q]
        indep = e["fam"] == 0
        new = ref is not None and t == 1 and frozenset([e["a"], e["b"]]) not in ref
        col = GREY if indep else (GOLD if new else INK)
        a.plot([x1, x2], [y1, y2], color=col, lw=0.9 if indep else 1.4, ls=":" if indep else "-", zorder=1)
        xm, ym = 0.5 * (x1 + x2), 0.5 * (y1 + y2)
        name = E.FAMN.get(e["fam"], str(e["fam"])) + (f" {e['rot']}°" if e["rot"] else "")
        tl = rf"$\tau = {e['tau']:.2f}$"
        if y1 == y2:
            a.text(xm, ym + 0.17, name, ha="center", va="bottom", fontsize=7.6, color=INK if not indep else MUTED)
            if not indep:
                a.text(xm, ym - 0.17, tl, ha="center", va="top", fontsize=7.4, color=MUTED)
        else:
            a.text(xm + 0.1, ym + 0.03, name, ha="left", va="bottom", fontsize=7.6, color=INK if not indep else MUTED)
            if not indep:
                a.text(xm + 0.1, ym - 0.03, tl, ha="left", va="top", fontsize=7.4, color=MUTED)
    for n_, (x, y) in pos.items():
        a.text(x, y, labels[n_], ha="center", va="center", fontsize=8 if t == 1 else 7.2, color=INK, zorder=3,
               bbox=dict(boxstyle="round,pad=0.25" if t > 1 else "circle,pad=0.25", fc="white", ec=BLUE, lw=0.9))
    a.set_xlim(-0.55, xspan + 0.55)
    a.set_ylim(-1.2 if any(y != 0 for _, y in pos.values()) else -0.75, 0.75)
    a.axis("off")


def vinecop():
    E.FAMN.update({12: "Plackett", 13: "FGM", 14: "AMH", 17: "BB2", 18: "BB3", 19: "BB5", 20: "BB9", 21: "BB10",
                   22: "BB4"})
    cv = cvine_example()
    u = cv.simulate(2000, seed=11)
    fixed = fv.Vinecop().set_cvine([4, 3, 2, 1]).select(u, num_threads=NT)
    free = fv.Vinecop.from_data(u, num_threads=NT)
    # 1. the trees of the hand-made C-vine and of the free (Dissmann) selection on its sample
    te, fe = vine_edges(cv), vine_edges(free)
    ref = {frozenset([e["a"], e["b"]]) for e in te if e["t"] == 1}
    f, ax = plt.subplots(3, 2, figsize=(7.2, 4.0), gridspec_kw=dict(hspace=0.3, wspace=0.08,
                                                                     height_ratios=[1.45, 1, 1]))
    for c, (edges, title) in enumerate(((te, rf"the C-vine of the example: $\ell = {cv.loglik(u):.1f}$"),
                                        (fe, rf"free selection on its sample: $\ell = {free.loglik(u):.1f}$"))):
        for t in range(1, 4):
            draw_tree(ax[t - 1][c], edges, t, ref=ref if c == 1 else None, xspan=3.0)
            ax[t - 1][c].text(-0.55, 0.62, f"$T_{t}$", fontsize=9, color=MUTED, ha="left", va="center")
        ax[0][c].set_title(title, fontsize=9.5, color=INK, pad=6)
    h = [plt.Line2D([], [], color=INK, lw=1.4), plt.Line2D([], [], color=GREY, lw=0.9, ls=":")]
    lab = ["pair copula", "independence"]
    if any(frozenset([e["a"], e["b"]]) not in ref for e in fe if e["t"] == 1):
        h.insert(1, plt.Line2D([], [], color=GOLD, lw=1.4)); lab.insert(1, "tree-1 edge not in the C-vine")
    f.legend(h, lab, loc="lower center", ncol=len(h), fontsize=8.5, bbox_to_anchor=(0.5, -0.03))
    f.subplots_adjust(left=0.01, right=0.99, top=0.92, bottom=0.08)
    finish(f, "vinecop_trees.png")

    # 2. Rosenblatt transform as a goodness-of-fit check
    w_fix, w_tr = fixed.rosenblatt(u), fv.Vinecop.from_data(u, trunc_lvl=1, num_threads=NT).rosenblatt(u)
    f, ax = plt.subplots(1, 3, figsize=(7.2, 2.65), gridspec_kw=dict(width_ratios=[1, 1, 1.25]))
    for a, x, ttl in ((ax[0], u, r"sample: $(u_2, u_3)$"), (ax[1], w_fix, r"Rosenblatt: $(w_2, w_3)$")):
        zz = norm.ppf(x[:, 1:3])
        a.scatter(zz[:, 0], zz[:, 1], s=1.6, color=BLUE, alpha=0.45, linewidths=0)
        zaxes(a, 3.5, labels=False)
        a.set_title(ttl + rf", $\hat\tau = {fv.kendall_tau(x[:, 1], x[:, 2]):+.3f}$", color=INK, fontsize=9)
        a.set_xlabel("normal score"); style(a, grid=False)
    a = ax[2]
    pairs = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    n = len(u)
    sd = np.sqrt(2 * (2 * n + 5) / (9 * n * (n - 1)))      # sd of Kendall's tau under independence
    a.axhspan(-2 * sd, 2 * sd, color="#E3EDF8", lw=0, zorder=0, label=r"$\pm 2$ sd under independence")
    xs = np.arange(len(pairs))
    a.plot(xs - 0.12, [fv.kendall_tau(w_fix[:, i], w_fix[:, j]) for i, j in pairs], "o", ms=5, color=BLUE,
           mec="white", label="C-vine, families selected")
    a.plot(xs + 0.12, [fv.kendall_tau(w_tr[:, i], w_tr[:, j]) for i, j in pairs], "D", ms=4.5, mfc="white", mec=GOLD,
           mew=1.1, label="vine truncated after tree 1")
    a.axhline(0, color=GREY, lw=0.6)
    a.set_xticks(xs); a.set_xticklabels([f"({i + 1},{j + 1})" for i, j in pairs], fontsize=8)
    a.set_ylabel(r"Kendall $\hat\tau$ of $(w_i, w_j)$")
    a.set_title("dependence left after the transform", color=INK, fontsize=9)
    a.legend(fontsize=7.2, loc="lower left", labelspacing=0.25, handlelength=1.2)
    style(a)
    for a, t in zip(ax, "abc"):
        tag(a, t)
    f.tight_layout(w_pad=1.0)
    finish(f, "vinecop_rosenblatt.png")

    # 3. pdf_grid: predictive densities, and what the shared edges save
    g = np.linspace(0.0025, 0.9975, 200)
    rows = [np.argmin(np.abs(u[:, 0] - q)) for q in (0.1, 0.5, 0.9)]
    ug = np.repeat(u[rows][:, None, :], len(g), axis=1)
    ug[:, :, 1] = g
    dens = fixed.pdf_grid(ug, target=2)
    f, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
    a = ax[0]
    for k, (r, col) in enumerate(zip(rows, (BLUE, GREEN, GOLD))):
        cd = dens[k] / np.trapezoid(dens[k], g)
        a.plot(g, cd, color=col, lw=1.4, label=rf"$u_1 = {u[r, 0]:.2f}$, $u_3 = {u[r, 2]:.2f}$, $u_4 = {u[r, 3]:.2f}$")
        a.axvline(u[r, 1], color=col, lw=0.8, ls=":")
    a.set_xlim(0, 1); a.set_ylim(0, None)
    a.set_xlabel("$u_2$"); a.set_ylabel(r"$c(u_2 \mid u_1, u_3, u_4)$")
    a.set_title("pdf_grid, target 2, normalised; dotted: observed $u_2$", color=INK, fontsize=9)
    a.legend(fontsize=7.2, loc="upper center", labelspacing=0.25, handlelength=1.4, frameon=True, facecolor="white",
             edgecolor="none", framealpha=0.95)
    style(a)
    a = ax[1]
    truth, x = decaying_dvine(10, 500, seed=3)
    bins = [5, 10, 20, 50, 100, 200, 400]
    for tgt, col, mk in ((10, BLUE, "o"), (5, GREEN, "s")):
        tg, tp = [], []
        for b in bins:
            gg = np.repeat(x[:, None, :], b, axis=1)
            gg[:, :, tgt - 1] = np.linspace(0.01, 0.99, b)
            flat = gg.reshape(-1, 10)
            tg.append(timeit(lambda: truth.pdf_grid(gg, tgt, num_threads=4), reps=3)[0])
            tp.append(timeit(lambda: truth.pdf(flat, num_threads=4), reps=3)[0])
        a.plot(bins, np.array(tp) / np.array(tg), "-" + mk, color=col, ms=4, mec="white",
               label=f"target {tgt}" + (" (end of the D-vine path)" if tgt == 10 else " (middle of the path)"))
    a.axhline(1, color=GREY, lw=0.7, ls=":")
    a.set_xscale("log")
    a.set_xlabel("bins per sample"); a.set_ylabel(r"time of pdf / time of pdf_grid")
    a.set_title("10-dim. D-vine, 500 samples, 4 threads", color=INK, fontsize=9)
    a.legend(fontsize=7.4, loc="lower right")
    a.set_ylim(0, None)
    style(a)
    tag(ax[0], "a"); tag(ax[1], "b")
    f.tight_layout(w_pad=1.6)
    finish(f, "vinecop_pdf_grid.png")

    # 4. relax: the non-simplified example of the Python page
    rng = np.random.default_rng(1)
    n = 2000
    u2, w1, w3 = rng.random(n), rng.random(n), rng.random(n)
    x1 = fv.Bicop("clayton", parameters=[2.0]).hinv1(np.column_stack([u2, w1]))
    rho = -0.6 + 1.4 * u2
    v3 = norm.cdf(rho * norm.ppf(w1) + np.sqrt(1 - rho ** 2) * norm.ppf(w3))
    x3 = fv.Bicop("frank", parameters=[5.0]).hinv1(np.column_stack([u2, v3]))
    x = np.column_stack([x1, u2, x3])
    vc = fv.Vinecop().set_dvine([1, 2, 3]).select(x, family_set=["indep", "gaussian", "clayton", "frank"],
                                                  num_threads=NT)
    simp = fv.Vinecop().set_dvine([1, 2, 3]).select(x, family_set=["indep", "gaussian", "clayton", "frank"],
                                                    num_threads=NT)
    vc.relax(x, target=0, mode="test", num_threads=NT)
    (t_rel, e_rel), = vc.relaxed_edges
    fam = vc.get_pair_copula(t_rel, e_rel).family
    wk, pk = vc.wk[0], vc.pk[0, :, 0]
    # pseudo-true parameter: the family fitted to 20 000 draws of the true conditional Gaussian copula
    zg = np.linspace(-1.9, 1.9, 13)
    pt = [fv.Bicop(fam).fit(fv.Bicop("gaussian", 0, [-0.6 + 1.4 * norm.cdf(z)]).simulate(20000, seed=2),
                            num_threads=NT).parameters[0] for z in zg]
    f, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
    a = ax[0]
    ww = np.linspace(-2.6, 2.6, 300)
    a.plot(ww, np.interp(ww, wk, pk), color=BLUE, lw=1.6, label="relaxed edge (linear between anchors)")
    a.plot(wk, pk, "o", ms=5, color=BLUE, mec="white", zorder=4)
    a.axhline(simp.par[t_rel, e_rel, 0], color=GREY, lw=1.1, ls="--",
              label=f"simplified vine: one {fam} parameter")
    a.plot(zg, pt, "D", ms=4, mfc="white", mec=GOLD, mew=1.1, label=f"{fam} closest to the true copula")
    a.set_xlabel(r"$w = \Phi^{-1}(u_2)$ (loading 1 on the single conditioning variable)")
    a.set_ylabel(f"{fam} parameter of edge (1, 3 | 2)")
    a.set_title(rf"$\ell$: simplified {simp.loglik(x):.1f}, non-simplified {vc.loglik(x):.1f}", color=INK, fontsize=9)
    a.legend(fontsize=7.2, loc="upper left", labelspacing=0.25)
    style(a)
    a = ax[1]
    gg = np.linspace(0.0025, 0.9975, 200)
    for q, col in zip((0.1, 0.5, 0.9), (BLUE, GREEN, GOLD)):
        base = np.array([0.5, q, 0.5])
        ug = np.repeat(base[None, None, :], len(gg), axis=1)
        ug[0, :, 2] = gg
        for m, ls in ((vc, "-"), (simp, "--")):
            dd = m.pdf_grid(ug, target=3)[0]
            a.plot(gg, dd / np.trapezoid(dd, gg), color=col, lw=1.3 if ls == "-" else 0.9, ls=ls)
        a.plot([], [], color=col, lw=1.3, label=f"$u_2 = {q}$")
    a.plot([], [], color=INK, lw=1.3, label="non-simplified")
    a.plot([], [], color=INK, lw=0.9, ls="--", label="simplified")
    a.set_xlim(0, 1); a.set_ylim(0, None)
    a.set_xlabel("$u_3$"); a.set_ylabel(r"$c(u_3 \mid u_1 = 0.5, u_2)$")
    a.set_title("predictive density of variable 3 (pdf_grid)", color=INK, fontsize=9)
    a.legend(fontsize=7.2, loc="upper center", ncol=2, labelspacing=0.25, columnspacing=1.0)
    style(a)
    tag(ax[0], "a"); tag(ax[1], "b")
    f.tight_layout(w_pad=1.6)
    finish(f, "vinecop_relax.png")


# ------------------------------------------------------------------------------------------ discrete
def ecdf_pair(k):
    """Empirical F(x) and F(x-) of integer data with the denominator n + 1 (as on the page)."""
    s = np.sort(k)
    n = len(k)
    return np.searchsorted(s, k, side="right") / (n + 1.0), np.searchsorted(s, k, side="left") / (n + 1.0)


def discrete():
    # 1. the input of a discrete column, and what ignoring the ties does to the estimate
    levels = [2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50]
    R, n = 30, 1000
    est = {"disc": [], "mid": []}
    for K in levels:
        dd, mm = [], []
        for r in range(R):
            v = fv.Bicop("gaussian", 0, [0.7]).simulate(n, seed=100 * K + r)
            k = np.minimum(np.floor(K * v[:, 1]), K - 1)
            up, um = ecdf_pair(k)
            vc = fv.Vinecop.from_data(np.column_stack([v[:, 0], up, um]), var_types=["c", "d"],
                                      family_set=["gaussian"], num_threads=1)
            dd.append(vc.par[0, 0, 0])
            mm.append(fv.Bicop("gaussian").fit(np.column_stack([v[:, 0], 0.5 * (up + um)]), num_threads=1).parameters[0])
        est["disc"].append(dd); est["mid"].append(mm)
    f, ax = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw=dict(width_ratios=[1, 1.3]))
    a = ax[0]
    v = fv.Bicop("gaussian", 0, [0.7]).simulate(n, seed=9)
    k = np.minimum(np.floor(4 * v[:, 1]), 3)
    up, um = ecdf_pair(k)
    m = 120
    a.hlines(v[:m, 0], um[:m], up[:m], color=BLUE, lw=0.7, alpha=0.75)
    a.plot(0.5 * (um[:m] + up[:m]), v[:m, 0], "o", ms=1.8, color=GOLD, mew=0)
    for b in np.unique(up):
        a.axvline(b, color=GRIDC, lw=0.6, zorder=0)
    a.plot([], [], color=BLUE, lw=1.0, label=r"$[u^-, u] = [F(x^-), F(x)]$")
    a.plot([], [], "o", ms=3, color=GOLD, label=r"mid-rank $(u + u^-)/2$")
    a.set_xlim(0, 1); a.set_ylim(0, 1)
    a.set_xlabel(r"$u_2$ of a count with 4 levels"); a.set_ylabel(r"$u_1$ (continuous)")
    a.set_title(f"input of a discrete column, first {m} rows", color=INK, fontsize=9)
    a.legend(fontsize=7.2, loc="upper left", frameon=True, facecolor="white", edgecolor="none", framealpha=0.95,
             handlelength=1.2)
    style(a, grid=False)
    a = ax[1]
    for key, col, mk, lab in (("disc", BLUE, "o", "discrete likelihood (var_types d)"),
                              ("mid", GOLD, "s", "mid-ranks treated as continuous")):
        e = np.array(est[key])
        a.fill_between(levels, np.percentile(e, 10, axis=1), np.percentile(e, 90, axis=1), color=col, alpha=0.15, lw=0)
        a.plot(levels, e.mean(axis=1), "-" + mk, color=col, ms=4, mec="white", label=lab)
    a.axhline(0.7, color=GREY, lw=0.8, ls="--", label=r"latent $\rho = 0.7$")
    a.set_xscale("log"); a.set_xticks(levels); a.set_xticklabels([str(x) for x in levels], fontsize=8)
    a.minorticks_off()
    a.set_xlabel("number of distinct values of the count"); a.set_ylabel(r"estimated Gaussian $\hat\rho$")
    a.set_title(rf"$n = {n}$, {R} samples per level (band: 10-90 %)", color=INK, fontsize=9)
    a.legend(fontsize=7.4, loc="lower right")
    style(a)
    tag(ax[0], "a"); tag(ax[1], "b")
    f.tight_layout(w_pad=1.6)
    finish(f, "discrete_ties.png")

    # 2. the page's example: randomised Rosenblatt residuals
    gen = fv.Vinecop().set_dvine([1, 2, 3])
    for e in range(2):
        gen.set_pair_copula(0, e, fv.Bicop("gaussian", 0, [0.7]))
    v = gen.simulate(2000, seed=2)
    k = np.minimum((4.0 * v[:, 2] ** 1.5).astype(int), 3)
    up, um = ecdf_pair(k)
    u = np.column_stack([v[:, 0], v[:, 1], up, um])
    vc = fv.Vinecop.from_data(u, var_types=["c", "c", "d"], num_threads=4)
    r = vc.rosenblatt(u, seed=3)
    rn = vc.rosenblatt(u, randomize_discrete=False)
    f, ax = plt.subplots(1, 3, figsize=(7.2, 2.6))
    a = ax[0]
    edges = np.linspace(0, 1, 21)
    a.hist(r[:, 2], bins=edges, density=True, color=BLUE, alpha=0.85, edgecolor="white", lw=0.4, label="randomised")
    a.hist(rn[:, 2], bins=edges, density=True, histtype="step", color=GOLD, lw=1.3,
           label=r"$F(x_3)$, not randomised")
    a.axhline(1, color=INK, lw=0.8)
    a.set_xlim(0, 1); a.set_xlabel("$w_3$ (the count)"); a.set_ylabel("density")
    a.legend(fontsize=7.0, loc="upper right", handlelength=1.2)
    a.set_ylim(0, a.get_ylim()[1] * 1.35)
    a.set_title("residual of the discrete column", color=INK, fontsize=9)
    style(a, grid=False)
    a = ax[1]
    a.plot([0, 1], [0, 1], color=GREY, lw=0.8, ls=":")
    for j, col in zip(range(3), (BLUE, GREEN, GOLD)):
        s = np.sort(r[:, j])
        a.plot(s, np.arange(1, len(s) + 1) / len(s), color=col, lw=1.2,
               label=f"$w_{j + 1}$: mean {r[:, j].mean():.3f}, var {r[:, j].var():.4f}")
    a.set_xlim(0, 1); a.set_ylim(0, 1); a.set_aspect("equal")
    a.set_xlabel("$w$"); a.set_ylabel("empirical distribution")
    a.legend(fontsize=6.8, loc="upper left", handlelength=1.2, labelspacing=0.2)
    a.set_title("all three columns: uniform", color=INK, fontsize=9)
    style(a)
    a = ax[2]
    zz = norm.ppf(r[:, [0, 2]])
    a.scatter(zz[:, 0], zz[:, 1], s=1.6, color=BLUE, alpha=0.45, linewidths=0)
    zaxes(a, 3.5, labels=False)
    a.set_xlabel(r"$\Phi^{-1}(w_1)$"); a.set_ylabel(r"$\Phi^{-1}(w_3)$")
    a.set_title(rf"independent: $\hat\tau(w_1, w_3) = {fv.kendall_tau(r[:, 0], r[:, 2]):+.3f}$", color=INK,
                fontsize=9)
    style(a, grid=False)
    for a, t in zip(ax, "abc"):
        tag(a, t)
    f.tight_layout(w_pad=0.9)
    finish(f, "discrete_rosenblatt.png")


# ------------------------------------------------------------------------------------------ nonparametric
def nonparametric():
    rng = np.random.default_rng(5)                           # the mixture of the Python page
    a_ = fv.Bicop("clayton", parameters=[4.0]).simulate(3000, seed=1)
    b_ = fv.Bicop("gumbel", 90, [3.0]).simulate(3000, seed=2)
    u = np.where(rng.random((3000, 1)) < 0.5, a_, b_)
    tll = fv.Bicop("tll").fit(u, num_threads=1)
    bern = fv.Bicop("bernstein").fit(u, num_threads=1)
    best = fv.Bicop().select(u, num_threads=NT)
    f, ax = plt.subplots(1, 4, figsize=(7.2, 2.3))
    a = ax[0]
    zz = norm.ppf(u)
    a.scatter(zz[:, 0], zz[:, 1], s=1.0, color=BLUE, alpha=0.4, linewidths=0)
    zaxes(a, 3.0, labels=False)
    a.set_xticks([-2, 0, 2]); a.set_yticks([-2, 0, 2]); a.tick_params(labelsize=8, length=0)
    for sp in a.spines.values():
        sp.set_visible(False)
    a.set_title("the sample, $n = 3000$\nClayton + 90° Gumbel", fontsize=8.6, color=INK)
    for a, c, ttl in ((ax[1], tll, f"TLL\n$\\ell = {tll.loglik():.1f}$, edf {tll.npars:.1f}"),
                      (ax[2], bern, f"Bernstein, $m = {int(bern.aux[2])}$\n$\\ell = {bern.loglik():.1f}$, "
                                    f"edf {bern.npars:.0f}"),
                      (ax[3], best, f"best parametric (BIC): {best.family}\n$\\ell = {best.loglik():.1f}$, "
                                    f"{best.npars} parameters")):
        F.panel(a, zdens(c))
        a.set_title(ttl, fontsize=8.6, color=INK)
    f.tight_layout(pad=0.3, w_pad=0.6)
    finish(f, "nonparametric_mixture.png")

    # 2. effective degrees of freedom and cost against n
    ns = [500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]
    big = fv.Bicop("gaussian", 0, [0.5]).simulate(100000, seed=7)
    res = {k: [] for k in ("tll", "bernstein", "bb1")}
    edf = {k: [] for k in ("tll", "bernstein")}
    for n in ns:
        for k in res:
            t, c = timeit(lambda: fv.Bicop(k).fit(big[:n], num_threads=1), reps=3 if n <= 20000 else 1)
            res[k].append(t)
            if k in edf:
                edf[k].append(c.npars)
    f, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
    a = ax[0]
    for k, col, mk in (("tll", BLUE, "o"), ("bernstein", GREEN, "s")):
        a.plot(ns, edf[k], "-" + mk, color=col, ms=4, mec="white", label=k.upper() if k == "tll" else "Bernstein")
    a.plot(ns, edf["tll"][0] * (np.array(ns) / ns[0]) ** (1 / 3), color=GREY, lw=0.8, ls=":", label=r"$\propto n^{1/3}$")
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel("$n$"); a.set_ylabel("effective degrees of freedom")
    a.set_title(r"Gaussian sample, $\rho = 0.5$", color=INK, fontsize=9)
    a.legend(fontsize=7.6, loc="upper left")
    style(a)
    a = ax[1]
    for k, col, mk, lab in (("tll", BLUE, "o", "TLL"), ("bernstein", GREEN, "s", "Bernstein"),
                            ("bb1", GOLD, "D", "BB1 (parametric, for comparison)")):
        a.plot(ns, res[k], "-" + mk, color=col, ms=4, mec="white", label=lab)
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel("$n$"); a.set_ylabel("time of one fit [s], 1 thread")
    a.set_title("cost of the fit", color=INK, fontsize=9)
    a.legend(fontsize=7.6, loc="upper left")
    style(a)
    tag(ax[0], "a"); tag(ax[1], "b")
    f.tight_layout(w_pad=1.6)
    finish(f, "nonparametric_cost.png")


# ------------------------------------------------------------------------------------------ generic
def generic():
    # 1. the code scheme: which (og, ie) combinations are new, duplicates or excluded
    new = {f.code for f in fv.families() if f.code >= 10000}
    ogs = list(range(0, 8)) + list(range(11, 27))
    f, ax = plt.subplots(2, 1, figsize=(7.2, 4.4), gridspec_kw=dict(height_ratios=[14, 5.6], hspace=0.55))
    stat_col = {2: BLUE, 1: "#AFCBEA", 0: "white"}
    count = {0: 0, 1: 0, 2: 0}
    for a, k, ies in ((ax[0], 1, range(14)), (ax[1], 2, range(5))):
        for i, ie in enumerate(ies):
            for j, og in enumerate(ogs):
                code = 10000 * k + 100 * og + ie
                s = 2 if code in new else (1 if fv._fam_name(code) != "invalid" else 0)
                count[s] += 1
                a.add_patch(plt.Rectangle((j - 0.45, i - 0.45), 0.9, 0.9, fc=stat_col[s],
                                          ec=GRIDC if s == 0 else "white", lw=0.6))
        a.set_xlim(-0.6, len(ogs) - 0.4); a.set_ylim(len(ies) - 0.4, -0.6)
        a.set_xticks(range(len(ogs))); a.set_xticklabels([str(o) for o in ogs], fontsize=7.5)
        a.set_yticks(range(len(ies))); a.set_yticklabels([str(x) for x in ies], fontsize=7.5)
        a.set_aspect("equal")
        a.tick_params(length=0)
        for sp in a.spines.values():
            sp.set_visible(False)
        a.set_ylabel("inner function ie", fontsize=8.5)
        a.axvline(7.5, color=GREY, lw=0.8)
    ax[0].set_title(r"$k = 1$ (plain): code $= 10000 + 100\,\mathrm{og} + \mathrm{ie}$", color=INK, fontsize=9.5)
    ax[1].set_title(r"$k = 2$ (Khoudraji-wrapped): code $= 20000 + 100\,\mathrm{og} + \mathrm{ie}$", color=INK,
                    fontsize=9.5)
    ax[1].set_xlabel("outer generator og (0-7 single generators | 11-26 compositions po.pin)", fontsize=8.5)
    h = [plt.Rectangle((0, 0), 1, 1, fc=BLUE), plt.Rectangle((0, 0), 1, 1, fc="#AFCBEA"),
         plt.Rectangle((0, 0), 1, 1, fc="white", ec=GREY, lw=0.6)]
    f.legend(h, [f"new family ({count[2]})", f"valid duplicate ({count[1]})", "excluded"], loc="lower center",
             ncol=3, fontsize=8.5, bbox_to_anchor=(0.5, -0.02))
    f.subplots_adjust(left=0.08, right=0.99, top=0.95, bottom=0.12)
    finish(f, "generic_codes.png")

    # 2. the 13 inner Pickands functions, A(w) = -log C(exp(-(1 - w)), exp(-w)) of the extreme-value codes 100xx
    w = np.linspace(0.0, 1.0, 201)[1:-1]
    uv = np.column_stack([np.exp(-(1 - w)), np.exp(-w)])
    groups = [("symmetric (ie 1-4)", [(1, [2.0], "logistic, $\\theta = 2$"), (2, [1.0], "negative logistic, $\\theta = 1$"),
                                      (3, [1.5], "Hüsler-Reiss, $\\lambda = 1.5$"), (4, [0.8], "symmetric mixed, $\\theta = 0.8$")]),
              ("asymmetric logistic (ie 5-10)", [(5, [0.4, 0.9, 3.0], "5: $(0.4, 0.9, 3)$"), (6, [3.0, 0.4], "6: $(3, 0.4)$"),
                                                 (7, [3.0, 0.4], "7: $(3, 0.4)$"), (8, [0.4, 0.9, 2.0], "8: $(0.4, 0.9, 2)$"),
                                                 (9, [2.0, 0.4], "9: $(2, 0.4)$"), (10, [2.0, 0.4], "10: $(2, 0.4)$")]),
              ("mixed and Khoudraji (ie 11-13)", [(11, [0.9, 0.9], "11 asym. mixed: $(0.9, 0.9)$"),
                                                  (12, [2.0, 0.5, 0.95], "12 Khoudraji-HR: $(2, 0.5, 0.95)$"),
                                                  (13, [1.0, 0.5, 0.95], "13 Khoudraji-mixed: $(1, 0.5, 0.95)$")])]
    cols = [BLUE, GREEN, GOLD, INK, RED, LBLUE]
    f, ax = plt.subplots(1, 3, figsize=(7.2, 3.0), sharey=True)
    for a, (ttl, items) in zip(ax, groups):
        a.fill_between(w, np.maximum(w, 1 - w), 1, color="#F1F3F8", lw=0, zorder=0)
        a.plot(w, np.maximum(w, 1 - w), color=GREY, lw=0.6)
        for (ie, p, lab), col in zip(items, cols):
            A = -np.log(fv.Bicop(10000 + ie, 0, p).cdf(uv))
            a.plot(w, A, color=col, lw=1.3, ls="--" if ie in (7, 10) else "-", label=lab)
        a.set_xlim(0, 1); a.set_ylim(0.48, 1.02); a.set_xlabel("$w$")
        a.set_title(ttl, color=INK, fontsize=9)
        a.legend(fontsize=6.6, loc="lower center", handlelength=1.6, labelspacing=0.2)
        style(a)
    ax[0].set_ylabel(r"$A(w) = -\log C(e^{-(1-w)}, e^{-w})$")
    for a, t in zip(ax, "abc"):
        tag(a, t)
    f.tight_layout(w_pad=0.6)
    finish(f, "generic_pickands_all.png")

    # 3. tail dependence of the new Archimedean compositions, against the gen_lt / gen_ut flags
    specs = [(11300, [1.5, 1.5]), (11700, [1.0, 1.5]), (11800, [1.0, 5.0]), (12300, [3.0, 1.5]),
             (12400, [3.0, 1.0]), (12500, [3.0, 1.5])]
    q = np.logspace(-5, np.log10(0.5), 80)
    f, ax = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
    for (code, p), col in zip(specs, cols):
        c = fv.Bicop(code, 0, p)
        lam_l = c.cdf(np.column_stack([q, q])) / q
        lam_u = (2 * q - 1 + c.cdf(np.column_stack([1 - q, 1 - q]))) / q
        tau = fv.kendall_tau(*c.simulate(20000, seed=1).T)
        lab = f"{c.family} {code}, ({p[0]:g}, {p[1]:g}), " + rf"$\tau \approx {tau:.2f}$"
        ax[0].plot(q, lam_l, color=col, lw=1.3, label=lab)
        ax[1].plot(q, lam_u, color=col, lw=1.3)
    ax[0].set_title(r"lower tail $\lambda_L(q) = C(q, q)/q$", color=INK, fontsize=9)
    ax[1].set_title(r"upper tail $\lambda_U(q) = \{2q - 1 + C(1-q, 1-q)\}/q$", color=INK, fontsize=9)
    for a in ax:
        a.set_xscale("log"); a.set_xlim(1e-5, 0.5); a.set_ylim(-0.02, 1.0)
        a.set_xlabel("$q$ (distance from the corner)")
        style(a)
    ax[0].set_ylabel(r"$\lambda(q)$")
    h, l = ax[0].get_legend_handles_labels()
    f.legend(h, l, loc="lower center", ncol=2, fontsize=7.6, bbox_to_anchor=(0.5, -0.13), handlelength=1.8)
    tag(ax[0], "a"); tag(ax[1], "b")
    f.tight_layout(w_pad=1.0)
    finish(f, "generic_tails.png")


# ------------------------------------------------------------------------------------------ performance
def performance():
    # 1. thread scaling of selection and evaluation (the page's example), at most 16 threads
    d, n = 8, 10000
    truth = fv.Vinecop().set_dvine(list(range(1, d + 1)))
    for e in range(d - 1):
        truth.set_pair_copula(0, e, fv.Bicop("gumbel", 0, [1.8]))
    for e in range(d - 2):
        truth.set_pair_copula(1, e, fv.Bicop("clayton", 0, [0.8]))
    u = truth.simulate(n, seed=9)
    big = truth.simulate(200000, seed=10)
    threads = [1, 2, 4, 8, MAXT]
    sel, lls, ev = [], [], {"pdf": [], "rosenblatt": [], "inverse_rosenblatt": []}
    vc = None
    for k in threads:
        t, vc = timeit(lambda: fv.Vinecop.from_data(u, num_threads=k), reps=1)
        sel.append(t); lls.append(vc.loglik(u))
        print(f"  select, {k:2d} threads: {t:7.2f} s  loglik {lls[-1]:.6f}")
        ev["pdf"].append(timeit(lambda: vc.pdf(big, num_threads=k), reps=3)[0])
        ev["rosenblatt"].append(timeit(lambda: vc.rosenblatt(big, num_threads=k), reps=3)[0])
        ev["inverse_rosenblatt"].append(timeit(lambda: vc.inverse_rosenblatt(big, num_threads=k), reps=3)[0])
    same = len(set(lls)) == 1
    f, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
    series = [("vine selection, $n = 10\\,000$", sel, BLUE, "o"),
              ("pdf, $n = 200\\,000$", ev["pdf"], GREEN, "s"),
              ("rosenblatt", ev["rosenblatt"], GOLD, "D"),
              ("inverse_rosenblatt", ev["inverse_rosenblatt"], INK, "^")]
    for lab, t, col, mk in series:
        ax[0].plot(threads, t, "-" + mk, color=col, ms=4.5, mec="white", label=lab)
        ax[1].plot(threads, t[0] / np.array(t), "-" + mk, color=col, ms=4.5, mec="white", label=lab)
    ax[1].plot(threads, threads, color=GREY, lw=0.8, ls=":", label="ideal")
    for a in ax:
        a.set_xscale("log", base=2); a.set_xticks(threads); a.set_xticklabels([str(k) for k in threads])
        a.minorticks_off()
        a.set_xlabel("threads (nthreads)")
        style(a)
    ax[0].set_yscale("log"); ax[0].set_ylabel("wall-clock time [s]")
    ax[1].set_yscale("log", base=2); ax[1].set_yticks(threads); ax[1].set_yticklabels([str(k) for k in threads])
    ax[1].set_ylabel("speed-up over 1 thread")
    ax[0].set_title(f"8-dim. vine; log-likelihood {lls[0]:.6f}" + (" for every count" if same else ""),
                    color=INK, fontsize=8.8)
    ax[1].set_title("speed-up", color=INK, fontsize=9)
    ax[1].legend(fontsize=7.2, loc="upper left", labelspacing=0.25)
    tag(ax[0], "a"); tag(ax[1], "b")
    f.tight_layout(w_pad=1.6)
    finish(f, "performance_threads.png")

    # 2. what costs time: family set against n, and truncation in 20 dimensions
    t5 = fv.Vinecop().set_dvine([1, 2, 3, 4, 5])
    for e in range(4):
        t5.set_pair_copula(0, e, fv.Bicop("bb1", 0, [0.6, 1.4]))
    for e in range(3):
        t5.set_pair_copula(1, e, fv.Bicop("student", 0, [0.3, 5.0]))
    x = t5.simulate(20000, seed=5)
    ns = [500, 1000, 2000, 5000, 10000, 20000]
    sets = [("one_par", "one-parameter (5 + indep.)", GREY, "v"), ("parametric", "parametric, codes 0-11", GREEN, "s"),
            ("default", "default, codes 0-22", BLUE, "o"), ("all+generic", "codes 0-23 + 85 generic", GOLD, "D")]
    tt = {}
    for key, *_ in sets:
        fs = {"one_par": ["indep"] + list(fv.family_group("one_par")), "parametric": "parametric",
              "default": None, "all+generic": ["all_parametric", "generic"]}[key]
        tt[key] = [timeit(lambda: fv.Vinecop.from_data(x[:m], family_set=fs, num_threads=NT), reps=1)[0] for m in ns]
        print(f"  {key}: " + " ".join(f"{t:.2f}" for t in tt[key]))
    d20, x20 = decaying_dvine(20, 2000, seed=4)
    ks = [1, 2, 3, 4, 6, 8, 12, 19]
    tk, llk = [], []
    for kk in ks:
        t, v = timeit(lambda: fv.Vinecop.from_data(x20, trunc_lvl=kk, num_threads=NT), reps=1)
        tk.append(t); llk.append(v.loglik(x20))
    auto = fv.Vinecop.from_data(x20, select_trunc_lvl=True, num_threads=NT)
    f, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
    a = ax[0]
    a.axvspan(1800, 25000, color="#F1F3F8", lw=0, zorder=0)
    a.text(2000, 0.84, "generic pre-screen\n($n > 1800$)", fontsize=7.4, color=MUTED, va="top",
           transform=a.get_xaxis_transform())
    for key, lab, col, mk in sets:
        a.plot(ns, tt[key], "-" + mk, color=col, ms=4, mec="white", label=lab)
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(400, 25000)
    a.set_ylim(None, 3000)
    a.set_xlabel("$n$"); a.set_ylabel(f"selection time [s], {NT} threads")
    a.set_title("family set, $d = 5$", color=INK, fontsize=9)
    a.legend(fontsize=7.0, loc="upper left", labelspacing=0.25, ncol=2, columnspacing=0.8, handlelength=1.6)
    style(a)
    a = ax[1]
    edges = [sum(20 - t for t in range(1, k + 1)) for k in ks]
    a.plot(ks, tk, "-o", color=BLUE, ms=4, mec="white", label="selection time")
    for kk, t, e in zip(ks, tk, edges):
        a.annotate(f"{e}", (kk, t), xytext=(0, 6), textcoords="offset points", ha="center", fontsize=7, color=MUTED)
    a.axvline(auto.trunc, color=GOLD, lw=1.1, ls="--", label=f"select_trunc_lvl: {auto.trunc}")
    a.set_xticks([1, 4, 8, 12, 16, 19])
    a.set_xlabel("trunc_lvl"); a.set_ylabel(f"selection time [s], {NT} threads")
    a.set_ylim(0, max(tk) * 1.18)
    a.set_title(r"truncation, $d = 20$, $n = 2000$ (grey: edges fitted)", color=INK, fontsize=9)
    a.legend(fontsize=7.4, loc="lower right")
    style(a)
    tag(ax[0], "a"); tag(ax[1], "b")
    f.tight_layout(w_pad=1.6)
    finish(f, "performance_cost.png")


JOBS = {"bicop": bicop, "rotations": rotations, "fitting": fitting, "controls": controls, "vinecop": vinecop,
        "discrete": discrete, "nonparametric": nonparametric, "generic": generic, "performance": performance}

if __name__ == "__main__":
    for job in (sys.argv[1:] or list(JOBS)):
        JOBS[job]()
