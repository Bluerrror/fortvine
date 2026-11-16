"""python_quickstart.py -- the quick-start program of the docs with the fortvine Python package.

    pip install ./python matplotlib
    python example/python/python_quickstart.py [figure.png]

Prints the worked example of docs/python.md and, given a path, draws its figure.
"""
import sys

import numpy as np

import fortvine as fv

n = 1000

# ---- one pair copula: Gumbel with theta = 1.5 (Kendall tau = 1 - 1/theta = 1/3)
cop = fv.Bicop("gumbel", parameters=[1.5])
g = np.array([0.1, 0.5, 0.95])
print("density c(u,u)   :", cop.pdf(np.column_stack([g, g])).round(4))
print("h1(u | u)        :", cop.hfunc1(np.column_stack([g, g])).round(4))

u12 = cop.simulate(n, seed=42)                       # u2 = h1^-1(w | u1)
print("sample Kendall tau: %.4f" % fv.kendall_tau(u12[:, 0], u12[:, 1]))

est = fv.Bicop("gumbel").fit(u12)
print("Gumbel MLE theta : %.4f   loglik %.3f" % (est.parameters[0], est.loglik()))
best = fv.Bicop().select(u12)
print("BIC selection    :", best)

# ---- a 4-dimensional vine: a C-vine with Clayton pairs as the truth
truth = fv.Vinecop().set_cvine([1, 2, 3, 4])
truth.fam[0, :3] = 3; truth.par[0, :3, 0] = 3.0      # tree 1: Clayton(3)
truth.fam[1, :2] = 3; truth.par[1, :2, 0] = 1.0      # tree 2: Clayton(1)
u = truth.simulate(n, seed=1)

vc = fv.Vinecop.from_data(u)                         # structure, families, parameters
print(vc)
print("loglik %.2f   (truth %.2f)   parameters %d" % (vc.loglik(u), truth.loglik(u), vc.npars))
print("five draws from the fitted vine:")
print(vc.simulate(5, seed=2).round(4))

if len(sys.argv) > 1:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman"], "mathtext.fontset": "stix",
                         "font.size": 10, "axes.facecolor": "white", "figure.facecolor": "white"})
    ramp = LinearSegmentedColormap.from_list(
        "fv", ["#FFFFFF", "#E3EDF8", "#AFCBEA", "#5F95D0", "#1F5FA6", "#0D2C52"])
    blue, green, gold = "#1F5FA6", "#2E8B57", "#C7922B"

    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.6), dpi=140, constrained_layout=True)
    m = 200
    t = (np.arange(m) + 0.5) / m
    grid = np.column_stack([np.repeat(t, m), np.tile(t, m)])
    dens = best.pdf(grid).reshape(m, m).T              # [u2, u1]
    im = ax[0].imshow(dens, origin="lower", extent=(0, 1, 0, 1), cmap=ramp, vmin=0, vmax=3, aspect="equal")
    fig.colorbar(im, ax=ax[0], shrink=0.8, label=r"$c(u_1, u_2)$")
    ax[0].scatter(u12[:200, 0], u12[:200, 1], s=3, color=gold, lw=0)
    ax[0].set(title=r"(a) %s, $\theta = %.3f$" % (best.family.capitalize(), best.parameters[0]),
              xlabel=r"$u_1$", ylabel=r"$u_2$")

    sim = vc.simulate(n, seed=3)
    ax[1].scatter(u[:, 0], u[:, 3], s=2, color=blue, lw=0, label="data")
    ax[1].scatter(sim[:, 0], sim[:, 3], s=2, color=green, lw=0, label="fitted vine")
    ax[1].set(title="(b) variables 1 and 4", xlabel=r"$u_1$", ylabel=r"$u_4$", aspect="equal",
              xlim=(0, 1), ylim=(0, 1))
    ax[1].legend(loc="lower right", fontsize=8, markerscale=4, frameon=False)

    w = vc.rosenblatt(u)
    ax[2].hist(w.ravel(), bins=20, range=(0, 1), density=True, color="#AFCBEA", edgecolor=blue, lw=0.6)
    ax[2].axhline(1.0, color=gold, lw=1.2)
    ax[2].set(title="(c) Rosenblatt transform", xlabel=r"$w$", ylabel="density", xlim=(0, 1), ylim=(0, 1.5))
    for a in ax:
        a.tick_params(direction="out", length=3)
    fig.savefig(sys.argv[1], dpi=140)
    print("figure written to", sys.argv[1])
