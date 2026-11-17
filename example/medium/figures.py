"""Figures of the Medium post, from bc_*.json/npz and ca_*.json/npz (run the two experiments first)."""
import json, os, warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from sklearn.datasets import load_breast_cancer
from vinepred import VineClassifier

warnings.filterwarnings("ignore")
os.makedirs("img", exist_ok=True)
plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"],
                     "mathtext.fontset": "stix", "font.size": 12, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": "#8a8a85", "axes.grid": True,
                     "grid.color": "#e6e5e0", "grid.linewidth": 0.8, "axes.axisbelow": True,
                     "savefig.dpi": 200, "savefig.bbox": "tight", "legend.frameon": False})
INK, MUTED = "#1f1f1e", "#6b6a66"
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
STY = {"Vine copula": (BLUE, "o", "-"), "Naive Bayes": (ORANGE, "s", "--"),
       "Gradient boosting": (AQUA, "^", "-."), "Quantile boosting": (AQUA, "^", "-."),
       "Logistic regression": (YELLOW, "v", ":"), "Linear (Gaussian)": (YELLOW, "v", ":")}
BOX = dict(boxstyle="round,pad=0.45", fc="white", ec="#d6d5d0", lw=0.8, alpha=0.95)
DIV = LinearSegmentedColormap.from_list("bo", [BLUE, "#9ec5f4", "#f2f1ee", "#f5b393", ORANGE])

# ---------------------------------------------------------------- breast cancer
bc = json.load(open("bc_results.json"))
bcp = np.load("bc_pred.npz")
yb = bcp["y"]
BCM = ["Vine copula", "Naive Bayes", "Logistic regression", "Gradient boosting"]

# 1. the probability map of a two-feature vine
d = load_breast_cancer()
j1, j2 = list(d.feature_names).index("worst concave points"), list(d.feature_names).index("worst area")
X2 = np.log1p(d.data[:, [j1, j2]])
y = (d.target == 0).astype(int)
clf = VineClassifier(num_threads=12).fit(X2, y)
g1 = np.linspace(X2[:, 0].min() - 0.01, X2[:, 0].max() + 0.01, 220)
g2 = np.linspace(X2[:, 1].min() - 0.15, X2[:, 1].max() + 0.15, 220)
G1, G2 = np.meshgrid(g1, g2)
Pm = clf.predict_proba(np.column_stack([G1.ravel(), G2.ravel()]))[:, 1].reshape(G1.shape)
fig, ax = plt.subplots(figsize=(7.6, 5.6))
cs = ax.contourf(np.expm1(G1), np.expm1(G2), Pm, levels=np.linspace(0, 1, 21), cmap=DIV)
ax.contour(np.expm1(G1), np.expm1(G2), Pm, levels=[0.1, 0.5, 0.9], colors=[INK], linewidths=[0.7, 1.8, 0.7],
           linestyles=["--", "-", "--"])
for k, mk, lab in [(0, "o", "benign"), (1, "^", "malignant")]:
    s = y == k
    ax.scatter(d.data[s, j1], d.data[s, j2], marker=mk, s=22, c="white" if k == 0 else INK, edgecolors=INK,
               linewidths=0.6, label=f"{lab} (n = {s.sum()})", zorder=3)
ax.set_yscale("log")
ax.set_xlabel("worst concave points"); ax.set_ylabel("worst area (log scale)")
ax.set_title("A vine copula turns two measurements into P(malignant)", loc="left", color=INK)
cb = fig.colorbar(cs, ax=ax, ticks=[0, 0.1, 0.5, 0.9, 1]); cb.set_label("P(malignant | x)")
ax.legend(loc="upper left")
ax.text(1.0, -0.12, "contours: 50 % (solid), 10 % and 90 % (dashed)", transform=ax.transAxes, ha="right",
        color=MUTED, fontsize=10)
ax.grid(False)
fig.savefig("img/bc_map.png"); plt.close(fig)

# 2. reliability diagram with the probability scores on the plot
fig, ax = plt.subplots(figsize=(7.2, 6.0))
ax.plot([0, 1], [0, 1], color="#b9b8b2", lw=1, zorder=1)
edges = np.linspace(0, 1, 11)
for m in BCM:
    p = bcp[m.replace(" ", "_")]
    b = np.minimum((p * 5).astype(int), 4)
    xs = [p[b == k].mean() for k in range(5) if np.sum(b == k) >= 3]
    ys = [yb[b == k].mean() for k in range(5) if np.sum(b == k) >= 3]
    c, mk, ls = STY[m]
    ax.plot(xs, ys, color=c, marker=mk, ls=ls, lw=2, ms=8, mec="white", mew=1.5, label=m)
txt = "\n".join(f"{m:<20s} ECE {bc[m]['ECE']:.3f}  Brier {bc[m]['Brier']:.3f}  LL {bc[m]['Log loss']:.3f}"
                for m in BCM)
ax.text(0.98, 0.03, txt, transform=ax.transAxes, ha="right", va="bottom", family="monospace", fontsize=8.6,
        bbox=BOX, color=INK)
ax.set_xlabel("predicted P(malignant)"); ax.set_ylabel("observed share malignant")
ax.set_title("Calibration: the vine has the smallest calibration error", loc="left", color=INK)
ax.legend(loc="upper left"); ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
fig.savefig("img/bc_reliability.png"); plt.close(fig)

# 3. scorecard: every metric, every model
hi = ["Accuracy", "Balanced accuracy", "Sensitivity", "Specificity", "Precision", "F1", "MCC", "ROC AUC", "PR AUC"]
lo = ["Log loss", "Brier", "ECE"]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5.4), gridspec_kw={"width_ratios": [1.5, 1]})
for i, m in enumerate(BCM):
    c, mk, _ = STY[m]
    off = (i - 1.5) * 0.17
    a1.scatter([bc[m][k] for k in hi], np.arange(len(hi)) + off, color=c, marker=mk, s=60, ec="white", lw=1.2,
               label=m, zorder=3)
    a2.scatter([bc[m][k] for k in lo], np.arange(len(lo)) + off, color=c, marker=mk, s=60, ec="white", lw=1.2,
               zorder=3)
a1.set_yticks(range(len(hi)), hi); a1.invert_yaxis(); a1.set_xlabel("score (higher is better)")
a2.set_yticks(range(len(lo)), lo); a2.invert_yaxis(); a2.set_xlabel("loss (lower is better)")
a1.set_title("Breast cancer, 10-fold cross-validation, 569 patients", loc="left", color=INK)
a1.legend(loc="lower left", fontsize=10)
cm = "\n".join(f"{m:<20s} missed {bc[m]['Missed cancers']:>2d}  false alarms {bc[m]['False alarms']:>2d}"
               for m in BCM)
a2.set_xlim(0, 0.31)
fig.tight_layout(rect=(0, 0.13, 1, 1))
fig.text(0.62, 0.01, cm, ha="left", va="bottom", family="monospace", fontsize=9.5, bbox=BOX, color=INK)
fig.text(0.03, 0.04, "212 malignant and 357 benign tumours; every prediction comes\nfrom a model that never saw that patient.",
         color=MUTED, fontsize=11)
fig.savefig("img/bc_scorecard.png"); plt.close(fig)

# ---------------------------------------------------------------- California
ca = json.load(open("ca_results.json"))
cp = np.load("ca_pred.npz")
obs, LV = cp["obs"], cp["levels"]
CAM = ["Vine copula", "Quantile boosting", "Linear (Gaussian)"]
QK = {m: "Q_" + m.replace(" ", "_").replace("(", "").replace(")", "") for m in CAM}


def metric_box(m):
    r = ca[m]
    return (f"NSE   {r['NSE']:.3f}\nKGE   {r['KGE']:.3f}\nr     {r['r']:.3f}\nRMSE  {r['RMSE']:.3f}\n"
            f"MAE   {r['MAE']:.3f}\nPBIAS {r['PBIAS %']:+.1f} %\nCRPS  {r['CRPS']:.3f}\n"
            f"PICP90 {r['PICP 90 %']:.1f} %\nMPIW90 {r['MPIW 90']:.2f}\nIS90  {r['Interval score']:.2f}")


# 4. observed vs predicted, three models, all metrics on each panel
fig, axs = plt.subplots(1, 3, figsize=(16, 5.9), sharex=True, sharey=True)
for ax, m in zip(axs, CAM):
    Q = cp[QK[m]]
    c, mk, _ = STY[m]
    ax.hexbin(Q[:, 9], obs, gridsize=50, extent=(0, 5.6, 0, 5.2), cmap="Greys", mincnt=1, bins="log",
              linewidths=0)
    rng = np.random.default_rng(0).choice(len(obs), 120, replace=False)
    ax.errorbar(Q[rng, 9], obs[rng], xerr=[Q[rng, 9] - Q[rng, 0], Q[rng, -1] - Q[rng, 9]], fmt=mk, ms=4,
                color=c, ecolor=c, elinewidth=0.8, alpha=0.75, capsize=0)
    ax.plot([0, 5.6], [0, 5.6], color=INK, lw=1)
    ax.text(0.03, 0.97, metric_box(m), transform=ax.transAxes, va="top", family="monospace", fontsize=9,
            bbox=BOX, color=INK)
    ax.set_title(m, loc="left", color=c, fontweight="bold")
    ax.set_xlim(0, 5.6); ax.set_ylim(0, 5.3)
    ax.set_xlabel("predicted median ($100k), bars: 90 % interval")
axs[0].set_ylabel("observed median value ($100k)")
fig.suptitle(f"California housing, {ca['_info']['n_test']} held-out districts", x=0.05, ha="left", y=0.995,
             color=INK, fontsize=14)
fig.tight_layout(); fig.savefig("img/ca_scatter.png"); plt.close(fig)

# 5. predictive densities of six districts
yg = np.exp(cp["yg"])
Qv, Qq = cp[QK["Vine copula"]], cp[QK["Quantile boosting"]]
width = (Qv[:, -1] - Qv[:, 0]) / Qv[:, 9]
order = np.argsort(Qv[:, 9])                       # districts at fixed percentiles of the predicted median
pick = order[(np.array([5, 25, 45, 65, 85, 97]) * len(order) // 100)]
fig, axs = plt.subplots(2, 3, figsize=(12, 6.6))
for ax, i in zip(axs.ravel(), pick):
    fv_ = cp["dens_vine"][i] / yg                  # density in price, not log price
    top = fv_.max()
    ax.fill_between(yg, fv_, where=(yg >= Qv[i, 0]) & (yg <= Qv[i, -1]), color=BLUE, alpha=0.18, lw=0)
    ax.plot(yg, fv_, color=BLUE, lw=2, label="Vine copula: density, 90 % shaded")
    ax.plot([Qq[i, 0], Qq[i, -1]], [1.12 * top] * 2, color=AQUA, lw=3, solid_capstyle="butt",
            label="Quantile boosting: 90 % interval")
    ax.plot(Qq[i, 9], 1.12 * top, "^", color=AQUA, ms=7, mec="white")
    ax.axvline(obs[i], color=INK, lw=1.6)
    ax.text(0.97, 0.80, f"observed {obs[i]:.2f}\nvine median {Qv[i, 9]:.2f}\nvine 90 %: {Qv[i, 0]:.2f}–{Qv[i, -1]:.2f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9.5, color=INK, bbox=BOX)
    if Qv[i, 9] > 3:                               # keep the box off the peak
        ax.texts[-1].set_position((0.03, 0.80)); ax.texts[-1].set_ha("left")
    ax.set_xlim(0, 6); ax.set_ylim(0, 1.25 * top); ax.set_yticks([])
    ax.set_title(f"({cp['lat'][i]:.2f}°N, {-cp['lon'][i]:.2f}°W)", loc="left", fontsize=11, color=MUTED)
for ax in axs[1]:
    ax.set_xlabel("median house value ($100k)")
h, l = axs[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc="upper right", ncol=2, fontsize=10.5, bbox_to_anchor=(0.99, 0.985))
fig.suptitle("Not one number but a distribution (black: the true value)", x=0.02, ha="left", color=INK,
             fontsize=14)
fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig("img/ca_densities.png"); plt.close(fig)

# 6. calibration of every quantile + PIT histograms
fig = plt.figure(figsize=(12, 6.2))
gs = fig.add_gridspec(3, 2, width_ratios=[1.5, 1])
ax = fig.add_subplot(gs[:, 0])
ax.plot([0, 1], [0, 1], color="#b9b8b2", lw=1)
for m in CAM:
    c, mk, ls = STY[m]
    emp = (obs[:, None] <= cp[QK[m]]).mean(0)
    ax.plot(LV, emp, color=c, marker=mk, ls=ls, lw=1.8, ms=6, mec="white", label=f"{m}")
ax.set_xlabel("nominal quantile level"); ax.set_ylabel("share of districts below the quantile")
ax.set_title("Are the quantiles honest?", loc="left", color=INK); ax.legend(loc="upper left", fontsize=10.5)
for k, m in enumerate(CAM):
    a = fig.add_subplot(gs[k, 1])
    c, _, _ = STY[m]
    pit = (obs[:, None] > cp[QK[m]]).sum(1)                       # bin of 20 between the 19 quantiles
    h = np.bincount(pit, minlength=20) / len(obs) * 20
    a.bar(np.arange(20) / 20 + 0.025, h, width=0.045, color=c)
    a.axhline(1, color=INK, lw=1)
    a.set_ylim(0, 2.0); a.set_xticks([0, 0.5, 1] if k == 2 else []); a.set_yticks([0, 1])
    a.text(0.5, 1.72, f"{m}: PICP90 {ca[m]['PICP 90 %']:.1f} %", ha="center", fontsize=10.5, color=INK)
    a.grid(False)
a.set_xlabel("PIT (flat = calibrated)")
fig.tight_layout(); fig.savefig("img/ca_calibration.png"); plt.close(fig)

# 7. the uncertainty map
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 6.2))
o = np.argsort(Qv[:, 9])
s1 = a1.scatter(cp["lon"][o], cp["lat"][o], c=Qv[o, 9], cmap="Blues", s=5, vmin=0, vmax=5)
fig.colorbar(s1, ax=a1, shrink=0.8).set_label("predicted median ($100k)")
o = np.argsort(width)
s2 = a2.scatter(cp["lon"][o], cp["lat"][o], c=np.minimum(width[o], 2.0), cmap="Oranges", s=5, vmin=0, vmax=2.0)
fig.colorbar(s2, ax=a2, shrink=0.8).set_label("90 % interval width / median")
for a, t in [(a1, "What a district is worth"), (a2, "How sure the vine is")]:
    a.set_title(t, loc="left", color=INK); a.set_xlabel("longitude"); a.set_aspect(1.2); a.grid(False)
a1.set_ylabel("latitude")
fig.tight_layout(); fig.savefig("img/ca_map.png"); plt.close(fig)
print("figures written")
