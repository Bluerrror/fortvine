"""California housing: the full predictive distribution of a district's median house value.
80/20 split. Vine copula, gradient-boosted quantile regression (19 models) and a Gaussian linear model;
all scored on the price scale ($100k)."""
import json, warnings
import numpy as np
from scipy.stats import norm
from sklearn.datasets import fetch_california_housing
from sklearn.model_selection import train_test_split
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import LinearRegression
from vinepred import VineRegressor

warnings.filterwarnings("ignore")
LEVELS = np.round(np.arange(0.05, 0.96, 0.05), 2)          # 19 quantile levels


def features():
    h = fetch_california_housing()
    X = h.data.copy()
    for j in (2, 3, 4, 5):                                  # AveRooms, AveBedrms, Population, AveOccup
        X[:, j] = np.log(X[:, j])
    lat, lon = X[:, 6], X[:, 7]
    km = lambda a, b: np.log1p(np.hypot((lat - a) * 111, (lon - b) * 111 * np.cos(np.radians(a))))
    X = np.column_stack([X, km(34.05, -118.24), km(37.77, -122.42), lat + lon, lat - lon])
    return X, np.log(h.target), h.data[:, 6], h.data[:, 7]


def scores(obs, Q):
    """obs: prices; Q: (n, 19) predictive quantiles at LEVELS (price scale); point = median."""
    p = Q[:, 9]
    r = np.corrcoef(obs, p)[0, 1]
    lo, hi = Q[:, 0], Q[:, -1]
    pin = np.mean([np.mean(np.maximum(a * (obs - Q[:, i]), (a - 1) * (obs - Q[:, i]))) for i, a in enumerate(LEVELS)])
    return {"NSE": 1 - np.sum((obs - p) ** 2) / np.sum((obs - obs.mean()) ** 2),
            "KGE": 1 - np.sqrt((r - 1) ** 2 + (p.std() / obs.std() - 1) ** 2 + (p.mean() / obs.mean() - 1) ** 2),
            "r": r, "RMSE": np.sqrt(np.mean((obs - p) ** 2)), "MAE": np.mean(np.abs(obs - p)),
            "PBIAS %": 100 * np.sum(p - obs) / np.sum(obs),
            "CRPS": 2 * pin,                                    # quantile-score approximation, 19 levels
            "PICP 90 %": 100 * np.mean((obs >= lo) & (obs <= hi)), "MPIW 90": np.mean(hi - lo),
            "Interval score": np.mean(hi - lo + 20 * (lo - obs) * (obs < lo) + 20 * (obs - hi) * (obs > hi))}


if __name__ == "__main__":
    X, y, LAT, LON = features()
    idx = np.arange(len(y))
    tr, te = train_test_split(idx, test_size=0.2, random_state=0)
    Xtr, Xte, ytr, yte = X[tr], X[te], y[tr], y[te]
    hgb = lambda **k: HistGradientBoostingRegressor(max_iter=500, random_state=0, **k)
    Q, dens = {}, {}

    lin = LinearRegression().fit(Xtr, ytr)
    s = np.std(ytr - lin.predict(Xtr))
    Q["Linear (Gaussian)"] = lin.predict(Xte)[:, None] + s * norm.ppf(LEVELS)[None, :]

    Q["Quantile boosting"] = np.sort(np.column_stack(
        [hgb(loss="quantile", quantile=a).fit(Xtr, ytr).predict(Xte) for a in LEVELS]), 1)

    v = VineRegressor(family_set="parametric", num_threads=12).fit(Xtr, ytr)
    f, F = v.predict_dist(Xte, num_threads=12)
    Q["Vine copula"] = v.quantiles(F, LEVELS)
    dens["vine"] = f

    obs = np.exp(yte)
    Qp = {k: np.exp(q) for k, q in Q.items()}
    res = {k: scores(obs, q) for k, q in Qp.items()}
    for k, r in res.items():
        print(f"{k:18s}", " ".join(f"{a} {b:.4f}" for a, b in r.items()))
    res["_info"] = {"n_train": len(tr), "n_test": len(te), "capped_share": float(np.mean(np.exp(y) > 5)),
                    "trunc_vine": v.vine.trunc_lvl}
    json.dump(res, open("ca_results.json", "w"), indent=1)
    np.savez("ca_pred.npz", obs=obs, lat=LAT[te], lon=LON[te], levels=LEVELS, yg=v.yg,
             dens_vine=dens["vine"].astype(np.float32), medinc=X[te, 0],
             **{"Q_" + k.replace(" ", "_").replace("(", "").replace(")", ""): q
                for k, q in Qp.items()})
