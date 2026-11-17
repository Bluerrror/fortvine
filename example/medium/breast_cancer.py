"""Breast cancer (Wisconsin diagnostic): P(malignant | x) from one vine copula with the diagnosis as a
discrete variable. 10-fold stratified CV against naive Bayes (same margins, no dependence),
logistic regression and gradient boosting."""
import json, warnings
import numpy as np
from sklearn.datasets import load_breast_cancer
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegressionCV
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn import metrics as M
from vinepred import VineClassifier, VineBayes

warnings.filterwarnings("ignore")
data = load_breast_cancer()
FEAT = slice(20, 30)                              # the ten "worst" (largest-nucleus) measurements
NAMES = [n.replace("worst ", "") for n in data.feature_names[FEAT]]
X = np.log1p(data.data[:, FEAT])
y = (data.target == 0).astype(int)                # 1 = malignant
MODELS = {
    "Vine copula": lambda: VineClassifier(num_threads=12),
    "Naive Bayes": lambda: VineBayes(copula=False),
    "Logistic regression": lambda: make_pipeline(StandardScaler(), LogisticRegressionCV(Cs=20, max_iter=5000)),
    "Gradient boosting": lambda: HistGradientBoostingClassifier(random_state=0),
}


def ece(y, p, bins=10):
    b = np.minimum((p * bins).astype(int), bins - 1)
    return sum(abs(y[b == k].mean() - p[b == k].mean()) * np.mean(b == k) for k in range(bins) if np.any(b == k))


def scores(y, p):
    yh = (p > 0.5).astype(int)
    tn, fp, fn, tp = M.confusion_matrix(y, yh).ravel()
    pc = np.clip(p, 1e-15, 1 - 1e-15)
    return {"Accuracy": M.accuracy_score(y, yh), "Balanced accuracy": M.balanced_accuracy_score(y, yh),
            "Sensitivity": tp / (tp + fn), "Specificity": tn / (tn + fp), "Precision": M.precision_score(y, yh),
            "F1": M.f1_score(y, yh), "MCC": M.matthews_corrcoef(y, yh), "ROC AUC": M.roc_auc_score(y, p),
            "PR AUC": M.average_precision_score(y, p), "Log loss": M.log_loss(y, pc),
            "Brier": M.brier_score_loss(y, p), "ECE": ece(y, p), "Missed cancers": int(fn), "False alarms": int(fp)}


if __name__ == "__main__":
    P = {k: np.zeros(len(y)) for k in MODELS}
    for tr, te in StratifiedKFold(10, shuffle=True, random_state=0).split(X, y):
        for k, mk in MODELS.items():
            P[k][te] = mk().fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
    res = {k: scores(y, p) for k, p in P.items()}
    for k, r in res.items():
        print(f"{k:20s}", " ".join(f"{a} {b:.4f}" for a, b in r.items()))
    full = VineClassifier(num_threads=12).fit(X, y)
    res["_vine"] = {"trunc_lvl": full.vine.trunc_lvl, "families": sorted({str(f) for t in full.vine.families for f in t}),
                    "features": NAMES, "n": len(y), "malignant": int(y.sum())}
    json.dump(res, open("bc_results.json", "w"), indent=1)
    np.savez("bc_pred.npz", y=y, **{k.replace(" ", "_"): p for k, p in P.items()})
