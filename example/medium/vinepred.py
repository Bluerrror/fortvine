"""Probabilistic prediction with vine copulas (fortvine): the two models of the Medium post."""
import numpy as np
from scipy.stats import gaussian_kde
from scipy.special import logsumexp
import fortvine as fv

EPS = 1e-6


class Margin:
    """One variable's distribution: a Gaussian KDE, tabulated on a grid (pdf and cdf by interpolation)."""

    def __init__(self, x, m=2048):
        x = np.asarray(x, float)
        kde = gaussian_kde(x)
        h = np.sqrt(kde.covariance[0, 0])
        self.g = np.linspace(x.min() - 5 * h, x.max() + 5 * h, m)
        f = kde(self.g)
        c = np.concatenate([[0], np.cumsum((f[1:] + f[:-1]) / 2 * np.diff(self.g))])
        self.f, self.c = f / c[-1], c / c[-1]

    def pdf(self, x):
        return np.interp(x, self.g, self.f, left=0, right=0)

    def cdf(self, x):
        return np.clip(np.interp(x, self.g, self.c), EPS, 1 - EPS)


def to_u(margins, X):
    return np.column_stack([m.cdf(X[:, j]) for j, m in enumerate(margins)])


class VineBayes:
    """Class-conditional densities f(x | k) = prod_j f_kj(x_j) * c_k(u); P(k | x) by Bayes' rule.
    copula=False gives naive Bayes with the same margins (the independence copula)."""

    def __init__(self, copula=True, **select):
        self.copula, self.select = copula, select

    def fit(self, X, y):
        self.classes = np.unique(y)
        self.parts = []
        for k in self.classes:
            Xk = X[y == k]
            mg = [Margin(Xk[:, j]) for j in range(X.shape[1])]
            vine = fv.Vinecop.from_data(to_u(mg, Xk), **self.select) if self.copula else None
            self.parts.append((np.log(len(Xk) / len(X)), mg, vine))
        return self

    def log_joint(self, X):
        out = []
        for logprior, mg, vine in self.parts:
            lf = sum(np.log(np.maximum(m.pdf(X[:, j]), 1e-300)) for j, m in enumerate(mg))
            if vine is not None:
                lf = lf + np.log(np.maximum(vine.pdf(to_u(mg, X)), 1e-300))
            out.append(logprior + lf)
        return np.column_stack(out)

    def predict_proba(self, X):
        lj = self.log_joint(X)
        return np.exp(lj - logsumexp(lj, axis=1, keepdims=True))


class VineClassifier:
    """One vine on (x, y) with the class label y as a discrete variable (labels 0..K-1).
    P(y = k | x) is proportional to P(y = k) c(u_x, F_y(k), F_y(k-)); the margins of x cancel."""

    def __init__(self, **select):
        self.select = select

    def fit(self, X, y):
        self.mx = [Margin(X[:, j]) for j in range(X.shape[1])]
        self.prior = np.bincount(y) / len(y)
        self.Fy = np.cumsum(self.prior)
        self.vine = fv.Vinecop.from_data(self._u(X, y), var_types=["c"] * X.shape[1] + ["d"], **self.select)
        return self

    def _u(self, X, y):
        Fm = np.where(y > 0, self.Fy[np.maximum(y - 1, 0)], 0.0)
        return np.column_stack([to_u(self.mx, X), self.Fy[y], Fm])   # F(x), F(y), then F(y-)

    def predict_proba(self, X):
        n = len(X)
        lp = np.column_stack([np.log(np.maximum(self.vine.pdf(self._u(X, np.full(n, k))), 1e-300))
                              + np.log(self.prior[k]) for k in range(len(self.prior))])
        return np.exp(lp - logsumexp(lp, axis=1, keepdims=True))


class VineRegressor:
    """A vine on (x, y); predictive density f(y | x) = f_y(y) c(u_x, F_y(y)) / c_x(u_x), tabulated on a
    grid of y values (the denominator is a normalising constant)."""

    def __init__(self, bins=400, **select):
        self.bins, self.select = bins, select

    def fit(self, X, y):
        self.mx = [Margin(X[:, j]) for j in range(X.shape[1])]
        self.my = Margin(y)
        U = np.column_stack([to_u(self.mx, X), self.my.cdf(y)])
        self.vine = fv.Vinecop.from_data(U, **self.select)
        lo, hi = np.quantile(y, [0.0, 1.0])
        pad = 0.05 * (hi - lo)
        self.yg = np.linspace(lo - pad, hi + pad, self.bins)
        return self

    def predict_dist(self, X, num_threads=None):
        """Predictive density on self.yg, one row per sample, and its cdf."""
        n, d = X.shape
        ug = np.empty((n, self.bins, d + 1))
        ug[:, :, :d] = to_u(self.mx, X)[:, None, :]
        ug[:, :, d] = self.my.cdf(self.yg)[None, :]
        f = self.vine.pdf_grid(ug, d + 1, num_threads=num_threads) * self.my.pdf(self.yg)[None, :]
        dy = np.gradient(self.yg)
        f /= (f * dy).sum(1, keepdims=True)
        F = np.cumsum(f * dy, 1)
        return f, F / F[:, -1:]

    def quantiles(self, F, q):
        return np.array([[np.interp(p, Fi, self.yg) for p in q] for Fi in F])


def crps_grid(F, yg, yobs):
    """CRPS of predictive cdfs F (n, bins) on grid yg at observations yobs."""
    step = (yg[None, :] >= yobs[:, None]).astype(float)
    return ((F - step) ** 2 * np.gradient(yg)[None, :]).sum(1)


if __name__ == "__main__":
    # self-check: the regressor recovers a known Gaussian conditional, y | x ~ N(0.8 x, 0.36)
    rng = np.random.default_rng(1)
    x = rng.standard_normal(4000)
    y = 0.8 * x + 0.6 * rng.standard_normal(4000)
    r = VineRegressor(family_set=["gaussian"]).fit(x[:, None], y)
    f, F = r.predict_dist(np.array([[0.0], [1.5]]))
    q = r.quantiles(F, [0.5, 0.8413])
    assert abs(q[1, 0] - 1.2) < 0.08 and abs((q[0, 1] - q[0, 0]) - 0.6) < 0.06, q
    print("vinepred self-check ok", q.round(3).tolist())
