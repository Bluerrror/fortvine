# ex05_weights_discrete.py -- the Python version of example/ex05_weights_discrete.f90:
# (1) observation weights: a time-decay weighted fit that follows a change in dependence, and
# integer weights = duplicated rows; (2) a vine with a discrete (Poisson count) variable.
import numpy as np
from scipy.stats import poisson
import fortvine as fv

n, LAM = 3000, 2.0
rng = np.random.default_rng(11)

# ---------------------------------------------------------------- (1) weights
u1, w = rng.random(n), rng.random(n)
u2 = np.concatenate([fv.Bicop("clayton", 0, [1.0]).hinv1(np.column_stack([u1, w])[: n // 2]),   # tau 1/3
                     fv.Bicop("clayton", 0, [4.0]).hinv1(np.column_stack([u1, w])[n // 2:])])   # tau 2/3
u = np.column_stack([u1, u2])
fit = fv.Bicop("clayton").fit(u)
print("Clayton fit, equal weights:          theta = %6.3f" % fit.parameters[0])
wt = 0.998 ** (n - 1 - np.arange(n))              # exponential time decay, half-life ~350 rows
fit = fv.Bicop("clayton").fit(u, weights=wt)
print("Clayton fit, time-decay weights:     theta = %6.3f   (effective n = %6.1f)"
      % (fit.parameters[0], wt.sum() ** 2 / (wt ** 2).sum()))

wt = np.where(rng.random(n) > 0.5, 2.0, 1.0)      # integer weights = repeated rows
fit = fv.Bicop("clayton").fit(u, weights=wt)
print("weights 1/2:     theta = %.8f  loglik = %.4f" % (fit.parameters[0], fit.loglik()))
fit = fv.Bicop("clayton").fit(np.repeat(u, wt.astype(int), axis=0))
print("duplicated rows: theta = %.8f  loglik = %.4f" % (fit.parameters[0], fit.loglik()))

# ---------------------------------------------------------------- (2) a discrete variable
# D-vine 1-2-3: gaussian(0.7) and clayton(2) in tree 1, frank(2) in tree 2
truth = fv.Vinecop().set_dvine([1, 2, 3])
truth.fam[0, :2] = [1, 3]; truth.par[0, 0, 0], truth.par[0, 1, 0] = 0.7, 2.0
truth.fam[1, 0] = 5; truth.par[1, 0, 0] = 2.0
z = truth.simulate(n, seed=5)

# variable 2 becomes a Poisson(2) count; the data are F(x) of all three, then F(x-) of the count
k = poisson.ppf(z[:, 1], LAM)
u = np.column_stack([z[:, 0], poisson.cdf(k, LAM), z[:, 2], poisson.cdf(k - 1, LAM)])
print("\ncount variable: %d zeros, maximum %d, mean %.3f" % ((k == 0).sum(), k.max(), k.mean()))

vd = fv.Vinecop().set_dvine([1, 2, 3]).select(u, var_types=["c", "d", "c"])   # discrete-aware
umid = u[:, :3].copy(); umid[:, 1] = 0.5 * (u[:, 1] + u[:, 3])                # naive: mid-point
vn = fv.Vinecop().set_dvine([1, 2, 3]).select(umid)
print("discrete-aware fit:\n%s\nnaive fit:\n%s" % (vd, vn))

truth.var_types = ["c", "d", "c"]
print("loglik of the mixed data: truth %9.2f   discrete-aware fit %9.2f" % (truth.loglik(u), vd.loglik(u)))
print("first three mixed densities:", np.array2string(vd.pdf(u[:3]), precision=3))

r = vd.rosenblatt(u, seed=1)                      # randomised Rosenblatt transform (Brockwell 2007)
print("Rosenblatt column means (1/2):      ", r.mean(axis=0).round(4))
print("Rosenblatt column variances (1/12): ", r.var(axis=0).round(4))
print("Kendall tau (1,2) (1,3) (2,3), ~0:  ", np.round([fv.kendall_tau(r[:, 0], r[:, 1]),
      fv.kendall_tau(r[:, 0], r[:, 2]), fv.kendall_tau(r[:, 1], r[:, 2])], 4))
