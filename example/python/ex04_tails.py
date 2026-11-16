"""ex04_tails -- tail dependence of five pair copulas, three ways: the limit from the family
formula, the exact coefficient at a finite level from the copula CDF, and the empirical
coefficient of a simulated sample.  Then: which family does BIC select from each sample?

    lower tail  lambda_L(q) = C(q, q) / q                      -> lambda_L as q -> 0
    upper tail  lambda_U(q) = (2 q - 1 + C(1-q, 1-q)) / q      -> lambda_U as q -> 0

The Python version of example/ex04_tails.f90 (pip install ./python).
"""
import numpy as np

import fortvine as fv

n, m = 200000, 2000

# the five copulas (parameters as in the Fortran program)
cops = [fv.Bicop("gumbel", 0, [2.0]),
        fv.Bicop("clayton", 0, [2.0]),
        fv.Bicop("clayton", 180, [2.0]),
        fv.Bicop("tawn", 0, [0.4, 0.9, 3.0]),            # psi1, psi2, theta
        fv.Bicop(20300, 0, [3.0, 0.5, 0.9])]             # Khoudraji-Joe: theta, a, b

# limits (lambda_L, lambda_U) from the family formulas
th = 2.0
lim = [(0.0, 2 - 2 ** (1 / th)),                         # Gumbel
       (2 ** (-1 / th), 0.0),                            # Clayton
       (0.0, 2 ** (-1 / th))]                            # survival Clayton
for a, b, th in [(0.4, 0.9, 3.0), (0.5, 0.9, 3.0)]:      # Tawn, Khoudraji-Joe
    lim.append((0.0, a + b - (a ** th + b ** th) ** (1 / th)))


def simulate(c, seed, size=n):
    """size draws of (u1, u2) from c by the inverse Rosenblatt transform, reproducible from seed"""
    rng = np.random.default_rng(seed)
    u1, w = rng.random(size), rng.random(size)
    return np.column_stack([u1, c.hinv1(np.column_stack([u1, w]))])   # u2 = h1^-1(w | u1)


def lam_lower(c, q):
    return c.cdf([[q, q]])[0] / q


def lam_upper(c, q):                                      # q is the distance from 1
    return (2 * q - 1 + c.cdf([[1 - q, 1 - q]])[0]) / q


f8 = lambda xs: "".join(f"{x:8.3f}" for x in xs)
print("                        |      lower tail lambda_L(q)         |      upper tail lambda_U(q)")
print("copula              rot |      limit  q=1e-2  q=1e-4  sample |      limit  q=1e-2  q=1e-4  sample")
for k, c in enumerate(cops):
    u = simulate(c, 1000 + k + 1)
    low = np.mean((u[:, 0] < 0.01) & (u[:, 1] < 0.01)) / 0.01
    upp = np.mean((u[:, 0] > 0.99) & (u[:, 1] > 0.99)) / 0.01
    print(f"{c.family:<20s}{c.rotation:3d} |   "
          + f8([lim[k][0], lam_lower(c, 1e-2), lam_lower(c, 1e-4), low]) + " |   "
          + f8([lim[k][1], lam_upper(c, 1e-2), lam_upper(c, 1e-4), upp]))

# BIC selection from m observations of each copula: default family set (hand-coded parametric
# families with rotations), then the same plus every generic family
allfam = list(range(24)) + [f.code for f in fv.families() if f.code >= 10000]
print("\nBIC selection from 2000 observations:")
for k, c in enumerate(cops):
    u = simulate(c, 1000 + k + 1)[:m]
    d = fv.Bicop().select(u, num_threads=4)
    g = fv.Bicop().select(u, family_set=allfam, num_threads=4)
    print(f"  {c.family:<14s} default set: {d.family:<10s} rot{d.rotation:3d}"
          f"   + generic families: {g.family:<14s} rot{g.rotation:3d}")
