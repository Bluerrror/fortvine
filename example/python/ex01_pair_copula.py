"""ex01_pair_copula -- one bivariate copula: evaluate, simulate, fit a known family, select a family.

The Python version of example/ex01_pair_copula.f90 (pip install ./python).
"""
import numpy as np

import fortvine as fv

n = 2000

# a Clayton copula with theta = 2 (Kendall's tau = theta / (theta + 2) = 0.5)
truth = fv.Bicop("clayton", 0, [2.0])

# density at three points
g = np.array([0.1, 0.5, 0.9])
dens = truth.pdf(np.column_stack([g, g]))
print("Clayton(2) density on the diagonal:" + "".join(f"{x:10.4f}" for x in dens))

# simulate by inverse Rosenblatt: u2 = h^-1(w | u1) with independent uniforms u1, w
rng = np.random.default_rng(20261116)                # reproducible sample
u1, w = rng.random(n), rng.random(n)
u2 = truth.hinv1(np.column_stack([u1, w]))
u = np.column_stack([u1, u2])
print(f"empirical Kendall tau:{fv.kendall_tau(u1, u2):7.3f}")

# maximum likelihood fit of a given family
fitted = fv.Bicop("clayton").fit(u)
print(f"fitted theta:{fitted.parameters[0]:8.4f}   log-likelihood:{fitted.loglik():10.2f}")

# family selection by BIC over the default family set (with rotations)
chosen = fv.Bicop().select(u)
par = np.zeros(3); par[:chosen.npars] = chosen.parameters
print(f"selected: {chosen.family}  rotation {chosen.rotation}  par" + "".join(f"{p:8.3f}" for p in par))
