"""ex02_vine -- a 5-dimensional vine copula: simulate from a known D-vine, select a new vine from
the sample (structure, families, parameters), then density, Rosenblatt transform and simulation.

The Python version of example/ex02_vine.f90 (pip install ./python).
"""
import numpy as np

import fortvine as fv

n, d = 3000, 5

# a D-vine 1-2-3-4-5 with Gumbel pairs in tree 1 and Frank pairs in tree 2, independence above
truth = fv.Vinecop().set_dvine([1, 2, 3, 4, 5])
for e in range(d - 1):
    truth.set_pair_copula(0, e, fv.Bicop("gumbel", 0, [2.0]))
for e in range(d - 2):
    truth.set_pair_copula(1, e, fv.Bicop("frank", 0, [3.0]))
u = truth.simulate(n, seed=42)

# select a regular vine: Dissmann's algorithm, all parametric families, BIC
vc = fv.Vinecop.from_data(u, num_threads=4)
print(vc)
print(f"log-likelihood:{vc.loglik(u):12.2f}   parameters:{vc.npars}")

# density, Rosenblatt transform (independent uniforms if the model is right) and new samples
dens = vc.pdf(u)
w = vc.rosenblatt(u)
print("Kendall tau of the first two Rosenblatt columns (should be ~0):"
      f"{fv.kendall_tau(w[:, 0], w[:, 1]):8.4f}")
sim = vc.simulate(10, seed=7)
print("first simulated row:" + "".join(f"{x:8.4f}" for x in sim[0]))
for t in range(2):
    print(f"tree {t + 1} families: " + " ".join(vc.families[t]))
