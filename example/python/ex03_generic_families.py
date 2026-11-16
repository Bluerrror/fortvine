"""ex03_generic_families -- the generic engine (codes 10000+): list the families, use an Archimax
copula and an asymmetric Khoudraji copula, and select among classic + generic families.

The Python version of example/ex03_generic_families.f90 (pip install ./python).
"""
import numpy as np

import fortvine as fv

n = 3000

codes = [f.code for f in fv.families() if f.code >= 10000]     # the 85 new families
print(f"generic families: {len(codes)}; the first five:")
for f in fv.families():
    if f.code in codes[:5]:
        print(f"{f.code:7d}  {f.name:<40s}  ({f.npars} parameters)")

# Archimax(gamma | Husler-Reiss) = code 10203: gamma (Clayton) generator, Husler-Reiss Pickands
rng = np.random.default_rng(20261116)                          # reproducible sample
u1, w = rng.random(n), rng.random(n)
truth = fv.Bicop(10203, 0, [1.0, 1.5])
u = np.column_stack([u1, truth.hinv1(np.column_stack([u1, w]))])
fitted = fv.Bicop(10203).fit(u)
print(f"fitted {fitted.family}:" + "".join(f"{p:8.3f}" for p in fitted.parameters))

# Khoudraji-Clayton = code 20200: asymmetric, C(u,v) = u^(1-a) v^(1-b) Clayton(u^a, v^b)
truth = fv.Bicop(20200, 0, [2.0, 0.6, 0.95])
u = np.column_stack([u1, truth.hinv1(np.column_stack([u1, w]))])
fitted = fv.Bicop(20200).fit(u)
print("fitted Khoudraji-Clayton (theta, a, b):" + "".join(f"{p:8.3f}" for p in fitted.parameters))

# selection among every parametric family + every generic family
chosen = fv.Bicop().select(u, family_set=list(range(24)) + codes, num_threads=4)
print(f"selected among all families: {chosen.family}  rotation {chosen.rotation}")
