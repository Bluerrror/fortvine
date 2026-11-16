# fortvine for Python

Python binding to [fortvine](https://github.com/Bluerrror/fortvine), a vine-copula library
written in Fortran: 111 pair-copula families, maximum-likelihood fitting with exact parameter
scores, family selection, Dissmann vine selection, vine densities, Rosenblatt transforms and
simulation, parallelised with OpenMP.

The package is a thin layer over the compiled library. Every density, h-function, fit,
selection and vine operation runs in fortvine's Fortran code; Python converts NumPy arrays and
family names for the C interface and wraps the results. Class and method names follow
[pyvinecopulib](https://github.com/vinecopulib/pyvinecopulib), so switching is easy.

## Install

```bash
pip install fortvine
```

Wheels for Linux (x86_64), macOS (arm64, x86_64) and Windows (x86_64) contain the compiled
library, so no Fortran compiler is needed. On other platforms pip builds from the sdist, which
needs `gfortran` (set `FC` to use another one). From a clone of the repository:

```bash
pip install ./python
```

## Example

```python
import numpy as np
import fortvine as fv

# a pair copula
cop = fv.Bicop("clayton", rotation=0, parameters=[2.0])
u = cop.simulate(2000, seed=1)
print(cop.pdf(u[:3]), cop.hfunc1(u[:3]))
best = fv.Bicop().select(u, family_set=["gaussian", "clayton", "gumbel", "frank"])
print(best)                       # Bicop(family='clayton', rotation=0, parameters=[...])

# a vine copula
rng = np.random.default_rng(0)
x = rng.multivariate_normal(np.zeros(4), 0.6 ** np.abs(np.subtract.outer(range(4), range(4))), 3000)
u = fv.to_pseudo_obs(x)
vc = fv.Vinecop.from_data(u, criterion="bic", num_threads=4)
print(vc)
print(vc.loglik(), vc.npars)
sim = vc.simulate(1000, seed=2)
```

`fv.families()` lists every family with its code and number of parameters. The nonparametric
families `tll` and `bernstein` need the Fortran API (their estimate is not passed through the
C interface).

## Citing

If you use fortvine in published work, please cite

> Shojaeezadeh, S. A. (2026). fortvine: a complete compendium to work on vine copula modeling.
> In preparation.

```bibtex
@unpublished{shojaeezadeh2026fortvine,
  author = {Shojaeezadeh, Shahab Aldin},
  title  = {fortvine: a complete compendium to work on vine copula modeling},
  year   = {2026},
  note   = {In preparation. Software: https://github.com/Bluerrror/fortvine}
}
```

`fortvine.cite()` prints the same reference.

Documentation: <https://bluerrror.github.io/fortvine/python/>. Licence: GPL-3.0-or-later.
