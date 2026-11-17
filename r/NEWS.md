# fortvine 1.0.0

* First release: R bindings to the fortvine Fortran library through its C interface
  (`fortvine_capi`, the C interface to the typed API).
* Pair copulas: `bicop_dist()`, `dbicop()`, `pbicop()`, `hbicop()`, `rbicop()`, `bicop()`
  (fit or select, with observation `weights`, `fast_mle` and `tau_mode`), `bicop_par_bounds()`,
  `gen_flip()`, `logLik()`, `print()` and `summary()` methods.
* Nonparametric families `tll` and `bernstein` (family set `"nonparametric"`): their estimate
  is stored in the model (`aux`) and used by every evaluation; `npars` is the effective degrees
  of freedom.
* Vines: `dvine_structure()`, `cvine_structure()`, `as_rvine_matrix()`, `vinecop_dist()`,
  `vinecop()` (Dissmann selection or a fixed structure; `trunc_lvl = NA` and `threshold = NA`
  select the truncation level and the threshold by mBICv; `weights`; discrete margins through
  `var_types`), `dvinecop()`, `dvinecop_grid()` (density on a grid of one variable),
  `rvinecop()`, `rosenblatt()` (randomised for discrete variables), `inverse_rosenblatt()`.
* Non-simplified vines: `vinecop_relax()`, `vinecop_rel_params()`, `vinecop_unrelax()`.
* Utilities: `pseudo_obs()`, `fortvine_families()` (the 111 families of the library) and the
  fit cache (`fitcache_enable()`, `fitcache_clear()`, `fitcache_disable()`,
  `fitcache_stats()`).
