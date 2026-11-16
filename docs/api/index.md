# API reference

This section lists every public entity of the Fortran module `fortvine_api` with its exact
signature. For explanations, defaults in context and worked examples see the
[User manual](../manual/index.md); for the `bind(C)` routines see the
[C interface](../manual/c-api.md).

```fortran
use fortvine_api
```

All reals are `real(kind(1.0d0))` (written `real(dp)` below). Array arguments are
assumed-shape. Data matrices are `u(n, d)` with one row per observation.

## Contents

| page | entities |
|---|---|
| [fv_bicop](bicop.md) | the pair-copula type: fields; `pdf`, `logpdf`, `cdf`, `hfunc1`, `hfunc2`, `hinv1`, `hinv2`, `fit`, `select`, `npars`, `npars_eff` |
| [fv_vinecop](vinecop.md) | the vine-copula type: fields; `init`, `set_dvine`, `set_cvine`, `select`, `pdf`, `pdf_grid`, `loglik`, `rosenblatt`, `inverse_rosenblatt`, `simulate`, `npars`, `npars_eff`, `matrix`, `print`, `relax`, `unrelax`, `rel_params` |
| [fv_controls](controls.md) | the options record |
| [Functions](functions.md) | `fv_threads`, `fv_kendall_tau`, `fv_dependence`, `fv_tau_mat`, `fv_pseudo_obs_nd`, `fv_par_bounds`, `fv_fam_name`, `fv_generic_families`, `fam_is_np`, `gen_is`, `gen_valid`, `gen_dup`, `gen_npar`, `gen_exch`, `gen_lt`, `gen_ut`, `gen_flip`, `fv_fitcache_enable`, `fv_fitcache_capacity`, `fv_fitcache_clear`, `fv_fitcache_disable`, `fv_fitcache_stats` |
| [Family codes and constants](constants.md) | `FV_*` family codes, family sets, `FV_NFAM`, `FV_NCLASSIC`, `NAUX`, selection criteria, methods, tree criteria, tree algorithms, tau modes |

## Reading the argument tables

Each routine is documented with its Fortran interface, then a table with one row per dummy
argument:

| column | content |
|---|---|
| name | dummy-argument name, usable as a keyword (`call vc%select(u, ctl, fixed=.false.)`) |
| type | declared type and rank, e.g. `real(dp)(:, :)` |
| intent | `in`, `out`, `inout`; *optional* arguments are marked |
| description | meaning, shape requirements, defaults |

For type-bound procedures the passed object (`self`) is omitted from the table: `call
cop%pdf(u1, u2, out)` passes `cop` as `self`.

## Error handling

The library does not return error codes from evaluation routines. Inputs outside $(0, 1)$ are
clamped to $[10^{-10}, 1 - 10^{-10}]$, parameters outside a family's domain are clipped to the
nearest admissible value, and failed likelihood terms are replaced by a very small
log-density. `fv_vinecop%select` has an optional `info` output, and
`fv_vinecop%rosenblatt` stops with `error stop` if a vine with discrete variables is called
without `u_minus`. Array extents are not checked beyond what the compiler does; pass output
arrays of the documented shape.

## Thread safety

The routines keep no hidden state, with two exceptions:

* the fit cache (`fv_fitcache_*`) is a single global table. Concurrent fits share it safely
  (lookups and stores are serialised), but enabling, clearing or disabling it while fits run
  elsewhere is not safe;
* `fv_vinecop%simulate` draws from the intrinsic `random_number` generator and, with `seed`,
  reseeds it, which affects the rest of the program and other threads using that generator.
