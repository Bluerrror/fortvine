# `fv_controls`

Options of fitting and selection. Manual: [Controls](../manual/controls.md).

```fortran
type, public :: fv_controls
  integer  :: family_set(FV_NFAM) = [FV_ALL, FV_NONPAR, spread(0, 1, FV_NFAM - FV_NCLASSIC)]
  integer  :: nfam = size(FV_ALL) - 1        ! 23: codes 0..22
  integer  :: selcrit = FV_BIC
  integer  :: method = FV_MLE
  integer  :: tree_crit = FV_TREE_TAU
  integer  :: tree_algo = FV_PRIM
  integer  :: tau_mode = FV_TAU_EXACT
  integer  :: trunc_lvl = -1                 ! -1 = d - 1
  integer  :: nthreads = 0                   ! 0 = FORTVINE_THREADS or half the logical CPUs
  logical  :: preselect = .true.
  logical  :: allow_rotations = .true.
  logical  :: fast_mle = .false.
  real(dp) :: threshold = 0.0_dp
  real(dp) :: psi0 = 0.9_dp
  logical  :: select_trunc = .false.
  logical  :: select_threshold = .false.
end type fv_controls
```

## Fields

| name | type | default | allowed values | description |
|---|---|---|---|---|
| `family_set` | `integer(FV_NFAM)` | `FV_ALL`, `FV_NONPAR`, zeros | any valid family codes | candidate families, first `nfam` used |
| `nfam` | `integer` | 23 | 1 to `FV_NFAM` (192) | number of candidates |
| `selcrit` | `integer` | `FV_BIC` | `FV_LOGLIK`, `FV_AIC`, `FV_BIC`, `FV_MBIC`, `FV_MBICV` | selection criterion |
| `method` | `integer` | `FV_MLE` | `FV_MLE`, `FV_ITAU` | estimation method |
| `tree_crit` | `integer` | `FV_TREE_TAU` | `FV_TREE_TAU` ... `FV_TREE_FERREIRA` | spanning-tree weight |
| `tree_algo` | `integer` | `FV_PRIM` | `FV_PRIM`, `FV_KRUSKAL` | spanning-tree algorithm |
| `tau_mode` | `integer` | `FV_TAU_EXACT` | `FV_TAU_AUTO`, `FV_TAU_EXACT`, `FV_TAU_APPROX` | Kendall's $\tau$ for starts and preselection |
| `trunc_lvl` | `integer` | $-1$ | $-1$ or $0 \dots d-1$ | last tree fitted |
| `nthreads` | `integer` | 0 | $\ge 0$ | threads; 0 = automatic |
| `preselect` | `logical` | `.true.` | | tail-based preselection of candidates |
| `allow_rotations` | `logical` | `.true.` | | try rotations 90, 180, 270 |
| `fast_mle` | `logical` | `.false.` | | $\tau$ inversion for Gaussian, Clayton, Gumbel, Frank, Joe |
| `threshold` | `real(dp)` | 0 | $[0, 1]$ | edges with a weaker tree criterion become independence |
| `psi0` | `real(dp)` | 0.9 | $(0, 1)$ | mBIC / mBICv prior of a dependent edge |
| `select_trunc` | `logical` | `.false.` | | truncation level by mBICv |
| `select_threshold` | `logical` | `.false.` | | threshold by mBICv |

The constants are listed in [Family codes and constants](constants.md).

## Example

```fortran
type(fv_controls) :: ctl
ctl%nfam = FV_NCLASSIC          ! all 26 hand-coded families
ctl%selcrit = FV_MBICV
ctl%select_trunc = .true.
ctl%nthreads = 8
call vc%select(u, ctl)
ctl = fv_controls()             ! back to the defaults
```
