# Family codes and constants

All constants are `integer, parameter, public` in `fortvine_api` unless stated otherwise.

## Hand-coded families

| constant | code | name (`FV_FAMILY_NAME`) | parameters | rotations |
|---|---:|---|---:|:-:|
| `FV_INDEP` | 0 | `indep` | 0 | no |
| `FV_GAUSSIAN` | 1 | `gaussian` | 1 | no |
| `FV_STUDENT` | 2 | `student` | 2 | no |
| `FV_CLAYTON` | 3 | `clayton` | 1 | yes |
| `FV_GUMBEL` | 4 | `gumbel` | 1 | yes |
| `FV_FRANK` | 5 | `frank` | 1 | no |
| `FV_JOE` | 6 | `joe` | 1 | yes |
| `FV_BB1` | 7 | `bb1` | 2 | yes |
| `FV_BB6` | 8 | `bb6` | 2 | yes |
| `FV_BB7` | 9 | `bb7` | 2 | yes |
| `FV_BB8` | 10 | `bb8` | 2 | yes |
| `FV_TAWN` | 11 | `tawn` | 3 | yes |
| `FV_PLACKETT` | 12 | `plackett` | 1 | no |
| `FV_FGM` | 13 | `fgm` | 1 | no |
| `FV_AMH` | 14 | `amh` | 1 | no |
| `FV_GALAMBOS` | 15 | `galambos` | 1 | yes |
| `FV_HUSLER_REISS` | 16 | `husler_reiss` | 1 | yes |
| `FV_BB2` | 17 | `bb2` | 2 | yes |
| `FV_BB3` | 18 | `bb3` | 2 | yes |
| `FV_BB5` | 19 | `bb5` | 2 | yes |
| `FV_BB9` | 20 | `bb9` | 2 | yes |
| `FV_BB10` | 21 | `bb10` | 2 | yes |
| `FV_BB4` | 22 | `bb4` | 2 | yes |
| `FV_TEV` | 23 | `tev` | 2 | yes |
| `FV_TLL` | 24 | `tll` | effective df | no |
| `FV_BERNSTEIN` | 25 | `bernstein` | effective df | no |

Codes 0 to 11 are the parametric families that vinecopulib also provides (`FV_PARAMETRIC`).
The parameter order and boxes are in the [catalogue](../families/catalogue.md) and in
[Pair copulas](../manual/bicop.md#parameter-bounds).

```fortran
character(len=12), parameter, public :: FV_FAMILY_NAME(0:25)
```

Lower-case names of codes 0 to 25. For any code, including generic ones, use
[`fv_fam_name`](functions.md#fv_fam_name).

## Generic families

Codes $\ge 10000$ are built by the [generic engine](../manual/generic.md):

$$
\texttt{fam} = 10000\,k + 100\,\text{og} + \text{ie},
$$

with $k \in \{1, 2\}$ (plain / Khoudraji), og the outer generator (0 to 7, 11 to 26) and ie the
inner Pickands function (0 to 13). There is no named constant per generic code;
[`fv_generic_families`](functions.md#fv_generic_families) returns the list.

## Family sets and sizes

| constant | value | meaning |
|---|---|---|
| `FV_ALL` | `[0, 1, ..., 23]` (24 codes) | every hand-coded parametric family |
| `FV_PARAMETRIC` | `[0, 1, ..., 11]` (12 codes) | vinecopulib's parametric set |
| `FV_NONPAR` | `[24, 25]` | TLL and Bernstein |
| `FV_ONE_PAR` | `[1, 3, 4, 5, 6]` | the one-parameter families with closed-form fits and tau inversion |
| `FV_NCLASSIC` | 26 | number of hand-coded codes (0 to 25) |
| `FV_NFAM` | 192 | capacity of `fv_controls%family_set` |
| `NAUX` | 933 | length of the `aux` vector of TLL / Bernstein ($3 + 30 + 30 \cdot 30$) |

The default family set of `fv_controls` is `family_set(1:23) = [0, ..., 22]`, that is `FV_ALL`
without t-EV.

## Selection criteria (`fv_controls%selcrit`)

| constant | value |
|---|---:|
| `FV_LOGLIK` | 0 |
| `FV_AIC` | 1 |
| `FV_BIC` | 2 |
| `FV_MBIC` | 3 |
| `FV_MBICV` | 4 |

## Estimation methods (`fv_controls%method`)

| constant | value | |
|---|---:|---|
| `FV_MLE` | 0 | maximum likelihood |
| `FV_ITAU` | 1 | Kendall's tau inversion (families 0 to 6) |

## Tree criteria (`fv_controls%tree_crit`, `fv_dependence`)

| constant | value | measure |
|---|---:|---|
| `FV_TREE_TAU` | 0 | Kendall's $\tau$ |
| `FV_TREE_RHO` | 1 | Spearman's $\rho$ |
| `FV_TREE_JOE` | 2 | $-\tfrac12\log(1 - r^2)$ of the normal scores |
| `FV_TREE_HOEFFD` | 3 | Hoeffding's $D$ |
| `FV_TREE_XI` | 4 | symmetrised Chatterjee $\xi$ |
| `FV_TREE_MI` | 5 | mutual information (KSG) |
| `FV_TREE_FERREIRA` | 6 | Ferreira's tail-dependence coefficient |

## Tree algorithms (`fv_controls%tree_algo`)

| constant | value |
|---|---:|
| `FV_PRIM` | 0 |
| `FV_KRUSKAL` | 1 |

## Tau modes (`fv_controls%tau_mode`)

| constant | value | |
|---|---:|---|
| `FV_TAU_AUTO` | 0 | exact for $n < 800$, approximate above |
| `FV_TAU_EXACT` | 1 | exact $\tau_b$ |
| `FV_TAU_APPROX` | 2 | $\tfrac2\pi\arcsin$ of the normal-score correlation |

## Rotations

Rotations are plain integers: 0, 90, 180, 270 (degrees, counter-clockwise). There are no named
constants. See [Rotations and asymmetry](../manual/rotations.md).

## Kernel module

The kernel module `fortvine` additionally exports `fv_abort`
(`logical, public, volatile`, default `.false.`), the host abort flag described in
[Fitting](../manual/fitting.md#aborting-a-running-fit). It is not re-exported by
`fortvine_api`; `use fortvine, only: fv_abort` to reach it.
