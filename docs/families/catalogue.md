# Complete family catalogue

Generated from the library (`tools/catalogue.f90`): every family code fortvine accepts, its
name, number of parameters and parameter box. *Exch.* = exchangeable (C(u,v) = C(v,u));
asymmetric families are transposed exactly when a vine edge is flipped. Every family except
the rotationless ones also has 90, 180 and 270 degree rotations.

## Hand-coded families (codes 0-25)

| code | name | parameters | lower bounds | upper bounds | rotations |
|---:|---|---:|---|---|---|
| 0 | indep | 0 | - | - | no |
| 1 | gaussian | 1 | -1 | 1 | no |
| 2 | student | 2 | -1 2 | 1 30 | no |
| 3 | clayton | 1 | 1.0E-10 | 28 | yes |
| 4 | gumbel | 1 | 1 | 50 | yes |
| 5 | frank | 1 | -35 | 35 | no |
| 6 | joe | 1 | 1 | 30 | yes |
| 7 | bb1 | 2 | 0 1 | 7 7 | yes |
| 8 | bb6 | 2 | 1 1 | 6 8 | yes |
| 9 | bb7 | 2 | 1 0.01 | 6 25 | yes |
| 10 | bb8 | 2 | 1 1.0E-04 | 8 1 | yes |
| 11 | tawn | 3 | 0 0 1 | 1 1 60 | yes |
| 12 | plackett | 1 | 0.001 | 1000 | no |
| 13 | fgm | 1 | -1 | 1 | no |
| 14 | amh | 1 | -1 | 0.9999 | no |
| 15 | galambos | 1 | 0.01 | 30 | yes |
| 16 | husler_reiss | 1 | 0.01 | 30 | yes |
| 17 | bb2 | 2 | 0.01 0.01 | 30 30 | yes |
| 18 | bb3 | 2 | 1 0.01 | 10 30 | yes |
| 19 | bb5 | 2 | 1 0.01 | 10 30 | yes |
| 20 | bb9 | 2 | 1 0.01 | 10 30 | yes |
| 21 | bb10 | 2 | 0.01 0 | 30 1 | yes |
| 22 | bb4 | 2 | 0.01 0.01 | 15 15 | yes |
| 23 | tev | 2 | -0.95 1 | 0.99 30 | yes |
| 24 | tll | (smoothing) | - | - | no |
| 25 | bernstein | (smoothing) | - | - | no |

## Generic families (85 codes, `fv_generic_families()`)

Code = 10000 k + 100 og + ie (see [the generic engine](../manual/generic.md)).

### Archimedean

| code | name | parameters | lower bounds | upper bounds | exch. | lower tail | upper tail |
|---:|---|---:|---|---|:-:|:-:|:-:|
| 11300 | Arch(PS.Sibuya) | 2 | 1 1 | 10 10 | yes | no  | yes |
| 11700 | Arch(gamma.Sibuya) | 2 | 0.001 1 | 15 10 | yes | yes | yes |
| 11800 | Arch(gamma.logseries) | 2 | 0.001 0.001 | 15 30 | yes | yes | no  |
| 12300 | Arch(logseries.PS) | 2 | 0.001 1 | 30 10 | yes | no  | yes |
| 12400 | Arch(logseries.gamma) | 2 | 0.001 0.001 | 30 15 | yes | yes | no  |
| 12500 | Arch(logseries.Sibuya) | 2 | 0.001 1 | 30 10 | yes | no  | yes |

### extreme-value

| code | name | parameters | lower bounds | upper bounds | exch. | lower tail | upper tail |
|---:|---|---:|---|---|:-:|:-:|:-:|
| 10004 | Tawn-mixed | 1 | 0 | 1 | yes | no  | yes |
| 10006 | Tawn-type1 | 2 | 1 0.01 | 20 1 | no  | no  | yes |
| 10007 | Tawn-type2 | 2 | 1 0.01 | 20 1 | no  | no  | yes |
| 10008 | Joe-asymneglogistic | 3 | 0.01 0.01 0.01 | 1 1 20 | no  | no  | yes |
| 10009 | asymneglogistic-type1 | 2 | 0.01 0.01 | 20 1 | no  | no  | yes |
| 10010 | asymneglogistic-type2 | 2 | 0.01 0.01 | 20 1 | no  | no  | yes |
| 10011 | Tawn-asymmixed | 2 | 0 0 | 1 1 | no  | no  | yes |
| 10012 | Khoudraji-HR | 3 | 0.01 0.01 0.01 | 20 1 1 | no  | no  | yes |
| 10013 | Khoudraji-mixed | 3 | 0 0.01 0.01 | 1 1 1 | no  | no  | yes |

### Archimax (1-par outer)

| code | name | parameters | lower bounds | upper bounds | exch. | lower tail | upper tail |
|---:|---|---:|---|---|:-:|:-:|:-:|
| 10103 | Archimax(PS\|HuslerReiss) | 2 | 1 0.01 | 20 20 | yes | no  | yes |
| 10104 | Archimax(PS\|mixed) | 2 | 1 0 | 20 1 | yes | no  | yes |
| 10106 | Archimax(PS\|asymlogistic1) | 3 | 1 1 0.01 | 20 20 1 | no  | no  | yes |
| 10107 | Archimax(PS\|asymlogistic2) | 3 | 1 1 0.01 | 20 20 1 | no  | no  | yes |
| 10109 | Archimax(PS\|asymneglogistic1) | 3 | 1 0.01 0.01 | 20 20 1 | no  | no  | yes |
| 10110 | Archimax(PS\|asymneglogistic2) | 3 | 1 0.01 0.01 | 20 20 1 | no  | no  | yes |
| 10111 | Archimax(PS\|asymmixed) | 3 | 1 0 0 | 20 1 1 | no  | no  | yes |
| 10203 | Archimax(gamma\|HuslerReiss) | 2 | 0.001 0.01 | 20 20 | yes | yes | yes |
| 10204 | Archimax(gamma\|mixed) | 2 | 0.001 0 | 20 1 | yes | yes | yes |
| 10206 | Archimax(gamma\|asymlogistic1) | 3 | 0.001 1 0.01 | 20 20 1 | no  | yes | yes |
| 10207 | Archimax(gamma\|asymlogistic2) | 3 | 0.001 1 0.01 | 20 20 1 | no  | yes | yes |
| 10209 | Archimax(gamma\|asymneglogistic1) | 3 | 0.001 0.01 0.01 | 20 20 1 | no  | yes | yes |
| 10210 | Archimax(gamma\|asymneglogistic2) | 3 | 0.001 0.01 0.01 | 20 20 1 | no  | yes | yes |
| 10211 | Archimax(gamma\|asymmixed) | 3 | 0.001 0 0 | 20 1 1 | no  | yes | yes |
| 10302 | Archimax(Sibuya\|neglogistic) | 2 | 1 0.01 | 20 20 | yes | no  | yes |
| 10303 | Archimax(Sibuya\|HuslerReiss) | 2 | 1 0.01 | 20 20 | yes | no  | yes |
| 10304 | Archimax(Sibuya\|mixed) | 2 | 1 0 | 20 1 | yes | no  | yes |
| 10306 | Archimax(Sibuya\|asymlogistic1) | 3 | 1 1 0.01 | 20 20 1 | no  | no  | yes |
| 10307 | Archimax(Sibuya\|asymlogistic2) | 3 | 1 1 0.01 | 20 20 1 | no  | no  | yes |
| 10309 | Archimax(Sibuya\|asymneglogistic1) | 3 | 1 0.01 0.01 | 20 20 1 | no  | no  | yes |
| 10310 | Archimax(Sibuya\|asymneglogistic2) | 3 | 1 0.01 0.01 | 20 20 1 | no  | no  | yes |
| 10311 | Archimax(Sibuya\|asymmixed) | 3 | 1 0 0 | 20 1 1 | no  | no  | yes |
| 10402 | Archimax(logseries\|neglogistic) | 2 | 0.001 0.01 | 35 20 | yes | no  | yes |
| 10403 | Archimax(logseries\|HuslerReiss) | 2 | 0.001 0.01 | 35 20 | yes | no  | yes |
| 10404 | Archimax(logseries\|mixed) | 2 | 0.001 0 | 35 1 | yes | no  | yes |
| 10406 | Archimax(logseries\|asymlogistic1) | 3 | 0.001 1 0.01 | 35 20 1 | no  | no  | yes |
| 10407 | Archimax(logseries\|asymlogistic2) | 3 | 0.001 1 0.01 | 35 20 1 | no  | no  | yes |
| 10409 | Archimax(logseries\|asymneglogistic1) | 3 | 0.001 0.01 0.01 | 35 20 1 | no  | no  | yes |
| 10410 | Archimax(logseries\|asymneglogistic2) | 3 | 0.001 0.01 0.01 | 35 20 1 | no  | no  | yes |
| 10411 | Archimax(logseries\|asymmixed) | 3 | 0.001 0 0 | 35 1 1 | no  | no  | yes |

### Archimax (2-par outer)

| code | name | parameters | lower bounds | upper bounds | exch. | lower tail | upper tail |
|---:|---|---:|---|---|:-:|:-:|:-:|
| 10501 | Archimax(tiltedPS\|logistic) | 3 | 1 0.01 1 | 10 25 20 | yes | no  | yes |
| 10502 | Archimax(tiltedPS\|neglogistic) | 3 | 1 0.01 0.01 | 10 25 20 | yes | no  | yes |
| 10503 | Archimax(tiltedPS\|HuslerReiss) | 3 | 1 0.01 0.01 | 10 25 20 | yes | no  | yes |
| 10504 | Archimax(tiltedPS\|mixed) | 3 | 1 0.01 0 | 10 25 1 | yes | no  | yes |
| 10601 | Archimax(genSibuya\|logistic) | 3 | 1 0.001 1 | 8 1 20 | yes | no  | yes |
| 10602 | Archimax(genSibuya\|neglogistic) | 3 | 1 0.001 0.01 | 8 1 20 | yes | no  | yes |
| 10603 | Archimax(genSibuya\|HuslerReiss) | 3 | 1 0.001 0.01 | 8 1 20 | yes | no  | yes |
| 10604 | Archimax(genSibuya\|mixed) | 3 | 1 0.001 0 | 8 1 1 | yes | no  | yes |
| 10701 | Archimax(negbin\|logistic) | 3 | 0.01 0 1 | 30 0.99 20 | yes | no  | yes |
| 10702 | Archimax(negbin\|neglogistic) | 3 | 0.01 0 0.01 | 30 0.99 20 | yes | no  | yes |
| 10703 | Archimax(negbin\|HuslerReiss) | 3 | 0.01 0 0.01 | 30 0.99 20 | yes | no  | yes |
| 10704 | Archimax(negbin\|mixed) | 3 | 0.01 0 0 | 30 0.99 1 | yes | no  | yes |
| 11301 | Archimax(PS.Sibuya\|logistic) | 3 | 1 1 1 | 10 10 20 | yes | no  | yes |
| 11302 | Archimax(PS.Sibuya\|neglogistic) | 3 | 1 1 0.01 | 10 10 20 | yes | no  | yes |
| 11303 | Archimax(PS.Sibuya\|HuslerReiss) | 3 | 1 1 0.01 | 10 10 20 | yes | no  | yes |
| 11304 | Archimax(PS.Sibuya\|mixed) | 3 | 1 1 0 | 10 10 1 | yes | no  | yes |
| 11502 | Archimax(gamma.PS\|neglogistic) | 3 | 0.001 1 0.01 | 15 10 20 | yes | yes | yes |
| 11503 | Archimax(gamma.PS\|HuslerReiss) | 3 | 0.001 1 0.01 | 15 10 20 | yes | yes | yes |
| 11504 | Archimax(gamma.PS\|mixed) | 3 | 0.001 1 0 | 15 10 1 | yes | yes | yes |
| 11701 | Archimax(gamma.Sibuya\|logistic) | 3 | 0.001 1 1 | 15 10 20 | yes | yes | yes |
| 11702 | Archimax(gamma.Sibuya\|neglogistic) | 3 | 0.001 1 0.01 | 15 10 20 | yes | yes | yes |
| 11703 | Archimax(gamma.Sibuya\|HuslerReiss) | 3 | 0.001 1 0.01 | 15 10 20 | yes | yes | yes |
| 11704 | Archimax(gamma.Sibuya\|mixed) | 3 | 0.001 1 0 | 15 10 1 | yes | yes | yes |
| 11801 | Archimax(gamma.logseries\|logistic) | 3 | 0.001 0.001 1 | 15 30 20 | yes | yes | yes |
| 11802 | Archimax(gamma.logseries\|neglogistic) | 3 | 0.001 0.001 0.01 | 15 30 20 | yes | yes | yes |
| 11803 | Archimax(gamma.logseries\|HuslerReiss) | 3 | 0.001 0.001 0.01 | 15 30 20 | yes | yes | yes |
| 11804 | Archimax(gamma.logseries\|mixed) | 3 | 0.001 0.001 0 | 15 30 1 | yes | yes | yes |
| 11902 | Archimax(Sibuya.PS\|neglogistic) | 3 | 1 1 0.01 | 10 10 20 | yes | no  | yes |
| 11903 | Archimax(Sibuya.PS\|HuslerReiss) | 3 | 1 1 0.01 | 10 10 20 | yes | no  | yes |
| 11904 | Archimax(Sibuya.PS\|mixed) | 3 | 1 1 0 | 10 10 1 | yes | no  | yes |
| 12302 | Archimax(logseries.PS\|neglogistic) | 3 | 0.001 1 0.01 | 30 10 20 | yes | no  | yes |
| 12303 | Archimax(logseries.PS\|HuslerReiss) | 3 | 0.001 1 0.01 | 30 10 20 | yes | no  | yes |
| 12304 | Archimax(logseries.PS\|mixed) | 3 | 0.001 1 0 | 30 10 1 | yes | no  | yes |
| 12501 | Archimax(logseries.Sibuya\|logistic) | 3 | 0.001 1 1 | 30 10 20 | yes | no  | yes |
| 12502 | Archimax(logseries.Sibuya\|neglogistic) | 3 | 0.001 1 0.01 | 30 10 20 | yes | no  | yes |
| 12503 | Archimax(logseries.Sibuya\|HuslerReiss) | 3 | 0.001 1 0.01 | 30 10 20 | yes | no  | yes |
| 12504 | Archimax(logseries.Sibuya\|mixed) | 3 | 0.001 1 0 | 30 10 1 | yes | no  | yes |

### Khoudraji

| code | name | parameters | lower bounds | upper bounds | exch. | lower tail | upper tail |
|---:|---|---:|---|---|:-:|:-:|:-:|
| 20200 | Khoudraji-Clayton | 3 | 0.001 0.01 0.01 | 20 1 1 | no  | yes | no  |
| 20300 | Khoudraji-Joe | 3 | 1 0.01 0.01 | 20 1 1 | no  | no  | yes |
| 20400 | Khoudraji-Frank | 3 | 0.001 0.01 0.01 | 35 1 1 | no  | no  | no  |
