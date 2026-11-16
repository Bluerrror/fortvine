### Pair-copula evaluation (ns per point, one thread, n = 10^6)

| family | fortvine pdf | fortvine h | fortvine h⁻¹ | pyvinecopulib pdf | pyvinecopulib h | pyvinecopulib h⁻¹ | torchvine pdf | torchvine h | torchvine h⁻¹ |
|---|---|---|---|---|---|---|---|---|---|
| gaussian | 21 | 62 | 31 | 122 | 79 | 98 | 66 | 61 | 57 |
| student | 314 | 462 | 555 | 332 | 352 | 812 | 157 | 1345 | 1432 |
| clayton | 44 | 55 | 43 | 117 | 166 | 123 | 65 | 104 | 86 |
| gumbel | 61 | 82 | 165 | 228 | 223 | 321 | 191 | 133 | 1477 |
| frank | 30 | 30 | 23 | 81 | 156 | 4954 | 43 | 42 | 15 |
| joe | 54 | 75 | 352 | 113 | 280 | 437 | 281 | 189 | 2626 |
| bb1 | 107 | 120 | 770 | 230 | 372 | 12540 | 261 | 265 | 3325 |
| bb6 | 125 | 130 | 978 | 415 | 471 | 16971 | 337 | 266 | 4340 |
| bb7 | 111 | 122 | 829 | 268 | 411 | 15391 | 513 | 348 | 4443 |
| bb8 | 92 | 100 | 612 | 285 | 389 | 11066 | 283 | 256 | 3729 |
| tawn | 158 | 162 | 1177 | 391 | 338 | 10053 | n/a | n/a | n/a |

### Pair-copula MLE fit of a known family (ms, one thread, n = 5000)

| family | fortvine ms | fortvine loglik | pyvinecopulib ms | pyvinecopulib loglik | torchvine ms | torchvine loglik |
|---|---|---|---|---|---|---|
| gaussian | 1.15 | 1683.5 | 8.09 | 1683.5 | 1.11 | 1683.5 |
| student | 22.12 | 1806.8 | 243.55 | 1806.8 | 34.41 | 1806.7 |
| clayton | 1.60 | 2112.8 | 9.16 | 2112.8 | 0.51 | 2104.3 |
| gumbel | 1.71 | 1835.8 | 5.82 | 1835.8 | 0.21 | 1835.3 |
| frank | 1.46 | 1539.6 | 5.60 | 1539.6 | 1.47 | 1515.8 |
| joe | 1.35 | 2202.3 | 11.67 | 2202.3 | 0.51 | 2189.5 |
| bb1 | 4.66 | 1554.9 | 66.24 | 1554.9 | 94.77 | 1554.9 |
| bb6 | 9.42 | 2048.4 | 103.97 | 2048.4 | 89.94 | 2048.4 |
| bb7 | 5.31 | 1928.7 | 36.30 | 1928.7 | 110.28 | 1928.7 |
| bb8 | 12.24 | 1347.7 | 81.05 | 1347.7 | 142.39 | 1347.7 |
| tawn | 8.31 | 1571.1 | 45.57 | 1571.1 | n/a | n/a |

### Family selection among the 12 common families + rotations (ms, one thread, n = 5000)

| truth | fortvine ms | fortvine chosen | pyvinecopulib ms | pyvinecopulib chosen | torchvine ms | torchvine chosen |
|---|---|---|---|---|---|---|
| gaussian | 50.5 | gaussian | 2085.4 | gaussian | 798.7 | gaussian |
| student | 76.9 | student | 1665.4 | student | 1032.0 | student |
| clayton | 76.9 | clayton | 2149.9 | clayton | 1140.5 | bb7 |
| gumbel | 87.9 | gumbel | 1157.7 | gumbel | 1392.0 | gumbel |
| frank | 50.2 | frank | 1508.7 | frank | 963.0 | frank |
| joe | 76.8 | joe | 1868.4 | joe | 2206.0 | bb6 |
| bb1 | 80.1 | bb1 | 1967.8 | bb1 | 1049.8 | bb1 |
| bb6 | 88.8 | bb6 | 1360.1 | bb6 | 1002.5 | bb6 |
| bb7 | 80.6 | bb7 | 1435.1 | bb7 | 1052.8 | bb7 |
| bb8 | 90.8 | bb8 | 1474.0 | bb8 | 1012.4 | bb8 |
| tawn | 107.9 | tawn | 969.2 | tawn | 995.0 | bb6 |

### Vine copula selection, density and simulation (n = 2000)

| d | threads | fortvine select s | fortvine loglik | fortvine pdf s | fortvine simulate s | pyvinecopulib select s | pyvinecopulib loglik | pyvinecopulib pdf s | pyvinecopulib simulate s | torchvine select s | torchvine loglik | torchvine pdf s | torchvine simulate s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 | 1 | 0.492 | 2213.9 | 0.0016 | 0.0074 | 6.941 | 2213.9 | 0.0047 | 0.1065 | 3.983 | 2184.6 | 0.0107 | 0.0466 |
| 5 | 8 | 0.083 | 2213.9 | 0.0004 | 0.0010 | 1.257 | 2213.9 | 0.0021 | 0.0169 | n/a | n/a | n/a | n/a |
| 10 | 1 | 2.036 | 4675.3 | 0.0056 | 0.0147 | 34.327 | 4675.3 | 0.0455 | 0.3952 | 18.465 | 4622.4 | 0.0716 | 0.1361 |
| 10 | 8 | 0.333 | 4675.3 | 0.0009 | 0.0021 | 6.702 | 4675.3 | 0.0089 | 0.0745 | n/a | n/a | n/a | n/a |
| 20 | 1 | 14.054 | 11702.1 | 0.0187 | 0.0414 | 158.899 | 11702.1 | 0.1654 | 1.3845 | n/a | n/a | n/a | n/a |
| 20 | 8 | 1.192 | 11702.1 | 0.0038 | 0.0057 | 27.569 | 11702.1 | 0.0220 | 0.2235 | n/a | n/a | n/a | n/a |
