**Table 1. Evaluation of the density $c(u,v)$: time per point [ns], one thread, $n = 10^6$.**

| | fortvine | pyvinecopulib | torchvine |
|---|---:|---:|---:|
| Gaussian | 21 | 122 | 66 |
| Student $t$ | 314 | 332 | 157 |
| Clayton | 44 | 117 | 65 |
| Gumbel | 61 | 228 | 191 |
| Frank | 30 | 81 | 43 |
| Joe | 54 | 113 | 281 |
| BB1 | 107 | 230 | 261 |
| BB6 | 125 | 415 | 337 |
| BB7 | 111 | 268 | 513 |
| BB8 | 92 | 285 | 283 |
| Tawn | 158 | 391 | — |

**Table 2. Evaluation of the $h$-function $h_1(v\mid u)$: time per point [ns], one thread, $n = 10^6$.**

| | fortvine | pyvinecopulib | torchvine |
|---|---:|---:|---:|
| Gaussian | 62 | 79 | 61 |
| Student $t$ | 462 | 352 | 1345 |
| Clayton | 55 | 166 | 104 |
| Gumbel | 82 | 223 | 133 |
| Frank | 30 | 156 | 42 |
| Joe | 75 | 280 | 189 |
| BB1 | 120 | 372 | 265 |
| BB6 | 130 | 471 | 266 |
| BB7 | 122 | 411 | 348 |
| BB8 | 100 | 389 | 256 |
| Tawn | 162 | 338 | — |

**Table 3. Evaluation of the inverse $h_1^{-1}(w\mid u)$: time per point [ns], one thread, $n = 10^6$.**

| | fortvine | pyvinecopulib | torchvine |
|---|---:|---:|---:|
| Gaussian | 31 | 98 | 57 |
| Student $t$ | 555 | 812 | 1432 |
| Clayton | 43 | 123 | 86 |
| Gumbel | 165 | 321 | 1477 |
| Frank | 23 | 4954 | 15 |
| Joe | 352 | 437 | 2626 |
| BB1 | 770 | 12540 | 3325 |
| BB6 | 978 | 16971 | 4340 |
| BB7 | 829 | 15391 | 4443 |
| BB8 | 612 | 11066 | 3729 |
| Tawn | 1177 | 10053 | — |

**Table 4. Maximum-likelihood fit of the true family: time [ms], one thread, $n = 5000$.**

| | fortvine | pyvinecopulib | torchvine |
|---|---:|---:|---:|
| Gaussian | 1.15 | 8.09 | 1.11 |
| Student $t$ | 22.12 | 243.55 | 34.41 |
| Clayton | 1.60 | 9.16 | 0.51 |
| Gumbel | 1.71 | 5.82 | 0.21 |
| Frank | 1.46 | 5.60 | 1.47 |
| Joe | 1.35 | 11.67 | 0.51 |
| BB1 | 4.66 | 66.24 | 94.77 |
| BB6 | 9.42 | 103.97 | 89.94 |
| BB7 | 5.31 | 36.30 | 110.28 |
| BB8 | 12.24 | 81.05 | 142.39 |
| Tawn | 8.31 | 45.57 | — |

**Table 5. Maximised log-likelihood $\hat\ell$ of the fit in Table 4, and the difference $\Delta\ell = \hat\ell_{\mathrm{lib}} - \hat\ell_{\mathrm{fortvine}}$.**

| | $\hat\ell$ fortvine | $\Delta\ell$ pyvinecopulib | $\Delta\ell$ torchvine |
|---|---:|---:|---:|
| Gaussian | 1683.53 | 0.00 | -0.01 |
| Student $t$ | 1806.78 | 0.00 | -0.11 |
| Clayton | 2112.84 | 0.00 | -8.51 |
| Gumbel | 1835.76 | 0.00 | -0.51 |
| Frank | 1539.56 | 0.00 | -23.71 |
| Joe | 2202.32 | 0.00 | -12.77 |
| BB1 | 1554.86 | 0.00 | 0.00 |
| BB6 | 2048.42 | 0.00 | 0.00 |
| BB7 | 1928.69 | 0.00 | 0.00 |
| BB8 | 1347.72 | 0.00 | 0.00 |
| Tawn | 1571.08 | 0.00 | — |

**Table 6. Selection among the 12 common families and their rotations (BIC): time [ms], one thread, $n = 5000$.**

| | fortvine | pyvinecopulib | torchvine |
|---|---:|---:|---:|
| Gaussian | 50.5 | 2085.4 | 798.7 |
| Student $t$ | 76.9 | 1665.4 | 1032.0 |
| Clayton | 76.9 | 2149.9 | 1140.5 |
| Gumbel | 87.9 | 1157.7 | 1392.0 |
| Frank | 50.2 | 1508.7 | 963.0 |
| Joe | 76.8 | 1868.4 | 2206.0 |
| BB1 | 80.1 | 1967.8 | 1049.8 |
| BB6 | 88.8 | 1360.1 | 1002.5 |
| BB7 | 80.6 | 1435.1 | 1052.8 |
| BB8 | 90.8 | 1474.0 | 1012.4 |
| Tawn | 107.9 | 969.2 | 995.0 |

**Table 7. Family chosen by BIC in Table 6 (true family in the first column).**

| | fortvine | pyvinecopulib | torchvine |
|---|---:|---:|---:|
| Gaussian | Gaussian | Gaussian | Gaussian |
| Student $t$ | Student $t$ | Student $t$ | Student $t$ |
| Clayton | Clayton | Clayton | BB7 |
| Gumbel | Gumbel | Gumbel | Gumbel |
| Frank | Frank | Frank | Frank |
| Joe | Joe | Joe | BB6 |
| BB1 | BB1 | BB1 | BB1 |
| BB6 | BB6 | BB6 | BB6 |
| BB7 | BB7 | BB7 | BB7 |
| BB8 | BB8 | BB8 | BB8 |
| Tawn | Tawn | Tawn | BB6 |

**Table 8. Regular-vine selection (Dissmann, BIC): time [s], $n = 2000$.**

| | fortvine | pyvinecopulib | torchvine |
|---|---:|---:|---:|
| $d = 5$, 1 thread | 0.492 | 6.941 | 3.983 |
| $d = 5$, 8 threads | 0.083 | 1.257 | — |
| $d = 10$, 1 thread | 2.036 | 34.327 | 18.465 |
| $d = 10$, 8 threads | 0.333 | 6.702 | — |
| $d = 20$, 1 thread | 14.054 | 158.899 | — |
| $d = 20$, 8 threads | 1.192 | 27.569 | — |

**Table 9. Log-likelihood $\hat\ell$ of the selected vine, $n = 2000$.**

| | fortvine | pyvinecopulib | torchvine |
|---|---:|---:|---:|
| $d = 5$, 1 thread | 2213.9 | 2213.9 | 2184.6 |
| $d = 5$, 8 threads | 2213.9 | 2213.9 | — |
| $d = 10$, 1 thread | 4675.3 | 4675.3 | 4622.4 |
| $d = 10$, 8 threads | 4675.3 | 4675.3 | — |
| $d = 20$, 1 thread | 11702.1 | 11702.1 | — |
| $d = 20$, 8 threads | 11702.1 | 11702.1 | — |

**Table 10. Vine density on the $n = 2000$ sample: time [ms].**

| | fortvine | pyvinecopulib | torchvine |
|---|---:|---:|---:|
| $d = 5$, 1 thread | 1.63 | 4.69 | 10.65 |
| $d = 5$, 8 threads | 0.40 | 2.11 | — |
| $d = 10$, 1 thread | 5.55 | 45.46 | 71.64 |
| $d = 10$, 8 threads | 0.93 | 8.92 | — |
| $d = 20$, 1 thread | 18.68 | 165.40 | — |
| $d = 20$, 8 threads | 3.75 | 22.05 | — |

**Table 11. Simulation of $n = 2000$ draws from the selected vine: time [ms].**

| | fortvine | pyvinecopulib | torchvine |
|---|---:|---:|---:|
| $d = 5$, 1 thread | 7.37 | 106.52 | 46.62 |
| $d = 5$, 8 threads | 1.05 | 16.87 | — |
| $d = 10$, 1 thread | 14.75 | 395.20 | 136.14 |
| $d = 10$, 8 threads | 2.09 | 74.53 | — |
| $d = 20$, 1 thread | 41.38 | 1384.54 | — |
| $d = 20$, 8 threads | 5.67 | 223.50 | — |

