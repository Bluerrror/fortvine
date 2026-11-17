**Table 1. Evaluation of the density $c(u,v)$: time per point [ns], one thread, $n = 10^6$.**

| | fortvine (Fortran) | fortvine (Python) | fortvine (R) | pyvinecopulib | rvinecopulib | torchvine |
|---|---:|---:|---:|---:|---:|---:|
| Gaussian | 21 | 20 | 66 | 85 | 155 | 35 |
| Student $t$ | 301 | 323 | 670 | 187 | 310 | 145 |
| Clayton | 43 | 39 | 160 | 66 | 110 | 39 |
| Gumbel | 60 | 57 | 200 | 114 | 275 | 162 |
| Frank | 29 | 26 | 143 | 55 | 180 | 17 |
| Joe | 51 | 50 | 195 | 90 | 280 | 127 |
| BB1 | 105 | 103 | 340 | 156 | 340 | 189 |
| BB6 | 122 | 116 | 430 | 274 | 970 | 234 |
| BB7 | 104 | 103 | 390 | 156 | 350 | 233 |
| BB8 | 87 | 85 | 360 | 123 | 127 | 195 |
| Tawn | 154 | 152 | 390 | 201 | 420 | — |

**Table 2. Evaluation of the $h$-function $h_1(v\mid u)$: time per point [ns], one thread, $n = 10^6$.**

| | fortvine (Fortran) | fortvine (Python) | fortvine (R) | pyvinecopulib | rvinecopulib | torchvine |
|---|---:|---:|---:|---:|---:|---:|
| Gaussian | 61 | 62 | 127 | 72 | 165 | 29 |
| Student $t$ | 441 | 479 | 940 | 251 | 400 | 690 |
| Clayton | 53 | 50 | 225 | 104 | 100 | 55 |
| Gumbel | 79 | 78 | 265 | 111 | 165 | 64 |
| Frank | 30 | 27 | 140 | 93 | 305 | 13 |
| Joe | 71 | 69 | 320 | 161 | 470 | 88 |
| BB1 | 119 | 113 | 420 | 245 | 560 | 125 |
| BB6 | 129 | 127 | 500 | 319 | 810 | 135 |
| BB7 | 114 | 114 | 470 | 246 | 710 | 116 |
| BB8 | 94 | 93 | 400 | 189 | 230 | 80 |
| Tawn | 156 | 154 | 430 | 154 | 265 | — |

**Table 3. Evaluation of the inverse $h_1^{-1}(w\mid u)$: time per point [ns], one thread, $n = 10^6$.**

| | fortvine (Fortran) | fortvine (Python) | fortvine (R) | pyvinecopulib | rvinecopulib | torchvine |
|---|---:|---:|---:|---:|---:|---:|
| Gaussian | 30 | 31 | 60 | 68 | 140 | 28 |
| Student $t$ | 352 | 376 | 770 | 463 | 710 | 701 |
| Clayton | 41 | 39 | 190 | 71 | 133 | 48 |
| Gumbel | 160 | 158 | 380 | 209 | 240 | 725 |
| Frank | 22 | 21 | 100 | 2917 | 10490 | 15 |
| Joe | 336 | 350 | 1040 | 281 | 600 | 1263 |
| BB1 | 752 | 756 | 2700 | 8282 | 19080 | 1691 |
| BB6 | 862 | 869 | 3250 | 11101 | 27650 | 2000 |
| BB7 | 767 | 784 | 3030 | 8461 | 24360 | 1613 |
| BB8 | 582 | 592 | 2520 | 6188 | 8050 | 1364 |
| Tawn | 1010 | 1011 | 2590 | 5070 | 9190 | — |

**Table 4. Maximum-likelihood fit of the true family: time [ms], one thread, $n = 5000$.**

| | fortvine (Fortran) | fortvine (Python) | fortvine (R) | pyvinecopulib | rvinecopulib | torchvine |
|---|---:|---:|---:|---:|---:|---:|
| Gaussian | 1.08 | 0.82 | 1.77 | 4.06 | 8.33 | 0.23 |
| Student $t$ | 21.12 | 22.78 | 48.00 | 154.94 | 309.00 | 19.03 |
| Clayton | 1.45 | 1.22 | 4.76 | 3.34 | 11.03 | 0.21 |
| Gumbel | 1.59 | 1.38 | 3.41 | 5.29 | 19.38 | 0.22 |
| Frank | 1.30 | 1.11 | 3.06 | 2.79 | 8.38 | 0.36 |
| Joe | 1.23 | 1.06 | 2.83 | 4.44 | 14.76 | 0.22 |
| BB1 | 4.49 | 4.63 | 11.85 | 32.75 | 84.00 | 33.56 |
| BB6 | 9.47 | 8.93 | 23.08 | 64.28 | 263.00 | 47.91 |
| BB7 | 5.00 | 5.19 | 14.55 | 29.68 | 127.00 | 36.07 |
| BB8 | 18.84 | 19.18 | 53.00 | 38.12 | 143.00 | 39.23 |
| Tawn | 8.14 | 8.73 | 20.00 | 38.56 | 333.00 | — |

**Table 5. Maximised log-likelihood $\hat\ell$ of the fit in Table 4, and the difference $\Delta\ell = \hat\ell_{\mathrm{lib}} - \hat\ell_{\mathrm{fortvine}}$.**

| | $\hat\ell$ fortvine | $\Delta\ell$ fortvine (Python) | $\Delta\ell$ fortvine (R) | $\Delta\ell$ pyvinecopulib | $\Delta\ell$ rvinecopulib | $\Delta\ell$ torchvine |
|---|---:|---:|---:|---:|---:|---:|
| Gaussian | 1683.53 | 0.00 | 0.00 | 0.00 | 0.00 | -0.01 |
| Student $t$ | 1806.78 | 0.00 | 0.00 | 0.00 | 0.00 | -0.11 |
| Clayton | 2112.84 | 0.00 | 0.00 | 0.00 | 0.00 | -8.51 |
| Gumbel | 1835.76 | 0.00 | 0.00 | 0.00 | 0.00 | -0.51 |
| Frank | 1539.56 | 0.00 | 0.00 | 0.00 | 0.00 | -23.71 |
| Joe | 2202.32 | 0.00 | 0.00 | 0.00 | 0.00 | -12.77 |
| BB1 | 1554.86 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| BB6 | 2048.42 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| BB7 | 1928.69 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| BB8 | 1347.72 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| Tawn | 1571.08 | 0.00 | 0.00 | 0.00 | 0.00 | — |

**Table 6. Selection among the 12 common families and their rotations (BIC): time [ms], one thread, $n = 5000$.**

| | fortvine (Fortran) | fortvine (Python) | fortvine (R) | pyvinecopulib | rvinecopulib | torchvine |
|---|---:|---:|---:|---:|---:|---:|
| Gaussian | 53.1 | 50.4 | 136.7 | 1006.1 | 2910.0 | 419.0 |
| Student $t$ | 78.6 | 75.6 | 200.0 | 827.3 | 2690.0 | 420.6 |
| Clayton | 75.6 | 73.8 | 195.0 | 978.1 | 2140.0 | 530.2 |
| Gumbel | 83.6 | 84.7 | 230.0 | 805.5 | 2270.0 | 419.6 |
| Frank | 51.0 | 48.8 | 136.7 | 892.7 | 3110.0 | 408.3 |
| Joe | 76.4 | 75.9 | 195.0 | 987.8 | 2270.0 | 811.9 |
| BB1 | 77.5 | 79.1 | 220.0 | 908.7 | 2230.0 | 388.2 |
| BB6 | 86.4 | 88.5 | 225.0 | 829.4 | 2190.0 | 474.8 |
| BB7 | 78.8 | 81.5 | 220.0 | 821.3 | 1940.0 | 407.8 |
| BB8 | 88.5 | 88.3 | 235.0 | 945.7 | 2400.0 | 415.5 |
| Tawn | 99.6 | 99.9 | 265.0 | 676.4 | 2250.0 | 424.0 |

**Table 7. Family chosen by BIC in Table 6 (true family in the first column).**

| | fortvine (Fortran) | fortvine (Python) | fortvine (R) | pyvinecopulib | rvinecopulib | torchvine |
|---|---:|---:|---:|---:|---:|---:|
| Gaussian | Gaussian | Gaussian | Gaussian | Gaussian | Gaussian | Gaussian |
| Student $t$ | Student $t$ | Student $t$ | Student $t$ | Student $t$ | Student $t$ | Student $t$ |
| Clayton | Clayton | Clayton | Clayton | Clayton | Clayton | BB7 |
| Gumbel | Gumbel | Gumbel | Gumbel | Gumbel | Gumbel | Gumbel |
| Frank | Frank | Frank | Frank | Frank | Frank | Frank |
| Joe | Joe | Joe | Joe | Joe | Joe | BB6 |
| BB1 | BB1 | BB1 | BB1 | BB1 | BB1 | BB1 |
| BB6 | BB6 | BB6 | BB6 | BB6 | BB6 | BB6 |
| BB7 | BB7 | BB7 | BB7 | BB7 | BB7 | BB7 |
| BB8 | BB8 | BB8 | BB8 | BB8 | BB8 | BB8 |
| Tawn | Tawn | Tawn | Tawn | Tawn | Tawn | BB6 |

**Table 8. Regular-vine selection (Dissmann, BIC): time [s], $n = 2000$.**

| | fortvine (Fortran) | fortvine (Python) | fortvine (R) | pyvinecopulib | rvinecopulib | torchvine |
|---|---:|---:|---:|---:|---:|---:|
| $d = 5$, 1 thread | 0.477 | 0.482 | 1.280 | 3.553 | 8.780 | 1.643 |
| $d = 5$, 8 threads | 0.117 | 0.082 | 0.235 | 0.666 | 1.990 | — |
| $d = 10$, 1 thread | 2.100 | 2.124 | 5.760 | 18.316 | 40.630 | 5.916 |
| $d = 10$, 8 threads | 0.344 | 0.326 | 0.910 | 3.607 | 8.210 | — |
| $d = 20$, 1 thread | 8.754 | 8.946 | 24.190 | 84.577 | 174.080 | — |
| $d = 20$, 8 threads | 1.228 | 1.229 | 3.330 | 14.916 | 29.840 | — |

**Table 9. Log-likelihood $\hat\ell$ of the selected vine, $n = 2000$.**

| | fortvine (Fortran) | fortvine (Python) | fortvine (R) | pyvinecopulib | rvinecopulib | torchvine |
|---|---:|---:|---:|---:|---:|---:|
| $d = 5$, 1 thread | 2213.9 | 2213.9 | 2213.9 | 2213.9 | 2213.9 | 2184.6 |
| $d = 5$, 8 threads | 2213.9 | 2213.9 | 2213.9 | 2213.9 | 2213.9 | — |
| $d = 10$, 1 thread | 4675.3 | 4675.3 | 4675.3 | 4675.3 | 4675.3 | 4622.4 |
| $d = 10$, 8 threads | 4675.3 | 4675.3 | 4675.3 | 4675.3 | 4675.3 | — |
| $d = 20$, 1 thread | 11702.1 | 11702.1 | 11702.1 | 11702.1 | 11702.1 | — |
| $d = 20$, 8 threads | 11702.1 | 11702.1 | 11702.1 | 11702.1 | 11702.1 | — |

**Table 10. Vine density on the $n = 2000$ sample: time [ms].**

| | fortvine (Fortran) | fortvine (Python) | fortvine (R) | pyvinecopulib | rvinecopulib | torchvine |
|---|---:|---:|---:|---:|---:|---:|
| $d = 5$, 1 thread | 1.54 | 1.57 | 5.77 | 4.40 | 15.50 | 3.36 |
| $d = 5$, 8 threads | 0.46 | 0.37 | 1.10 | 1.97 | 2.77 | — |
| $d = 10$, 1 thread | 5.57 | 5.75 | 20.67 | 21.47 | 56.67 | 20.40 |
| $d = 10$, 8 threads | 0.97 | 0.91 | 3.26 | 5.74 | 9.68 | — |
| $d = 20$, 1 thread | 17.38 | 17.40 | 66.00 | 64.56 | 195.00 | — |
| $d = 20$, 8 threads | 2.86 | 2.52 | 9.38 | 13.20 | 32.00 | — |

**Table 11. Simulation of $n = 2000$ draws from the selected vine: time [ms].**

| | fortvine (Fortran) | fortvine (Python) | fortvine (R) | pyvinecopulib | rvinecopulib | torchvine |
|---|---:|---:|---:|---:|---:|---:|
| $d = 5$, 1 thread | 7.15 | 7.28 | 25.38 | 56.46 | 195.00 | 14.03 |
| $d = 5$, 8 threads | 1.18 | 1.07 | 3.83 | 10.36 | 28.18 | — |
| $d = 10$, 1 thread | 14.69 | 14.86 | 51.43 | 236.93 | 810.00 | 38.47 |
| $d = 10$, 8 threads | 2.13 | 2.14 | 8.11 | 38.98 | 116.67 | — |
| $d = 20$, 1 thread | 38.67 | 39.22 | 130.00 | 633.95 | 2170.00 | — |
| $d = 20$, 8 threads | 5.86 | 5.39 | 19.38 | 102.25 | 330.00 | — |

