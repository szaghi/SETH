# SETH

**S**moothness-weighted **E**NO, **T**ENO and **H**ybrids: high-order non-oscillatory reconstruction and
interpolation for HPC, in modern Fortran.

[![CI](https://github.com/szaghi/SETH/actions/workflows/ci.yml/badge.svg)](https://github.com/szaghi/SETH/actions/workflows/ci.yml)
[![Matrix](https://github.com/szaghi/SETH/actions/workflows/matrix.yml/badge.svg)](https://github.com/szaghi/SETH/actions/workflows/matrix.yml)

SETH is extracted from the WENO machinery of ADAM and succeeds WenOOF: one `pure`, device-callable kernel per scheme,
coefficients as compile-time constants, no dynamic dispatch where the cells are swept.

| Family | Kernels |
|---|---|
| WENO-JS | `seth_weno{1,3,5,7,9}_js` (+ `_interp`, `_wr`): ADAM's schemes, bitwise |
| WENO-Z | `seth_wenoz{3,5,7,9}` (+ `_interp`) |
| TENO | `seth_teno5`, `seth_teno6` (+ `_interp`), `seth_teno6_opt`, `seth_teno6_neno6` |
| Positivity | `seth_positivity_euler`: the flux limiter of Hu, Adams & Shu (2013), `seth_euler_admissible` |

Twenty-six verification programs, which are also the test suite, check every kernel with gfortran, nvfortran and
ifx and on the GPU. They cover the coefficients against their derivations and ADAM, the orders at a point, at
critical points and over whole domains, invariances and robustness, golden values, Sod, Lax, 123, Shu-Osher and
Woodward-Colella in a 1-D Euler harness, Burgers against its exact solution, long-time advection, and the spectral and
dispersion properties. A mutation run catches 21 of 21 planted defects. On an RTX 4070 the WENO-5/7 kernels run
1.26-1.31x faster than ADAM's.

```bash
fobis fetch && fobis build --mode tests-gnu && bash scripts/run_tests.sh
```

Documentation: https://szaghi.github.io/SETH. License: GPL-3.0.
