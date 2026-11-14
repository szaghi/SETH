# Verification

`fobis build --mode tests-gnu && bash scripts/run_tests.sh` runs eleven test programs; every one sweeps every kernel
it applies to, through one registry (`src/tests/seth_test_schemes.F90`). They pass with gfortran 16.2, nvfortran 26.9
and ifx 2026.2, debug (all checks and floating-point traps on) and optimised; the OpenACC build compiles every kernel
as `acc routine seq`. Where a test checks a number from the literature, the source is the one recorded in the
[llm-wiki](https://szaghi.github.io/llm-wiki/) bundles `weno-teno` and `high-order-capturing-schemes`.

## The coefficients

The tables are generated, not typed. The generators check, in exact rational arithmetic, before writing anything:

- every candidate is exact for polynomials up to its degree; the TENO candidates equal Li, Fu & Adams' Table 1;
- the derived Jiang–Shu indicators equal ADAM's WENO-5/7/9 tables (up to their scales 1, 240, 5040);
- the derived linear weights equal ADAM's (WENO, $S = 2..5$) and Fu et al.'s (TENO6), and are exact to the order;
- the order of every $\tau$ on smooth data, measured exactly: 3, 5, 7, 7, 5, 6;
- the $LDL^T$ factors reproduce every indicator exactly.

## The tests

| Test | What it proves |
|---|---|
| `seth_test_weno_js_adam` | the WENO-JS kernels against ADAM's code copied verbatim: 1.4 M results bitwise equal without fused multiply-add (gfortran), within 9 ulp of the data with it |
| `seth_test_weno_js_order`, `seth_test_wenoz`, `seth_test_teno` | design order on smooth data, both states, reconstruction and interpolation; boundedness across a jump; exact mirroring; TENO returns its linear scheme on smooth data; a $10^{30}$ jump stays finite |
| `seth_test_registry` | the test registry: every kernel dispatched, `_wr` kernels equal the plain ones bitwise, weights ratios in $[0,1]$ |
| `seth_test_properties` | on every scheme: polynomial exactness, constants, translation, mirror (bitwise for Z and TENO), power-of-2 scaling (bitwise where the theory says so; plain WENO-JS asserted *not* invariant), finite results from $10^{-200}$ to each scheme's limit, the zero stencil |
| `seth_test_critical_points` | the orders the literature measures at critical points (below) |
| `seth_test_advection` | the derivative's order inside a finite-difference scheme; SSP-RK3 advection: TENO equals its linear scheme over a period of a Gaussian (Fu et al. 2016); the Jiang–Shu composite profile bounded, WENO5-Z less dissipative than WENO5-JS (Borges et al. 2008) |
| `seth_test_spectral` | modified wavenumber of every linear background: no anti-dissipation, error slope $o+1$ at low wavenumbers, the resolution efficiencies of Martín et al. 2006 |

## Critical points

The setting of Henrick et al. (2005) and Borges et al. (2008): the derivative $(F_{1/2} - F_{-1/2})/h$ at a node
where $f' = 0$, $\varepsilon = 10^{-40}$, Tang et al.'s functions $e^x - \sum_{k \le n} x^k/k!$:

| Scheme | first-order point: literature | SETH | second-order point: literature | SETH |
|---|---|---|---|---|
| WENO-JS, $r = 2..5$ | $\max(2r-2-n_{cp}, n_{cp})$ = 1, 3, 5, 7 | 1.0, 3.0, 5.0, 6.6–7.2 | 2, 2, 4, 6 | 2.0, 2.0, 4.0, 5.9–6.1 |
| WENO5-Z, $p = 1$ | 4 (3.93–3.99) | 4.00 | 2 | 2.06 |
| WENO5-Z, $p = 2$ | 5 (4.80–4.98) | 5.0 | — | 2.02 |
| TENO5 | 5 | 5.0 | contested | 2.0 |
| TENO6 / TENO6-opt | 6 / 5 | 6.0 / 5.0 | contested (Fu 2023: none; Zhang 2024b: 2) | 2.0 |

At second-order critical points SETH's TENO measures 2, which sides with Zhang et al. (2024b); the test asserts only
a floor of 1.7 there, since the literature disagrees.

## What the tests found

- The WENO-Z and TENO kernels first evaluated the indicators as expanded quadratic forms, which round to about
  $-10^{-14}$ on constant data; with $\varepsilon = 10^{-40}$ TENO then returned NaN on a constant stencil (a
  freestream). The property test caught it; the indicators are now sums of squares.
- The weights the first TENO6-opt used were the NENO6 paper's (4th order), not Fu et al.'s (5th order).
- Two published orders of $\tau_7$ are one too high (measured exactly, confirmed by Li et al. 2022).
