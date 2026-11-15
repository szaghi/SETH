# Verification

`fobis build --mode tests-gnu && bash scripts/run_tests.sh` runs fifteen test programs; every one sweeps every kernel
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
| `seth_test_convergence` | whole-domain $L_\infty$ convergence on Henrick et al.'s field $\sin(\pi x - \sin(\pi x)/\pi)$ (critical points that are not symmetry points): the spatial operator and SSP-RK3 advection with $\Delta t \sim h^{o/3}$ |
| `seth_test_tables` | the compiled coefficient tables re-verified in Fortran: candidates and linear schemes exact, indicators 0 on constants and 1 (or the Balsara–Shu scale) on unit slope, $\tau$ of its order |
| `seth_test_golden` | regression: 5632 values of every kernel on 64 stencils built in exact binary arithmetic, against `src/tests/golden/seth_golden_v1.txt` |
| `seth_test_device` | every kernel inside an OpenACC parallel loop against the host (in the `tests-nvf-oac` build; on the host otherwise) |

Random stencils come from a fixed, printed seed (`SETH_TEST_SEED` overrides it), so a failure reproduces.

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

## Verification and regression

Most tests verify: they compare the kernels with something they do not share code with (exact rationals, analytic
functions, theorems, published numbers), and assert orders, properties and bounds. Such assertions let a small
change through: a few ulp, or an error that grows inside its bound. The **golden test** pins the values. Its inputs
are exact in binary (no libm), so the same build reproduces them bitwise; across compilers the states agree within
1.1e-15 of the data (nvfortran's fused multiply-adds), the weights ratios, which read the weights themselves, within
3.5e-10. A golden file is never edited: a deliberate change of results gets `seth_golden_v2.txt`.

## On the device

`seth_test_device`, built with `fobis build --mode tests-nvf-oac` and run on an RTX 4070: all 29 cases (the 26
kernels and the descaler) agree with the host; 23 bitwise on every stencil, the others within 8.7e-16 of the data.

## Whole-domain convergence

On Henrick et al.'s field, $L_\infty$ over every node, $\varepsilon = 10^{-40}$:

| Scheme | operator | SSP-RK3 advection |
|---|---|---|
| WENO-JS 5 / 7 / 9 | 3.3 / 5.0 / 6.9 (Henrick's critical-point rates 3, 5, 7) | — |
| WENO5-Z, $p$ = 1 / 2 | 4.6 / 5.2 | — / 5.0 |
| WENO-Z7 / Z9 | 7.0 / 8.9 | — |
| TENO5 / TENO6 / TENO6-opt / TENO6-NENO6 | 5.0 / 6.0 / 5.2 / 4.0 | 5.0 / 6.0 / 5.1 / 4.0 |

On $\sin(\pi x)$ itself WENO5-JS shows 5: its critical points are symmetry points, where the degraded errors cancel;
Henrick et al. chose their field for this reason.

## Testing the tests: mutation testing

`python3 scripts/mutation_test.py` corrupts one thing at a time in a scratch copy of the tree (a coefficient of each
generated table, a linear weight, $\tau$, an indicator pivot, the TENO cutoff and power, the WENO-JS exponent, the
mirroring, the use of $\varepsilon$, the descaler, the weights ratio, the global indicator), rebuilds from scratch and
runs the suite; a mutant is killed when a test fails. **All 18 mutants are killed** (`notes/mutation-run-2026-11-15.txt`);
the runs before that found the gaps:

- **Two TENO mutants survived** every test, a $\tau_6$ coefficient and an indicator pivot: they change only which
  candidates TENO cuts, and none of the test stencils sits near the cutoff, so the outputs stayed bitwise equal. The
  compiled tables are now re-verified (`seth_test_tables`), which kills both.
- **Thin spots** (killed by a single test): ignoring $\varepsilon$ (only the critical-point test passes a non-default
  one), the descaler (only the comparison with ADAM calls it), and changes that keep a valid scheme (the WENO-JS
  exponent, WENO-Z's $p = 2$ branch), which only the golden values and the ADAM comparison can see.
- **FoBiS does not rebuild the files that `#include` a changed `.INC`**: the first mutation runs kept a stale mutant
  compiled in. The script rebuilds every mutant with `--force-compile`; the same applies to development, so rebuild
  with `--force-compile` after editing a kernel body.

## What the tests found

- The WENO-Z and TENO kernels first evaluated the indicators as expanded quadratic forms, which round to about
  $-10^{-14}$ on constant data; with $\varepsilon = 10^{-40}$ TENO then returned NaN on a constant stencil (a
  freestream). The property test caught it; the indicators are now sums of squares.
- The weights the first TENO6-opt used were the NENO6 paper's (4th order), not Fu et al.'s (5th order).
- Two published orders of $\tau_7$ are one too high (measured exactly, confirmed by Li et al. 2022).
