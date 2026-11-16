# Test programs

`fobis build --mode tests-gnu && bash scripts/run_tests.sh` runs 35 programs (`fpm test` the same set). Every one
exits non-zero on a violated assertion. They pass with gfortran 16.2, nvfortran 26.9 and ifx 2026.2, debug (all
checks and floating-point traps on) and optimised, and in the OpenACC build on the GPU. Random stencils come from a
fixed, printed seed (`SETH_TEST_SEED` overrides it), so a failure reproduces.

## Single precision

The tests compute in the library's kind `R_P`. Their kind-sensitive literals (tolerances, regularisers, data
magnitudes) live in `src/tests/seth_test_parameters.F90`, one set per kind, and `scripts/check_literals.py` fails on
any left in a program. The R4P build (`fobis build --mode tests-gnu-r4 && bash scripts/run_tests.sh --exe-dir exe-r4`)
leaves out the two comparisons with ADAM, which are R8P only, has its own golden file
(`src/tests/golden/seth_golden_r4_v1.txt`), and passes its 33 programs on gfortran, nvfortran and ifx and in the
OpenACC build on the GPU. Its expectations are measured in R4P, the worst of the three compilers plus 5 %, each with
its derivation next to it; the R8P set keeps every test as it was. What R4P changes:

- **Orders.** An order is taken on a pair of grids whose finer error is above R4P's round-off, which for a derivative
  sits near $10^{-5}$ on these grids and, in a time-accurate run, grows with the number of steps. Schemes of order 6
  or more reach it between the first two grids: no order can be observed, and the test accepts that they converged
  as far as R4P allows (a degraded scheme would stay above it). The pairs that remain are coarse, with the same
  errors as in R8P there, so their orders are pre-asymptotic: the loss is recorded per test (0.2 to 0.8 of an order;
  2.65 at the critical points of the interpolations).
- **Data.** Test data beyond R4P's range or resolution are per kind: magnitudes whose squares overflow, a step of
  $10^{-7}$ on 1 (one ulp), densities of $10^{-8}$ next to $O(1)$ cells (their one-sided states round to negative
  pressures, which no limiter can prevent), cell averages that cancel on fine grids.
- **One R4P property.** The affine-invariant descaler alone does not keep WENO-JS ENO on a small step over a large
  offset (subtract the mean, [API](/guide/api)). The cost hybrid of TENO6 with LCR on Shu–Osher, once 3 to 9 times
  pure TENO6 in R4P (2.1 in R8P), was the harness's transmissive inflow drifting under the central background, in both
  kinds; with the inflow held fixed it is within 1.05 in R4P and 1.11 in R8P on every compiler
  ([H3](./hybrids/harness)).

| Program | What it asserts | Page |
|---|---|---|
| `seth_test_kinds` | the working kinds of PENF | — |
| `seth_test_registry` | the test registry: every kernel dispatched, names unique, `_wr` kernels equal the plain ones bitwise | [C2](./contracts/properties) |
| `seth_test_weno_js_adam` | WENO-JS against ADAM's kernels copied verbatim, 1.4 M results (R8P only) | [C1](./contracts/adam-equivalence) |
| `seth_test_central_adam` | `seth_keep_pe6` bitwise equal to ADAM's `compute_kep_flux`, `seth_ducros` to `compute_ducros` (R8P only) | [H2](./hybrids/central) |
| `seth_test_weno_js_order` | WENO-JS design order, reconstruction and interpolation | [A1](./accuracy/order-of-accuracy) |
| `seth_test_weno_js_eno` | WENO-JS across a unit step at every cell boundary | [C2](./contracts/properties) |
| `seth_test_wenoz` | WENO-Z order, ENO, mirror; reconstruction and interpolation | [A1](./accuracy/order-of-accuracy), [C2](./contracts/properties) |
| `seth_test_teno` | TENO order, linear scheme on smooth data, cutoff, large jump, mirror; reconstruction and interpolation | [A1](./accuracy/order-of-accuracy), [C2](./contracts/properties) |
| `seth_test_properties` | exactness, constants, translation, mirror, scaling, robustness, on every scheme | [C2](./contracts/properties) |
| `seth_vv_hybrid` | hybrids in the 1-D Euler harness: the sensor quiet on a smooth wave, the central/TENO5 junction at fifth order, conservation with random switching, Sod/Lax/Shu–Osher against the pure schemes | [H3](./hybrids/harness) |
| `seth_vv_roe` | the Roe flux path with entropy fix of the 1-D Euler harness: design order on a smooth wave (the fix silent), Sod and Lax, the transonic rarefaction with and without the fix, the 123 problem; the blast waves at Fu et al.'s cutoffs without the limiter ([issue #5](https://github.com/szaghi/SETH/issues/5)) | [S5](./shocks/blast#the-flux-issue-5) |
| `seth_test_central` | the KEEP-PE central fluxes of order 2 to 8: order, consistency, the kinetic-energy identity, pressure equilibrium, the device; the Ducros sensors | [H2](./hybrids/central) |
| `seth_test_sensors` | the four sensors on Li, Chen & Ren's static tests: their wavenumbers reproduced, the jump of test (b) flagged alone by theirs, constants smooth, steps not; critical points under refinement; invariance under the field-range normalisation | [H1](./hybrids/sensors) |
| `seth_vv_sensors` | the four sensors on Shu–Osher and the blast waves: positivity, each sensor's cost hybrid of TENO5 and TENO6 within 3 times the pure scheme, the flagged fractions ordered as the static tests predict | [H1](./hybrids/sensors) |
| `seth_test_sensor_fields` | the stencil sensors against the Ducros sensors on analytic 2-D fields: Taylor–Green unflagged by Ducros, every shock cell flagged, Ducros blind to a shock in a strong vortex and to a contact, the compression form blind to an expansion | [H1](./hybrids/sensors) |
| `seth_test_linear` | the linear kernels: design order, equal to TENO5/6/8 and the TENO5/6 interpolations on smooth data (8 eps), mirror bitwise | [Schemes](/guide/schemes#hybrids) |
| `seth_test_mp` | the monotonicity-preserving bound: monotone data inside $[v_0, v_1]$, steps of WENO9-JS, WENO-Z9 and TENO6 inside $[0, 1]$, smooth data unchanged, $\sin^4$ orders with and without it, mirror and scaling bitwise | [Schemes](/guide/schemes#monotonicity-preserving-bound) |
| `seth_test_critical_points` | the literature's orders at first- and second-order critical points (reconstructions) | [A2](./accuracy/critical-points) |
| `seth_test_interpolation` | the interpolations over a whole domain and at critical points | [A2](./accuracy/critical-points), [A3](./accuracy/convergence) |
| `seth_test_convergence` | whole-domain $L_\infty$ of the operator and of SSP-RK3 advection on Henrick's field | [A3](./accuracy/convergence) |
| `seth_test_advection` | the derivative's order in a finite-difference scheme; TENO equal to its linear scheme over a period; the composite profile bounded | [A1](./accuracy/order-of-accuracy), [S7](./shocks/advection) |
| `seth_test_spectral` | modified wavenumber of every linear background | [F1](./spectral/linear-spectra) |
| `seth_test_tables` | the compiled coefficient tables re-verified | [C3](./contracts/tables) |
| `seth_test_golden` | 8192 golden values | [C4](./contracts/golden) |
| `seth_test_device` | every kernel on the device against the host, the positivity limiter included | [C5](./contracts/device), [C7](./contracts/positivity) |
| `seth_test_positivity` | the positivity limiter: admissible, maximal, untouched where not needed | [C7](./contracts/positivity) |
| `seth_vv_sod`, `seth_vv_lax`, `seth_vv_test123` | shock tubes against the exact solution; positivity with and without the limiter | [S1](./shocks/sod), [S2](./shocks/lax), [S3](./shocks/test123) |
| `seth_vv_shu_osher`, `seth_vv_blast` | against a finer reference | [S4](./shocks/shu-osher), [S5](./shocks/blast) |
| `seth_vv_burgers` | Burgers against the exact solution, before and after the shock | [S6](./shocks/burgers) |
| `seth_vv_euler_order` | the design order through the complete Euler solver | [A4](./accuracy/euler-order) |
| `seth_vv_advection` | the composite profile over four periods | [S7](./shocks/advection) |
| `seth_vv_adr` | the approximate dispersion relation of every nonlinear scheme | [F2](./spectral/adr) |

The mutation run (`scripts/mutation_test.py`, [C6](./contracts/mutation)) checks the suite itself.

## Verification and regression

Most programs verify: they compare the kernels with something they share no code with (exact rationals, analytic and
exact solutions, theorems, published numbers) and assert orders, properties and bounds. The L1 bounds of the
nonlinear cases are recorded measurements plus 5%, a regression guard on top of the physics checks (positivity,
low-dissipation orderings, convergence). The golden test pins the values themselves.
