# Precision

SETH computes in one real kind per build, PENF's `R_P`: R8P (double precision) by default, R4P (single precision)
with the `*-r4` build modes. This page says what R4P changes, what it is verified against, where it stops, and what
it buys, measured.

## Building in R4P

```bash
fobis build --mode static-gnu-r4      # the library, gfortran (static-nvf-r4, static-intel-r4, shared-*-r4)
fobis build --mode tests-gnu-r4 && bash scripts/run_tests.sh --exe-dir exe-r4
```

The modes add `-D_R_P_IS_R4P`, which PENF and SETH must both see; a build where they disagree stops at compile time.
There are no generic interfaces: every kernel is compiled once, for that kind, and a program uses one kind (the two
builds export the same module and procedure names). The kind-dependent literals (regularisers, back-off steps) live in `seth_kind_parameters`, one set
per kind ([Kinds and constants](/reference/kinds)).

## What changes in R4P

- `SETH_TENO_EPS` is $10^{-30}$ (the papers' $10^{-40}$ is below the smallest normal R4P number, and a GPU flushes
  it to zero); pass it as `eps` to WENO-Z and TENO in either kind.
- The WENO-JS, WENO-Z and TENO weights are computed as ratios to the smallest indicator: the same weights, with no
  factor above 1, where the R8P forms overflow R4P's range ($3.4\cdot10^{38}$).
- The positivity limiter backs off by $5\cdot10^{-7}$ (about 4 ulps) instead of $10^{-15}$.

## Where R4P stops

- **Range.** The smoothness indicators are quadratic in the data: data above about $10^{19}$ overflow them, about
  $3\cdot10^{15}$ for WENO9-JS (its integer coefficients reach $2.5\cdot10^6$, and the indicators are computed on the
  raw data, before any descaler).
- **Cancellation.** The expanded WENO-JS indicator cancels to about $10^{-7} v^2$: data with a large mean need the
  descaler `rmu`; the affine-invariant descaler `seth_weno_descaling_ai` alone does not keep WENO-JS ENO on a small
  step over a large offset; subtract the stencil's mean first (Wang & Don).
- **Order.** A scheme of order 6 or more reaches R4P's round-off within the first refinements of a smooth problem
  (a derivative's error near $10^{-5}$ on 20 to 640 points): its design order cannot be observed, only that it
  converged as far as R4P allows. On the coarser grids, above that floor, the errors equal R8P's, so the orders
  measured there are pre-asymptotic.

## What is verified

The test suite runs in R4P (34 programs; the comparisons with ADAM are R8P only) on gfortran, nvfortran, ifx and the
GPU, against expectations measured in R4P: the worst of the three compilers plus 5 %, each with its derivation in
`src/tests/seth_test_parameters.F90` ([test programs](/verification/unit-tests#single-precision)). Every case is held
to its bound. On Shu–Osher the R4P error is R8P's: TENO6 0.1142 against 0.1129 (400 cells, $L_1(\rho)$ against a
4× reference), and the cost hybrids of TENO5 and TENO6, driven by each of the four
sensors, within 0.88 to 1.19 of their pure schemes on Shu–Osher and the blast, in both kinds (gfortran).

## What it buys

`src/bench/seth_bench_kinds.F90` times SETH's kernels in `R_P` (the kernels of `seth_bench_weno` without ADAM's,
whose code is R8P), on the same field and with the same timing: 2²⁴ faces, two warm-up passes, the median of 20
(GPU) or 10 (CPU) timed passes. The same source is built in both kinds (`bench-kinds-nvf-oac` and `bench-nvf-oac-r4`,
`bench-kinds-gnu` and `bench-gnu-r4`) and run in five interleaved rounds, the order alternating. $S$ is the R8P time
over the R4P time per kernel; the decision rule was fixed before measuring: the median $S$ over the 18 kernels at
2 or more recommends R4P on that platform, between 1.2 and 2 leaves it an option, below 1.2 it is not worth it.

| ns per face (median of 5 runs) | RTX 4070, R8P | R4P | $S$ | i9-13900KF core, R8P | R4P | $S$ |
|---|---|---|---|---|---|---|
| `seth_weno5_js` | 1.067 | 0.114 | 9.36 | 31.0 | 46.2 | 0.67 |
| `seth_weno7_js` | 1.814 | 0.132 | 13.74 | 38.4 | 44.9 | 0.85 |
| `seth_wenoz5` ($p = 1$) | 1.043 | 0.351 | 2.97 | 55.7 | 55.7 | 1.00 |
| `seth_teno6` | 2.450 | 0.467 | 5.25 | 111.6 | 90.8 | 1.23 |
| `seth_teno5` | 1.504 | 0.330 | 4.56 | 89.3 | 63.8 | 1.40 |
| `seth_teno5_opt` | 1.499 | 0.327 | 4.58 | 89.7 | 64.8 | 1.38 |
| `seth_teno5_az` | 1.532 | 0.334 | 4.59 | 90.4 | 63.9 | 1.41 |
| `seth_teno8` | 3.988 | 0.781 | 5.11 | 153.0 | 143.6 | 1.07 |
| `seth_teno5_a` | 2.090 | 0.399 | 5.24 | 102.4 | 85.1 | 1.20 |
| `seth_teno8_a` | 4.606 | 0.810 | 5.69 | 157.7 | 151.6 | 1.04 |
| `seth_weno5_js` + `seth_weno_descaling` | 1.163 | 0.136 | 8.55 | 33.5 | 48.9 | 0.69 |
| `seth_weno5_js` + `seth_weno_descaling_ai` | 1.415 | 0.153 | 9.25 | 38.9 | 53.4 | 0.73 |
| `seth_weno5_js` + `seth_mp_bound` | 1.489 | 0.129 | 11.54 | 38.8 | 54.2 | 0.72 |
| `seth_linear5` | 0.206 | 0.105 | 1.96 | 11.4 | 12.6 | 0.90 |
| `seth_linear7` | 0.207 | 0.105 | 1.97 | 7.1 | 6.1 | 1.17 |
| `seth_sensor_lcr` | 0.384 | 0.112 | 3.43 | 23.1 | 25.6 | 0.91 |
| cost hybrid of WENO-Z5 (`seth_sensor_lcr`, `seth_linear5` / `seth_wenoz5`) | 0.731 | 0.130 | 5.62 | 26.6 | 30.2 | 0.88 |
| `seth_keep_pe6` | 0.849 | 0.293 | 2.90 | 44.8 | 45.3 | 0.99 |
| **median over the 18 kernels** | | | **5.17** | | | **0.99** |

**RTX 4070: median $S$ = 5.17, R4P recommended** where its limits above are acceptable. Its FP64 runs at 1/64 of its
FP32 rate: the kernels that compute (WENO-JS, TENO) gain 3 to 14 times, those that only move data (the linear
kernels) the factor 2 of the halved traffic. On a data-centre GPU (A30, H100: FP64 at half the FP32 rate) expect
about 2 at most, from the traffic and the rate alike: not measured here.

**One CPU core (gfortran): median $S$ = 0.99, R4P not worth it.** The face loop is not vectorised (gfortran
`-O3` packs only the stencil copies; each kernel is a call into the library), so halving the width buys nothing per
operation: TENO5 and TENO6 gain 1.2 to 1.4, TENO8, WENO-Z and KEEP-PE nothing, the WENO-JS
family loses (0.67 to 0.85). Part of that loss is the R4P weight form, the ratios to the smallest indicator that
keep the weights inside R4P's range: built in R4P with the R8P form instead, WENO-JS runs 15 to 18 % faster (WENO5-JS
37.9 against 44.6 ns per face, three interleaved runs on 2²² faces), every other kernel unchanged; the rest of the
loss is not attributed. A vectorised caller (several faces per SIMD
register, which needs the kernels inlined into the caller's loop) could change this: not measured. The box was shared with other jobs during these
runs: a few kernels spread widely over the five runs (TENO8 in R4P by 150 ns per face on a median of 144), but
the median $S$ is far from the 1.2 threshold.

Raw outputs and the analysis: `notes/r4p-benchmark.md` and `notes/bench-kinds/`.
