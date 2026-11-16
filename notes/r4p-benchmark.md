# R4P against R8P: what single precision buys, kernel by kernel

Night 2026-11-16 (phase 5 of issue #4, 6.2). Program: `src/bench/seth_bench_kinds.F90`, the SETH variants of
`seth_bench_weno` in `R_P` (ADAM's variants are R8P by construction, and the old `bench-*-r4` modes, which pointed at
`seth_bench_weno`, did not build: REAL(8) passed to REAL(4) dummies). Same field (2²⁴ faces, a smooth wave with a jump
every 2048 cells), same timing (2 warm-up passes, then the median of `reps` timed passes, `!$acc wait` +
`system_clock`), and a device-vs-host check on 65536 faces (64 ulps; every run passed, with 0 faces beyond it).
Raw outputs: `notes/bench-kinds/<platform>/{r8p,r4p}/run-<k>.txt`; analysis: `python3 notes/bench-kinds/speedup.py
notes/bench-kinds/<platform>`.

## Decision rule (fixed before measuring, night LOG 13:36)

S = t(R8P)/t(R4P) per kernel, the median over the kernels, per platform, from interleaved runs:
S ≥ 2 recommend R4P there; 1.2 ≤ S < 2 an option, documented, not advertised; S < 1.2 not worth it.

## Runs

| Platform | Builds | Runs |
|---|---|---|
| RTX 4070 (device 1, idle), driver 610.88 | `bench-kinds-nvf-oac` / `bench-nvf-oac-r4`, nvfortran 26.9, `-O4 -fast -acc -gpu=cc89` (release profile) | 5 rounds, order alternating, 20 timed passes, 13:51–13:52 |
| one i9-13900KF core (`taskset -c 2`), WSL2 | `bench-kinds-gnu` / `bench-gnu-r4`, gfortran 16.2.0, `-O3 -funroll-loops` | 5 rounds, order alternating, 10 timed passes, 13:53–14:35, the box shared with other jobs (load 2–11) |

## Results (ns per face, median of the 5 runs; spread max − min in brackets)

| Kernel | GPU R8P | GPU R4P | S | CPU R8P | CPU R4P | S |
|---|---|---|---|---|---|---|
| weno5_js | 1.067 (0.004) | 0.114 (0.001) | 9.36 | 31.0 (3.4) | 46.2 (1.9) | 0.67 |
| weno7_js | 1.814 (0.011) | 0.132 (0.001) | 13.74 | 38.4 (1.9) | 44.9 (1.4) | 0.85 |
| wenoz5 (p = 1) | 1.043 (0.007) | 0.351 (0.004) | 2.97 | 55.7 (3.7) | 55.7 (5.2) | 1.00 |
| teno6 | 2.450 (0.005) | 0.467 (0.001) | 5.25 | 111.6 (8.2) | 90.8 (8.7) | 1.23 |
| teno5 | 1.504 (0.004) | 0.330 (0.003) | 4.56 | 89.3 (4.3) | 63.8 (20.6) | 1.40 |
| teno5_opt | 1.499 (0.004) | 0.327 (0.001) | 4.58 | 89.7 (5.8) | 64.8 (16.2) | 1.38 |
| teno5_az | 1.532 (0.003) | 0.334 (0.001) | 4.59 | 90.4 (2.0) | 63.9 (14.3) | 1.41 |
| teno8 | 3.988 (0.007) | 0.781 (0.002) | 5.11 | 153.0 (13.1) | 143.6 (150.2) | 1.07 |
| teno5_a | 2.090 (0.004) | 0.399 (0.003) | 5.24 | 102.4 (4.1) | 85.1 (21.3) | 1.20 |
| teno8_a | 4.606 (0.010) | 0.810 (0.002) | 5.69 | 157.7 (25.8) | 151.6 (42.7) | 1.04 |
| weno5_js + descaling | 1.163 (0.004) | 0.136 (0.006) | 8.55 | 33.5 (15.8) | 48.9 (3.4) | 0.69 |
| weno5_js + descaling_ai | 1.415 (0.003) | 0.153 (0.003) | 9.25 | 38.9 (16.0) | 53.4 (6.6) | 0.73 |
| weno5_js + mp_bound | 1.489 (0.001) | 0.129 (0.002) | 11.54 | 38.8 (23.4) | 54.2 (5.5) | 0.72 |
| linear5 | 0.206 (0.001) | 0.105 (0.003) | 1.96 | 11.4 (4.5) | 12.6 (1.3) | 0.90 |
| linear7 | 0.207 (0.001) | 0.105 (0.001) | 1.97 | 7.1 (4.3) | 6.1 (1.3) | 1.17 |
| sensor_lcr | 0.384 (0.003) | 0.112 (0.002) | 3.43 | 23.1 (10.5) | 25.6 (1.7) | 0.91 |
| hybrid lcr + linear5 / wenoz5 | 0.731 (0.001) | 0.130 (0.002) | 5.62 | 26.6 (5.9) | 30.2 (1.5) | 0.88 |
| keep_pe6 | 0.849 (0.004) | 0.293 (0.003) | 2.90 | 44.8 (1.5) | 45.3 (1.8) | 0.99 |
| **median S** | | | **5.17** | | | **0.99** |

## Decisions

- **RTX 4070: S = 5.17 → R4P recommended** (within the R4P limits of `docs/guide/precision.md`). The linear kernels
  gain the factor 2 of the halved traffic (they move data); the weighted ones 3 to 14 (FP64 at 1/64 of FP32 on this
  card). A data-centre GPU (FP64 at half rate) was not measured: about 2 at most is the expectation, from the traffic
  and the rate, not a measurement.
- **One CPU core: S = 0.99 → not worth it.** gfortran `-O3` does not vectorise the face loops (`-fopt-info-vec`:
  only the stencil copies are packed; each kernel is a call into the separately compiled library), so the narrower
  kind buys nothing per operation. The run-to-run spread is large for a few kernels (the box was shared), but the
  median is far from 1.2.

## The WENO-JS loss on the CPU

The WENO-JS family is slower in R4P on one core (S 0.67–0.85). Test: an R4P build with the R8P weight form (the
`#ifdef _R_P_IS_R4P` of `seth_weno_js_body.INC` turned off, in a scratch copy only), against the R4P build, three
interleaved runs on 2²² faces, core 2, 14:36–14:41: the R4P form costs WENO-JS 15–18 % (weno5_js 44.6 against
37.9 ns per face; descaling, descaling_ai, mp_bound 0.89–0.90), every other kernel within 0.98–1.05. The rest of the
loss is not attributed (the R8P runs were on 2²⁴ faces at another time: no split across runs is claimed). On the GPU
the form's cost was not isolated; weno5_js gains 9.4 with it.
