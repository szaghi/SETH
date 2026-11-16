# Performance

`src/bench/seth_bench_weno.F90` compares SETH's kernels with ADAM's on the same stencils and memory traffic
(2²⁴ faces, `!$acc wait` + `system_clock`, median of 20 passes): `fobis build --mode bench-nvf-oac`, then run with
`CUDA_VISIBLE_DEVICES` set.

| ns per face | SETH | ADAM per-width kernel | ADAM generic kernel |
|---|---|---|---|
| WENO-5, RTX 4070 | **0.97** | 1.27 (SETH 1.31x faster) | 1.38 (1.42x) |
| WENO-7, RTX 4070 | **1.76** | 2.22 (1.26x) | 2.36 (1.34x) |
| WENO-5, one i9-13900KF core, gfortran | 23.1 | 21.3 (SETH 0.92x) | 49.8 |
| WENO-7, one core, gfortran | 34.6 | 38.2 (1.10x) | 78.3 |

![ns per face of SETH's kernels and ADAM's, on the RTX 4070 and on one CPU core with gfortran and nvfortran](/guide/perf-kernels.png)

*GPU: the median of four runs of the program, the error bar their spread (the wide bar of SETH's WENO-5 is the first
kernel the program launches, 6 % slower in run 1 only). WENO-Z5 and TENO6 have no ADAM counterpart and are shown for
the cost of their weights: WENO-Z5 costs what WENO5-JS costs, TENO6 (four candidates, a six-point global
indicator) 1.3 times what WENO7-JS costs. CPU: one run per compiler.*

Against ADAM's host routines (FLUME's CPU path) SETH is about 2x faster on one core. The results are bitwise equal to
ADAM's on the device.

**Why.** The coefficients are compile-time constants: no table loads, no array descriptors, constant trip counts the
compiler unrolls, nothing to keep on the device.

**Caveat.** The RTX 4070 runs FP64 at 1/64 of its FP32 rate, which favours whichever kernel issues fewer instructions;
the production margin is the one of a data-centre GPU (A30, H100). Details and raw outputs: `notes/m3-gpu-benchmark.md`.

## Regressions

Every kernel added after v0.1.0 is benchmarked against a recorded baseline before it is merged:
`notes/bench-baseline/` holds three raw runs per build (the RTX 4070 and one core of the i9-13900KF), and
`scripts/bench_compare.py` fails when a variant's median exceeds the baseline's by more than 3 % plus the spread of
the baseline's runs. On the GPU that spread is below 1 %; on one WSL2 CPU core it is 5 to 20 %, so the CPU check only
catches large regressions.

```bash
python3 scripts/bench_compare.py notes/bench-baseline/rtx4070-nvf-oac run-1.txt run-2.txt run-3.txt
```

## The kernels added after v0.1.0

ns per face, the median of three runs of `seth_bench_weno` (raw outputs in `notes/m4-bench-data/`): the GPU column on
the RTX 4070 (`bench-nvf-oac`), the CPU column on one pinned core of the i9-13900KF (`bench-gnu`). The ratio is to the
kernel the cost budget of the roadmap (issue #4) was set against, on the GPU.

| Kernel | GPU | CPU | GPU ratio | against | budget |
|---|---|---|---|---|---|
| `seth_teno5` | 1.435 | 79.5 | — | — | — |
| `seth_teno5_opt` | 1.435 | 79.1 | 1.00 | TENO5 | 1.02 |
| `seth_teno5_az` | 1.461 | 80.8 | 1.02 | TENO5 | 1.05 |
| `seth_teno5_a` | 2.025 | 80.0 | **1.41** | TENO5 | 1.25 |
| `seth_teno8` | 3.977 | 131.1 | **1.66** | TENO6 | 1.6 |
| `seth_teno8_a` | 4.624 | 136.7 | 1.16 | TENO8 | 1.25 |
| `seth_weno5_js` + `seth_weno_descaling_ai` | 1.336 | 23.6 | (descaler) **~3** | the mean-\|v\| descaler | 1.15 |
| `seth_weno5_js` + `seth_mp_bound` | 1.413 | 29.8 | (bound) **0.43** | WENO5-JS | 0.35 |
| `seth_linear5` | 0.095 | 2.3 | 0.10 | WENO5-JS | 0.25 |
| `seth_linear7` | 0.120 | 2.7 | 0.07 | WENO7-JS | 0.25 |
| `seth_sensor_lcr` | 0.298 | 15.9 | — | — | — |
| sensor + `seth_linear5` | 0.393 | 18.2 | 0.40 | WENO-Z5 | 0.6 |
| `seth_keep_pe6` | 0.777 | 31.8 | 1.00 | ADAM's `compute_kep_flux` | 1.0 |
| `seth_weno4_c` | 0.783 | 13.9 | **1.66** | WENO3-JS (0.473, 8.1) | 1.4 |
| `seth_weno6_c` | 1.316 | 27.5 | **1.33** | WENO5-JS | 1.3 |
| `seth_weno8_c` | 2.151 | 39.7 | 1.20 | WENO7-JS | 1.25 |

Bold: over the budget set before coding, a finding rather than an edit of the budget. On one CPU core the picture
differs: the adaptive cutoff of TENO5-A costs nothing measurable (80.0 against 79.5) while it costs 41% on the GPU,
whose FP64 rate (1/64 of FP32 on this card) makes every division count; the data-centre GPU was not measured. The
old kernels did not slow down (`bench_compare.py`: ok on both machines).

The centred WENO (2026-11-16, four GPU runs and three CPU runs, `notes/m4-bench-data/30-weno-c-*`): its blend adds a
near-constant 0.31–0.35 ns per face on the GPU (a central dot product and two FP64 divisions per state), so the
budget set as a ratio is missed by the two cheapest kernels, WENO4-C most; on one core the ratios are 1.71 (WENO4-C),
1.13 and 1.09 with gfortran, 1.09, 1.12 and 1.04 with nvfortran. The old kernels did not slow down
(`bench_compare.py`: ok on both machines).

## Single precision

`src/bench/seth_bench_kinds.F90` times the same kernels in the working kind, built in R8P and in R4P
(`bench-kinds-nvf-oac` against `bench-nvf-oac-r4`, `bench-kinds-gnu` against `bench-gnu-r4`). On the RTX 4070 R4P is
5.17 times faster (the median over 18 kernels: 2 for the linear ones, 3 to 14 for the weighted ones); on one CPU core
0.99 times, nothing gained. What R4P changes, its limits and the per-kernel table: [Precision](/guide/precision).

## Device-kernel rules: arrays and loop extents

Two rules of ADAM's device kernels, measured on SETH by `src/bench/seth_bench_device_arrays.F90`
(`fobis build --mode bench-arrays-nvf-oac`; RTX 4070, ns per face):

| Rule | Measured on | Follows it | Breaks it | Cost |
|---|---|---|---|---|
| No automatic local arrays (bounds from a dummy argument) in a device routine | WENO-JS 5 | 0.97 | 553 | 569× |
| | positivity limiter, 1-D | 2.70 | 899 | 333× |
| Loop extents known at compile time | WENO-JS 5 (over-dimensioned arrays, run-time `S`) | 0.97 | 1.04 | 7 % |
| | positivity limiter, loops to `nv` | 2.67 | 2.70 | 1 % |

![ns per face, log scale, of WENO5-JS and of the positivity limiter following or breaking the two rules, on the GPU and on one CPU core](/guide/perf-device-rules.png)

*Blue: follows both rules (the library); light blue: arrays of constant bounds but run-time loop extents; umber:
automatic arrays. On the GPU automatic arrays cost two to three orders of magnitude and run-time extents 7 %; on a
CPU core both cost a factor of 1.3 to 1.8, together up to 4. The limiter in 1-D (`nv` = 3) and 3-D (`nv` = 5); the
library's loops to the dummy `nv` cost nothing on the GPU against `nv` constant.*

An automatic array in a device routine is allocated from the device heap at every call. Besides the time, a launch can
fail with the default heap (`DEVICE FORTRAN AUTO ALLOCATION FAILED`, `NVCOMPILER_ACC_CUDA_HEAPSIZE`).

SETH's WENO-JS, WENO-Z and TENO kernels follow both rules by construction: every bound and extent is a named constant
of the kernel. The first version of `seth_positivity_euler` had automatic arrays of size `nv`; its locals are now
dimensioned `SETH_POSITIVITY_NV_MAX = 5` and it runs 333× faster on the device with the same results. The audit of
every device routine, routine by routine, and the raw outputs: `notes/gpu-kernel-audit.md`.

::: tip Writing a kernel of your own
Declare no local array whose bound depends on a dummy argument, and pass no array expression of such a size as an
argument. Over-dimension with a named constant and loop to the dummy.
:::

---

*The figures are drawn by `docs/scripts/make_performance_figures.py` from the raw outputs of the two benchmarks
(`notes/m3-bench-data/`, `notes/gpu-audit-data/`), light and dark.*
