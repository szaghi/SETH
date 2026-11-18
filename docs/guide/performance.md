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
| `seth_weno4_c` | 0.740 | 12.5 | **1.56** | WENO3-JS (0.475, 7.8) | 1.4 |
| `seth_weno6_c` | 1.271 | 27.3 | 1.28 | WENO5-JS | 1.3 |
| `seth_weno8_c` | 2.085 | 39.4 | 1.16 | WENO7-JS | 1.25 |

Bold: over the budget set before coding, a finding rather than an edit of the budget. On one CPU core the picture
differs: the adaptive cutoff of TENO5-A costs nothing measurable (80.0 against 79.5) while it costs 41% on the GPU,
whose FP64 rate (1/64 of FP32 on this card) makes every division count; the data-centre GPU was not measured. The
old kernels did not slow down (`bench_compare.py`: ok on both machines).

The centred WENO (2026-11-16, four GPU runs and three CPU runs, `notes/m4-bench-data/30-weno-c-*`): its blend adds a
near-constant 0.31–0.35 ns per face on the GPU (a central dot product and two FP64 divisions per state), so the
budget set as a ratio is missed by the two cheapest kernels, WENO4-C most; on one core the ratios are 1.71 (WENO4-C),
1.13 and 1.09 with gfortran, 1.09, 1.12 and 1.04 with nvfortran. The old kernels did not slow down
(`bench_compare.py`: ok on both machines).

The same evening the blend went to one division per state in R8P (the larger ratio chosen by cross products,
[Schemes](/guide/schemes#centred-weno-weno4-6-8-c)); the table shows it. Measured against the two-division build in
the same session, interleaved (four GPU rounds, three per compiler on one core, `notes/m4-bench-data/31-weno-c-1div-*`):
on the GPU it saves 0.045 ns per face (WENO4-C, WENO6-C) and 0.065 (WENO8-C), and WENO6-C's ratio falls from 1.33 to
1.28, within its budget; WENO4-C's falls from 1.66 to 1.56, still over. The divisions were about 30% of the blend's
cost, not most of it, as the change had assumed: the rest is its other FP64 operations (the central sum, the
indicators), presumably counted at this card's FP64 rate like any other (not measured separately). On one core the
gain is within the spread of the runs: gfortran 1.60, 1.14, 1.08 (WENO4-C still over), nvfortran 1.08, 1.04, 1.04.
WENO4-C's miss is accepted (2026-11-17): the budget stays as frozen and recorded as missed, and the kernel as it is.

## The non-polynomial composites (issue #6)

ns per face against `seth_teno5` on the same device and in the same run, the median of the runs of
`seth_bench_weno` (raw outputs in `notes/m5-bench-data/42-*` and `42b-*`, 2026-11-18): the GPU column on the RTX 4070
(`bench-nvf-oac`, $2^{24}$ faces, 20 passes, series 42b: four runs on a quiet device, every ratio within 1 %; series
42's three GPU runs, taken next to other jobs, spread up to 8 % and agree within it), the CPU columns on one pinned core of the i9-13900KF ($2^{22}$ faces, 11
passes, three runs; `bench-gnu`, `bench-nvf`). The field is the benchmark's smooth wave with a jump every 2048 cells:
mostly smooth faces, where the early exit acts. The budgets are issue #6's. The pieces are timed alone: the
detectors on stored masks, the overrides forced on every face (in a scheme they act on a few faces only), the mask on
one bias, the unit of its budget (two biases cost twice).

| Kernel | GPU | gfortran | nvfortran | budget |
|---|---|---|---|---|
| `seth_teno5` (ns per face) | 1.50 | 84.0 | 38.9 | — |
| `seth_teno5_mask`, one bias | 0.41 | 0.44 | 0.28 | 0.70 |
| `seth_dr_detect` / `seth_thinc_detect` | 0.05 / 0.05 | 0.11 / 0.11 | 0.12 / 0.11 | 0.10 |
| `seth_dr_override`, every face | **0.19** | 0.05 | **0.32** | 0.15 |
| `seth_thinc_override`, every face | **0.48** | **0.32** | **0.73** | 0.30 |
| `seth_teno5_st` | 0.88 | 0.87 | 0.98 | 1.05 |
| `seth_teno5_rbf` ($s_1 = s_2 = 0.09$) | **1.27** | 0.91 | 1.06 | 1.25 |
| `seth_teno5_dr` | 2.74 | 2.85 | **3.05** | 3.0 |
| `seth_teno5_thinc` | **3.49** | **3.72** | **3.49** | 3.0 |
| `seth_teno5_rbf_dr` | 2.97 | 2.76 | **3.14** | 3.0 |
| DR / THINC / RBF+DR, no early exit | 5.96 / 6.77 / 6.16 | 6.29 / 7.14 / 6.20 | 4.44 / 5.01 / 4.67 | — |
| DR / THINC / RBF+DR, stored masks (two passes) | 1.82 / 1.83 / 2.09 | 1.95 / 1.95 / 1.86 | 1.64 / 1.64 / 1.70 | — |

![The non-polynomial kernels against their budgets, and the composites against the density of jumps](/guide/perf-nonpoly.png)

Experiment E1 of the plan: the self-contained kernels recompute the masks of the neighbouring cells at every face,
in that face's characteristic basis (the default, owner decision D3). TENO5-DR fits the budget of 3 times TENO5 on the
GPU and with gfortran and misses it narrowly with nvfortran; TENO5-THINC misses it everywhere, because its detector
needs three downwind masks per state where DR's needs two; the complete scheme of Liang & Fu, RBF+DR, fits it on the
GPU and with gfortran and misses it with nvfortran. The stored-mask variant (each cell's two masks once, then the faces), which a caller builds from the exported
pieces, costs 1.6–2.1 times TENO5 for any detector; on scalar data it is equal bit for bit, on systems it computes
each mask in its own face's basis and the results differ ([NP4](/verification/nonpolynomial/dr#results)). The pieces
fit their budgets except the overrides forced on every face, the cost a scheme pays only where a detector fires; the
Gaussian-RBF kernel misses its budget on the GPU only (the sensor's divisions). Published whole-solver costs, for
orientation only: TENO5-THINC 1.17 times TENO5 per time step (Takagi et al. 2022, Table 3), the complete RBF+DR scheme
1.38–1.42 (Liang & Fu 2024, Table 1).

**The cost against the density of jumps** (PERF7; the same programs, the field's jumps every $\infty$, 2048, 64, 8
cells, `42-fields-*`): on one CPU core the per-face kernels go from 2.9 (DR) and 3.8 (THINC) times TENO5 on a smooth
field to 3.8 and 4.5 with a jump every 8 cells. On the GPU they go from 2.7 and 3.4 to 5.3 and 5.8: above the kernels
without the early exit (5.2), since neighbouring faces then take different paths through the detector and the
threads of a warp wait for each other. The stored-mask variant stays at 1.7–2.0 everywhere. A solver whose faces are
mostly smooth gains from the early exit; one with many discontinuities, on the GPU, does better storing the masks.

The old kernels did not slow down. `bench_compare.py` is ok on the gfortran core (series 42). On the GPU every old
variant, ADAM's copies included, runs 4–5 % slower than the 2026-11-15 baseline, and so does the benchmark of the
unchanged base commit 538a471 run in the same hour (`42c-*`, 4.4–4.7 %): the device is slower than three days
before, the kernels are not. The nvfortran core had shown `keep_pe6` and ADAM's `compute_kep_flux` at 71 ns per face in series 41
against 39 before: the bench program of that commit (fa76aa0) built with the same library runs them at 73.7, the one
before at 41.2, tonight's at 40.7. A code-generation effect of the benchmark program, not of `seth_kep`, which did not
change.

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
