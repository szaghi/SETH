# M3 phase A: SETH's kernels against ADAM's, on the GPU and the CPU

Night 2026-11-14 (owner's go-ahead for the WSL GPU at 22:10). Program: `src/bench/seth_bench_weno.F90`, modes
`bench-nvf-oac`, `bench-nvf`, `bench-gnu`. Raw outputs: `notes/m3-bench-data/`.

## What is compared

2^24 faces of a 1-D field (a smooth wave with a jump every 2048 cells, two biased fields), every variant reading the
same stencils and writing the two states of every face:

| Variant | Kernel |
|---|---|
| SETH weno5_js / weno7_js | SETH: coefficients compile-time constants, explicit-shape stencils |
| ADAM generic S=3 / S=4 | ADAM `weno_reconstruct_upwind_dev`: run-time S, tables as assumed-shape device arrays |
| ADAM width s3 / s4 | ADAM `weno_reconstruct_upwind_s<N>_dev`: S a constant, tables still arrays (FLUME's FNL fast path) |
| ADAM host S=3 / S=4 | ADAM's host `weno_reconstruct_upwind` (FLUME's CPU path; CPU builds only) |
| SETH wenoz5, teno6 | for reference (more work per face than WENO-JS) |

ADAM's kernels are byte-for-byte ADAM's (`src/bench/adam_fnl_*`); ADAM's variants pack the stencils as FLUME does.
Timing: 2 warm-up passes, then 20 (GPU) or 10 (CPU) timed passes, `!$acc wait` + `system_clock`, median reported.
**Correctness, the first device execution of SETH's kernels:** on GPU 1 the device results of SETH equal the host
SETH results, and ADAM's FNL kernels equal SETH's, **bitwise**, in all three runs.

## GPU: RTX 4070 (GPU 1, idle), nvfortran 26.9 `-O4 -fast -acc -gpu=cc89`

ns per face, median of 20 passes; three runs of the whole program (run 1 / 2 / 3):

| Variant | run 1 | run 2 | run 3 | SETH speed-up |
|---|---|---|---|---|
| **SETH weno5_js** | 1.028 | 0.966 | 0.966 | |
| ADAM generic S=3 | 1.371 | 1.376 | 1.376 | **1.42x** |
| ADAM width s3 | 1.264 | 1.268 | 1.268 | **1.31x** |
| **SETH weno7_js** | 1.755 | 1.761 | 1.761 | |
| ADAM generic S=4 | 2.353 | 2.353 | 2.357 | **1.34x** |
| ADAM width s4 | 2.214 | 2.216 | 2.218 | **1.26x** |
| SETH wenoz5 (p=1) | 1.018 | 1.019 | 1.018 | |
| SETH teno6 (C_T=1e-7) | 2.365 | 2.366 | 2.366 | |

Run-to-run spread below 0.5 %, except SETH weno5_js in run 1 (6 %, the first kernel the program launches).

Run 4 (`bench-gpu1-run4-sumsquares.txt`), after the WENO-Z/TENO indicators became sums of squares (night step 10):
WENO-JS and ADAM unchanged (0.967, 1.269, 1.762, 2.215 ns), SETH wenoz5 0.961 ns (was 1.018), teno6 2.330 ns (was
2.366).

## CPU: one core of an i9-13900KF, ns per face, median of 10 passes

| Variant | gfortran 16.2 -O3 | nvfortran 26.9 -O4 -fast |
|---|---|---|
| **SETH weno5_js** | 23.1 | 28.2 |
| ADAM host S=3 (FLUME CPU) | 47.1 (SETH 2.04x) | 47.8 (SETH 1.69x) |
| ADAM width s3 | 21.3 (**SETH 0.92x**) | 29.8 (SETH 1.06x) |
| ADAM generic S=3 | 49.8 | 51.1 |
| **SETH weno7_js** | 34.6 | 38.1 |
| ADAM host S=4 (FLUME CPU) | 71.5 (SETH 2.06x) | 86.7 (SETH 2.28x) |
| ADAM width s4 | 38.2 (SETH 1.10x) | 39.4 (SETH 1.04x) |
| ADAM generic S=4 | 78.3 | 85.4 |
| SETH wenoz5 / teno6 | 48.2 / 75.2 | 25.6 / 71.9 |

## Reading

- **On the GPU SETH wins every comparison**, 1.26-1.31x against ADAM's fastest kernels (the per-width ones FLUME uses)
  and 1.34-1.42x against the generic ones. The design bet of the kernel API (compile-time coefficients: no table
  loads, no descriptors, constant trip counts) is confirmed on this device.
- **On the CPU SETH is about 2x faster than what FLUME's CPU backend calls today** (ADAM's host routines, run-time S),
  and level with ADAM's per-width kernels, with one exception: gfortran runs ADAM's WENO-5 width kernel 8 % faster
  than SETH's. Not investigated tonight; SETH copies the stencils into ADAM's row layout before the arithmetic (kept
  for the bitwise agreement with ADAM), a copy the width kernel does not make.
- WENO-Z5 costs what WENO-JS costs on the GPU (1.02 vs 0.97 ns); TENO6 costs about what ADAM's generic WENO-7 does
  (its 6-point global indicator adds a 21-term quadratic form).

## Caveat: a consumer GPU

The RTX 4070 runs FP64 at 1/64 of its FP32 rate. These kernels are FP64 arithmetic with little memory traffic per
flop, so on this card they are bound by FP64 throughput, which favours whichever kernel issues fewer instructions.
On a data-centre GPU (A30, FP64 at 1/2 of FP32: medooza) the balance shifts toward memory and latency, and the margin
can shrink or grow. The ranking is unlikely to invert, because SETH does strictly less work per face (no table loads,
no index arithmetic), but the number to quote for production is the A30 one.

## Verdict for M3 phase B (FLUME calls SETH)

Proceed. SETH is at least as fast as ADAM's best kernel on every platform measured except one CPU case (-8 %), bitwise
equal to ADAM's results, and removes the table mirroring. Before porting: one run of this benchmark on a medooza A30
(left to the owner: tonight's permission covered the WSL GPU only).

Outcome (2026-11-22): ADAM ported without that run (ADAM issue 139) and measured the kernels inside FLUME on an A30
instead (ADAM commit 4210fc8f): the WENO-9 face kernel 13.6 % faster than ADAM's former per-width kernel, WENO-5
1.7 % slower (ADAM issue 141). This benchmark on an A30 is item E0 of #7, below.

## GPU: A30 (item E0 of #7)

medooza `gnode23`, one A30 (driver 545.23, max SM clock 1440 MHz), nvfortran 26.9 `-O4 -fast -acc -gpu=cc80,cuda12.9`
(`bench-nvf-oac --varset iac`), the default problem (2^24 faces, 20 timed passes, median), three runs of the whole
program in each Slurm job. Job 1369294 ran this program as on main (raw: `bench-a30-run{1,2,3}.txt`); job 1369295 ran
it with variants 47-58 added: every width, with the `_wr` and `_interp_wr` variants (raw:
`bench-a30-widths-run{1,2,3}.txt`; GPU state `bench-a30-state.txt`). Every check passes in all six runs. SETH on the
device equals SETH on the host bitwise for every WENO-JS variant and its weights ratios, and ADAM's kernels equal
SETH's bitwise.

ns per face, the median of runs 2 and 3 of job 1369294. Run 1 starts with SETH weno5_js at 0.226 ns (0.235 in job
1369295), as in the RTX 4070's first run: it is the first kernel the program launches.

| Variant | A30 | SETH speed-up | RTX 4070 (run 4) |
|---|---|---|---|
| **SETH weno5_js** | 0.184 | | 0.967 |
| ADAM generic S=3 | 0.79 | **4.3x** | 1.376 |
| ADAM width s3 | 0.496 | **2.70x** | 1.269 |
| **SETH weno7_js** | 0.291 | | 1.762 |
| ADAM generic S=4 | 1.13 | **3.9x** | 2.353 |
| ADAM width s4 | 0.629 | **2.16x** | 2.215 |
| SETH wenoz5 (p=1) | 0.472 | | 0.961 |
| SETH teno5 / teno6 | 0.501 / 0.996 | | — / 2.330 |
| SETH keep_pe6 / ADAM compute_kep_flux | 0.575 / 0.574 | 1.00x | |

Every width, job 1369295, runs 2 and 3:

| Width | `_js` | `_js_wr` | `_js_interp_wr` |
|---|---|---|---|
| WENO-1 | 0.102 | 0.135 | 0.135 |
| WENO-3 | 0.136 | 0.167 | 0.167 |
| WENO-5 | 0.180 | 0.22 | 0.22 |
| WENO-7 | 0.291 | 0.29 | 0.29 |
| WENO-9 | 0.408 | 0.47 | 0.47 |

- **The A30 margin is larger than the RTX 4070's**: 2.2-2.7x against ADAM's per-width kernels, against 1.26-1.31x.
  This is the direction the caveat above expected: with FP64 at half the FP32 rate, the table loads and index
  arithmetic SETH avoids weigh more.
- **Inside FLUME the margin does not carry over** (ADAM issue 141, closed): there the face kernel adds the
  characteristic projection, the per-face local arrays and the flux assembly around the reconstruction. FLUME's
  WENO-9 face kernel is 13.9 % faster and its WENO-5 kernel 1.6 % slower. ADAM measured that the WENO-5 cost is
  register allocation of the whole kernel, at the same occupancy, not the reconstruction itself.
- The weights ratio (`_wr`) costs 0.03-0.06 ns per face, except at WENO-7, where it costs nothing measurable. The
  interpolation tables cost what the reconstruction ones do.
- On the A30, WENO-Z5 costs 2.6x WENO-JS5 (0.472 vs 0.184 ns), whereas on the RTX 4070 the two cost the same. Its
  global indicator and the power with an integer exponent are not hidden behind FP64 throughput here. Not
  investigated.
