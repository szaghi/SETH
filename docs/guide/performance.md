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

Against ADAM's host routines (FLUME's CPU path) SETH is about 2x faster on one core. The results are bitwise equal to
ADAM's on the device.

**Why.** The coefficients are compile-time constants: no table loads, no array descriptors, constant trip counts the
compiler unrolls, nothing to keep on the device.

**Caveat.** The RTX 4070 runs FP64 at 1/64 of its FP32 rate, which favours whichever kernel issues fewer instructions;
the production margin is the one of a data-centre GPU (A30, H100). Details and raw outputs: `notes/m3-gpu-benchmark.md`.
