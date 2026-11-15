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

## Device-kernel rules: arrays and loop extents

Two rules of ADAM's device kernels, measured on SETH by `src/bench/seth_bench_device_arrays.F90`
(`fobis build --mode bench-arrays-nvf-oac`; RTX 4070, ns per face):

| Rule | Measured on | Follows it | Breaks it | Cost |
|---|---|---|---|---|
| No automatic local arrays (bounds from a dummy argument) in a device routine | WENO-JS 5 | 0.97 | 553 | 569× |
| | positivity limiter, 1-D | 2.70 | 899 | 333× |
| Loop extents known at compile time | WENO-JS 5 (over-dimensioned arrays, run-time `S`) | 0.97 | 1.04 | 7 % |
| | positivity limiter, loops to `nv` | 2.67 | 2.70 | 1 % |

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
