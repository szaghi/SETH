---
title: Time a kernel on a GPU
---

# Time a kernel on a GPU

A reference: how SETH's benchmarks time the kernels, so that your numbers compare with
[Performance](/guide/performance).

- **Enough work.** $2^{24}$ faces per launch on the GPU: a launch of a few thousand faces measures the launch, not the
  kernel.
- **Synchronise before reading the clock.** OpenACC launches are asynchronous with `async`, and a timer read before
  the device finished times nothing: `!$acc wait`, then `system_clock`.
- **Warm up, then take the median.** The first launch pays for the module load and the memory pools; SETH's benchmark
  reports the median of 20 passes, and repeats whole runs to see their spread.
- **Data on the device.** Copy the fields once, outside the timed loop (`!$acc data`), as in
  [tutorial 11](/tutorial/11-gpu); a copy per launch measures PCIe.
- **Pick the GPU.** `CUDA_VISIBLE_DEVICES` set, so the run is on the card you mean.
- **Know the card's FP64 rate.** The RTX 4070 runs FP64 at 1/64 of its FP32 rate: a division costs far more there
  than on a data-centre GPU, and the ratios between kernels change.
- **No automatic arrays in device code.** A local array sized by a dummy argument is allocated from the device heap
  at every call: two to three orders of magnitude slower
  ([measured](/guide/performance#device-kernel-rules-arrays-and-loop-extents)).

The benchmark itself: `src/bench/seth_bench_weno.F90`, built with `fobis build --mode bench-nvf-oac` (GPU) or
`bench-gnu` (one CPU core); `scripts/bench_compare.py` compares two runs against regressions.
