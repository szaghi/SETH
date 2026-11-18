---
title: 11. On the GPU
---

# 11. On the GPU

Every kernel is `pure`, `!$acc routine seq` and `!$omp declare target`, and needs no table, object or copy on the
device: its coefficients are compile-time constants. So the solver of chapter 3 runs on a GPU by putting its loops
there, and calling the kernel inside the face loop, one thread per face.

## The data on the device

One data region for the whole run: the fields are copied in once and out once, and every loop of the time step works
on the device copies.

<<< @/examples/snippets/t11_gpu-data.F90{fortran}

## The face loop

<<< @/examples/snippets/t11_gpu-residual.F90{fortran}

- Each thread copies its face's six values into `vp` and `vm`, private arrays of **constant bounds**. A local array
  whose bounds depend on a dummy argument would be allocated from the device heap at every call, hundreds of times
  slower ([Performance](/guide/performance#device-kernel-rules-arrays-and-loop-extents)).
- The kernel call is the same as on the host.

## Building it for the GPU

The device code of the kernels must be compiled by the same compiler, with OpenACC: declare SETH and PENF as FoBiS
dependencies of your project ([Installation](/guide/installation#as-a-fobis-dependency)), so that their sources are
compiled with your flags, `-acc -gpu=cc89` for nvfortran on an RTX 4070.

## Running it

Built with gfortran, the directives are comments and the program runs on the host; that is the run the documentation
build checks:

<<< @/examples/output/t11_gpu.txt{text}

Built with nvfortran and OpenACC, on the GPU, captured once and kept in `docs/examples/device/`:

<<< @/examples/device/t11_gpu.txt{text}

The same numbers to the printed digits. The two builds need not agree bitwise: nvfortran contracts multiplications
and additions into fused multiply-adds and gfortran does not; how far apart they may be, kernel by kernel, is
[C5](/verification/contracts/device).

::: tip What you learned
To run on the GPU, offload your loops and keep your fields on the device; the kernels need nothing. Reference:
[Performance](/guide/performance).
:::

<details><summary>The whole program</summary>

<<< @/examples/snippets/t11_gpu.F90{fortran}

</details>
