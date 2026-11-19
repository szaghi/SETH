---
title: 11. On the GPU
---

# 11. On the GPU

Every kernel is `pure`, `!$acc routine seq` and `!$omp declare target`, and needs no table, object or copy on the
device: its coefficients are compile-time constants. So the solver of chapter 3 runs on a GPU by putting its loops
there, and calling the kernel inside the face loop, one thread per face.

## The data on the device

The field lives in device memory for the whole run, copied in once and out once. The time step is
[SELAH](https://szaghi.github.io/SELAH/)'s, as in chapters 6 to 10, and on the GPU SELAH works on device memory: the
field is a FUNDAL device array (`dev_alloc`: on the GPU with OpenACC, on the host without), and the stepper's
registers are allocated there too. The work arrays of the residual live in a data region, so a loop names the device
arrays it reads with `deviceptr` (the field, SELAH's register) and the others with `present`.

<<< @/examples/snippets/t11_gpu-data.F90{fortran}

## The face loop

<<< @/examples/snippets/t11_gpu-residual.F90{fortran}

- Each thread copies its face's six values into `vp` and `vm`, private arrays of **constant bounds**. A local array
  whose bounds depend on a dummy argument would be allocated from the device heap at every call, hundreds of times
  slower ([Performance](/guide/performance#device-kernel-rules-arrays-and-loop-extents)).
- The kernel call is the same as on the host.

## Building it for the GPU

The device code of the kernels must be compiled by the same compiler, with OpenACC: declare SETH, PENF, FUNDAL and
SELAH as FoBiS dependencies of your project ([Installation](/guide/installation#as-a-fobis-dependency)), so that their
sources are compiled with your flags: `-acc -gpu=cc89` for nvfortran on an RTX 4070, and SELAH's device build,
`-DDEV_OAC -DCOMPILER_NVF -DFUNDAL_NO_MPI` with FUNDAL's `src/lib` among the include directories. SETH's
`fobis build --mode example-t11-nvf-oac` builds this chapter that way.

## Running it

Built with gfortran, the directives are comments and the program runs on the host; that is the run the documentation
build checks:

<<< @/examples/output/t11_gpu.txt{text}

Built with nvfortran and OpenACC (`example-t11-nvf-oac`), on the GPU, captured and kept in `docs/examples/device/`:

<<< @/examples/device/t11_gpu.txt{text}

The same numbers to the printed digits. The two builds need not agree bitwise: nvfortran contracts multiplications
and additions into fused multiply-adds and gfortran does not; how far apart they may be, kernel by kernel, is
[C5](/verification/contracts/device).

::: tip What you learned
To run on the GPU, offload your loops and keep your fields on the device; the kernels need nothing, and SELAH's
stepper takes the device field as it is. Reference:
[Performance](/guide/performance).
:::

<details><summary>The whole program</summary>

<<< @/examples/snippets/t11_gpu.F90{fortran}

</details>
