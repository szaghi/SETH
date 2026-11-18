# About

SETH (**S**moothness-weighted **E**NO, **T**ENO and **H**ybrids) reconstructs and interpolates with the weighted
essentially non-oscillatory family of schemes: WENO, TENO and the hybrids built on them. Each scheme is one `pure`
kernel that computes the two states at one face; you call it inside your own loops, on the host or on a GPU, with no
setup and nothing to copy to the device.

SETH succeeds WenOOF, whose run-time polymorphism in the inner loops cost too much on HPC hardware. In the Genesis
story Seth is the son of Adam who carries the line on after Abel: here WenOOF is Abel.

## What it has

- **WENO-JS** 1/3/5/7/9, **WENO-Z** 3/5/7/9, **centred WENO** 4/6/8, and **TENO** 5 to 8 with the optimised, the
  AZ and the adaptive-cutoff variants: reconstruction from cell averages and interpolation from point values.
- **Non-polynomial TENO**: the DR and THINC overrides, stencil targeting and the Gaussian-RBF scheme.
- **The pieces of a hybrid**: the linear schemes, four smoothness sensors, the KEEP-PE central fluxes and the Ducros
  shock sensors.
- **Limiters**: the monotonicity-preserving bound, and the positivity-preserving flux limiter for the Euler equations.
- **Double or single precision**, chosen at build time; gfortran, nvfortran and ifx; CPU and GPU (OpenACC, verified on
  the device; the kernels also carry `!$omp declare target`).

Every public name, family by family: the [feature map](/reference/).

## How to read this documentation

| If you want to | Read |
|---|---|
| build it and see one kernel work | [Installation](./installation), then the [Quick start](./quickstart) |
| learn to use the kernels, from one face to the Euler equations on a GPU | the [Tutorial](/tutorial/), in order |
| solve one task | the [Cookbook](/cookbook/) |
| look up a kernel, its arguments, its constants | the [Reference](/reference/), and formal's [API](/api/) |
| understand a scheme | [Schemes](./schemes) |
| know how far to trust it | [Verification](/verification/): every claim measured, every case a program |

The pages follow one another in that order: the *Next* link at the foot of each page reads the documentation
straight through.
