---
title: Quick start
---

# Quick start

One program, compiled against the library, that reconstructs the two states at every face of a step with WENO5-JS
and with the linear scheme WENO5 falls back to on smooth data.

## Build the library

```bash
git clone https://github.com/szaghi/SETH && cd SETH
fobis fetch                       # PENF, into src/third_party/
fobis build --mode static-gnu     # lib/gnu/static/libseth.a and its modules in lib/gnu/static/mod
```

Other compilers and the single-precision build: [Installation](./installation).

## The program

`use seth` gives every kernel, its constants and the kinds (`R_P` is the working real kind):

<<< @/examples/snippets/quickstart-use.F90{fortran}

A kernel computes **one face** from the cells around it, indexed by their offset from the face: cell 0 is left of
it, cell 1 right of it. WENO5 reads six cells, `-2:3`. The same array goes in twice, once for the state left of the
face (`rp`) and once for the state right of it (`rm`):

<<< @/examples/snippets/quickstart-face.F90{fortran}

- `SETH_WENO_JS_EPS` is the regulariser of the weights; `1._R_P` is `rmu`, the descaler of the smoothness
  indicators (1: Jiang and Shu's weights).
- There is no setup and no object: the coefficients are compile-time constants, and every kernel is `pure`, so the
  call goes inside your own loops, on the host or on a GPU.

<details><summary>The whole program</summary>

<<< @/examples/snippets/quickstart.F90{fortran}

</details>

## Compile and run

```bash
gfortran -I lib/gnu/static/mod quickstart.F90 lib/gnu/static/libseth.a -o quickstart
```

<<< @/examples/output/quickstart.txt{text}

The linear scheme overshoots to 1.18 before the step and undershoots to −0.18 after it. WENO5 puts its weight on the
candidate stencils that do not cross the step and keeps both states within $[0, 1]$; at the face of the step, each
state takes the value of its own side.

## Where next

| To | Read |
|---|---|
| learn the kernels one at a time, from this face to a finite-volume solver | [Tutorial](/tutorial/) |
| solve one task: choose a scheme, bound a reconstruction, a central flux | [Cookbook](/cookbook/) |
| find the kernel, its arguments and its constants | [Reference](/reference/) |
| see the conventions every kernel shares | [Conventions](./conventions) |
| see what each scheme is | [Schemes](./schemes) |
| see how SETH is checked | [Verification](/verification/) |
