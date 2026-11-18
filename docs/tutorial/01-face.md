---
title: 1. One face by hand
---

# 1. One face by hand

A SETH kernel computes the two states at **one face** $x_{i+1/2}$ from the cells around it. The cells are indexed by
their offset from the face: cell 0 is left of it, cell 1 right of it. WENO5 reads six cells, `-2:3`.

## Smooth data

The program fills the six cells with the cell averages of $\sin 2\pi x$, cells of width 0.02, the face at
$x = 0.1$, and calls WENO5-JS and its linear scheme:

<<< @/examples/snippets/t01_face-smooth.F90{fortran}

- The same array goes in twice: as `vp`, read with a **left** bias, it gives the state left of the face, `rp`; as
  `vm`, read with a **right** bias, the state right of it, `rm`.
- `SETH_WENO_JS_EPS` regularises the smoothness indicators; `1._R_P` is `rmu`, their descaler (chapter 4).
- `seth_weno5_js_wr` also returns the **weights ratios** `wrp`, `wrm`: $\min_k \omega_k/d_k$, the smallest nonlinear
  weight over its linear one. Near 1, the scheme is close to its linear scheme; near 0, it has dropped a candidate.
- `seth_linear5` is the linear scheme WENO5 tends to on smooth data, one row of coefficients.

## Across a jump

A step between cells 1 and 2, which both states' stencils cross:

<<< @/examples/snippets/t01_face-jump.F90{fortran}

## The mirror

`vm` is not shifted: the kernel mirrors it, reading cell $m$ of the right-biased stencil as cell $1-m$ of a
left-biased one. So the right state of the mirrored data is the left state of the data:

<<< @/examples/snippets/t01_face-mirror.F90{fortran}

## Running it

<<< @/examples/output/t01_face.txt{text}

- On smooth data both schemes are within $10^{-5}$ of the exact value, and the weights ratios are near 1.
- Across the jump, the linear scheme overshoots to 1.05 and 1.18. WENO5-JS puts all the weight on the candidate
  stencils that stay on the side of the 1s, and returns 1: the weights ratios are 0.
- The mirrored call returns the same state twice.

::: tip What you learned
A kernel is one face: two stencils in, two states out, no setup. Reference: [Conventions](/guide/conventions),
[WENO-JS](/reference/weno-js).
:::

<details><summary>The whole program</summary>

<<< @/examples/snippets/t01_face.F90{fortran}

</details>

Next: [2. A line of faces](./02-line).
