---
title: 10. Non-polynomial TENO
---

# 10. Non-polynomial TENO

A polynomial reconstruction smears a contact discontinuity over more cells at every step; TENO5 keeps six cells
across the jump of a square wave after one period. The non-polynomial schemes of the TENO line keep TENO5's cutoff
and, where its flags show a discontinuity **inside** a cell, replace the state by the face value of a steep
non-polynomial profile through the cell: a hyperbolic tangent (THINC, Takagi et al. 2022) or a two-steepness
interpolation (DR, Liang et al. 2022).

## The kernels

The composite kernels read the flags of the neighbouring cells, each from the five values centred on it: both states
together read the twelve cells $-5..6$, so the solver needs **six ghost cells**.

<<< @/examples/snippets/t10_nonpoly-kernels.F90{fortran}

- `SETH_DR_T1`, `SETH_DR_T2` and `SETH_THINC_BETA` are the papers' steepnesses; the `1._R_P` is the descaler of the
  overrides' first-order guard (1: the papers' form).
- The RBF kernels take their shape parameters in the grid's units, $s = 3h$: their results depend on the length unit
  ([NP8](/verification/nonpolynomial/units)).

## Running it

<<< @/examples/output/t10_nonpoly.txt{text}

DR keeps each jump within 2 cells and THINC within 4, against TENO5's 6, with no new extrema, and the $L_1$ error
falls by a factor of 2 to 2.8. The two schemes with the DR override, DR alone and the complete scheme of Liang & Fu
(2024) with its Gaussian-RBF rows, have the smallest errors, 4% apart: their order on this problem is within the reach
of round-off, which can flip a detector's decision at a cell or two (the step's arithmetic changed with SELAH's
integrator, and so did these last digits).

The price, measured in two dimensions ([NP9](/verification/nonpolynomial/two-d)): the overrides keep a contact
sharp, but leave small ripples in the smooth flow behind a shock that TENO5 does not.

::: tip What you learned
The non-polynomial kernels are drop-in replacements for TENO5 at the call, with a wider stencil: twelve cells, six
ghost cells. Reference: [Non-polynomial TENO](/reference/nonpolynomial).
:::

<details><summary>The whole program</summary>

<<< @/examples/snippets/t10_nonpoly.F90{fortran}

</details>

Next: [11. On the GPU](./11-gpu).
