---
title: 5. TENO and its cutoff
---

# 5. TENO and its cutoff

TENO does not blend its candidates smoothly as WENO does. It computes a normalised smoothness measure $\chi_k$ per
candidate and compares it with a cutoff $C_T$: a candidate above it keeps its optimal linear weight $d_k$, a candidate
below it gets weight 0. So on smooth data TENO **is** its linear scheme, to round-off, and near a discontinuity it is
a linear scheme on the candidates that do not cross it.

## When TENO is linear

Four fields on 64 periodic cells, and the faces where TENO5 differs from `seth_linear5`, for three cutoffs:

<<< @/examples/snippets/t05_teno-linear.F90{fortran}

## At a critical point

Where $f' = f'' = 0$ the candidates' indicators all vanish together, their ratios are $O(1)$, and the cutoff drops
candidates on perfectly smooth data: TENO5 falls to second order there. TENO5-AZ (Tang 2024) uses a cutoff that
vanishes with the data's own roughness, $C_T = \min(\lambda_d, \tau)$, and keeps fifth order. The program measures
both in the setting of Henrick et al. (2005), finite differences at the critical point:

<<< @/examples/snippets/t05_teno-critical.F90{fortran}

## Running it

<<< @/examples/output/t05_teno.txt{text}

- On the sines, even at 4 cells per wave, TENO5 is its linear scheme at every face, for every cutoff.
- At the two steps, 10 faces (the five whose stencils cross each step) drop candidates, whatever the cutoff.
- The Gaussian of 2.6 cells is where the cutoff matters: a larger $C_T$ drops candidates at more faces. TENO's
  measure is relative to the data, so the Gaussian's steep flanks count as rough at this resolution.
- At the critical point TENO5 converges at second order; TENO5-AZ at fifth, once $h$ is small enough that $\tau$
  falls below the ceiling $\lambda_d$.

TENO5-AZ's cutoff is dimensional: `rmu = 1` keeps the property, a descaler makes it scale invariant and loses it
([TENO](/reference/teno)). The adaptive TENO5-A and TENO8-A take the cutoff from a sensor instead.

::: tip What you learned
TENO is linear where the cutoff judges the data smooth; the cutoff decides, and its value matters only in the
under-resolved cases. Critical points are its weak spot, and TENO5-AZ the remedy. Reference: [TENO](/reference/teno).
:::

<details><summary>The whole program</summary>

<<< @/examples/snippets/t05_teno.F90{fortran}

</details>

Next: [6. A nonlinear flux](./06-burgers).
