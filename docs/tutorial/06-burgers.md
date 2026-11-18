---
title: 6. A nonlinear flux
---

# 6. A nonlinear flux

Burgers' equation, $u_t + (u^2/2)_x = 0$, from $u_0 = 0.5 + \sin 2\pi x$ on the periodic unit interval. The sine
steepens and breaks into a shock at $t = 1/(2\pi)$; at $t = 0.5$ the exact solution has one maximum, one minimum and
a shock between them.

## The splitting

With $f' = u$ both split fluxes are nonzero, $f^\pm = \tfrac12(u^2/2 \pm \alpha u)$. $\alpha$ must bound $|f'|$ over
the stencils: the program takes 10% more than $\max|u|$ over the field, at every stage. With $\alpha = \max|f'|$
exactly, $f^-$ has a second-order critical point where $|u|$ peaks, and TENO falls to second order there (chapter 5;
[S6](/verification/shocks/burgers)).

<<< @/examples/snippets/t06_burgers-residual.F90{fortran}

WENO-Z7 reads cells $-3..4$: the program has four ghost cells, enough for every scheme it calls.

## Running it

<<< @/examples/output/t06_burgers.txt{text}

- The three nonlinear schemes keep the exact solution's two extrema: no oscillation at the shock. Their total
  variation decreases from the initial 4.0, as the exact solution's does.
- The linear scheme makes six new extrema, oscillations around the shock, and its total variation barely decreases.
- All conserve the mass to round-off: the face fluxes are differenced, whatever the scheme.

::: tip What you learned
A nonlinear flux changes only the splitting: $\alpha$ a little above $\max|f'|$, recomputed at each stage. The
kernels do not know the equation.
:::

<details><summary>The whole program</summary>

<<< @/examples/snippets/t06_burgers.F90{fortran}

</details>

Next: [7. A system](./07-euler).
