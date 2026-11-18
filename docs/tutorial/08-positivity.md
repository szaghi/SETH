---
title: 8. Positivity
---

# 8. Positivity

High-order schemes can produce a negative density or pressure next to strong shocks and near vacuum, and the
solver stops: the sound speed is no longer real. The positivity-preserving flux limiter of Hu, Adams & Shu (2013)
blends each face flux with a first-order one just enough to keep the neighbouring cells admissible:

$$ \hat F = \theta\, \hat F_{\rm high} + (1 - \theta)\, \hat F_{\rm low}, \qquad \theta \in [0, 1], $$

$\theta$ the largest value for which the one-sided states $q_i - \lambda\hat F$ and $q_{i+1} + \lambda\hat F$,
$\lambda = 2\Delta t/h$, keep $\rho$ and $p$ above a floor.

## Limiting the face flux

The solver of chapter 7, with reflecting walls, and `seth_positivity_euler` called on each face flux. The first-order
flux is Lax–Friedrichs with $\alpha = \max(|u| + c)$, positive for $\lambda\alpha \le 1$, which CFL 0.4 satisfies:

<<< @/examples/snippets/t08_positivity-residual.F90{fortran}

`SETH_POSITIVITY_EPS` ($10^{-13}$) is the floor; `seth_euler_admissible` lets the program check every cell after each
stage.

## Running it

The blast waves of Woodward and Colella: pressures of 1000, 0.01 and 100 in three regions, two strong shocks that
reflect off the walls and collide.

<<< @/examples/output/t08_positivity.txt{text}

Without the limiter TENO5 produces a negative pressure in the first third of the run. With it, the run reaches the
final time; the limiter acted on 225 face fluxes over the whole run, and left every other flux unchanged.

::: tip What you learned
The limiter is one call per face, after the high-order flux; it costs nothing where it is not needed ($\theta = 1$)
and keeps the run alive where it is. Reference: [Limiters](/reference/limiters); measured on the stress cases:
[C7](/verification/contracts/positivity).
:::

<details><summary>The whole program</summary>

<<< @/examples/snippets/t08_positivity.F90{fortran}

</details>

Next: [9. A hybrid](./09-hybrid).
