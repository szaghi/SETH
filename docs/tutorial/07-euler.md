---
title: 7. A system
---

# 7. A system

The Euler equations of an ideal gas, $q = (\rho, \rho u, E)$, on Sod's shock tube: a rarefaction, a contact and a
shock. Applied to each conservative variable separately, a nonlinear reconstruction mixes waves that travel at
different speeds and oscillates. The standard remedy is to reconstruct in **characteristic** variables, where the
system decouples locally into three scalar waves.

## Characteristic variables

At each face, the eigenvectors of the flux Jacobian at the Roe average of the two neighbouring cells:

<<< @/examples/snippets/t07_euler-eigen.F90{fortran}

## The face flux

For each face: project the stencil's states and fluxes, $w = Lq$ and $g = Lf$; split each field with its own
$\alpha_k = 1.1\max|\lambda_k|$; reconstruct each field with the kernel; and go back with $\hat F = R\,(r_p + r_m)$.
The component-wise variant skips the projection and splits all three variables with one $\alpha$:

<<< @/examples/snippets/t07_euler-residual.F90{fortran}

This is finite differences: the cells hold point values and the face fluxes are reconstructed from the point fluxes
(Shu & Osher 1989). The kernels are the same as for finite volumes; only the meaning of the data changes.

## Running it

<<< @/examples/output/t07_euler.txt{text}

- The exact solution is computed by the program (Toro's Newton iteration on the pressure between the waves).
- The characteristic splitting has the smaller $L_1$ error, and its velocity overshoots the exact plateau $u^*$ by
  $6\cdot10^{-4}$.
- The component-wise splitting overshoots by $6.5\cdot10^{-3}$, and the excess of total variation of the velocity
  over the exact $2u^*$, a measure of the oscillations, is 13 times larger.

::: tip What you learned
For a system, the kernel is called once per characteristic field, on the projected stencil. The projection is the
caller's: SETH reconstructs scalars. The cost of the projection is what a [hybrid](./09-hybrid) saves where the flow is smooth.
:::

<details><summary>The whole program</summary>

<<< @/examples/snippets/t07_euler.F90{fortran}

</details>

Next: [8. Positivity](./08-positivity).
