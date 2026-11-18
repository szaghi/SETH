---
title: A central branch for the Euler equations
---

# A central branch for the Euler equations

The non-dissipative branch of a hybrid for compressible turbulence: a kinetic-energy- and pressure-equilibrium-
preserving flux where the flow is smooth, switched off at shocks by a sensor on the velocity gradient.

## The KEEP-PE flux

The kernel takes the primitive state $(\rho, u, v, w, p)$ of the cells around the face and returns the flux vector:

<<< @/examples/snippets/cookbook_euler-keep.F90{fortran}

A density wave in uniform velocity and pressure must stay so: the property the flux is named for.

<<< @/examples/output/cookbook-keep.txt{text}

## The shock sensors

`seth_ducros` is $\theta^2/(\theta^2 + |\boldsymbol\omega|^2 + \varepsilon)$, $\theta = \nabla\cdot\mathbf u$,
$\boldsymbol\omega = \nabla\times\mathbf u$: near 1 at shocks, near 0 in vortices. The compression-only variant of De
Vanna et al. is zero in expansions too. Both take the velocity gradient the caller computes:

<<< @/examples/snippets/cookbook_euler-ducros.F90{fortran}

<<< @/examples/output/cookbook-ducros.txt{text}

A hybrid takes the central flux where the sensor is below a threshold and the characteristic WENO or TENO flux of
[tutorial 7](/tutorial/07-euler) elsewhere; the measurements on the shock cases: [H2–H3](/verification/hybrids/).

<details><summary>The whole program</summary>

<<< @/examples/snippets/cookbook_euler.F90{fortran}

</details>
