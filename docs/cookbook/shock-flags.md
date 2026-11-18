---
title: Flag the shocks
---

# Flag the shocks

The `_wr` kernels of WENO-JS also return the weights ratios $\min_k \omega_k/d_k$ of the two states: near 1 where the
nonlinear weights are close to the linear ones, toward 0 where the scheme dropped a candidate. They cost nothing more
than the reconstruction, and flag the faces near a discontinuity: for a limiter, a hybrid's switch, or a refinement
criterion.

<<< @/examples/snippets/cookbook_schemes-flags.F90{fortran}

<<< @/examples/output/cookbook-flags.txt{text}

On a sine of 40 cells with a step, the faces whose stencils cross one of the two jumps are flagged, and no face of
the smooth sine. The threshold is yours; the ratios are continuous, so a blend can use them directly. The dedicated
sensors ([Sensors](/reference/sensors)) are cheaper when the nonlinear kernel is not called anyway.
