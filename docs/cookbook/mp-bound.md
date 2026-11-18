---
title: Bound a reconstruction
---

# Bound a reconstruction

`seth_mp_bound` bounds the two states of a face by the monotonicity-preserving interval of Suresh & Huynh (1997):
call it after any reconstruction kernel, on the five cells around each upwind cell. On monotone data the interval is
the two neighbouring cell values; near a smooth extremum its curvature terms widen it, so an accurate state is left
alone.

<<< @/examples/snippets/cookbook_schemes-mp.F90{fortran}

The ninth-order linear scheme across a step, the step moved through the stencil, and at the maximum of a sine:

<<< @/examples/output/cookbook-mp.txt{text}

The bound removes the oscillations of the linear scheme, up to 19% of the jump, and does not touch the smooth
maximum. `SETH_MP_ALPHA` = 4 is Suresh & Huynh's value for Runge–Kutta stepping, whose monotonicity proof needs CFL
$\le 0.2$. Cost: 0.42 ns per face on the RTX 4070 ([Limiters](/reference/limiters)).
