# Lax shock tube

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| S2 | Shock capturing | P | CPU (gfortran, nvfortran, ifx; debug and optimised) | yes | `seth_vv_lax` |

## Case and scope

The shock tube of Lax (1954): a strong contact follows the shock, and low-dissipation schemes oscillate there first.
The exact solution is known, so as for [Sod](./sod) the error is measured, under recorded bounds, and the
low-dissipation claims of WENO-Z5 and TENO5 are asserted against WENO5-JS on each path.

## Models and numerics

The harness of [Sod](./sod#models-and-numerics): 1-D Euler, the splitting path (characteristic Lax–Friedrichs
splitting with the reconstructions) and the interpolation path (HLLC with the interpolations and the flux
corrections), SSP-RK3 at CFL 0.5, the positivity limiter run off and on.

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain, grid | $[0, 1]$, $N = 200$; report mode 100 … 1600 |
| Left state, $x < 0.5$ | $(\rho, u, p) = (0.445, 0.698, 3.528)$ |
| Right state | $(0.5, 0, 0.571)$ |
| Boundaries | transmissive |
| Final time | $t = 0.13$ (Jiang & Shu's $t = 1.3$ on $[-5, 5]$) |

## Program

```bash
./exe/seth_vv_lax                                  # test run
scripts/run_vv.sh --only seth_vv_lax               # report mode, figures, tables
```

## Results

![Lax shock tube, 200 cells, t = 0.13](/verification/lax.png)

*Density, velocity and pressure at $t = 0.13$, 200 cells; inset: the contact and the shock.*

<!--@include: ../_tables/lax.md-->

## Findings

- WENO-Z5 and TENO5 are closer to the exact solution than WENO5-JS on both paths; the interpolation path is again
  the more accurate.
- The positivity limiter limits no face.
