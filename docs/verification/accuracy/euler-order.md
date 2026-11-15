# Euler: order of accuracy

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| A4 | Accuracy | P | CPU (gfortran, nvfortran, ifx; debug and optimised) | yes | `seth_vv_euler_order` |

## Case and scope

The design order through a complete solver: a density wave advected by the 1-D Euler equations, with characteristic
projections, flux splitting or a Riemann solver, flux corrections and SSP-RK3 all in the loop. An error in any of
those, or in how the kernels are called, shows as a lost order.

## Models and numerics

The harness of [Sod](../shocks/sod#models-and-numerics), both paths, the positivity limiter on (it must never act).

| Item | Choice |
|---|---|
| Expected order | the reconstruction's on the splitting path; $\min(o, 6)$ on the interpolation path, whose flux corrections are $O(h^6)$ |
| Time step | $\Delta t = 0.4h(h/h_0)^{(o-3)/3}$ times the CFL step, $h_0 = 2/20$: RK3's time error falls with the spatial one |
| Order | the finest pair whose finer $L_1$ is above $10^{-12}$ |
| Asserted | every WENO-Z (but Z3) and TENO at its expected order $- 0.5$; every scheme at least 2; the limiter never acts. WENO-JS is printed: the wave's critical points ($\rho' = 0$ at $x = \pm 1/2$) degrade its weights (Henrick et al. 2005). |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[-1, 1]$, periodic |
| State | $\rho = 1 + 0.2\sin\pi x$, $u = 1$, $p = 1$ |
| Final time | $t = 2$ (one period); the exact solution is the initial one |
| Grids | $N = 20 … 320$ (report mode), $20 … 160$ (test run) |

## Program

```bash
./exe/seth_vv_euler_order
scripts/run_vv.sh --only seth_vv_euler_order
```

## Results

![Euler order of accuracy](/verification/euler-order.png)

*$L_1(\rho)$ of the density wave after one period: splitting path (left), interpolation path (right).*

<!--@include: ../_tables/euler-order.md-->

## Findings

- Both paths reach the design order through the complete solver: WENO-Z5 and TENO5 at 5.00, TENO6 at 6.00, on the
  splitting and the interpolation path (test run).
- The positivity limiter never acts on the smooth wave.
- For the 5th-order schemes the interpolation path is about three times more accurate than the splitting path (the
  HLLC flux carries the contact wave without the splitting's added dissipation); for TENO6 the two coincide.
