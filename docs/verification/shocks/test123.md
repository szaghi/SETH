# 123 problem

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| S3 | Shock capturing | P | CPU (gfortran, nvfortran, ifx; debug and optimised) | yes | `seth_vv_test123` |

## Case and scope

The 123 problem of Einfeldt et al. (1991) (Toro's test 2): two strong rarefactions leave a near-vacuum region
($\rho, p \approx 0.02$) between them, the standard positivity stress test. Every scheme of both paths runs twice:
without any positivity safeguard (recording whether $\rho$ and $p$ stay positive) and with the positivity-preserving
flux limiter of Hu, Adams & Shu (2013) (counting the faces it limits).

What the case catches: a scheme or path that goes negative near vacuum; a limiter that fails to restore positivity.

## Models and numerics

The harness of [Sod](./sod#models-and-numerics). A run that produces a cell with $\rho \le 0$ or $p \le 0$ at any
Runge–Kutta stage is stopped and recorded as not positive (the next stage would take the square root of a negative
pressure, a floating-point trap in the debug builds).

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain, grid | $[0, 1]$, $N = 200$; report mode 100 … 1600 |
| Left state, $x < 0.5$ | $(\rho, u, p) = (1, -2, 0.4)$ |
| Right state | $(1, 2, 0.4)$ |
| Boundaries | transmissive |
| Final time | $t = 0.15$ |

## Program

```bash
./exe/seth_vv_test123
scripts/run_vv.sh --only seth_vv_test123
```

## Results

![123 problem, 200 cells, t = 0.15](/verification/123.png)

*Density, velocity and pressure at $t = 0.15$, 200 cells, positivity limiter on.*

<!--@include: ../_tables/123.md-->

## Findings

- Without the limiter, every reconstruction up to 7th order stays positive on the splitting path, and so do the
  WENO3-JS, WENO5-JS and WENO-Z3 interpolations (asserted). WENO-Z9 on the splitting path does not.
- **The higher-order interpolations lose positivity without the limiter**: WENO7/9-JS, WENO-Z5/7/9, TENO5 and
  TENO6. WENO-Z5 and TENO5 need CFL ≤ 0.1, and TENO6 fails even there. The flux corrections are not the cause: without them the second-order flux fails too.
  The cause is high-order, low-dissipation face states fed to HLLC next to vacuum. The splitting path, whose global
  Lax–Friedrichs splitting adds dissipation, survives at CFL 0.5. This is a property of how a client uses the
  kernels near vacuum, not a kernel defect: a code using the interpolation path needs a positivity safeguard there
  (to become a SETH feature: [issue #2](https://github.com/szaghi/SETH/issues/2)).
- With the limiter every scheme stays positive (asserted). The limiter also acts on schemes that would have stayed
  positive, because it checks the stricter one-sided states $q \mp 2\lambda F$. Their $L_1$ is then slightly larger
  than without it.
