---
title: Linear schemes
---

# Linear schemes

The linear background of a scheme, $\sum_k d_k q_k$ collapsed to one row of face coefficients: the smooth branch of
a hybrid. They equal TENO5, TENO6, TENO8 and the TENO5/6 interpolations on smooth data to $4.4\cdot10^{-16}$ of the data, and cost one dot
product per state (0.095 ns per face for `seth_linear5` on the RTX 4070, a tenth of WENO5-JS). Module `seth_linear`.

```fortran
call seth_linear<N>       (vp, vm, rp, rm)   ! N = 3, 5, 7, 9, 6, 8
call seth_linear<N>_interp(vp, vm, rp, rm)   ! N = 3, 5, 7, 9, 6: from point values at the cell centres
```

| Kernel | Bounds | Order | The linear scheme of |
|---|---|---|---|
| `seth_linear3`, `seth_linear3_interp` | `-1:2` | 3 | WENO3 (upwind) |
| `seth_linear5`, `seth_linear5_interp` | `-2:3` | 5 | WENO5, TENO5 |
| `seth_linear7`, `seth_linear7_interp` | `-3:4` | 7 | WENO7 |
| `seth_linear9`, `seth_linear9_interp` | `-4:5` | 9 | WENO9 |
| `seth_linear6`, `seth_linear6_interp` | `-2:3` | 6 | TENO6 |
| `seth_linear8` | `-3:4` | 6 | TENO8 (optimised dispersion, explicit eighth-order dissipation) |

A linear background is not a central scheme: the upwind ones keep their dissipation. A non-dissipative branch for
the Euler equations is [Central fluxes](./central). What a hybrid needs besides: [Sensors](./sensors);
how to build one: [Schemes](/guide/schemes#hybrids). API: [`seth_linear`](/api/src/lib/seth_linear).
