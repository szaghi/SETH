---
title: WENO-Z
---

# WENO-Z

The candidates and linear weights of WENO-JS with the weights of Borges et al. (2008), orders 3 to 9. Module
`seth_wenoz`; what the scheme is and which $p$ to use: [Schemes](/guide/schemes#weno-z-weno-z3-5-7-9).

## Calls

```fortran
call seth_wenoz<N>       (vp, vm, eps, p, rp, rm)   ! N = 3, 5, 7, 9
call seth_wenoz<N>_interp(vp, vm, eps, p, rp, rm)   ! the same, from point values at the cell centres
```

| Argument | Declaration | Meaning |
|---|---|---|
| `vp`, `vm` | `real(R_P), intent(in) :: v(1-S:S)`, $S = (N+1)/2$ | left- and right-biased values |
| `eps` | `real(R_P), intent(in)` | regulariser: `SETH_TENO_EPS` ($10^{-40}$ in R8P, $10^{-30}$ in R4P) |
| `p` | `integer(I4P), intent(in)` | exponent of $\tau/(\beta_k + \varepsilon)$: 1 or 2 |
| `rp`, `rm` | `real(R_P), intent(out)` | left and right states at the face |

## The kernels

| Kernel | Bounds | Order | $\tau$ | Recommended $p$ |
|---|---|---|---|---|
| `seth_wenoz3`, `seth_wenoz3_interp` | `-1:2` | 3 | $\lvert\beta_0-\beta_1\rvert$ | 1 |
| `seth_wenoz5`, `seth_wenoz5_interp` | `-2:3` | 5 | $\lvert\beta_0-\beta_2\rvert$ | 1 |
| `seth_wenoz7`, `seth_wenoz7_interp` | `-3:4` | 7 | $\lvert\beta_0+3\beta_1-3\beta_2-\beta_3\rvert$ | 1 |
| `seth_wenoz9`, `seth_wenoz9_interp` | `-4:5` | 9 | $\lvert\beta_0-\beta_4\rvert$ | **2**: with $p = 1$ it is not robust at shocks |

The interpolation kernels use the face coefficients of the Lagrange interpolant of each candidate and the linear
weights of WENO-JS interpolation (equal to ADAM's); the indicators, $\tau$ and the weights are shared with the
reconstruction. Measured orders: [A1](/verification/accuracy/order-of-accuracy); $p$ on the shock cases:
[Schemes](/guide/schemes#weno-z-weno-z3-5-7-9). API: [`seth_wenoz`](/api/src/lib/seth_wenoz).

WENO-Z takes no descaler: $\tau/\beta_k$ is a ratio of indicators, so the weights depend on the data's scale only
through `eps`.
