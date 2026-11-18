---
title: Limiters
---

# Limiters

Two safeguards applied after a reconstruction: a bound on the states of a face, and a limiter on the flux of an
Euler solver. Modules `seth_mp` and `seth_positivity`.

## Monotonicity-preserving bound

```fortran
call seth_mp_bound(vp, vm, alpha, rp, rm)   ! vp, vm(-2:3); rp, rm intent(inout): bounded in place
```

| Argument | Declaration | Meaning |
|---|---|---|
| `vp`, `vm` | `real(R_P), intent(in) :: v(-2:3)` | the data, as passed to the reconstruction (of a wider stencil: its `-2:3` part) |
| `alpha` | `real(R_P), intent(in)` | `SETH_MP_ALPHA` = 4, Suresh & Huynh's value for Runge–Kutta stepping (CFL $\le 0.2$) |
| `rp`, `rm` | `real(R_P), intent(inout)` | the states of any reconstruction kernel, bounded |

Each state becomes the median of itself and Suresh & Huynh's (1997) interval, from the five values around its upwind
cell: on monotone data the interval is $[v_0, v_1]$, near a smooth extremum the curvature terms widen it so that an
accurate state is left alone. Cost on the RTX 4070: 0.42 ns per face. What it is and what it was measured on:
[Schemes](/guide/schemes#monotonicity-preserving-bound). API: [`seth_mp`](/api/src/lib/seth_mp).

## Positivity limiter

```fortran
call seth_positivity_euler(nv, gamma, lam, eps_rho, eps_p, ql, qr, f_high, f_low, f, theta)
ok = seth_euler_admissible(nv, gamma, q, eps_rho, eps_p)
```

The flux limiter of Hu, Adams & Shu (2013) at one face of an ideal-gas Euler solver:
`f = theta f_high + (1 - theta) f_low`, `theta` the largest value in $[0, 1]$ for which the one-sided states
`ql - lam f` and `qr + lam f` keep $\rho \ge$ `eps_rho` and $p \ge$ `eps_p`.

| Argument | Meaning |
|---|---|
| `nv` | 3, 4 or 5: the state is $(\rho, \rho u_1, \dots, \rho u_d, E)$ in $d = $ `nv` $- 2$ dimensions |
| `gamma` | ratio of specific heats |
| `lam` | $2\Delta t/h$ in 1-D; in $d$ dimensions $2\Delta t/(w_k h_k)$ for the faces of direction $k$, the partition weights $w_k$ summing to 1 |
| `eps_rho`, `eps_p` | the floors; Hu, Adams & Shu take $\min(10^{-13}$, the initial minimum$)$; `SETH_POSITIVITY_EPS` $= 10^{-13}$ |
| `ql`, `qr` | the conservative states of the cells left and right of the face |
| `f_high`, `f_low` | the high-order face flux and a first-order positivity-preserving one (Lax–Friedrichs with $\alpha \ge \max(\vert u_n\vert + c)$ is, for `lam` $\alpha \le 1$) |
| `f`, `theta` | the limited flux and its factor (1: `f_high` unchanged) |
| `q` | the conservative state `seth_euler_admissible` checks; `ok` true when $\rho \ge$ `eps_rho` and $p \ge$ `eps_p` |

- If even `f_low` leaves a one-sided state below the floors (`lam` too large, or an inadmissible cell), the kernel
  returns `theta = 0` and `f = f_low`; `seth_euler_admissible` lets the caller check its cells.
- **`nv` at most `SETH_POSITIVITY_NV_MAX` = 5.** The kernel's local arrays have that constant dimension, since an
  array sized by a dummy argument is allocated from the device heap at every call
  ([measured](/guide/performance#device-kernel-rules-arrays-and-loop-extents): 333 times slower). For a larger `nv` it
  returns `theta = 0` and `f = f_low`.
- **Round-off.** The floors hold for the one-sided states as the kernel evaluates them. A caller that evaluates
  `ql - lam f` differently (a fused multiply-add) gets states that differ by $\epsilon(|q| + |\lambda f|)$: choose the
  floors above it.

Verification: [C7](/verification/contracts/positivity). API: [`seth_positivity`](/api/src/lib/seth_positivity).
