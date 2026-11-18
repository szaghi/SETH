---
title: Central fluxes
---

# Central fluxes and shock sensors

The non-dissipative branch of a hybrid for the Euler equations, and the shock sensors that switch it off. Modules
`seth_kep` and `seth_shock_sensors`.

## KEEP-PE fluxes

A kinetic-energy-preserving face flux above second order is not a reconstruction followed by a flux function: it is
a sum of two-point fluxes of the full state. These kernels therefore take the primitive state of the cells around the
face and return the flux vector: Pirozzoli's (2010) conservative form of order $2L$ with the two-point flux of Shima
et al. (2021, KEEP-PE).

```fortran
call seth_keep_pe6(d, gamma, w, f)
```

| Argument | Declaration | Meaning |
|---|---|---|
| `d` | `integer(I4P), intent(in)` | normal direction: 1 = x, 2 = y, 3 = z |
| `gamma` | `real(R_P), intent(in)` | ratio of specific heats |
| `w` | `real(R_P), intent(in) :: w(5, 1-L:L)` | $\rho, u, v, w, p$ of the cells around the face, cell 0 on its left |
| `f` | `real(R_P), intent(out) :: f(5)` | the flux of $\rho, \rho u, \rho v, \rho w, E$ |

| Kernel | Bounds of `w` | Order |
|---|---|---|
| `seth_keep_pe2` | `(5, 0:1)` | 2 (the two-point flux itself) |
| `seth_keep_pe4` | `(5, -1:2)` | 4 |
| `seth_keep_pe6` | `(5, -2:3)` | 6; ADAM's `compute_kep_flux`, bitwise |
| `seth_keep_pe8` | `(5, -3:4)` | 8 |
| `seth_keep_pe8_opt` | `(5, -3:4)` | 6, the dispersion-optimised coefficients of Liang & Fu (2024b), TENO-S's central branch |

A uniform velocity and pressure stay so across any density profile; on a periodic domain of uniform pressure the
discrete kinetic energy is conserved by the convective terms ([H2](/verification/hybrids/central)). API:
[`seth_kep`](/api/src/lib/seth_kep).

## Shock sensors

```fortran
call seth_ducros(eps, g, phi)              ! phi = theta^2/(theta^2 + |omega|^2 + eps)
call seth_ducros_compression(eps, g, psi)  ! psi = max(-theta/sqrt(theta^2 + |omega|^2 + eps), 0)
```

| Argument | Declaration | Meaning |
|---|---|---|
| `eps` | `real(R_P), intent(in)` | floor of the denominator |
| `g` | `real(R_P), intent(in) :: g(3,3)` | the velocity gradient, `g(a,b)` $= \partial u_a/\partial x_b$, computed by the caller |
| `phi`, `psi` | `real(R_P), intent(out)` | the sensor in $[0, 1]$: near 1 at a shock |

`seth_ducros` is Ducros et al.'s (1999) sensor, $\theta = \nabla\cdot\mathbf u$, $\boldsymbol\omega = \nabla\times\mathbf u$;
`seth_ducros_compression` De Vanna et al.'s (2023a) variant, zero in expansions. API:
[`seth_shock_sensors`](/api/src/lib/seth_shock_sensors).
