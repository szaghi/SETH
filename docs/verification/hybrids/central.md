# Central fluxes

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| H2 | Hybrids | P, contract | CPU (gfortran, nvfortran, ifx), device | — | `seth_test_central` |

## Case and scope

The non-dissipative branch of a hybrid: the kinetic-energy- and pressure-equilibrium-preserving (KEEP-PE) flux, in
Pirozzoli's conservative form of order 2, 4, 6 and 8, and the eight-point form with Liang & Fu's (2024b)
dispersion-optimised coefficients; and the two Ducros sensors that decide where it must give way. A KEEP-PE flux is a
sum of two-point products of the full state, so these kernels take the primitive state of the cells around a face,
not a scalar stencil ([Schemes](/guide/schemes#hybrids)).

What the case catches: a wrong coefficient (the order), a wrong two-point flux (the kinetic-energy identity or the
pressure equilibrium), a difference from ADAM's flux, from which SETH's sixth-order one was extracted. ADAM's FLUME
calls `seth_keep_pe6` and `seth_ducros` since ADAM commit 0e085764 (ADAM issue 139), and ADAM commit ba37c4cb deleted
its own: the copy the comparison reads is frozen history.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Order | $(F_{i+1/2} - F_{i-1/2})/h$ of a smooth periodic state against the exact derivative of the Euler flux (complex step), $N = 32 \dots 256$, each direction | [`seth_test_central.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_test_central.F90) |
| Kinetic energy | random $\rho$, $\mathbf u$ and uniform $p$ on 64 periodic cells: $\sum_i \mathbf u_i \cdot R(\rho\mathbf u)_i - \tfrac12\lvert\mathbf u_i\rvert^2 R(\rho)_i$ against its scale | |
| Pressure equilibrium | random $\rho \in [0.01, 1]$, uniform $\mathbf u$ and $p$: $R(\rho\mathbf u) - \mathbf u R(\rho)$, $R(E) - \tfrac12\lvert\mathbf u\rvert^2 R(\rho)$ | |
| ADAM | `seth_keep_pe6` against `compute_kep_flux` copied verbatim (`src/tests/seth_adam_reference.F90`), 20 000 random states (R8P only) | [`seth_test_central_adam.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_test_central_adam.F90) |
| Ducros | the classical $\phi = \theta^2/(\theta^2 + \lvert\omega\rvert^2 + \varepsilon)$ and De Vanna et al.'s compression-only $\psi$ on dilatation, expansion, rotation, and $\psi = \sqrt\phi$ in compression | |

## Program

```bash
./exe/seth_test_central
```

## Results

| Kernel | observed order ($N = 128 \to 256$) | kinetic-energy identity | pressure equilibrium |
|---|---|---|---|
| `seth_keep_pe2` | 2.00 | $6\cdot10^{-18}$ | $2\cdot10^{-15}$ |
| `seth_keep_pe4` | 4.00 | $2\cdot10^{-17}$ | $2\cdot10^{-15}$ |
| `seth_keep_pe6` | 6.00 | $2\cdot10^{-17}$ | $2\cdot10^{-15}$ |
| `seth_keep_pe8` | 7.97 ($64 \to 128$; round-off below) | $5\cdot10^{-17}$ | $4\cdot10^{-15}$ |
| `seth_keep_pe8_opt` | 5.99 | $1\cdot10^{-17}$ | $3\cdot10^{-15}$ |

`seth_keep_pe6` equals ADAM's flux on 100 000 of 100 000 values, bitwise (gfortran 16.2); on the device it equals the
host. Cost on the RTX 4070: 0.777 ns per face against ADAM's 0.774.

## Findings

- Every order is the design one, the kinetic energy of the convective terms is conserved to round-off, and a uniform
  velocity and pressure stay uniform across density ratios of 100.
- The eight-point optimised coefficients of Liang & Fu, printed to 16 digits, give a sixth-order scheme, as they
  state; its error at a given grid is a third larger than the plain sixth-order one on this smooth state, the price
  of a wider resolved band.
- What is not verified here: the junction of a central and a shock-capturing flux (roadmap step 4.3), and any
  turbulent flow (ADAM's issue 119).
