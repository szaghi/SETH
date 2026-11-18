---
title: TENO
---

# TENO

Targeted ENO (Fu, Hu & Adams 2016 and its successors): a candidate either keeps its optimal linear weight or is
dropped, by a cutoff $C_T$ on its normalised smoothness measure. On smooth data TENO is its linear scheme exactly.
Module `seth_teno`; what each variant is: [Schemes](/guide/schemes#teno-teno5-teno5-opt-teno5-az-teno5-a-teno6-teno6-opt-teno6-neno6-teno8-teno8-a).

## Calls

```fortran
call seth_teno5       (vp, vm, eps, ct, rp, rm)                   ! vp, vm(-2:3)
call seth_teno5_opt   (vp, vm, eps, ct, rp, rm)
call seth_teno5_az    (vp, vm, eps, ct, rmu, rp, rm)              ! ct: the ceiling lambda_d
call seth_teno5_a     (vp, vm, eps, cr, xi, a1, a2, rmu, rp, rm)  ! the cutoff from a sensor
call seth_teno6       (vp, vm, eps, ct, rp, rm)
call seth_teno6_opt   (vp, vm, eps, ct, rp, rm)
call seth_teno6_neno6 (vp, vm, eps, ct, rp, rm)
call seth_teno8       (vp, vm, eps, ct, rp, rm)                   ! vp, vm(-3:4)
call seth_teno8_a     (vp, vm, eps, cr, xi, a1, a2, rmu, rp, rm)  ! vp, vm(-3:4)
call seth_teno5_interp(vp, vm, eps, ct, rp, rm)                   ! from point values at the cell centres
call seth_teno6_interp(vp, vm, eps, ct, rp, rm)
```

| Argument | Declaration | Meaning |
|---|---|---|
| `vp`, `vm` | `real(R_P), intent(in)` | left- and right-biased values, bounds as in the table below |
| `eps` | `real(R_P), intent(in)` | regulariser: `SETH_TENO_EPS` |
| `ct` | `real(R_P), intent(in)` | cutoff $C_T$: the kernel's `SETH_*_CT` constant |
| `rmu` | `real(R_P), intent(in)` | inverse descaler: 1 (the papers' schemes), or `seth_weno_descaling(...)` |
| `cr`, `xi`, `a1`, `a2` | `real(R_P), intent(in)` | the adaptive cutoff's $C_r$, $\xi$, $a_1$, $a_2$: the `SETH_TENO{5,8}_A_*` constants |
| `rp`, `rm` | `real(R_P), intent(out)` | left and right states at the face |

## The kernels

| Kernel | Bounds | Order | Cutoff constant | Value |
|---|---|---|---|---|
| `seth_teno5` | `-2:3` | 5 | `SETH_TENO5_CT` | $10^{-5}$ |
| `seth_teno5_interp` | `-2:3` | 5 | `SETH_TENO5_CT` | $10^{-5}$ |
| `seth_teno5_opt` | `-2:3` | 4 | `SETH_TENO5_OPT_CT` | $10^{-5}$ |
| `seth_teno5_az` | `-2:3` | 5, also at critical points | `SETH_TENO5_AZ_CT` (the ceiling $\lambda_d$) | $10^{-5}$ |
| `seth_teno5_a` | `-2:3` | 5 | `SETH_TENO5_A_CR`, `SETH_TENO5_A_XI`, `SETH_TENO5_A_A1`, `SETH_TENO5_A_A2` | 0.24, $10^{-3}$, 10, 5: $C_T$ from $10^{-10}$ to $10^{-5}$ |
| `seth_teno6` | `-2:3` | 6 | `SETH_TENO6_CT` | $10^{-7}$ |
| `seth_teno6_interp` | `-2:3` | 6 | `SETH_TENO6_CT` | $10^{-7}$ |
| `seth_teno6_opt` | `-2:3` | 5 | `SETH_TENO6_OPT_CT` | $10^{-7}$ |
| `seth_teno6_neno6` | `-2:3` | 4 | `SETH_TENO6_NENO6_CT` | $10^{-6}$ |
| `seth_teno8` | `-3:4` | 6 | `SETH_TENO8_CT` | $10^{-7}$ |
| `seth_teno8_a` | `-3:4` | 6 | `SETH_TENO8_A_CR`, `SETH_TENO8_A_XI`, `SETH_TENO8_A_A1`, `SETH_TENO8_A_A2` | 0.23, $10^{-3}$, 10.5, 3.5: $C_T$ from $10^{-10}$ to $10^{-7}$ |

`SETH_TENO_EPS` is $10^{-40}$ in R8P (the papers' value) and $10^{-30}$ in R4P, whose smallest normal number is
$1.2\cdot10^{-38}$ ([Kinds and constants](./kinds)). TENO5-opt, TENO6-opt and TENO6-NENO6 have no interpolation
twin: their weights are tuned for the reconstruction.

- **Critical points.** TENO5 and TENO8 fall to second order where $u'' = u' = 0$; TENO5-AZ keeps fifth order there
  with `rmu = 1` (Tang's scheme) and loses it with a descaler, which makes it scale invariant instead: no `rmu` keeps
  both ([A2](/verification/accuracy/critical-points)).
- **Spectral properties.** TENO6 and TENO8 are anti-dissipative over a band of wavenumbers where the cutoff starts
  to drop candidates ([F2](/verification/spectral/adr)).

API: [`seth_teno`](/api/src/lib/seth_teno).
