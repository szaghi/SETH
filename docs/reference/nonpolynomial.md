---
title: Non-polynomial TENO
---

# Non-polynomial TENO

Four schemes of the TENO line that keep TENO5's cutoff and use its flags ([issue
#6](https://github.com/szaghi/SETH/issues/6)): two replace TENO5's state where a discontinuity sits inside the
cell (DR, THINC), one selects the widest stencil free of flagged cells (stencil targeting), one adds a Gaussian-RBF
correction to it. Module `seth_nonpoly`; what each scheme is:
[Schemes](/guide/schemes#non-polynomial-teno-teno5-dr-teno5-thinc-stencil-targeting-gaussian-rbf).

## The composite kernels

```fortran
call seth_teno5_dr    (vp, vm, eps, ct, t1, t2, rmu, rp, rm)              ! vp, vm(-5:6)
call seth_teno5_thinc (vp, vm, eps, ct, beta, rmu, rp, rm)                ! vp, vm(-5:6)
call seth_teno5_st    (vp, vm, eps, ct, rp, rm)                           ! vp, vm(-2:3)
call seth_teno5_rbf   (vp, vm, eps, ct, cr, s1, s2, rp, rm)               ! vp, vm(-2:3)
call seth_teno5_rbf_dr(vp, vm, eps, ct, cr, s1, s2, t1, t2, rmu, rp, rm)  ! vp, vm(-5:6)
```

| Kernel | Bounds | Scheme | Source |
|---|---|---|---|
| `seth_teno5_dr` | `-5:6` | TENO5 with the DR override | Liang, Xiao, Shyy & Fu 2022 |
| `seth_teno5_thinc` | `-5:6` | TENO5 with the THINC override | Takagi, Fu, Wakimura & Xiao 2022 |
| `seth_teno5_st` | `-2:3` | stencil targeting: fifth order, fourth beside a jump (TENO5: third) | Liang & Fu 2024, §3.1 |
| `seth_teno5_rbf` | `-2:3` | stencil targeting plus the Gaussian-RBF correction | Liang & Fu 2024, Eqs. 25, 30–32 |
| `seth_teno5_rbf_dr` | `-5:6` | the RBF rows with the DR override: Liang & Fu 2024's complete scheme | Liang & Fu 2024 |

| Argument | Value to pass | Meaning |
|---|---|---|
| `eps`, `ct` | `SETH_TENO_EPS`, `SETH_TENO5_CT` | TENO5's: the detectors are calibrated on its flags |
| `t1`, `t2` | `SETH_DR_T1` = $\tanh 0.91$, `SETH_DR_T2` = $\tanh 1.23$ | DR's two steepnesses |
| `beta` | `SETH_THINC_BETA` = 1.6 | THINC's steepness |
| `rmu` | 1 | scales the overrides' first-order guard: the state falls to the cell value where $(f_{i+1}-f_i)(f_i-f_{i-1})\,r_\mu^2 \le$ `SETH_NONPOLY_GUARD` ($10^{-15}$); 1 is the papers' form |
| `cr` | `SETH_RBF_CR` = 0.4 | $c_r$ of the RBF scheme's sensor |
| `s1`, `s2` | `SETH_RBF_ALPHA1*dx`, `SETH_RBF_ALPHA2*dx` (both $3\,\Delta x$) | the RBF shape parameters, in the units of the grid |

- **Twelve cells.** DR, THINC and the complete scheme read the flags of the neighbouring cells, each from the five
  values centred on it: both states together read `-5..6`. A solver needs **6 ghost cells** for them and projects all
  twelve cells into the face's characteristic basis ([Conventions](/guide/conventions)).
- **Length units.** The RBF kernels' results depend on `dx`: the correction is $s^2$ times a row, $s$ in the grid's
  units. The validated range is $s \le 0.12$ ([NP6](/verification/nonpolynomial/rbf), [NP8](/verification/nonpolynomial/units)).

Every kernel on one stencil, a unit jump on a slope (a compiled example, `docs/examples/src/reference_nonpoly.F90`):

<<< @/examples/snippets/reference_nonpoly-composites.F90{fortran}

## The stored-mask pair

The same three schemes in two passes ([issue #6](https://github.com/szaghi/SETH/issues/6), D12): a first pass
computes each cell's masks once and the caller stores them; a second pass reconstructs every face from TENO5's six
cells and the stored masks. The self-contained kernels above recompute up to seven masks per face instead.

```fortran
call seth_teno5_face_masks(vp, vm, eps, ct, ml, mr)                                ! vp, vm(-2:3); ml, mr integers
call seth_teno5_dr_stored    (vp, vm, mp, mm, eps, ct, t1, t2, rmu, rp, rm)          ! vp, vm(-2:3), mp, mm(-3:3)
call seth_teno5_thinc_stored (vp, vm, mp, mm, eps, ct, beta, rmu, rp, rm)
call seth_teno5_rbf_dr_stored(vp, vm, mp, mm, eps, ct, cr, s1, s2, t1, t2, rmu, rp, rm)
```

| Kernel | Pass | What it does |
|---|---|---|
| `seth_teno5_face_masks` | 1, at every face $i+1/2$ | `ml` the mask of cell $i$ (left state), `mr` the mask of cell $i+1$ on the mirrored data (right state), in the face's basis; the caller stores `msl(i) = ml`, `msr(i+1) = mr` |
| `seth_teno5_dr_stored` | 2 | TENO5's states, the DR detector on the stored masks, the DR override |
| `seth_teno5_thinc_stored` | 2 | TENO5's states, the THINC detector, the THINC override |
| `seth_teno5_rbf_dr_stored` | 2 | the RBF states, the DR detector, the DR override: the complete scheme of Liang & Fu 2024 |

The second pass at face $i+1/2$ reads `mp(c) = msl(i+c)` and `mm(c) = msr(i+1-c)`, $c = -3..3$; the other arguments
are the self-contained kernels'. On a field of $N$ cells the first pass runs over the faces $-3..N+3$, so the stored
arrays need three cells of halo beyond the faces reconstructed.

- **When to use it.** It costs 1.8–2.1 times `seth_teno5` per face, against about 3 (DR, RBF+DR) and 3.5 (THINC)
  for the self-contained kernels, and it does not slow down where jumps are dense: on the GPU, with a jump every 8
  cells, the self-contained kernels reach 5.3 and 5.8 times TENO5 and the stored-mask variant stays at 1.8
  ([Performance](/guide/performance#the-non-polynomial-composites-issue-6)). The caller keeps two integer arrays per
  characteristic field and exchanges their halo.
- **What changes on a system.** Each mask is computed in the basis of the face where it was stored, not of the face
  being reconstructed: on the Euler equations the results differ from the self-contained kernels' (both positive on
  the shock tubes, Shu–Osher and the blast waves, [NP4](/verification/nonpolynomial/dr#results)). On scalar data the
  two are equal bit for bit.

The pair on a field of 16 cells, equal to `seth_teno5_dr` on every face (the same example):

<<< @/examples/snippets/reference_nonpoly-stored.F90{fortran}

## The pieces

For a caller that builds its own composite (another detector, another override, another basis for the masks):

```fortran
call seth_teno5_mask(v, eps, ct, m)                         ! v(-2:2) centred on a cell: m = delta_0 + 2 delta_1 + 4 delta_2
z = seth_teno5_zeta(m)                                      ! zeta of the mask (Liang et al. 2022, Eq. 17): 2, 1, -1 or 0
call seth_dr_detect(mp, mm, lp, lm)                         ! mp, mm(-3:3): the masks around the face; lp, lm in 0..2
call seth_thinc_detect(mp, mm, lp, lm)                      ! lp, lm in 0..1
call seth_dr_override(vp, vm, lp, lm, t1, t2, rmu, rp, rm)  ! vp, vm(-1:2); rp, rm intent(inout)
call seth_thinc_override(vp, vm, lp, lm, beta, rmu, rp, rm)
```

| Name | What it does |
|---|---|
| `seth_teno5_mask` | the 3-bit mask of TENO5's flags on the stencil centred on one cell |
| `seth_teno5_zeta` | $\zeta$ of a mask: 2 where both one-sided candidates are cut, 1 the downwind, $-1$ the upwind, 0 none |
| `seth_dr_detect` | the DR levels of the two states (0 none, 1 or 2 the steepness) from the masks of cells $-3..3$ |
| `seth_thinc_detect` | the THINC levels of the two states (0 or 1) |
| `seth_dr_override` | where a level is 1 or 2, replaces the state by DR's interpolation |
| `seth_thinc_override` | where a level is 1, replaces the state by the face value of the THINC function |

TENO5-DR built from its pieces, equal to `seth_teno5_dr` bit for bit (the same example):

<<< @/examples/snippets/reference_nonpoly-pieces.F90{fortran}

<<< @/examples/output/reference-nonpoly.txt{text}

`mp(c)` is the mask of cell `c` of the left-biased data, `mm(c)` the mask of cell `c` of the mirrored right-biased
data `vm(1-m)`. The constants `SETH_DR_T1` and `SETH_DR_T2` live in `seth_nonpoly_tables`, `SETH_NONPOLY_GUARD` in
`seth_kind_parameters`; `seth` re-exports them. Verification: [Non-polynomial](/verification/nonpolynomial/). API:
[`seth_nonpoly`](/api/src/lib/seth_nonpoly).
