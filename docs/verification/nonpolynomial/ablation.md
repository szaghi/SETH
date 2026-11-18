# Ablation

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| NP7 | Non-polynomial | P | CPU (gfortran) | — | `seth_vv_nonpoly` (report mode), `seth_vv_burgers`, `seth_vv_advection` |

## Case and scope

Liang & Fu (2024) combine three pieces in one scheme: stencil targeting (ST, [NP5](./targeting)), the Gaussian-RBF
correction ([NP6](./rbf)) and the DR override ([NP4](./dr)). They print no table of results, only figures of the
complete scheme. Gate N5 of [issue #6](https://github.com/szaghi/SETH/issues/6) asks for the table they lack: every
piece alone and in combination, on the same cases, so that what each buys is measured. The complete scheme is
`seth_teno5_rbf_dr` (step 3.4): `seth_teno5_rbf`'s states, the DR detector on the masks of plain TENO5's cutoff (the
paper's flags), the DR override; with $s_1 = s_2 = 0$ it is ST with the DR override.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Schemes | TENO5; TENO5-ST; TENO5-RBF; TENO5-THINC; TENO5-DR; ST+DR; RBF+DR (Liang & Fu's complete scheme); the RBF kernels with $s = 3\,\Delta x$ ($\alpha_1 = \alpha_2 = 3$), the grid's $\Delta x$ passed by the harness (the registry's `rbf_dx`) | [`seth_test_schemes.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_test_schemes.F90) (the unswept extras) |
| Solver | the 1-D Euler harness, characteristic Lax–Friedrichs splitting, SSP-RK3, CFL 0.5; one scheme at a time, the detector activations counted per characteristic field (`count_activations`) | [`seth_vv_euler.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_euler.F90) |
| Cases | Sod, Lax and the 123 problem (with the positivity limiter) at $N = 200$ against the exact solution: $L_1(\rho)$, the overshoot beyond the exact extrema, the cells across the contact (Sod, Lax); Shu–Osher at $N = 200$ on $[-5, 5]$, $t = 1.8$, against WENO5-JS on 3200 cells; the blast waves at $N = 400$, reflective ends, the positivity limiter on, against WENO5-JS on 6400 cells | [`seth_vv_nonpoly.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_nonpoly.F90) |

| Scalar cases | Burgers after the shock ($N = 640$, $t = 1.5/\pi$) and the Jiang–Shu composite profile after four periods ($N = 200$), the seven schemes in `seth_vv_burgers` and `seth_vv_advection`, which also hold the six non-polynomial ones to TENO5's bounds in the test run | [`seth_vv_burgers.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_burgers.F90), [`seth_vv_advection.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_advection.F90) |

## Program

```bash
scripts/run_vv.sh --only seth_vv_nonpoly
```

The ablation runs in report mode only (its references take minutes); the test suite runs the rest of the program.

## Results

<!--@include: ../_tables/nonpoly-ablation.md-->

<!--@include: ../_tables/nonpoly-ablation-scalar.md-->

## Findings

- **On the shock tubes the gain is the DR override's**: every scheme with it halves TENO5's error on Sod and cuts
  Lax's by a third to a half; stencil targeting and the RBF term alone change TENO5's by less than 1 %.
- **On Shu–Osher the gain is the RBF term's**: TENO5-RBF halves TENO5's error and the complete scheme does a little
  better, while DR and ST alone change it by a few %. It is the same term that grows long waves at this $s$ (0.15,
  outside the validated range $s \le 0.12$, [NP6](./rbf)): the entropy waves behind the shock are sharpened by an
  anti-dissipative term, not by a more accurate one, and the result depends on the length unit ([NP8](./units)).
- **On the blast waves the overrides cut the faces the positivity limiter must correct** from about 500 to 70–110,
  with the error unchanged within 5 %. **On the 123 problem they do the opposite**: 120–154 limited faces against
  TENO5's 76, the overshoot up from $3\cdot10^{-3}$ to $5\text{–}6\cdot10^{-3}$; THINC lowers the error by a fifth,
  DR by 3 %.
- **The overrides cost overshoot on Lax**: $5.5\cdot10^{-3}$ (DR) and $8.1\cdot10^{-3}$ (THINC) against TENO5's
  $4.0\cdot10^{-3}$, while cutting the error by 28–41 %; on Sod no scheme overshoots beyond round-off. The contact is 2–3 cells
  with an override against TENO5's 4–5.
- **On the scalar cases** DR cuts the composite profile's error by 36 % (THINC 33 %) and Burgers' whole-domain error by
  a fifth; nothing changes the smooth region behind the Burgers shock, except the RBF term, which nearly triples it
  ($7\cdot10^{-11}$ against $2.6\cdot10^{-11}$).
