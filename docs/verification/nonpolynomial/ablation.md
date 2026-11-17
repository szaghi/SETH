# Ablation

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| NP7 | Non-polynomial | P | CPU (gfortran) | — | `seth_vv_nonpoly` (report mode) |

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
| Solver | the 1-D Euler harness, characteristic Lax–Friedrichs splitting, SSP-RK3, CFL 0.5 | [`seth_vv_euler.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_euler.F90) |
| Cases | Sod and Lax at $N = 200$ against the exact solution; Shu–Osher at $N = 200$ on $[-5, 5]$, $t = 1.8$, against WENO5-JS on 3200 cells; the blast waves at $N = 400$, reflective ends, the positivity limiter on, against WENO5-JS on 6400 cells | [`seth_vv_nonpoly.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_nonpoly.F90) |

## Program

```bash
scripts/run_vv.sh --only seth_vv_nonpoly
```

The ablation runs in report mode only (its references take minutes); the test suite runs the rest of the program.

## Results

<!--@include: ../_tables/nonpoly-ablation.md-->

## Findings

- **On the shock tubes the gain is the DR override's**: every scheme with it halves TENO5's error on Sod and cuts
  Lax's by a third to a half; stencil targeting and the RBF term alone change TENO5's by less than 1 %.
- **On Shu–Osher the gain is the RBF term's**: TENO5-RBF halves TENO5's error and the complete scheme does a little
  better, while DR and ST alone change it by a few %. It is the same term that grows long waves at this $s$ (0.15,
  outside the validated range $s \le 0.12$, [NP6](./rbf)): the entropy waves behind the shock are sharpened by an
  anti-dissipative term, not by a more accurate one, and the result depends on the length unit ([NP8](./units)).
- **On the blast waves the overrides cut the faces the positivity limiter must correct** from about 500 to 70–110,
  with the error unchanged within 5 %.
