# DR override

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| NP4 | Non-polynomial | P, C | CPU (gfortran, nvfortran, ifx), device | — | `seth_test_nonpoly`, `seth_vv_nonpoly`, the registry's cases |

## Case and scope

TENO5-DR (Liang, Xiao, Shyy & Fu 2022) replaces TENO5's state where its detector finds a discontinuity inside a cell
by a two-steepness interpolation through the three cell values, sharper than any polynomial. This case verifies the
override on its own (steps 2.1–2.2 of [issue #6](https://github.com/szaghi/SETH/issues/6)), the self-contained kernel
`seth_teno5_dr` (step 2.3) and its cost (experiment E1, step 2.4); the halo-width proof (E4), contact spreading in time
(E2) and the steepening of smooth gradients (E3) are planned.

$$\hat f_{i+1/2} = \frac{f_{i+1}+f_{i-1}}{2} + \frac{f_{i+1}-f_{i-1}}{2}\,\frac{\tanh A + B}{1 + B\tanh A},\qquad
B = \frac{2f_i - f_{i+1} - f_{i-1}}{f_{i+1}-f_{i-1}}$$

(Eq. 31), with $A = 0.91$ where the detector's Eq. 24 holds and $A = 1.23$ where Eq. 32 holds too; the cell value
$f_i$ where $(f_{i+1}-f_i)(f_i-f_{i-1})\,r_\mu^2 \le 10^{-15}$ (Eq. 33). The detector reads $\zeta$ of the cells
$i-3..i+2$ ([Detectors](./detectors)).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Override | `seth_dr_override`: $\tanh A$ passed as `SETH_DR_T1`, `SETH_DR_T2` (generated, the nearest doubles) | [`seth_nonpoly.F90`](https://github.com/szaghi/SETH/blob/main/src/lib/seth_nonpoly.F90) |
| Kernel | `seth_teno5_dr(vp, vm(-5:6), …)`: TENO5, the masks of the neighbouring cells recomputed per face, the detector (downwind cells first: an early exit with the same result), the override | |
| Specs | on $10^6$ strictly monotone triples: level 0 bitwise inert; the state in $[f_i, f_{i+1}]$; mirror and $2^k$ scaling bitwise; equal to $m + h\tanh(A + \operatorname{atanh} B)$ to $10^{-13}$; linear data give $f_i + d\tanh A$; finite on $10^6$ triples of any sign; the guard | [`seth_test_nonpoly.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_test_nonpoly.F90) |
| N1, N2 | the kernel equals a two-pass, stored-mask reference bitwise on nine periodic scalar fields; inert on the smooth ones | |
| Solver | the 1-D Euler harness, characteristic Lax–Friedrichs splitting, SSP-RK3, CFL 0.5; Sod and Lax at $N = 100$ (step 2.2, spec 7); the cases of the registry ([shock capturing](../shocks/)) | [`seth_vv_nonpoly.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_nonpoly.F90) |
| E1 | the same faces with the masks stored once per cell, each in the characteristic basis of its own face (the harness's `stored_masks`), against the default (every mask in the basis of the face being reconstructed) | [`seth_vv_euler.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_euler.F90) |

## Program

```bash
./exe/seth_test_nonpoly
./exe/seth_vv_nonpoly
scripts/run_vv.sh --only seth_vv_nonpoly
```

## Results

Sod and Lax at $N = 100$, against the exact solution; the cells across the contact are those whose density lies
between 5 % and 95 % of the exact contact jump, within 0.1 of the exact contact (asserted: TENO5-DR's $L_1$ not above
TENO5's, strictly fewer cells):

<!--@include: ../_tables/nonpoly-contact.md-->

Experiment E1: the default kernel against the stored-mask variant, $N = 200$ (Sod, Lax), Shu–Osher on $[-5, 5]$ and the
blast waves at $N = 400$:

<!--@include: ../_tables/nonpoly-e1.md-->

## Findings

- **TENO5-DR halves TENO5's error on the shock tubes**, and its contact is one to two cells narrower; on the cases
  of the registry it is at least as good as TENO5 everywhere it is measured (Sod, Lax, 123, Shu–Osher, the composite
  profile, the blast waves with a quarter of TENO5's limited faces), and keeps fifth order on smooth data and at
  TENO5's critical points.
- **The gain is a splitting effect on Sod.** On the Roe path (one upwind state per field) TENO5-DR's $L_1(\rho)$ on
  Sod is 0.94 of TENO5's, against 0.51 with flux splitting; on Lax the gain survives on both paths
  ([the Roe path](../shocks/blast#the-flux-issue-5), `seth_vv_roe`).
- **The basis of a neighbour's mask changes the result on systems.** On scalar data the two variants are equal bit
  for bit; on Euler they differ by up to 0.03 in density on the shock tubes and 1.9 on the blast waves, with neither
  better on every case. The papers do not say which they meant.
