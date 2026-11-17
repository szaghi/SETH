# DR override

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| NP4 | Non-polynomial | P, C | CPU (gfortran, nvfortran, ifx), device | — | `seth_test_nonpoly`, `seth_vv_nonpoly`, the registry's cases |

## Case and scope

TENO5-DR (Liang, Xiao, Shyy & Fu 2022) replaces TENO5's state where its detector finds a discontinuity inside a cell
by a two-steepness interpolation through the three cell values, sharper than any polynomial. This case verifies the
override on its own (steps 2.1–2.2 of [issue #6](https://github.com/szaghi/SETH/issues/6)), the self-contained kernel
`seth_teno5_dr` (step 2.3), its cost (experiment E1, step 2.4), conservation (N4), the halo width (N6, E4), contact
spreading in time (E2) and the steepening of smooth gradients (E3).

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
| N4, N6 | linear advection, periodic, a square wave on a sine, $N = 200$, one period: the sum of $u$ conserved to $10^{-13}$; the domain split in two blocks with $G$ ghost cells refilled from each other at every stage: $G = 6$ bitwise equal to one block, $G = 5$ (the sixth cell unfilled) different | [`seth_vv_nonpoly.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_nonpoly.F90) |
| E2 | a square wave advected 10 periods, $N = 100$: the cells across its jumps after each period | |
| E3 | a smooth periodic plateau whose two gradients are $w$ cells wide, one period, $N = 100$: $L_1$ against TENO5's | |

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

Conservation and the halo (N4, N6):

<!--@include: ../_tables/nonpoly-halo.md-->

Experiment E2, the contact in time:

![The width of a square wave over 10 periods, TENO5, TENO5-DR, TENO5-THINC](/verification/nonpoly-e2.png)

<!--@include: ../_tables/nonpoly-e2.md-->

Experiment E3, smooth gradients:

<!--@include: ../_tables/nonpoly-e3.md-->

## Findings

- **TENO5-DR halves TENO5's error on the shock tubes**, and its contact is one to two cells narrower; on the cases
  of the registry it is at least as good as TENO5 everywhere it is measured (Sod, Lax, 123, Shu–Osher, the composite
  profile, the blast waves with a quarter of TENO5's limited faces), and keeps fifth order on smooth data and at
  TENO5's critical points.
- **The gain is a splitting effect on Sod.** On the Roe path (one upwind state per field) TENO5-DR's $L_1(\rho)$ on
  Sod is 0.94 of TENO5's, against 0.51 with flux splitting; on Lax the gain survives on both paths
  ([the Roe path](../shocks/blast#the-flux-issue-5), `seth_vv_roe`).
- **The contact stops spreading.** Over ten periods a square wave keeps 4–5 cells across its jumps with TENO5-DR
  (4–8 with TENO5-THINC), where TENO5's widen to 12–15 (E2).
- **No smooth gradient is steepened**: a gradient two cells wide or wider is reconstructed bitwise as by TENO5 (the
  detector never fires); at one cell, a jump in all but name, both overrides act and cut the error 3 to 9 times (E3).
- **The halo is 6 cells, proved**: with 6 ghost cells a two-block domain equals one block bit for bit, with the sixth
  cell unfilled it does not; the scheme conserves to $2\cdot10^{-14}$ with thousands of overrides active (N4, N6).
  (A five-cell halo padded by constant continuation happened to give the same result in this run: the sixth cell
  matters only when the detector reaches the far upwind mask at a block boundary.)
- **The basis of a neighbour's mask changes the result on systems.** On scalar data the two variants are equal bit
  for bit; on Euler they differ by up to 0.03 in density on the shock tubes and 1.9 on the blast waves, with neither
  better on every case. The papers do not say which they meant.
