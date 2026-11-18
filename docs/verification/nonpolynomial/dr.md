# DR override

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| NP4 | Non-polynomial | P, C | CPU (gfortran, nvfortran, ifx), device | contact zooms, E2 (square wave, Sod), E3; animations | `seth_test_nonpoly`, `seth_vv_nonpoly`, `seth_vv_euler_order`, the registry's cases |

## Case and scope

TENO5-DR (Liang, Xiao, Shyy & Fu 2022) replaces TENO5's state where its detector finds a discontinuity inside a cell
by a two-steepness interpolation through the three cell values, sharper than any polynomial. This case verifies the
override on its own (steps 2.1–2.2 of [issue #6](https://github.com/szaghi/SETH/issues/6)), the self-contained kernel
`seth_teno5_dr` (step 2.3), its cost (experiment E1, step 2.4), conservation (N4), the halo width (N6, E4), contact
spreading in time (E2), the steepening of smooth gradients (E3), and the paper's own accuracy test, the Gaussian pulse
of its Table 4 (step 2.2, spec 6).

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
| Specs | on $10^6$ strictly monotone triples: level 0 bitwise inert; the state in $[f_i, f_{i+1}]$; mirror and $2^k$ scaling bitwise; equal to $m + h\tanh(A + \operatorname{atanh} B)$ to $10^{-13}$; the denominator $1 + B\tanh A$ positive; an offset $\mu$ of the data moves the state by $\mu$ to $10^{-12}|\mu|$; linear data give $f_i + d\tanh A$; finite on $10^6$ triples of any sign; the guard; `SETH_DR_T1`, `SETH_DR_T2` within one ulp of $\tanh 0.91$, $\tanh 1.23$ (`seth_test_tables`) | [`seth_test_nonpoly.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_test_nonpoly.F90) |
| N1, N2 | the kernel equals a two-pass, stored-mask reference bitwise on nine periodic scalar fields; inert on the smooth ones; the same for `seth_teno5_rbf_dr` against TENO5-RBF; the library's stored-mask pair (`seth_teno5_face_masks` and the three `*_stored` kernels, step 2.6) equal to both, and on $10^5$ stencils with independent left and right data equal to the self-contained kernels; on the Euler density wave (`seth_vv_euler_order`, $N = 20..160$) the density of DR, THINC and RBF+DR against their base schemes' | [`seth_vv_euler_order.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_euler_order.F90) |
| Spec 6 | Liang et al. 2022, §4.1: $e^{-300(x-0.5)^2}$ on $[0, 1]$, periodic, one period, exact cell averages, SSP-RK3 at CFL $0.1\,(80/N)^{2/3}$ (the paper's time step is not stated), $N = 80, 160, 320$ in the suite, 640 and 1280 in report mode | [`seth_vv_nonpoly.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_nonpoly.F90) |
| Solver | the 1-D Euler harness, characteristic Lax–Friedrichs splitting, SSP-RK3, CFL 0.5; Sod and Lax at $N = 100$ (step 2.2, spec 7); the cases of the registry ([shock capturing](../shocks/)) | [`seth_vv_nonpoly.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_nonpoly.F90) |
| E1 | the same faces with the masks stored once per cell, each in the characteristic basis of its own face (the harness's `stored_masks`, through the library's stored-mask pair), against the default; TENO5-DR, TENO5-THINC and the complete RBF+DR scheme (every mask in the basis of the face being reconstructed) | [`seth_vv_euler.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_euler.F90) |
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

![Sod and Lax at N = 100, the density across the contact: exact, TENO5, TENO5-DR, TENO5-THINC, TENO5-RBF+DR](/verification/nonpoly-contact.png)

<div class="seth-casts">
<figure class="seth-cast seth-cast-wide" style="--ar: 2">
  <video class="light-only" autoplay loop muted playsinline poster="/verification/casts/nonpoly-sod-light.png"><source src="/verification/casts/nonpoly-sod-light.webm" type="video/webm"><source src="/verification/casts/nonpoly-sod-light.mp4" type="video/mp4"></video>
  <video class="dark-only" autoplay loop muted playsinline poster="/verification/casts/nonpoly-sod-dark.png"><source src="/verification/casts/nonpoly-sod-dark.webm" type="video/webm"><source src="/verification/casts/nonpoly-sod-dark.mp4" type="video/mp4"></video>
  <figcaption>Sod's density in time, N = 100: TENO5, TENO5-DR, TENO5-THINC.</figcaption>
</figure>
</div>

The Gaussian pulse of Liang et al. 2022, Table 4 (spec 6):

<!--@include: ../_tables/nonpoly-pulse.md-->

Experiment E1: the default kernel against the stored-mask pair, $N = 200$ (Sod, Lax), Shu–Osher on $[-5, 5]$ and the
blast waves at $N = 400$:

<!--@include: ../_tables/nonpoly-e1.md-->

Conservation and the halo (N4, N6):

<!--@include: ../_tables/nonpoly-halo.md-->

Experiment E2, the contact in time:

![The width of a square wave over 10 periods, TENO5, TENO5-DR, TENO5-THINC](/verification/nonpoly-e2.png)

![The cells across Sod's contact in time, N = 100, TENO5, TENO5-DR, TENO5-THINC](/verification/nonpoly-e2-sod.png)

<div class="seth-casts">
<figure class="seth-cast seth-cast-wide" style="--ar: 2">
  <video class="light-only" autoplay loop muted playsinline poster="/verification/casts/nonpoly-square-light.png"><source src="/verification/casts/nonpoly-square-light.webm" type="video/webm"><source src="/verification/casts/nonpoly-square-light.mp4" type="video/mp4"></video>
  <video class="dark-only" autoplay loop muted playsinline poster="/verification/casts/nonpoly-square-dark.png"><source src="/verification/casts/nonpoly-square-dark.webm" type="video/webm"><source src="/verification/casts/nonpoly-square-dark.mp4" type="video/mp4"></video>
  <figcaption>The square wave over 10 periods, N = 100; under it, the faces where TENO5-DR and TENO5-THINC override TENO5's state.</figcaption>
</figure>
</div>

<!--@include: ../_tables/nonpoly-e2.md-->

Experiment E3, smooth gradients:

![L1 after one period against the width of the gradients, TENO5, TENO5-DR, TENO5-THINC](/verification/nonpoly-e3.png)

<!--@include: ../_tables/nonpoly-e3.md-->

<div class="seth-casts">
<figure class="seth-cast seth-cast-wide" style="--ar: 2">
  <video class="light-only" autoplay loop muted playsinline poster="/verification/casts/nonpoly-blast-light.png"><source src="/verification/casts/nonpoly-blast-light.webm" type="video/webm"><source src="/verification/casts/nonpoly-blast-light.mp4" type="video/mp4"></video>
  <video class="dark-only" autoplay loop muted playsinline poster="/verification/casts/nonpoly-blast-dark.png"><source src="/verification/casts/nonpoly-blast-dark.webm" type="video/webm"><source src="/verification/casts/nonpoly-blast-dark.mp4" type="video/mp4"></video>
  <figcaption>The blast waves under TENO5-DR, N = 400, the limiter on; under them, where the DR detector fires on the density (the scheme detects on each characteristic field).</figcaption>
</figure>
</div>

## Findings

- **The paper's accuracy test holds to round-off, not bit for bit.** On the Gaussian pulse TENO5-DR is TENO5 bit for
  bit at $N = 80$ and 160; at $N = 320$ and beyond the detector fires on the tails, where $|u| \approx 10^{-15}$, the
  scheme's round-off read as structure by TENO5's cutoff (its indicators, $\sim 10^{-30}$, lie far above
  $\varepsilon = 10^{-40}$), and the solution moves by $8\cdot10^{-15}$. The frozen specification (zero overrides,
  bitwise TENO5) is asserted to a noise floor of $10^{-12}$, the form the owner accepted (issue #6, D9 and D10, 2026-11-18). The $L_\infty$ errors match
  the paper's to 0.3–1 % from $N = 160$ (15 % at $N = 80$, where the discretisations differ most). The same holds on
  the Euler density wave: with $u = p = 1$ the acoustic characteristic fields are constant to round-off, and DR,
  THINC and RBF+DR differ from their base schemes by $2\cdot10^{-14}$ at most.

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
  for bit; on Euler they differ by up to 0.04 in density on the shock tubes and 1.9 on the blast waves, with neither
  better on every case (RBF+DR included). The papers do not say which they meant. Both ship: the self-contained
  kernels are the default, the stored-mask pair costs under two times TENO5 against three
  ([issue #6](https://github.com/szaghi/SETH/issues/6), D12).
