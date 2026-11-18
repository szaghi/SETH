# THINC override

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| NP3 | Non-polynomial | P, C | CPU (gfortran, nvfortran, ifx), device | E2 (with DR); the square-wave animation | `seth_test_nonpoly`, `seth_vv_nonpoly`, `seth_test_device`, the registry's cases |

## Case and scope

TENO5-THINC (Takagi, Fu, Wakimura & Xiao 2022) replaces TENO5's state, where its detector (Eq. 20) finds a
discontinuity inside a cell, by the face value of a hyperbolic-tangent step through the cell that conserves the cell
value. This case verifies the override (step 2.1 of [issue #6](https://github.com/szaghi/SETH/issues/6)), the
self-contained kernel `seth_teno5_thinc` (step 2.3) and its cost (E1, step 2.4).

$$h(x) = f_+ + f_-\tanh\big(\beta(X - d)\big),\qquad \hat f_{i+1/2} = f_+ + f_-\,\frac{A + T}{1 + AT},\quad
T = \tanh\frac{\beta}{2},\ A = \frac{\tanh(\alpha\beta/2)}{T},\ \alpha = \frac{f_i - f_+}{f_-}$$

with $f_\pm$ the mean and half difference of $f_{i\pm1}$, $\beta = 1.6$ (Eqs. 15–16; the form with one $\tanh$ and no
exponential is derived and checked by `scripts/gen_nonpoly_tables.py`), and Takagi's first-order guard, Eq. 17: the
cell value where $(f_{i+1}-f_i)(f_i-f_{i-1})\,r_\mu^2 \le 10^{-15}$.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Override | `seth_thinc_override`, `SETH_THINC_BETA` | [`seth_nonpoly.F90`](https://github.com/szaghi/SETH/blob/main/src/lib/seth_nonpoly.F90) |
| Kernel | `seth_teno5_thinc(vp, vm(-5:6), …)`: TENO5, the masks of the cells $-3..3$ recomputed per face (downwind first), Takagi's detector, the override | |
| Specs | on $10^6$ strictly monotone triples: level 0 bitwise inert; the state in $[f_i, f_{i+1}]$; mirror and $2^k$ scaling bitwise; the cell average of the step function through the face value equals $f_i$ to $10^{-13}$, by quadrature as the plan froze it (ten-point Gauss–Legendre on each quarter of the cell; measured $3.5\cdot10^{-15}$, the closed-form integral of $\tanh$ printed beside it); finite on any triple; the guard; on the device (case 59) bitwise the host's | [`seth_test_nonpoly.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_test_nonpoly.F90) |
| N1, N2 | as for [DR](./dr) | |
| Solver | as for [DR](./dr) | [`seth_vv_nonpoly.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_nonpoly.F90) |

## Program

```bash
./exe/seth_test_nonpoly
./exe/seth_vv_nonpoly
```

## Results

The tables of [DR](./dr#results) carry TENO5-THINC beside TENO5-DR (the contact at $N = 100$, experiment E1); the
cases of the registry carry both ([shock capturing](../shocks/)).

![The width of a square wave over 10 periods, TENO5, TENO5-DR, TENO5-THINC](/verification/nonpoly-e2.png)

<div class="seth-casts">
<figure class="seth-cast seth-cast-wide" style="--ar: 2">
<video class="light-only" autoplay loop muted playsinline poster="/verification/casts/nonpoly-square-light.png" aria-label="Square wave with the overrides"><source src="/verification/casts/nonpoly-square-light.webm" type="video/webm"><source src="/verification/casts/nonpoly-square-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/verification/casts/nonpoly-square-dark.png" aria-label="Square wave with the overrides"><source src="/verification/casts/nonpoly-square-dark.webm" type="video/webm"><source src="/verification/casts/nonpoly-square-dark.mp4" type="video/mp4"></video>
<figcaption>The square wave over 10 periods, N = 100; under it, the faces where TENO5-DR and TENO5-THINC override TENO5's state.</figcaption>
</figure>
</div>

## Findings

- **Ten Gauss–Legendre points over the whole cell are not enough to verify conservation to $10^{-13}$**: $\tanh$ of
  slope $\beta = 1.6$ has poles $\pi/(2\beta)$ off the real axis, which limit that rule to about $4\cdot10^{-13}$;
  on quarters of the cell it is exact to round-off. In R4P a face can round to $\bar f \pm \Delta$ (the step
  saturated), which leaves the step's position undetermined: those triples (2 in $10^6$) are counted and skipped.

- **TENO5-THINC improves on TENO5 less than TENO5-DR** on the shock tubes ($L_1(\rho)$ on Sod at $N = 200$: 1.86e-3
  against 2.73e-3 and DR's 1.39e-3), and not at all on Shu–Osher (0.171 for both). Its detector misses the jumps of a
  square wave and of Liang et al.'s profiles 2 and 3 (the defect Liang et al. report and correct with the DR detector).
- **It is over the cost budget.** The per-face kernel recomputes seven masks per state on non-smooth faces and three
  on smooth ones: 3.5 times TENO5 on the GPU against the budget of 3 ([Performance](/guide/performance)); the
  stored-mask variant costs 1.8 times.
