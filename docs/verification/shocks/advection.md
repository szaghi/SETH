# Composite profile: long-time advection

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| S7 | Shock capturing | P | CPU (gfortran, nvfortran, ifx; debug and optimised) | yes, animation | `seth_vv_advection` |

## Case and scope

Jiang & Shu (1996), example 2: a Gaussian pulse, a square wave, a triangle and a half ellipse, advected over four
periods. One problem tests dissipation (the square wave's corners, the ellipse), dispersion (the pulse) and the
non-oscillatory property at once.

What the case catches: a growing oscillation (the solution must stay within $[-0.05, 1.05]$); a broken
low-dissipation claim (WENO-Z5 and TENO5 must have a smaller $L_1$ error than WENO5-JS: Borges et al. 2008; Fu, Hu
& Adams 2016).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equation | $u_t + u_x = 0$, periodic, upwind (the splitting with $\alpha = 1$), SSP-RK3, CFL 0.4 | [`seth_vv_scalar.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_scalar.F90) |
| Profile | constants as llm-wiki high-order-capturing-schemes `jiang-shu-composite-advection.md` | [`seth_vv_advection.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_advection.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain, grid | $[-1, 1]$, periodic, 200 cells |
| Final time | $t = 8$ (four periods) |
| Schemes | test run: WENO5-JS, WENO-Z5, TENO5, TENO6; report mode: every reconstruction of order $\ge 3$ |

## Program

```bash
./exe/seth_vv_advection
scripts/run_vv.sh --only seth_vv_advection
```

## Results

![Composite profile after four periods](/verification/advection.png)

*The composite profile after four periods, 200 cells; right: the square wave and the half ellipse.*

<!--@include: ../_tables/advection.md-->

<div class="seth-casts">
<figure class="seth-cast seth-cast-wide" style="--ar: 2.25">
<video class="light-only" autoplay loop muted playsinline poster="/verification/casts/advection-light.png" aria-label="Composite profile advection"><source src="/verification/casts/advection-light.webm" type="video/webm"><source src="/verification/casts/advection-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/verification/casts/advection-dark.png" aria-label="Composite profile advection"><source src="/verification/casts/advection-dark.webm" type="video/webm"><source src="/verification/casts/advection-dark.mp4" type="video/mp4"></video>
<figcaption>Four periods, 200 cells: the exact profile, WENO5-JS, WENO-Z5, TENO5, TENO6.</figcaption>
</figure>
</div>

## Findings

- WENO-Z5 and TENO5 dissipate less than WENO5-JS (asserted); TENO6 least of the four.
- No scheme leaves $[-0.05, 1.05]$; the largest overshoot of the quick set is TENO5's, about 1% at the square wave.
