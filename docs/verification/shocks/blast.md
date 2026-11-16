# Woodward–Colella interacting blast waves

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| S5 | Shock capturing | P | CPU (gfortran, nvfortran, ifx; debug and optimised) | yes, animation | `seth_vv_blast` |

## Case and scope

Woodward & Colella (1984): between reflecting walls, pressures of 1000, 0.01 and 100 launch two strong blast waves
that reflect and collide, leaving narrow density peaks and contacts. Just before the collision the low-pressure gap
shrinks to one or two cells, where any polynomial reconstruction can go negative; Harten et al. (1987, §7) needed a
local positivity condition there. The reference is WENO5-JS on 16 times the cells (4 times in the test run).

Every scheme runs twice: without any positivity safeguard, recording whether $\rho$ and $p$ stay positive, and with
the positivity-preserving flux limiter of Hu, Adams & Shu (2013), counting the faces it limits.

What the case catches: the robustness of each scheme and path on strong shocks; a limiter that fails; a loss of
resolution (L1 under recorded bounds).

## Models and numerics

The harness of [Sod](./sod#models-and-numerics); reflecting walls (mirror ghosts, normal momentum reversed).

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain, grid | $[0, 1]$; $N = 400$ against $N = 6400$ (report mode); $N = 200$ against $800$ (test run) |
| State | $\rho = 1$, $u = 0$; $p = 1000$ for $x < 0.1$, $0.01$ in the middle, $100$ for $x > 0.9$ |
| Boundaries | reflecting walls |
| Final time | $t = 0.038$ |
| Reference | WENO5-JS, splitting path (it needs no limiter) |

## Program

```bash
./exe/seth_vv_blast
scripts/run_vv.sh --only seth_vv_blast
```

## Results

![Woodward–Colella blast, 400 cells, t = 0.038](/verification/blast.png)

*Density at $t = 0.038$, 400 cells, positivity limiter on, against the 6400-cell reference; right: the contacts and
peaks of the collision.*

<!--@include: ../_tables/blast.md-->

<div class="seth-casts">
<figure class="seth-cast seth-cast-wide" style="--ar: 2.25">
<video class="light-only" autoplay loop muted playsinline poster="/verification/casts/blast-light.png" aria-label="Woodward–Colella blast waves, density"><source src="/verification/casts/blast-light.webm" type="video/webm"><source src="/verification/casts/blast-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/verification/casts/blast-dark.png" aria-label="Woodward–Colella blast waves, density"><source src="/verification/casts/blast-dark.webm" type="video/webm"><source src="/verification/casts/blast-dark.mp4" type="video/mp4"></video>
<figcaption>Density from $t = 0$ to 0.038, 400 cells, positivity limiter on: WENO5-JS and TENO6 (splitting path).</figcaption>
</figure>
</div>

## Findings

- **Without the limiter only WENO3/5/7-JS and WENO-Z3 on the splitting path stay positive** (WENO5-JS asserted);
  WENO9-JS, the higher-order WENO-Z, every TENO and every interpolation lose positivity. WENO-Z5, TENO5 and TENO6
  lose it in the low-pressure gap before the collision. WENO-Z7 fails within the first three
  steps, next to the initial 1000 : 0.01 jump. WENO-Z9 runs with $p = 2$ ([Burgers](./burgers)) and loses
  positivity too; with $p = 1$ it undershot that jump by 0.6% of its height after one step, enough to drive the
  0.01 side negative. Lowering the CFL number does not
  help (measured at 0.2–0.5; for WENO-Z9 with $p = 1$ only), nor does refining to 800 cells.
- With the limiter every scheme completes the run (asserted). On the splitting path WENO-Z5, TENO5 and TENO6 need
  about 100 to 500 limited faces, WENO-Z9 ($p = 2$) 585, WENO-Z7 ($p = 1$) 1826; the interpolation path needs
  5500 to 7100.
- Borges et al. (2008) and Fu, Hu & Adams (2016) report WENO-Z5 and TENO results on this problem at 400 cells, and
  neither paper mentions a positivity safeguard ([issue #3](https://github.com/szaghi/SETH/issues/3)). Borges et al.
  (§5.2) project the Lax–Friedrichs-split fluxes on the Roe eigenvectors, at CFL 0.4; Fu et al. (§5.3) evaluate the
  flux with "a Roe scheme with entropy-fix" and write in their conclusions that the limiter of Hu, Adams & Shu "can
  also be applied to the TENO schemes to further improve robustness". SETH's harness splits the fluxes globally with
  $\alpha = 1.1\max|\lambda|$, so it differs from Fu et al.'s in the flux, and reproduces their results only with
  the limiter.
- **The cutoff decides much of it.** Later TENO papers run this case with $C_T = 10^{-4}$ (llm-wiki weno-teno
  `benchmarks/woodward-colella-blast-wave.md`). In this harness, without the limiter:

  | Scheme | $C_T$ | positive | limited faces with the limiter | $L_1(\rho)$ |
  |---|---|---|---|---|
  | TENO5 | $10^{-5}$ (Fu et al.) | no | 515 | 0.0884 |
  | TENO5 | $10^{-4}$ | no | 78 | 0.0916 |
  | TENO6 | $10^{-7}$ (Fu et al.) | no | 197 | 0.0584 |
  | TENO6 | $10^{-4}$ | **yes** | 13 | 0.0875 |

  A larger cutoff buys robustness with dissipation: TENO6 at $10^{-4}$ needs no safeguard and loses half its accuracy
  advantage. Measured on 2026-11-15 by running the case with the cutoff of the registry changed (not a registry
  default).
- The positivity-preserving flux limiter is a SETH kernel, `seth_positivity_euler` ([C7](../contracts/positivity),
  [API](/guide/api#positivity-limiter)); the harness uses it.
