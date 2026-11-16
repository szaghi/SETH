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

- **Without the limiter only WENO3/5/7-JS, WENO-Z3 and WENO4/6-C on the splitting path stay positive** (WENO5-JS
  asserted); WENO9-JS, the higher-order WENO-Z, WENO8-C, every TENO and every interpolation lose positivity. WENO-Z5, TENO5 and TENO6
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
  also be applied to the TENO schemes to further improve robustness". With that flux TENO6 at their cutoff needs no
  safeguard, TENO5 at theirs still does ([below](#the-flux-issue-5)).
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
## The flux (issue #5)

[Issue #5](https://github.com/szaghi/SETH/issues/5) asked whether the flux explains the difference with Fu, Hu &
Adams (2016). Their "Roe scheme with entropy-fix" is not defined in the paper; in the lineage it is Jiang & Shu's
(1996) WENO-RF, which ran this very case: Roe's flux splitting with entropy fix of Shu & Osher (1989, Algorithm 2.3),
upwind by the sign of each characteristic field's Roe eigenvalue, local Lax–Friedrichs where that eigenvalue changes
sign (`notes/issue5-roe-sources.md`). The harness has it as `PATH_ROE` ([models](./sod#models-and-numerics)), checked
on its own by `seth_vv_roe`: design order on the density wave, Sod and Lax within 0.71–1.18 of the splitting path's
$L_1$, the transonic rarefaction of Toro's test 1 resolved with the fix and broken by an expansion shock without it
(error 4–13 times larger).

Same case and grid as above (N = 400 against 6400 cells), every reconstruction; positive without the limiter at CFL
0.5 (the Roe path also at Fu et al.'s 0.6, with the same outcome for every scheme but WENO6-C, positive at 0.5
only), then with the limiter:

| Scheme | positive, GLF | positive, Roe | limited faces, GLF | limited faces, Roe | $L_1(\rho)$, GLF | $L_1(\rho)$, Roe |
|---|---|---|---|---|---|---|
| WENO3-JS | yes | yes | 0 | 0 | 0.2148 | 0.1364 |
| WENO5-JS | yes | yes | 0 | 0 | 0.1312 | 0.0843 |
| WENO7-JS | yes | **no** | 2 | 41 | 0.1084 | 0.0688 |
| WENO9-JS | **no** | **no** | 67 | 72 | 0.0895 | 0.0612 |
| WENO-Z3 | yes | yes | 0 | 0 | 0.1832 | 0.1096 |
| WENO-Z5 | **no** | yes | 104 | 0 | 0.0977 | 0.0603 |
| WENO-Z7 | **no** | **no** | 1826 | 2023 | 0.0743 | 0.0534 |
| WENO-Z9 | **no** | **no** | 585 | 892 | 0.0774 | 0.0522 |
| TENO5 ($C_T = 10^{-5}$) | **no** | **no** | 515 | 36 | 0.0884 | 0.0570 |
| TENO5-opt | **no** | yes | 111 | 0 | 0.0798 | 0.0567 |
| TENO5-AZ | **no** | **no** | 824 | 1701 | 0.0913 | 0.0552 |
| TENO6 ($C_T = 10^{-7}$) | **no** | yes | 197 | 5 | 0.0584 | 0.0502 |
| TENO6-opt | **no** | yes | 251 | 8 | 0.0644 | 0.0504 |
| TENO6-NENO6 | **no** | yes | 84 | 0 | 0.0710 | 0.0535 |
| TENO8 | **no** | yes | 418 | 2 | 0.0568 | 0.0463 |
| TENO5-A | **no** | **no** | 1930 | 21 | 0.0865 | 0.0527 |
| TENO8-A | **no** | **no** | 760 | 3 | 0.0557 | 0.0429 |
| WENO4-C | yes | yes | 0 | 0 | 0.1904 | 0.1178 |
| WENO6-C | yes | yes | 1 | 3 | 0.1190 | 0.0754 |
| WENO8-C | **no** | **no** | 4 | 1623 | 0.1013 | 0.0703 |

- **The answer is: largely, not entirely.** With the Roe flux, WENO-Z5, TENO6 at Fu et al.'s $C_T = 10^{-7}$, TENO5-opt,
  TENO6-opt, TENO6-NENO6 and TENO8 need no safeguard; TENO5 at their $C_T = 10^{-5}$ still loses positivity (36
  limited faces instead of 515), and so do TENO5-AZ, TENO5-A, TENO8-A, WENO-Z7/9 and WENO9-JS. The Roe flux is not
  uniformly more robust: WENO7-JS, positive with GLF, is not with it, and WENO8-C needs 1623 limited faces with it
  against 4 with GLF ([centred WENO](../accuracy/centred-weno)). $L_1$ falls by 14–40% on every scheme: the
  global splitting dissipates everywhere at the rate of the largest wave speed of the domain (that of the 1000 : 0.01
  jump), the Roe flux at the local one.
- **What does it** (one-off variants of the harness on the quick schemes, recorded in the note): the entropy fix plays
  no part (Roe without it: the same outcomes); local Lax–Friedrichs in every field, with the local $\alpha$ of the
  fix, makes WENO-Z5 positive but not TENO6; global Lax–Friedrichs without the 10% margin changes nothing. The local
  wave speed rescues WENO-Z5, and TENO6 needs the unsplit upwind reconstruction. At $C_T = 10^{-4}$ TENO5 is
  positive on the Roe path, not on GLF.
- The harness keeps global Lax–Friedrichs as its splitting path (it needs no sonic switch and is what the other TENO
  papers use); the Roe path exists to answer this question, and runs in `seth_vv_roe` only.

## Notes

- The positivity-preserving flux limiter is a SETH kernel, `seth_positivity_euler` ([C7](../contracts/positivity),
  [API](/guide/api#positivity-limiter)); the harness uses it.
