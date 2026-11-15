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
  lose it in the low-pressure gap before the collision. WENO-Z7 and WENO-Z9 fail within the first three
  steps, next to the initial 1000 : 0.01 jump. After one step WENO-Z9 undershoots that jump by 0.6% of its height,
  enough to drive the 0.01 side negative; WENO9-JS, with its heavier weights, stays positive at that step (it fails
  later). Lowering the CFL number does not
  help (measured at 0.2–0.5), nor does refining to 800 cells.
- With the limiter every scheme completes the run (asserted). On the splitting path WENO-Z5, TENO5 and TENO6 need
  about 100 to 500 limited faces, WENO-Z7 and WENO-Z9 thousands; the interpolation path needs 5000 to 15000.
- Borges et al. (2008) and Fu, Hu & Adams (2016) report WENO-Z5 and TENO results on this problem at 400 cells; Fu
  et al. use a Roe flux with entropy fix, not Lax–Friedrichs splitting. SETH's harness reproduces those results only
  with the limiter. Whether the published runs needed a safeguard is not stated in the summaries of those papers in
  the llm-wiki: an open point, under investigation in [issue #3](https://github.com/szaghi/SETH/issues/3).
- The positivity-preserving flux limiter will become a SETH feature ([issue #2](https://github.com/szaghi/SETH/issues/2));
  today it lives in the test harness.
