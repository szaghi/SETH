# Shu–Osher shock / entropy-wave interaction

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| S4 | Shock capturing | P | CPU (gfortran, nvfortran, ifx; debug and optimised) | yes, animation | `seth_vv_shu_osher` |

## Case and scope

Shu & Osher (1989), example 8: a Mach 3 shock runs into the density waves $\rho = 1 + 0.2\sin 5x$. Behind it the
waves are compressed into a band of short wavelength that a dissipative scheme flattens, so this is the standard
test of resolution next to a shock. There is no exact solution: the reference is WENO5-JS on 16 times the cells
(4 times in the test run). The splitting path runs every reconstruction, the interpolation path every interpolation.

What the case catches: a loss of resolution (L1 against the reference, under recorded bounds); a broken
low-dissipation claim (WENO-Z5 and TENO5 must be closer to the reference than WENO5-JS, on each path); a loss of
positivity.

## Models and numerics

The harness of [Sod](./sod#models-and-numerics), positivity limiter off. $L_1(\rho)$ is measured at the coarse cell
centres against the mean of the two reference cells that share each centre (the refinement ratio is even).

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain, grid | $[-5, 5]$; $N = 400$ against $N = 6400$ (report mode); $N = 200$ against $800$ (test run) |
| $x < -4$ | $(\rho, u, p) = (3.857143, 2.629369, 10.33333)$ |
| $x \ge -4$ | $(1 + 0.2\sin 5x, 0, 1)$ |
| Boundaries | transmissive |
| Final time | $t = 1.8$ |
| Reference | WENO5-JS, splitting path |

## Program

```bash
./exe/seth_vv_shu_osher
scripts/run_vv.sh --only seth_vv_shu_osher
```

## Results

![Shu–Osher, 400 cells, t = 1.8](/verification/shu-osher.png)

*Density at $t = 1.8$, 400 cells, against the 6400-cell reference; right: the compressed entropy waves.*

<!--@include: ../_tables/shu-osher.md-->

<div class="seth-casts">
<figure class="seth-cast seth-cast-wide" style="--ar: 2.25">
<video class="light-only" autoplay loop muted playsinline poster="/verification/casts/shu-osher-light.png" aria-label="Shu–Osher, density"><source src="/verification/casts/shu-osher-light.webm" type="video/webm"><source src="/verification/casts/shu-osher-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/verification/casts/shu-osher-dark.png" aria-label="Shu–Osher, density"><source src="/verification/casts/shu-osher-dark.webm" type="video/webm"><source src="/verification/casts/shu-osher-dark.mp4" type="video/mp4"></video>
<figcaption>Density from $t = 0$ to 1.8, 400 cells: WENO5-JS and TENO6 (splitting path).</figcaption>
</figure>
</div>

## Findings

- The ordering the literature reports holds on both paths: TENO6 resolves the entropy band best, then TENO5 and
  WENO-Z5, then WENO5-JS (asserted for WENO-Z5 and TENO5 against WENO5-JS).
- The interpolation path is closer to the reference than the splitting path at equal order.
