# Verification

Every capability of SETH is backed by a verification case: a Fortran program in `src/tests/` that runs the kernels
and asserts the result against something it shares no code with, such as an exact solution, a fine reference, a
published number, a theorem or ADAM's own code. The same programs are the test suite (`bash scripts/run_tests.sh`
runs all 31, on six compiler builds and on the GPU), and with `SETH_VV_DIR` set they run every scheme on the report
grids and write the data from which these pages are drawn (`scripts/run_vv.sh`). The figures are real SETH results
and can be regenerated after any change of the kernels.

Each card names what its page verifies. Every page has the same sections: case and scope, models and numerics with
links to the code, initial and boundary conditions, the program and how to run it, results, findings. The tag of a
card says what kind of check it is:

- **physics**: a run against an exact or reference solution, an order of accuracy, a spectral property;
- **contract**: a bitwise identity, an equivalence with ADAM, a coefficient table re-derived, a regression baseline,
  host–device agreement, the suite's own power to detect defects.

<div class="seth-casts">
<figure class="seth-cast" style="--ar: 2.25">
<video class="light-only" autoplay loop muted playsinline poster="/verification/casts/shu-osher-light.png" aria-label="Shu–Osher problem, density"><source src="/verification/casts/shu-osher-light.webm" type="video/webm"><source src="/verification/casts/shu-osher-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/verification/casts/shu-osher-dark.png" aria-label="Shu–Osher problem, density"><source src="/verification/casts/shu-osher-dark.webm" type="video/webm"><source src="/verification/casts/shu-osher-dark.mp4" type="video/mp4"></video>
<figcaption>Shu–Osher: a Mach 3 shock runs into entropy waves, 400 cells, WENO5-JS against TENO6 (splitting, characteristic variables). <a href="./shocks/shu-osher">Case S4</a>.</figcaption>
</figure>
<figure class="seth-cast" style="--ar: 2.25">
<video class="light-only" autoplay loop muted playsinline poster="/verification/casts/blast-light.png" aria-label="Woodward–Colella blast waves, density"><source src="/verification/casts/blast-light.webm" type="video/webm"><source src="/verification/casts/blast-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/verification/casts/blast-dark.png" aria-label="Woodward–Colella blast waves, density"><source src="/verification/casts/blast-dark.webm" type="video/webm"><source src="/verification/casts/blast-dark.mp4" type="video/mp4"></video>
<figcaption>Woodward–Colella: two blast waves between walls collide, 400 cells, positivity limiter on. <a href="./shocks/blast">Case S5</a>.</figcaption>
</figure>
<figure class="seth-cast" style="--ar: 2.25">
<video class="light-only" autoplay loop muted playsinline poster="/verification/casts/burgers-light.png" aria-label="Burgers equation, shock formation"><source src="/verification/casts/burgers-light.webm" type="video/webm"><source src="/verification/casts/burgers-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/verification/casts/burgers-dark.png" aria-label="Burgers equation, shock formation"><source src="/verification/casts/burgers-dark.webm" type="video/webm"><source src="/verification/casts/burgers-dark.mp4" type="video/mp4"></video>
<figcaption>Burgers: the shock forms from a sine and moves, against the exact solution, 80 cells. <a href="./shocks/burgers">Case S6</a>.</figcaption>
</figure>
<figure class="seth-cast" style="--ar: 2.25">
<video class="light-only" autoplay loop muted playsinline poster="/verification/casts/advection-light.png" aria-label="Composite profile advection"><source src="/verification/casts/advection-light.webm" type="video/webm"><source src="/verification/casts/advection-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/verification/casts/advection-dark.png" aria-label="Composite profile advection"><source src="/verification/casts/advection-dark.webm" type="video/webm"><source src="/verification/casts/advection-dark.mp4" type="video/mp4"></video>
<figcaption>The Jiang–Shu composite profile over four periods, 200 cells: WENO5-JS, WENO-Z5, TENO5, TENO6. <a href="./shocks/advection">Case S7</a>.</figcaption>
</figure>
</div>

## What a user of the kernels should know

The findings that change how the kernels should be used, each documented on its page:

- The higher-order interpolations need a positivity safeguard near vacuum, and every interpolation on strong blast
  waves; on the Woodward–Colella blast so do WENO9-JS and the WENO-Z (but Z3) and TENO reconstructions
  ([S3](./shocks/test123), [S5](./shocks/blast)). SETH provides one: `seth_positivity_euler`, the flux limiter of
  Hu, Adams & Shu ([API](/guide/api#positivity-limiter), [C7](./contracts/positivity)).
- WENO-Z9 needs the exponent $p = 2$: with $p = 1$ it overshoots by 28% behind the Burgers shock and oscillates
  on the composite profile; with $p = 2$ it stays within the bounds of the other schemes, on the one-dimensional
  cases run here ([S6](./shocks/burgers); [issue #1](https://github.com/szaghi/SETH/issues/1)).
- With Lax–Friedrichs splitting, use $\alpha$ a little above $\max|f'|$. With $\alpha = \max|f'|$ exactly, the split
  flux has a second-order critical point where $|f'|$ peaks, and TENO falls to second order there
  ([S6](./shocks/burgers)).
- TENO6 is anti-dissipative in a narrow band near $\vartheta \approx 0.47\pi$ and leaves grid-scale packets behind a
  shock undamped ([F2](./spectral/adr), [S6](./shocks/burgers)).

## Accuracy

The design order on smooth data, at critical points and over whole domains, for the reconstructions and the
interpolations, and inside a complete Euler solver. [Group page](./accuracy/).

<div class="seth-gallery seth-catalog">

[<strong>Order of accuracy at a point</strong><span>every scheme, both states, reconstruction and interpolation, on a generic smooth function (not exp)</span><em>physics</em>](./accuracy/order-of-accuracy)

[<strong>Critical points</strong><span>the orders the literature measures where $f' = 0$ or $f' = f'' = 0$: Henrick, Borges, Fu, Zhang</span><em>physics</em>](./accuracy/critical-points)

[![Whole-domain convergence](/verification/convergence.png)<strong>Whole-domain convergence</strong><span>$L_\infty$ over every node of Henrick's field, operator and interpolation, critical points included</span><em>physics</em>](./accuracy/convergence)

[![Euler order of accuracy](/verification/euler-order.png)<strong>Euler: order of accuracy</strong><span>a density wave through the complete 1-D Euler solver, both paths, the limiter inactive</span><em>physics</em>](./accuracy/euler-order)

</div>

## Shock capturing

Exact Riemann solutions, fine references, an exact Burgers solution with a moving shock, long-time advection of
discontinuous profiles, positivity. [Group page](./shocks/).

<div class="seth-gallery seth-catalog">

[![Sod shock tube](/verification/sod.png)<strong>Sod shock tube</strong><span>$L_1(\rho)$ against the exact solution, convergence to 1600 cells, the low-dissipation claims</span><em>physics</em>](./shocks/sod)

[![Lax shock tube](/verification/lax.png)<strong>Lax shock tube</strong><span>a strong contact behind a shock, against the exact solution</span><em>physics</em>](./shocks/lax)

[![123 problem](/verification/123.png)<strong>123 problem</strong><span>near vacuum: positivity with and without the limiter, every scheme and path</span><em>physics</em>](./shocks/test123)

[![Shu–Osher](/verification/shu-osher.png)<strong>Shu–Osher</strong><span>entropy waves behind a Mach 3 shock against a 16× finer reference</span><em>physics</em>](./shocks/shu-osher)

[![Woodward–Colella blast](/verification/blast.png)<strong>Woodward–Colella blast</strong><span>colliding blast waves: which schemes need the positivity limiter, and how often</span><em>physics</em>](./shocks/blast)

[![Burgers](/verification/burgers.png)<strong>Burgers equation</strong><span>exact solution before and after the shock: design order, the shock, the smooth region behind it</span><em>physics</em>](./shocks/burgers)

[![Composite advection](/verification/advection.png)<strong>Composite profile</strong><span>four periods of Gaussian, square, triangle and ellipse: dissipation and boundedness</span><em>physics</em>](./shocks/advection)

</div>

## Spectral properties

The modified wavenumber of every linear background and the approximate dispersion relation of every nonlinear
scheme. [Group page](./spectral/).

<div class="seth-gallery seth-catalog">

[![Linear spectra](/verification/spectral.png)<strong>Linear backgrounds</strong><span>no anti-dissipation, error slope $o + 1$, the resolution of Martín et al. (2006)</span><em>physics</em>](./spectral/linear-spectra)

[![ADR](/verification/adr.png)<strong>Approximate dispersion relation</strong><span>Pirozzoli's procedure on every nonlinear scheme: TENO equals its linear scheme where resolved</span><em>physics</em>](./spectral/adr)

</div>

## Hybrids

The pieces of a hybrid scheme: the sensors that choose the linear branch. [Group page](./hybrids/).

<div class="seth-gallery seth-catalog">

[<strong>Central fluxes</strong><span>KEEP-PE of order 2 to 8: kinetic energy and pressure equilibrium to round-off, ADAM's flux bitwise</span><em>physics, contract</em>](./hybrids/central)

[![Sensors](/verification/sensors.png)<strong>Sensors</strong><span>four published sensors on Li, Chen & Ren's static tests: their wavenumbers 0.64 and 2.06 reproduced</span><em>physics</em>](./hybrids/sensors)

</div>

## Contracts

What must hold exactly, or bit for bit: equivalence with ADAM, invariances, the coefficient tables, the regression
baseline, the device, and the suite's own ability to catch defects. [Group page](./contracts/).

<div class="seth-gallery seth-catalog">

[<strong>Equivalence with ADAM</strong><span>WENO-JS bitwise equal to ADAM's kernels on 1.4 M results</span><em>contract</em>](./contracts/adam-equivalence)

[<strong>Properties and invariances</strong><span>exactness, constants, translation, mirror, power-of-2 scaling, robustness, ENO across jumps</span><em>contract</em>](./contracts/properties)

[<strong>Coefficient tables</strong><span>derived in exact arithmetic, cross-checked against papers and ADAM, re-verified when compiled</span><em>contract</em>](./contracts/tables)

[<strong>Golden regression</strong><span>8192 values of every kernel on exact binary inputs; v1 kept bitwise inside v2</span><em>contract</em>](./contracts/golden)

[<strong>Device</strong><span>every kernel in an OpenACC loop on the GPU against the host</span><em>contract</em>](./contracts/device)

[<strong>Mutation testing</strong><span>25 deliberate defects, each caught by at least one test</span><em>contract</em>](./contracts/mutation)

[<strong>Positivity limiter</strong><span>Hu, Adams & Shu's limiter as a SETH kernel: admissible, maximal, untouched where not needed</span><em>contract</em>](./contracts/positivity)

</div>

## Running it

The suite runs on gfortran 16.2, nvfortran 26.9 and ifx 2026.2, debug and optimised, and on the GPU:

```bash
fobis build --mode tests-gnu && bash scripts/run_tests.sh   # the suite (the verification cases on their quick set)
scripts/run_vv.sh                                           # every scheme on the report grids: data, tables, figures
python3 scripts/mutation_test.py --mode tests-gnu           # does every planted defect make a test fail?
```

The list of all test programs, with what each asserts and the page that documents it, is on the
[unit tests](./unit-tests) page.
