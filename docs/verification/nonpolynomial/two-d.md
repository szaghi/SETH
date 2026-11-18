# Two dimensions

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| NP9 | Non-polynomial | P, C | CPU (gfortran, nvfortran, ifx) | yes | `seth_vv_2d` |

## Case and scope

A minimal 2-D harness in SETH's test code (step 4.3 of [issue #6](https://github.com/szaghi/SETH/issues/6), fenced by
owner decision D8: a uniform Cartesian grid, no AMR, no MPI, no curvilinear grid, no viscous terms, no further case
without a new decision). It runs the composite kernels where the papers show their gains, on three cases of
$[0, 1]^2$: an oblique contact, Kelvin–Helmholtz (Liang & Fu 2024, Eq. 46) and the 2-D Riemann problem of
configuration 3 (their Eq. 47). The double Mach reflection, Rayleigh–Taylor, 3-D cases and turbulence are ADAM's.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Harness | 2-D Euler, $q = (\rho, \rho u, \rho v, E)$, $\gamma = 1.4$, 6 ghost cells; the residual is the sum of the 1-D characteristic splitting of the 1-D harness along $x$ and the same code on the data with the momenta swapped along $y$; SSP-RK3 at CFL 0.4; the registry's `reconstruct`, so every scheme of the [ablation](./ablation) runs (the RBF kernels with $s = 3\,\Delta x$) | [`seth_vv_euler2d.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_euler2d.F90) |
| Contact | $\rho = 1.5$ on the band $0.25 \le \operatorname{frac}(x + y) < 0.75$, 1 elsewhere, $u = v = p = 1$, periodic, one period | [`seth_vv_2d.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_2d.F90) |
| Kelvin–Helmholtz | $(2, -0.5, 0.01\sin 2\pi x, 2.5)$ for $0.25 \le y \le 0.75$, $(1, 0.5, 0.01\sin 2\pi x, 2.5)$ elsewhere, periodic, $t = 0.5$ | |
| Configuration 3 | the four states of Liang & Fu's Eq. 47, transmissive, $t = 0.3$ | |
| Asserted | $LR = I$ for the harness's eigenvectors; positivity; conservation to $10^{-13}$ on the periodic cases; symmetry about the diagonal on the contact (configuration 3's reported) | |
| Suite / report | the suite: $N = 48$, TENO5, TENO5-DR, RBF+DR (KH to $t = 0.1$); report mode, on request: $N = 400$, the seven schemes, the density fields | |

## Program

```bash
./exe/seth_vv_2d
scripts/run_vv.sh --only seth_vv_2d      # N = 400, seven schemes: hours on one core
```

## Results

In the suite every run is positive, the periodic cases conserve to about $10^{-14}$, and the solutions of the two
symmetric cases are symmetric bit for bit with gfortran and ifx, in both kinds; with nvfortran (`-fast`, and the
OpenACC builds) the contact stays symmetric to $3\cdot10^{-14}$ and configuration 3 does not: $2.6\cdot10^{-8}$ with
RBF+DR at $N = 48$.

The report run ($N = 400$, the seven schemes of the [ablation](./ablation), gfortran -O3, about 9.5 hours on one
core): every run positive, the periodic cases conserved to $8\cdot10^{-13}$ of their sums, and both symmetric cases
symmetric bit for bit, configuration 3 included.

<!--@include: ../_tables/nonpoly-2d.md-->

The density fields, one colour scale per case so that the panels compare by eye:

![The oblique contact after one period, density at N = 400, seven schemes](/verification/nonpoly-2d-contact.png)

![Kelvin–Helmholtz at t = 0.5, density at N = 400, seven schemes](/verification/nonpoly-2d-kh.png)

![Configuration 3 at t = 0.3, density at N = 400 on the window [0.15, 0.65]^2, seven schemes](/verification/nonpoly-2d-config3.png)

Measured on the fields:

<!--@include: ../_tables/nonpoly-2d-fields.md-->

## Findings

- **The symmetry rests on identical $x$ and $y$ passes** and holds bit for bit where the compiler treats them alike
  (gfortran, ifx). Where it does not (nvfortran), a round-off asymmetry appears, and on configuration 3 its slip lines
  amplify it, to $2.6\cdot10^{-8}$ in density by $t = 0.3$ on the coarse grid of the suite with the low-dissipation
  RBF+DR scheme. This is the mechanism for which the papers run their 2-D cases with a symmetry-preserving technique
  (Fleischmann, Adami & Adams 2019); the harness does not have it, and asserts symmetry only where nothing amplifies
  the round-off (the contact).
- **The overrides keep the 2-D contact sharp** as they keep the 1-D one ([E2](./dr)): after one period on the
  diagonal, a row crosses each jump in 2 cells (TENO5-DR), 2.5 (with ST or RBF) or 3 (THINC) against TENO5's 6, with
  no overshoot. Stencil targeting and the RBF term alone change nothing here (6 cells). The sharp edges follow the
  grid: in the figure the DR contact is visibly stepped where TENO5's is smooth.
- **The overrides leave ripples in smooth flow** (configuration 3): behind the shock, in $[0.55, 0.80]^2$, the
  density has 22 local extrema with TENO5 against 382 (THINC), 396 (DR), 404 (RBF+DR) and 608 (ST+DR), and its rms
  departure from a $7\times 7$ running mean is 1.7 (THINC) to 2.6 (ST+DR) times TENO5's (up to $5\cdot10^{-4}$ on a
  density of 1.5).
  The ripples are visible in the figure as faint wavefronts around the jet. Stencil targeting and the RBF term alone
  add few (44, 68 extrema, the rms unchanged). In 1-D the overrides leave smooth gradients bitwise TENO5's ([E3](./dr));
  in 2-D the likely source is the waves the sharpened slip lines and shock emit into the smooth region (a
  hypothesis: the run does not isolate it). The papers' 2-D figures (Liang
  et al. 2022, Liang & Fu 2024) are not measured this way; these ripples are a cost to weigh against the sharper
  interfaces, not a failure of an asserted property.
- **Kelvin–Helmholtz is not a measure of resolution here**: without the symmetry-preserving technique and with a
  perturbation of $10^{-2}$, the roll-ups at $t = 0.5$ differ between schemes in position and number, so the figure
  is shown for completeness, without a ranking.
