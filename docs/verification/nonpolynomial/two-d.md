# Two dimensions

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| NP9 | Non-polynomial | P, C | CPU (gfortran, nvfortran, ifx) | — | `seth_vv_2d` |

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
RBF+DR at $N = 48$. The density fields at $N = 400$ are
written by the report run; their figures are not yet drawn.

## Findings

- **The symmetry rests on identical $x$ and $y$ passes** and holds bit for bit where the compiler treats them alike
  (gfortran, ifx). Where it does not (nvfortran), a round-off asymmetry appears, and on configuration 3 its slip lines
  amplify it, to $2.6\cdot10^{-8}$ in density by $t = 0.3$ on the coarse grid of the suite with the low-dissipation
  RBF+DR scheme. This is the mechanism for which the papers run their 2-D cases with a symmetry-preserving technique
  (Fleischmann, Adami & Adams 2019); the harness does not have it, and asserts symmetry only where nothing amplifies
  the round-off (the contact).
