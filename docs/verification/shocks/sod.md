# Sod shock tube

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| S1 | Shock capturing | P | CPU (gfortran, nvfortran, ifx; debug and optimised) | yes | `seth_vv_sod` |

## Case and scope

The shock tube of Sod (1978): two gases at rest separated at $x = 0.5$, high pressure on the left; at $t > 0$ a
left-running rarefaction, a contact and a right-running shock. The exact solution is known (Toro 2009, §4.2–4.5), so
the error is measured, not judged by eye.

What the case catches:

- **a loss of accuracy** of a scheme in a complete Euler solver: $L_1(\rho)$ against the exact solution, under a
  recorded bound;
- **a broken low-dissipation claim**: WENO-Z5 (Borges et al. 2008) and TENO5 (Fu, Hu & Adams 2016) must be closer to
  the exact solution than WENO5-JS, on each path;
- **a limiter that acts where it should not**: on Sod the positivity limiter must limit no face;
- **non-convergence**: report mode measures $L_1$ from 100 to 1600 cells.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | 1-D Euler, ideal gas, $\gamma = 1.4$; point values at the cell centres, conservative finite differences, 5 ghost cells | [`seth_vv_euler.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_euler.F90) |
| Splitting path | Roe-averaged eigenvectors $L$, $R$ at each face; $w = Lq$, $g = Lf$; $g^\pm = (g \pm \alpha_k w)/2$ with $\alpha_k = 1.1\max|\lambda_k|$ over the domain; $g^+$ by the left-biased kernel (`vp`), $g^-$ by the right-biased one (`vm`); $F = R(r_p + r_m)$ | `residual`, PATH_SPLIT |
| Why $1.1\max$ | with $\alpha = \max|\lambda|$ exactly, $\partial_x g^- = (\lambda - \alpha)\partial_x w$ vanishes to second order where $|\lambda|$ peaks: a second-order critical point, where TENO falls to second order ([A2](../accuracy/critical-points)) | |
| Interpolation path | $w = Lq$ interpolated to both sides of the face by the `_interp` kernels; HLLC (Toro, Spruce & Speares 1994; Davis–Einfeldt speeds); $F = F_0 - \tfrac{h^2}{24}f_{xx} + \tfrac{7h^4}{5760}f_{xxxx}$, the derivatives from the six node fluxes around the face (exact coefficients, error $O(h^6)$); a non-positive face state falls back to first order | PATH_INTERP |
| Positivity limiter | $F = \theta F_{\rm high} + (1-\theta)F_{\rm LF}$ per face, $\theta$ the largest keeping $\rho$, $p \ge 10^{-13}$ in $q_i - 2\lambda F$ and $q_{i+1} + 2\lambda F$ (Hu, Adams & Shu 2013), by bisection; positive under CFL $\le 0.5$ | `limit` |
| Time | SSP-RK3, CFL 0.5; a stage with a non-positive cell stops the run (recorded as not positive) | `euler_step` |
| Exact solution | exact Riemann solver: Newton on the pressure functions from the two-rarefaction guess, sampling at $x/t$ (Toro 2009, ch. 4) | [`seth_vv_riemann.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_riemann.F90) |
| Error | $L_1(\rho) = h\sum|\rho_i - \rho(x_i, t)|$, point values at the cell centres, pairwise summation | [`seth_vv_cases.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_cases.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain, grid | $[0, 1]$, $N = 200$ (test run and figure); report mode also 100 … 1600 |
| Left state, $x < 0.5$ | $(\rho, u, p) = (1, 0, 1)$ |
| Right state | $(0.125, 0, 0.1)$ |
| Boundaries | transmissive |
| Final time | $t = 0.2$ |
| Schemes | test run: WENO5-JS, WENO-Z5, TENO5, TENO6 and their interpolations; report mode: every scheme of order $\ge 3$ |

## Program

```bash
fobis build --mode tests-gnu && ./exe/seth_vv_sod                  # the test run
SETH_VV_DIR=$PWD/exe/vv-data ./exe/seth_vv_sod                     # report mode
scripts/run_vv.sh --only seth_vv_sod                               # report mode, figures and tables
```

## Results

![Sod shock tube, 200 cells, t = 0.2](/verification/sod.png)

*Density, velocity and pressure at $t = 0.2$, 200 cells: the exact solution and four schemes (the last on the
interpolation path); inset: the contact.*

<!--@include: ../_tables/sod.md-->

![L1 convergence of the shock tubes](/verification/riemann-convergence.png)

*$L_1(\rho)$ from 100 to 1600 cells on the three shock tubes (dashed: interpolation path).*

## Findings

- On both paths the low-dissipation claims hold: WENO-Z5 and TENO5 have a smaller $L_1(\rho)$ than WENO5-JS, and
  TENO6 the smallest of the quick set.
- The interpolation path is more accurate than the splitting path at equal order. Its HLLC flux resolves the contact
  that the global Lax–Friedrichs splitting smears.
- The positivity limiter limits no face on Sod, for any scheme.
- $L_1$ converges at first order or slightly below (0.8 to 1.0 from 800 to 1600 cells, every scheme and path), as
  expected with a contact and a shock: the order of the scheme sets the error constant, not the rate.
