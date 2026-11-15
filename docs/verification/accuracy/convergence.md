# Whole-domain convergence

| ID | Group | Type | Backends | Figure | Programs |
|---|---|---|---|---|---|
| A3 | Accuracy | P | CPU (gfortran, nvfortran, ifx) | yes | `seth_test_convergence` (reconstructions), `seth_test_interpolation` (interpolations) |

## Case and scope

The pointwise tests measure one face. Here every node of a periodic domain counts: the error is the maximum over the
whole grid, critical points included, so a scheme whose weights lose order at critical points shows it in $L_\infty$.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Field | $u = \sin(\pi x - \sin(\pi x)/\pi)$ on $[-1, 1]$, periodic: Henrick et al.'s (2005) field, with first-order critical points that are not symmetry points (with $\sin \pi x$ the degraded errors would cancel) | [`seth_test_convergence.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_test_convergence.F90) |
| Operator | $D_i = (F_{i+1/2} - F_{i-1/2})/h$, $L_\infty(D - u')$, $N = 20 … 640$, $\varepsilon = 10^{-40}$ | |
| Advection | SSP-RK3 over one period, $\Delta t = 0.4h(h/h_0)^{(o-3)/3}$, $N = 20 … 160$, $L_\infty(u - u_0)$ | |
| Interpolation | point values at the cell centres, both face states, $L_\infty$ over every face, $N = 20 … 640$ | [`seth_test_interpolation.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_test_interpolation.F90) |
| Asserted | WENO-Z ($p = 2$) and TENO at their order; WENO-JS at least Henrick et al.'s critical-point rate $\max(2r-3, 1)$ (the operator) or $r + 1$ (the interpolation, which does not divide by $h$); WENO-Z5 ($p = 1$) at least 4; WENO-Z3 printed only | |

## Program

```bash
./exe/seth_test_convergence ; ./exe/seth_test_interpolation
scripts/run_vv.sh --only seth_test_convergence --only seth_test_interpolation
```

## Results

![Whole-domain operator convergence](/verification/convergence.png)

*$L_\infty$ of the finite-difference derivative over every node of Henrick's field (dotted: WENO-Z with $p = 1$).*

![Whole-domain interpolation convergence](/verification/interpolation.png)

*$L_\infty$ of the interpolated face values over every face.*

Observed orders on the finest pair above round-off (the last report run):

<!--@include: ../_tables/convergence-orders.md-->

<!--@include: ../_tables/interpolation-orders.md-->

## Findings

- WENO-Z with $p = 2$ and TENO keep their order over the whole domain, critical points included. WENO-JS falls to
  its critical-point rate, as Henrick et al. predict.
- WENO-Z3's worst face sits on the field's critical point, and its $L_\infty$ ratio jumps with the position of the
  grid relative to that point (measured on the interpolation: 4.2, 2.0, 4.3, 4.9 from 160 to 2560 cells). It is
  printed, not asserted, on both the operator and the interpolation.
