# Stencil targeting

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| NP5 | Non-polynomial | P, C | CPU (gfortran, nvfortran, ifx), device | — | `seth_test_nonpoly` |

## Case and scope

Liang & Fu (2024, §3.1) keep TENO5's cutoff but replace its weighted candidates by one polynomial row chosen by the
mask: the widest stencil that does not cross the flagged region ("as many smooth data as possible"). `seth_teno5_st`
is that scheme (step 3.1 of [issue #6](https://github.com/szaghi/SETH/issues/6)); its rows are the $\lambda = 0$ rows of
the Gaussian-RBF scheme ([NP6](./rbf)), generated and checked by `scripts/gen_nonpoly_tables.py`.

| Mask $(\delta_0\delta_1\delta_2)$ | Cells | Row | Source |
|---|---|---|---|
| (1,1,1) | $-2..2$ | upwind-5, $(2, -13, 47, 27, -3)/60$ | §3.1 |
| (1,0,1) | $-2..1$ | $3q_0/4 + q_2/4$ | §3.1, Fig. 2 (the paper's $R'_3$: its labels $R'_2$, $R'_3$ are swapped against §3.1) |
| (1,1,0) | $-1..2$ | $q_0/2 + q_1/2$ | §3.1, Fig. 2 (the paper's $R'_2$) |
| one candidate kept | its three cells | that candidate | §3.1 |
| (0,1,1) | $-2..2$ | $3q_1/4 + q_2/4$, TENO5's renormalised weights | not in the paper: owner decision D4 |

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Kernel | the body of `seth_teno5` with `SETH_TENO_ROWS`: TENO5's cutoff, then the row of each state's mask (the preprocessed code of every other kernel unchanged) | [`seth_nonpoly.F90`](https://github.com/szaghi/SETH/blob/main/src/lib/seth_nonpoly.F90) |
| Spec 14 | order 5.0 ± 0.1 at a point (the generic $g$, cell averages, face at 0.3), both states; equal to `seth_linear5` on smooth data; beside a jump between the cells $-2$ and $-1$ (mask (1,1,0)) the order of the smooth side's error 4.0 ± 0.2 where TENO5's is 3.0 ± 0.2; no overshoot across a unit jump | [`seth_test_nonpoly.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_test_nonpoly.F90) |
| Spec 16 | the rows: masks 1, 2, 4 are TENO5's candidates, 7 is `seth_linear5`'s row, 6 TENO5's renormalised weights, exactly | |
| E7 | the masks by data, $10^7$ stencils each: smooth cubics, a jump of random place and height ($10^{-6}..10^2$) on a linear function, a kink, independent values | |

## Program

```bash
./exe/seth_test_nonpoly
scripts/run_vv.sh --only seth_test_nonpoly
```

## Results

The orders (spec 14; the last row is the RBF kernel of [NP6](./rbf)):

<!--@include: ../_tables/nonpoly-st.md-->

Experiment E7, how often each mask occurs:

<!--@include: ../_tables/nonpoly-census.md-->

## Findings

- **Stencil targeting buys one order beside a jump**: with the upwind candidate cut, TENO5 falls to third order
  (its weights 2/3, 1/3 of the two remaining candidates) where the four-point row of stencil targeting keeps fourth.
- **The mask (0,1,1) occurs, but only on noise-like data** (0.05 % of independent values; never on smooth data, a
  single jump or a kink): the choice of decision D4 matters there only.
- **TENO5 never cuts a candidate at a kink** of these data ($C_T = 10^{-5}$): stencil targeting, like TENO5, treats a
  kink as smooth. And smooth cubics close to constant give masks other than (1,1,1) in 0.4 % of the draws, the
  rounding effect of [NP1](./flags).
