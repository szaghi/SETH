# Centred WENO

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| A5 | Accuracy | P, C | CPU (gfortran, nvfortran, ifx; debug and optimised), GPU | — | `seth_test_weno_c`, `seth_test_weno_c_adam` |

## Case and scope

The centred WENO (`seth_weno4_c`, `seth_weno6_c`, `seth_weno8_c`, [Schemes](/guide/schemes#centred-weno-weno4-6-8-c)):
ADAM's `weno_reconstruct_centered`, which blends the upwind WENO-JS of order $2S-1$ with the central reconstruction
of order $2S$, ported with ADAM's defect D1 fixed. Three questions: does it reach its order, is it as ENO as its
upwind part, and is it ADAM's scheme with exactly D1 corrected? The bounds were frozen before the kernel was
written; the two changed after measuring are named below, with the reason.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Order | cell averages of the generic function $g$ ([A1](./order-of-accuracy)), face at $x = 0.3$, $h = 0.4 \dots 0.0125$; the finest pair above round-off | [`seth_test_weno_c.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_test_weno_c.F90) |
| Properties | mirror ($v_m(m) = v_p(1-m)$), constant data, a jump of $10^{30}$ ($10^{16}$ in R4P) anywhere | |
| ENO | a unit step at every cell boundary of the big stencil, both orientations: the worst value outside $[0, 1]$, against the upwind WENO-JS of the same $S$ | |
| ADAM | 20000 random stencils per $S$, smooth and rough, magnitudes $10^{-6} \dots 10^{6}$, $r_\mu = 1$; ADAM's routine copied verbatim, and with its one D1 line corrected | [`seth_test_weno_c_adam.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_test_weno_c_adam.F90), [`seth_adam_reference.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_adam_reference.F90) |
| Asserted | order $\ge 2S - 0.3$ (R4P: less `CENTRED_ORDER_LOSS`); mirror to WENO-JS's tolerance; constant data to 4 ulps; finite at the jump; overshoot at most the upwind scheme's $+ 10^{-3}$ ($+ \tfrac1{280}$ for $S = 4$, below); SETH equal to ADAM's corrected routine to $64\varepsilon$ of the data, and ADAM's routine as it is different on the left state | |

## Program

```bash
./exe/seth_test_weno_c && ./exe/seth_test_weno_c_adam
scripts/run_vv.sh --only seth_test_weno_c --only seth_test_weno_c_adam
```

The same registry entries put the three kernels in every sweeping test: critical points ([A2](./critical-points)),
whole-domain convergence ([A3](./convergence)), properties ([C2](../contracts/properties)), the golden files
([C4](../contracts/golden), v7 and r4_v2), host against device ([C5](../contracts/device)), and the shock cases in
report mode.

## Results

Order at a point:

<!--@include: ../_tables/weno-c-orders.md-->

A unit step at every cell boundary of the big stencil:

<!--@include: ../_tables/weno-c-step.md-->

Against ADAM:

<!--@include: ../_tables/weno-c-adam.md-->

## Findings

- **The design orders hold**: 4, 6 and 8 at a point, both states. At critical points the orders are those of the
  upwind WENO-JS ([A2](./critical-points)): where a difference of the blend's cells $i-1..i+2$ vanishes, $\phi$ is
  $O(1)$ and the state is the upwind one. The whole-domain convergence on Henrick's field is that of WENO-JS for the
  same reason ([A3](./convergence)).
- **D1 is real, and fixed.** SETH's states equal ADAM's routine with the left weights computed on the shifted cells to
  round-off; ADAM's routine as it is gives the same right state and a different left one on most random stencils,
  by up to a few times the data's magnitude. SETH's fix needs no code of its own: the upwind part is
  `seth_weno<2S-1>_js`, which reads the face-offset stencil.
- **WENO8-C is not ENO for a jump in its outer cells.** The blend's detector reads cells $i-1..i+2$ for every $S$;
  a step in WENO8-C's cells $i-3$ or $i+4$ leaves $\phi = 0$, and the central row's edge coefficient passes: an
  overshoot of $\tfrac1{280} = 3.6\cdot10^{-3}$ of the jump (the frozen bound, the upwind scheme's overshoot plus
  $10^{-3}$, was raised by $\tfrac1{280}$ for $S = 4$ when this was measured). WENO6-C's edge coefficient is
  $+\tfrac1{60}$: a precursor of the jump inside the data's range. Widening the detector would change ADAM's scheme:
  left to the owner.
- **In the flows** (report mode, the same harness as every scheme): on Sod, Lax, Shu–Osher and the blast wave the
  $L_1(\rho)$ of each centred WENO lies between that of its upwind WENO-JS and that of the next WENO-JS order; on the
  123 problem too, except WENO8-C, 1% above WENO7-JS. WENO4-C and WENO6-C stay positive on the blast
  wave without the limiter; WENO8-C does not, and with the Roe flux it needs 1623 limited faces against 4 with the
  splitting ([S5](../shocks/blast#the-flux-issue-5)). After the Burgers shock WENO6-C and WENO8-C leave a
  grid-scale packet upstream of it, and their smooth-region error falls to $\sim3\cdot10^{-6}$ at 640 cells, where
  their upwind parts reach $10^{-10}$ and below ([S6](../shocks/burgers)). The cause is the blend: a small packet on
  a smooth background is damped at $\phi$ times the upwind rate, with the background's $\phi = O(h)$, a damping that
  does not grow with the resolution; the ADR, one wave of unit amplitude with no background, cannot show it.
- **The mirror is WENO-JS's**, not bitwise: the central sums and the blend mirror exactly, the upwind WENO-JS sums its
  mirrored tables in another order ([C2](../contracts/properties)). The frozen bound of 4 ulps was raised to WENO-JS's
  tolerance for that reason.
- **In R4P** the pairs above the round-off floor are the coarse, pre-asymptotic ones (the same errors as R8P there):
  WENO6-C's right state measures 4.69, hence `CENTRED_ORDER_LOSS` = 1.06 in R4P.
