---
title: Choose a scheme
---

# Choose a scheme

A reference, from the verification pages and the benchmarks; every number links to its measurement.

## By the data

| Your data | Kernel family | Why |
|---|---|---|
| cell averages (finite volumes), or point values of a split flux (finite differences, Shu & Osher) | the reconstruction kernels: `seth_weno<N>_js`, `seth_wenoz<N>`, `seth_teno*` | they reconstruct a point value from averages |
| point values, to interpolate to the face (WCNS, the alternative finite-difference formulation) | the `_interp` kernels | a reconstruction kernel on point values is second order ([recipe](./point-values)) |

## By what matters

| Need | Choice | Evidence |
|---|---|---|
| a default at fifth order | `seth_teno5` | the smallest $L_1(\rho)$ of the fifth-order schemes on Sod, Lax, 123, Shu–Osher and the blast waves ([S1–S5](/verification/shocks/)); on the blast waves only with the positivity limiter, which WENO5-JS does without |
| fifth order at critical points too | `seth_teno5_az` (`rmu = 1`), at first- and second-order ones; `seth_wenoz5` with $p = 2$ at first-order ones only | [A2](/verification/accuracy/critical-points) |
| higher order on smooth flows | `seth_teno6`, `seth_wenoz7` ($p = 1$), `seth_wenoz9` ($p = 2$ only) | [A1](/verification/accuracy/order-of-accuracy), [S6](/verification/shocks/burgers) |
| low dissipation for turbulence | `seth_teno6`, `seth_teno8`, or a hybrid with a linear or central branch | [F1–F2](/verification/spectral/), [Hybrids](/verification/hybrids/) |
| sharp contacts | `seth_teno5_dr`, `seth_teno5_thinc`, `seth_teno5_rbf_dr` (six ghost cells) | [Non-polynomial](/verification/nonpolynomial/); ripples in 2-D: [NP9](/verification/nonpolynomial/two-d) |
| data of any magnitude | WENO-JS with a [descaler](/reference/weno-js#descalers); WENO-Z and TENO need none | [Tutorial 4](/tutorial/04-weights) |
| positivity near vacuum and strong shocks | any kernel, with `seth_positivity_euler` on the face flux | [C7](/verification/contracts/positivity) |
| the least cost | a hybrid: a sensor, `seth_linear<N>` where smooth | [Tutorial 9](/tutorial/09-hybrid) |

## By cost

ns per face on an RTX 4070 in double precision ([Performance](/guide/performance)):

| Kernel | ns per face |
|---|---|
| `seth_linear5` | 0.095 |
| `seth_sensor_lcr` + `seth_linear5` | 0.39 |
| `seth_weno5_js` | 0.97 |
| `seth_teno5` | 1.44 |
| `seth_weno7_js` | 1.76 |
| `seth_teno8` | 3.98 |
| `seth_teno5_dr`, `seth_teno5_rbf_dr` | 3.9, 4.3 |
| `seth_teno5_thinc` | 5.0 |

This card runs FP64 at 1/64 of its FP32 rate, which favours the kernels with fewer operations, divisions above all;
on a data-centre GPU the ratios are smaller (not measured).

## Ghost cells

The stencil's half-width: 1 for WENO1, 2 for WENO3, WENO-Z3, WENO4-C and linear3; 3 for the five-point and six-point
schemes (WENO5, WENO-Z5, TENO5 and TENO6, WENO6-C); 4 for WENO7, WENO-Z7, TENO8 and WENO8-C; 5 for WENO9 and WENO-Z9;
**6** for the composite non-polynomial kernels ([Reference](/reference/#every-kernel-at-a-glance)).
