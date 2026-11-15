# Properties and invariances

| ID | Group | Type | Backends | Figure | Programs |
|---|---|---|---|---|---|
| C2 | Contracts | C | CPU (gfortran, nvfortran, ifx; debug and optimised) | no | `seth_test_properties`, `seth_test_weno_js_eno`, `seth_test_wenoz`, `seth_test_teno`, `seth_test_registry` |

## Case and scope

The invariances every scheme must have whatever its data, checked on every scheme of the registry (24: the
reconstructions and the interpolations), on random stencils from a fixed, printed seed (`SETH_TEST_SEED` overrides
it).

## Setup and results

| Property | Check | Tolerance |
|---|---|---|
| Exactness | polynomials of degree below the narrowest candidate are reproduced whatever the weights (cell averages for reconstructions, point values for interpolations) | 64 ε of the data |
| Constants | a constant stencil returns the constant | 16 ε |
| Translation | $R(v + c) = R(v) + c$ (the indicators are blind to constants) | $10^{-11}(|c| + \max|v|)$ |
| Mirror | $v_m(m) = v_p(1-m)$ gives $r_m = r_p$ | bitwise for WENO-Z and TENO; $10^{-11}$ for WENO-JS (mirrored tables summed in another order) |
| Power-of-2 scaling | $R(2^k v) = 2^k R(v)$ | bitwise for WENO-Z, TENO and scale-invariant WENO-JS; plain WENO-JS is asserted **not** scale invariant (its absolute ε) |
| Robustness | the zero stencil, and magnitudes from $10^{-200}$ to each scheme's limit, give finite states within the data's range ($10^{140}$; plain WENO-JS $10^{30}$, as $(\varepsilon + IS)^{w}$ overflows near $10^{38}$–$10^{77}$) | debug builds trap overflow, invalid, zero divide |
| ENO | a unit step at every cell boundary of the stencil, both orientations, stays in $[0, 1]$ | $10^{-6}$ (WENO-JS, WENO-Z: the crossing candidates keep weights of $O(\varepsilon^{w})$); $10^{-12}$ (TENO: they are cut) |
| Large jump | a step of $10^{30}$ stays finite and within $[0, 10^{30}]$ (TENO's normalised $\chi$ avoids $\gamma^6$ overflow) | round-off |
| Linear scheme | TENO on smooth data returns $\sum_k d_k q_k$ | 4 ε |
| Registry | every kernel dispatched, names unique, `_wr` kernels give the plain kernels' states bitwise, weights ratios in $[0, 1]$ | bitwise |

All hold for every scheme, the six WENO-Z and TENO interpolations included.
