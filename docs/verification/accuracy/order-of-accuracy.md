# Order of accuracy at a point

| ID | Group | Type | Backends | Figure | Programs |
|---|---|---|---|---|---|
| A1 | Accuracy | P | CPU (gfortran, nvfortran, ifx), GPU build compiled | no | `seth_test_weno_js_order`, `seth_test_wenoz`, `seth_test_teno` |

## Case and scope

At a single face, data of a smooth function sampled on cells of width $h$ halved six times: cell averages for the
reconstructions, point values at the cell centres for the interpolations. The error of each state against the exact
face value must fall at the design order: $2S-1$ for WENO-JS and WENO-Z, 5 and 6 for TENO5 and TENO6, and 5 and 4
for the optimised TENO6 weights (Fu et al.'s set loses one order, NENO6's two).

What the case catches:

- **a wrong coefficient** in a face polynomial or a linear weight, which drops the order;
- **a mismatch between data kind and kernel**, a reconstruction fed point values or the reverse;
- **an order that only holds by accident**: the test function is generic (below), so no cancellation hides an
  error term.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Test function | $g(x) = \sin(0.9x + 0.4) + 0.5\,e^{0.6x} + 0.3\cos(1.7x - 1.1)$, face at $x = 0.3$ ($g' > 1$ there: no critical point) | [`seth_test_functions.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_test_functions.F90) |
| Cell averages | 10-point Gauss–Legendre on each cell: no cancellation, error far below round-off | `gavg` |
| Grids | $h_0 = 0.4$ ($0.8$ for $S = 5$), six halvings | the three programs |
| Order | the finest pair whose finer error is above $100\,\varepsilon\,\vert g\vert$; asserted $\ge$ design $- 0.5$ | idem |
| Kernels | every reconstruction and interpolation, `eps` = 1e-6 (WENO-JS), 1e-40 (WENO-Z, TENO), WENO-Z with $p = 1$ and 2 | [`seth_weno_js.F90`](https://github.com/szaghi/SETH/blob/main/src/lib/seth_weno_js.F90), [`seth_wenoz.F90`](https://github.com/szaghi/SETH/blob/main/src/lib/seth_wenoz.F90), [`seth_teno.F90`](https://github.com/szaghi/SETH/blob/main/src/lib/seth_teno.F90) |

Why $g$ and not $e^x$: the cell averages of $e^x$ are $\sinh(h/2)/(h/2)$ times its point values, so the two data
kinds cannot be distinguished; and $e^x$ satisfies $f'f'''' = f''f'''$, which cancels leading terms of the
nonlinear weights' errors. The wavenumbers of $g$ stay below 2. With a mode of wavenumber 2.9, the 9th-order schemes,
whose stencils span 10 cells, were still pre-asymptotic when their error reached round-off, with ratios swinging
between 7.5 and 12.7.

## Programs

```bash
fobis build --mode tests-gnu
./exe/seth_test_weno_js_order ; ./exe/seth_test_wenoz ; ./exe/seth_test_teno
```

## Results

Observed orders (gfortran 16.2, `tests-gnu`, the last report run):

<!--@include: ../_tables/order-of-accuracy.md-->

On smooth data every TENO kernel also returns its linear scheme, $\sum_k d_k q_k$, to 0 ulp (asserted within
$4\varepsilon$).

## Findings

- Every kernel reaches its design order on both states, the six WENO-Z and TENO interpolations included.
- The 9th-order kernels needed the generic function with low wavenumbers. With the first choice of $g$ (a mode of
  wavenumber 2.9), WENO9-JS and WENO-Z9 measured 8.2–8.7 above round-off, and $e^x$ had let them pass. The error was
  in the test, not the kernels: their error reached round-off before the asymptotic range.
- The round-off floor is $100\varepsilon|g|$. It was $1000\varepsilon$, which discarded pairs whose errors, near
  $10^{-13}$, were still 500 times above round-off.
