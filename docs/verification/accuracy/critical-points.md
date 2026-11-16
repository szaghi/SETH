# Critical points

| ID | Group | Type | Backends | Figure | Programs |
|---|---|---|---|---|---|
| A2 | Accuracy | P | CPU (gfortran, nvfortran, ifx) | no | `seth_test_critical_points` (reconstructions), `seth_test_interpolation` (interpolations) |

## Case and scope

Where $f' = 0$ the smoothness indicators of the candidates tend to zero together, and the nonlinear weights of
WENO-JS deviate from the optimal ones by $O(1)$ relative amounts (Henrick et al. 2005): the scheme loses order. The
Z-type and TENO weights were designed against that. The case measures the order at first-order ($f' = 0$,
$f'' \ne 0$) and second-order ($f' = f'' = 0$) critical points and compares it with the orders the literature
reports.

What the case catches: a weight formula that is right on generic data but wrong at critical points, a wrong $\tau$
(its order fixes the order at critical points), a wrong $p$ or $\varepsilon$.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Reconstructions | the setting of Henrick et al. (2005) and Borges et al. (2008): point values $f_i = f(ih)$, the critical point at node 0, the derivative $(F_{1/2} - F_{-1/2})/h$ from the kernels' face fluxes | [`seth_test_critical_points.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_test_critical_points.F90) |
| Functions | Tang et al.'s $e^x - \sum_{k \le n} x^k/k!$, critical point of order $n$ at 0 | idem |
| Interpolations | the face value at a critical point placed on the face: $f = g - g'(0)x$ (first order), $f = g - g'(0)x - g''(0)x^2/2$ (second order), $g$ the generic function of [A1](./order-of-accuracy) | [`seth_test_interpolation.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_test_interpolation.F90) |
| $\varepsilon$ | $10^{-40}$ (WENO-Z, TENO), as Borges et al. | |
| Order | finest pair above round-off | |

## Programs

```bash
fobis build --mode tests-gnu
./exe/seth_test_critical_points ; ./exe/seth_test_interpolation
```

## Results

Reconstructions, the derivative at the critical point, against the literature (asserted within 0.5 where the literature is settled):

| Scheme | first order: literature | SETH | second order: literature | SETH |
|---|---|---|---|---|
| WENO-JS, $r = 2..5$ | $\max(2r-2-n_{cp}, n_{cp})$ = 1, 3, 5, 7 | 1.0, 3.0, 5.0, 6.6–7.2 | 2, 2, 4, 6 | 2.0, 2.0, 4.0, 5.9–6.1 |
| WENO5-Z, $p = 1$ | 4 (3.93–3.99) | 4.00 | 2 | 2.06 |
| WENO5-Z, $p = 2$ | 5 (4.80–4.98) | 5.0 | — | 2.02 |
| TENO5 | 5 | 5.0 | contested | 2.0 |
| TENO6 / TENO6-opt | 6 / 5 | 6.0 / 5.0 | contested (Fu 2023: none; Zhang 2024b: 2) | 2.0 |

SETH's measured orders for every reconstruction (the last report run):

<!--@include: ../_tables/critical-points.md-->

Interpolations, the face value (asserted at first-order points for WENO-Z with $p = 2$ and TENO at their design
order):

<!--@include: ../_tables/interpolation-orders.md-->

## Findings

- Every reconstruction order agrees with the literature where the literature is settled.
- At second-order critical points SETH's TENO measures 2, which sides with Zhang et al. (2024b) against Fu (2023);
  the test asserts only a floor there, since the literature disagrees.
- No published analysis of interpolation at critical points was found: the interpolation orders are measured, and
  asserted only for WENO-Z ($p = 2$) and TENO at first-order points, where they keep the design order as the
  reconstructions do. TENO5 and TENO6 interpolations keep 5 and 6 even at second-order points, unlike their
  reconstructions (the face value does not divide by $h$).
- Second-order critical points of the *split fluxes* matter in practice: with $\alpha = \max|f'|$ exactly, the
  Lax–Friedrichs split flux $f^- = (f - \alpha u)/2$ has one where $|f'|$ peaks, and TENO falls to second order
  there. The nonlinear harness splits with $\alpha = 1.1\max|f'|$ ([Burgers](../shocks/burgers)).
- In R4P the asymptotic orders at critical points cannot be observed: the errors fall below R4P's round-off on the
  grids where R8P measures them. The coarse grids R4P can use carry R8P's errors exactly, and their orders are
  pre-asymptotic (WENO-Z7, $p = 2$: 4.0 at $h = 0.4 \to 0.2$ against 7). TENO5-AZ's error at the second-order point
  drops from $4\cdot10^{-4}$ to $2\cdot10^{-10}$ as its cutoff switches, below R4P's round-off: Tang's 5 is verified
  in R8P only ([test programs](../unit-tests#single-precision)).
