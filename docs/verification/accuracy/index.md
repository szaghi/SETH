# Accuracy

Does every kernel reach its design order, where the literature says it should, and keep the order the literature
says it keeps at critical points? Four cases answer this, from a single face to a complete Euler solver. All of them
use test functions without special structure. Until 2026-11-15 the pointwise tests used $e^x$, which is degenerate
for them: its cell averages are a constant multiple of its point values, so a reconstruction and an interpolation
cannot be told apart, and it satisfies $f'f'''' = f''f'''$ and similar identities that cancel error terms of the
nonlinear weights.

| Case | Type | What it verifies |
|---|---|---|
| [A1 Order of accuracy at a point](./order-of-accuracy) | P | design order $2S-1$ (TENO: 5, 6) of every reconstruction and interpolation, both states |
| [A2 Critical points](./critical-points) | P | the orders measured by Henrick, Borges, Fu and Zhang where $f' = 0$ and $f' = f'' = 0$ |
| [A3 Whole-domain convergence](./convergence) | P | $L_\infty$ over every node of Henrick's field: the operator, SSP-RK3 advection, the interpolations |
| [A4 Euler order of accuracy](./euler-order) | P | the design order through the complete 1-D Euler solver, splitting and interpolation paths |
