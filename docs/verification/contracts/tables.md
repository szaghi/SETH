# Coefficient tables

| ID | Group | Type | Backends | Figure | Programs |
|---|---|---|---|---|---|
| C3 | Contracts | C | Python (exact rationals); CPU | no | `scripts/gen_weno_js_tables.py`, `scripts/gen_z_teno_tables.py`, `seth_test_tables` |

## Case and scope

Every coefficient SETH uses is generated, not typed, and checked twice: in exact rational arithmetic by the generator
before it writes anything, and again in double precision on the compiled constants, so that a hand edit of a
generated file, an emission or rounding bug, or a stale table cannot ship unnoticed.

## Setup

The generators check, in exact arithmetic:

| Check | Source of truth |
|---|---|
| every candidate exact for polynomials up to its degree, for cell averages (reconstruction) and for point values at the cell centres (interpolation) | the definition |
| the TENO candidates' face coefficients | Li, Fu & Adams (2020), Table 1 |
| the Jiang–Shu indicators of the 3-, 4-, 5-point stencils | ADAM's WENO-5/7/9 tables (scales 1, 240, 5040) |
| the max-order linear weights, exact to the design order | ADAM's WENO tables ($S = 2 … 5$); Fu et al.'s TENO6 (0.45, 0.3, 0.05, 0.2) |
| the interpolation coefficients and weights of $S = 2 … 5$ | ADAM's `initialize_interpolation_S*` |
| the TENO5 interpolation weights 5/8, 5/16, 1/16 | Hiejima (2022, WCNS-T) and Chu, Kurganov & Na (2022) |
| the TENO6 interpolation weights 15/32, 5/16, 1/32, 3/16 | derived: positive, 6th order (no published set) |
| the order of every $\tau$, on cell averages and on point values of a generic polynomial (not exp) | Z3 3, Z5 5, Z7 7, Z9 7, TENO5 5, TENO6 6 on both |
| the $LDL^T$ factors reproduce every indicator exactly (indicators are sums of squares: $\beta \ge 0$ by construction) | the indicator |

`seth_test_tables` re-checks the compiled constants: candidates and linear schemes exact, weights summing to 1 and
positive, indicators 0 on constants and 1 (or the Balsara–Shu scale) on unit slope, $\tau$ of its order on data of
the generic function, for the reconstruction and the interpolation tables.

## Findings

- Using the reconstruction indicators for the interpolations keeps every $\tau$ order. The alternative, integrating
  the derivatives of the point-value interpolant, drops WENO-Z7's $\tau$ from 7 to 5 and TENO6's from 6 to 4
  (`notes/interpolation-plan.md`).
- With $e^x$ as the test function the $\tau$ orders cannot distinguish averages from point values. The generator now
  measures them on a generic polynomial.
- Two TENO mutants (a $\tau$ coefficient, an indicator pivot) once passed every other test. They change only TENO's
  cut decisions, which no test stencil sits near, and `seth_test_tables` exists because of them
  ([mutation testing](./mutation)).
