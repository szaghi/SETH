# Golden regression

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| C4 | Contracts | C | CPU (gfortran, nvfortran, ifx; debug and optimised) | no | `seth_test_golden` |

## Case and scope

The verification cases assert orders, properties and bounds. A change of a few ulp, or an error that grows inside its
bound, passes them. The golden test pins the values themselves.

## Setup

| Item | Value |
|---|---|
| Kernels | every scheme of the registry (WENO-Z also with $p = 2$), both states, and the weights ratios of WENO-JS |
| Stencils | 64, exact in binary (no libm, no random-number generator): 16 smooth cubics with dyadic coefficients, 20 jumps (every cell boundary, two shapes), 12 scaled ($2^{-498} … 2^{96}$), 16 rough from an integer LCG |
| File | `src/tests/golden/seth_golden_v6.txt`, 9472 values, written once by the gfortran debug build (its compiler and options in the header) |
| Tolerance | states: $10^{-13}$ of the stencil's magnitude; weights ratios: $10^{-9}$ (they read the weights themselves, which move by about $w\,\delta IS/(\varepsilon + IS)$ when the indicators round differently) |
| Rule | a golden file is never edited; a deliberate change of results gets a new file and the test's default path |

## Results

Measured on 2026-11-15 (the final gate of the night):

| Build | Bitwise equal | Worst state difference (of the data) | Worst weights-ratio difference |
|---|---|---|---|
| gfortran debug (the writer), gfortran optimised | 8192 / 8192 | 0 | 0 |
| ifx debug | 8192 / 8192 | 0 | 0 |
| ifx optimised | 7359 / 8192 | $4.4\cdot10^{-16}$ | $2.2\cdot10^{-16}$ |
| nvfortran debug / optimised | 6862 / 6853 of 8192 | $1.1\cdot10^{-15}$ | $3.5\cdot10^{-10}$ |

## Findings

- v2 (2026-11-15) adds the six WENO-Z and TENO interpolations. Its lines for the schemes of v1 are identical to v1's,
  character for character, checked when v2 was written with the same build; v1 is kept.
- The inputs are exact in binary: an earlier version built them with libm and differed between gfortran -O0 and -O3
  before any kernel ran.
