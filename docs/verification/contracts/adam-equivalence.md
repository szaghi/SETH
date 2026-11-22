# Equivalence with ADAM

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| C1 | Contracts | C | CPU (gfortran, nvfortran, ifx; debug and optimised) | no | `seth_test_weno_js_adam` |

## Case and scope

SETH's WENO-JS kernels are extracted from ADAM, and ADAM calls them (milestone M3, ADAM issue 139): FLUME on the CPU
and with OpenACC, PRISM's CPU backend and RIPPLE (ADAM commits dede5114, 73ac233a, b1683349, b4398d8c), with every
ADAM regression case bitwise the same on both backends. They had to give ADAM's results, so the test compares them
with ADAM's own routines, copied verbatim into `src/tests/seth_adam_reference.F90`, on the same data. ADAM has since
deleted those routines (ADAM commit ba37c4cb): the copy is frozen history, kept as the fixed reference of this
contract.

## Setup

| Item | Value |
|---|---|
| Schemes | WENO-JS, $S = 1 … 5$, reconstruction and interpolation tables |
| Weights | Jiang–Shu ($\mu = 1$) and the scale-invariant descaler (Don et al. 2022), `seth_weno_descaling` |
| Data | smooth and rough stencils at magnitudes $10^{-6} … 10^{6}$, fixed seed |
| Mapping | ADAM's rows: `v(2,m) = vp(m)`, `v(1,m) = vm(m+1)`; `rp = r(2)`, `rm = r(1)` |
| Compared | both states, both weights ratios, the descaler; a state's difference against the magnitude of its data |
| Results | about 1.4 million |

## Results

| Build | Bitwise equal | Worst difference |
|---|---|---|
| gfortran (debug and -O3), nvfortran debug | 100% | 0 |
| nvfortran -fast | 91% | $2.0\cdot10^{-15}$ of the data |
| ifx | 98.5% | idem |

## Findings

- Without fused multiply-adds SETH is ADAM, bit for bit: the kernel body keeps ADAM's operation order on purpose.
- With them, the compilers contract SETH's and ADAM's copies differently; the differences stay within 9 ulp of the
  data.
