# Equivalence with ADAM

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| C1 | Contracts | C | CPU (gfortran, nvfortran, ifx; debug and optimised) | no | `seth_test_weno_js_adam`, `seth_test_weno6_point_adam` |

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

## The seam WENO-6 transfer

`seth_weno6_point` and `seth_weno6_jumps` ([WENO-6 point interpolation](/reference/weno6-point)) are ADAM's 2:1 AMR
seam transfer (ADAM issue 68, C3 and D2), moved into SETH by issue #7, E2; ADAM's step S7 (ADAM issue 139) replaces
its own routines with them. `seth_test_weno6_point_adam` compares them with ADAM's, copied verbatim from ADAM commit
2ea71524 into `src/tests/seth_adam_seam_reference.F90`: built by gfortran every result must be bitwise equal.

| Item | Value |
|---|---|
| Tables | the 169 constants SETH compiles at the targets ADAM uses (linear and cubic weights at the midpoints and the quarters, the regulariser, the jump ratio) against ADAM's `SEAM_D_WENO6`, `SEAM_W_CUBIC_MID`, `SEAM_D_FILL_WENO6`, `SEAM_W_CUBIC_Q`, `SEAM_WENO6_EPS`, `SEAM_WENO6_JUMP_RATIO` |
| Mapping | `seam_restriction_weno6_1d(u, code)` = `seth_weno6_point(u, 4*code + 6)`, $code = 0..4$; `seam_fill_weno6_1d(u, sub, p)` = `seth_weno6_point(u, 4*p + 2*sub - 3)`, $p = 1..6$, $sub = 1, 2$; `seam_weno6_jumps` = `seth_weno6_jumps` |
| Data | 100000 lines of six values, seven kinds (smooth samples, independent random values, steps with noise, constants, a flat sub-stencil next to a ramp, ADAM's `rnd**8`, small variations on a large offset), magnitudes $10^{-12} … 10^{12}$, fixed seed |
| Results | 1.7 million values (17 targets per line) and the 100000 lines' jump flags |

| Build | Tables | Values bitwise equal | Jump flags |
|---|---|---|---|
| gfortran 16.2, -O3 and debug | 169 of 169 | 100% | all equal |
| nvfortran 26.9, -O4 -fast (CPU) | 169 of 169 | 100% | all equal |

The kernel keeps ADAM's operation order, and the generator writes every weight as ADAM writes it (an integer over 16,
128 or 320), so each constant folds to ADAM's double; it derives the weights from their definition and asserts they
equal ADAM's four tables (`scripts/gen_weno6_tables.py`). On the device: [C5](./device), case 64.
