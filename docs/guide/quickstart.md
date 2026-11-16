# Quick Start

```fortran
use seth
real(R8P) :: vp(-2:3), vm(-2:3), rp, rm

! the values around the face i+1/2, cells by offset from it (cell 0 left, cell 1 right):
! vp reconstructed with a left bias (e.g. f+), vm with a right bias (e.g. f-); the same array for both to get the
! two states of one field
call seth_weno5_js(vp, vm, SETH_WENO_JS_EPS, 1._R8P, rp, rm)   ! rp: left state, rm: right state
```

The WENO-JS kernels take `(vp, vm, eps, rmu, rp, rm)`, with bounds `1-S:S` for WENO-(2S-1); the other families
replace `rmu` by their own parameter:

| Kernels | Scheme |
|---|---|
| `seth_weno1_js`, `seth_weno3_js`, `seth_weno5_js`, `seth_weno7_js`, `seth_weno9_js` | WENO-JS reconstruction from cell averages |
| `..._interp` | WENO interpolation from point values (Jiang, Shu & Zhang 2013) |
| `..._wr`, `..._interp_wr` | the same, also returning the weights ratios `wrp`, `wrm` (1 on smooth data, toward 0 near a discontinuity) |
| `seth_wenoz3`, `seth_wenoz5`, `seth_wenoz7`, `seth_wenoz9` | WENO-Z (Borges et al. weights; Don & Borges' tau_7, Castro's tau_9): `(vp, vm, eps, p, rp, rm)`, `p` = 1 (Z3, Z5, Z7) or 2 (Z9): [why](/guide/schemes#weno-z-weno-z3-5-7-9) |
| `seth_teno5`, `seth_teno5_opt`, `seth_teno6`, `seth_teno6_opt`, `seth_teno6_neno6` | TENO (Fu, Hu & Adams): `(vp, vm, eps, ct, rp, rm)`, bounds `-2:3`, `eps = SETH_TENO_EPS`, `ct = SETH_TENO5_CT`, `SETH_TENO5_OPT_CT`, `SETH_TENO6_CT`, `SETH_TENO6_OPT_CT`, `SETH_TENO6_NENO6_CT` |

- **No tables, no setup:** the coefficients are compile-time constants, so a kernel needs nothing on the device.
  Every kernel is `pure`, `!$acc routine seq` and `!$omp declare target`: call it from inside your own loops.
- **`rmu`** is the inverse descaler of the smoothness indicators: `1` for the Jiang-Shu weights, or
  `seth_weno_descaling(S, 1._R8P, vp, vm)` for scale-invariant weights (Don et al. 2022), or a magnitude of your own.
- **`eps`**: `SETH_WENO_JS_EPS` is ADAM's value. WENO-7 and WENO-9 keep Balsara & Shu's integer smoothness
  indicators (`WENO7_JS_D_SCALE = 240`, `WENO9_JS_D_SCALE = 5040`), so the same `eps` is relatively smaller there.

## Tests

```bash
fobis fetch
fobis build --mode tests-gnu && bash scripts/run_tests.sh
fobis rule --ex makecoverage-analysis   # coverage into docs/public/coverage.json
```

The tests check the WENO-JS kernels against ADAM's original code (bitwise without fused multiply-add), every kernel's
order on smooth data and boundedness across a jump, TENO's exact return to its linear scheme on smooth data, and the
exact mirroring of the right state. Where each WENO-Z and TENO formula comes from: `notes/m2-sources.md`.
