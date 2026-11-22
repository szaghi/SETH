---
title: WENO-JS
---

# WENO-JS

Jiang and Shu's weighted ENO of order $2S-1$, for $S = 1$ to 5, extracted from ADAM and bitwise equal to ADAM's
former kernels ([C1](/verification/contracts/adam-equivalence)), which ADAM replaced with these (ADAM issue 139).
Module `seth_weno_js`; what the scheme is: [Schemes](/guide/schemes#weno-js-weno-1-3-5-7-9).

## Calls

```fortran
call seth_weno<N>_js          (vp, vm, eps, rmu, rp, rm)             ! N = 1, 3, 5, 7, 9
call seth_weno<N>_js_wr       (vp, vm, eps, rmu, rp, rm, wrp, wrm)
call seth_weno<N>_js_interp   (vp, vm, eps, rmu, rp, rm)
call seth_weno<N>_js_interp_wr(vp, vm, eps, rmu, rp, rm, wrp, wrm)
```

| Argument | Declaration | Meaning |
|---|---|---|
| `vp`, `vm` | `real(R_P), intent(in) :: v(1-S:S)`, $S = (N+1)/2$ | left- and right-biased values (cells $1-S..S-1$ read by `vp`, $2-S..S$ by `vm`) |
| `eps` | `real(R_P), intent(in)` | regulariser of the smoothness indicators: `SETH_WENO_JS_EPS` |
| `rmu` | `real(R_P), intent(in)` | inverse descaler of the indicators: `1` (Jiang–Shu), or a [descaler](#descalers) |
| `rp`, `rm` | `real(R_P), intent(out)` | left and right states at the face |
| `wrp`, `wrm` | `real(R_P), intent(out)` | weights ratios $\min_k \omega_k/d_k \in [0, 1]$: 1 on smooth data, toward 0 near a jump |

## The kernels

| Kernel | $S$ | Bounds | Order | Input |
|---|---|---|---|---|
| `seth_weno1_js` | 1 | `0:1` | 1 | cell averages |
| `seth_weno3_js` | 2 | `-1:2` | 3 | cell averages |
| `seth_weno5_js` | 3 | `-2:3` | 5 | cell averages |
| `seth_weno7_js` | 4 | `-3:4` | 7 | cell averages |
| `seth_weno9_js` | 5 | `-4:5` | 9 | cell averages |
| `seth_weno1_js_interp`, `seth_weno3_js_interp`, `seth_weno5_js_interp`, `seth_weno7_js_interp`, `seth_weno9_js_interp` | 1 to 5 | `1-S:S` | $2S-1$ | point values at the cell centres (Jiang, Shu & Zhang 2013) |
| `seth_weno1_js_wr`, `seth_weno3_js_wr`, `seth_weno5_js_wr`, `seth_weno7_js_wr`, `seth_weno9_js_wr` | 1 to 5 | `1-S:S` | $2S-1$ | as the plain kernel, also returning `wrp`, `wrm` |
| `seth_weno1_js_interp_wr`, `seth_weno3_js_interp_wr`, `seth_weno5_js_interp_wr`, `seth_weno7_js_interp_wr`, `seth_weno9_js_interp_wr` | 1 to 5 | `1-S:S` | $2S-1$ | as the interpolation kernel, also returning `wrp`, `wrm` |

WENO1 is the first-order upwind value (one candidate, weight 1): it is there to complete the family. API:
[`seth_weno_js`](/api/src/lib/seth_weno_js).

## Constants

| Constant | Value | Meaning |
|---|---|---|
| `SETH_WENO_JS_EPS` | $10^{-6}$ | ADAM's regulariser, absolute, on the indicators as the tables scale them |
| `WENO1_JS_D_SCALE`, `WENO3_JS_D_SCALE`, `WENO5_JS_D_SCALE` | 1 | the indicator of unit-slope data: these indicators are normalised |
| `WENO7_JS_D_SCALE` | 240 | WENO7 keeps Balsara & Shu's integer indicators |
| `WENO9_JS_D_SCALE` | 5040 | WENO9 likewise |

With the integer indicators of WENO7 and WENO9, the same `eps` is relatively smaller there: divide by the
`D_SCALE` to compare. The `D_SCALE` constants live in `seth_weno_js_tables` and are re-exported by `seth`.

## Descalers

```fortran
rmu = seth_weno_descaling(S, sigma, vp, vm)  ! sigma = 0: rmu = 1 exactly; sigma = 1: 1/mean|v| (Don et al. 2022)
rmu = seth_weno_descaling_ai(S, vp, vm)      ! affine invariant: 1/(mean deviation from the mean) (Wang & Don 2022)
```

| Argument | Declaration | Meaning |
|---|---|---|
| `S` | `integer(I4P), intent(in)` | candidates of the scheme, WENO-$(2S-1)$ |
| `sigma` | `real(R_P), intent(in)` | 0: Jiang–Shu weights; 1: scale-invariant weights |
| `vp`, `vm` | `real(R_P), intent(in) :: v(1-S:S)` | the face's two stencils, as passed to the kernel |
| result | `real(R_P)` | `rmu`, to pass to the kernel |

The indicators are quadratic in the data, and `eps` is absolute: without a descaler the weights depend on the
magnitude of the data. With $r_\mu = 1/\mu$ the indicators are those of the data divided by $\mu$:

- `seth_weno_descaling` ($\mu$ the mean $|v|$) makes the result scale bitwise under $v \to 2^k v$, but fails for a
  small fluctuation on a large offset (a pressure, a density): $\mu$ is then the offset and the scheme oscillates;
- `seth_weno_descaling_ai` ($\mu$ the mean deviation from the stencil's mean) ignores the offset too, up to the
  round-off of the expanded indicators, which $\mu^{-2}$ amplifies: not enough for WENO7/9-JS;
- the recommended use for data with an offset: subtract the stencil's mean from `vp` and `vm`, call the kernel with
  `seth_weno_descaling_ai`, add the mean back to the states (Wang & Don's Remark 2). No overshoot with any compiler,
  for one more pass over the stencil.

The measurements: [Schemes](/guide/schemes#weno-js-weno-1-3-5-7-9). API: [`seth_descaling`](/api/src/lib/seth_descaling).
