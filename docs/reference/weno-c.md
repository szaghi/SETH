---
title: Centred WENO
---

# Centred WENO

ADAM's centred WENO (`WENO_C_4/6/8`), ported with its defect D1 fixed: WENO-JS of order $2S-1$ blended with the
central reconstruction of order $2S$ on all $2S$ cells of the stencil. Module `seth_weno_c`; the blend:
[Schemes](/guide/schemes#centred-weno-weno4-6-8-c).

## Calls

```fortran
call seth_weno4_c(vp, vm, eps, rmu, rp, rm)   ! vp, vm(-1:2)
call seth_weno6_c(vp, vm, eps, rmu, rp, rm)   ! vp, vm(-2:3)
call seth_weno8_c(vp, vm, eps, rmu, rp, rm)   ! vp, vm(-3:4)
```

The arguments are those of [WENO-JS](./weno-js): `eps` > 0 (`SETH_WENO_JS_EPS`) and `rmu` (1, or a descaler) go to
the upwind weights and to the blend's indicators.

| Kernel | $S$ | Bounds | Upwind part | Order |
|---|---|---|---|---|
| `seth_weno4_c` | 2 | `-1:2` | WENO3-JS | 4 |
| `seth_weno6_c` | 3 | `-2:3` | WENO5-JS | 6 |
| `seth_weno8_c` | 4 | `-3:4` | WENO7-JS | 8 |

The design orders hold at a point; at critical points the orders are the upwind WENO-JS's. WENO8-C is not ENO for a
jump in its two outer cells (an overshoot of $\tfrac1{280}$ of the jump): use it on smooth data, or with the
positivity limiter ([A5](/verification/accuracy/centred-weno)). API: [`seth_weno_c`](/api/src/lib/seth_weno_c).
