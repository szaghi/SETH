---
title: Kinds and constants
---

# Kinds and constants

## Kinds

`use seth` gives PENF's kinds: `I4P`, `I8P` (integers of 4 and 8 bytes), `R4P`, `R8P` (reals of 4 and 8 bytes) and
`R_P`, the **working real kind**. Every real argument of every kernel is `real(R_P)`.

`R_P` is fixed when the library is built: R8P by default, R4P with the `*-r4` modes (`-D_R_P_IS_R4P`, which PENF and
SETH must both see; a build where they disagree stops with a division by zero in a constant expression). One
precision per build: no generic interfaces, the kernels are compiled once, for that kind.

```bash
fobis build --mode static-gnu       # R_P = R8P
fobis build --mode static-gnu-r4    # R_P = R4P (also static-nvf-r4, static-intel-r4, shared-*-r4)
```

What R4P changes, what it is verified against, where it stops, and what it costs and buys: [Precision](/guide/precision).

## Constants that depend on the kind

The literals whose meaning depends on the kind live in `seth_kind_parameters`, one set per kind:

| Constant | R8P | R4P | Why |
|---|---|---|---|
| `SETH_TENO_EPS` | $10^{-40}$ | $10^{-30}$ | the papers' $10^{-40}$ is below R4P's smallest normal number ($1.2\cdot10^{-38}$), and a GPU flushes it to zero; pass it as `eps` to WENO-Z and TENO in either kind |
| `SETH_NONPOLY_GUARD` | $10^{-15}$ | $10^{-15}$ | the papers' first-order guard of the DR and THINC overrides: a product of two differences of the data, not a step against 1 |

`seth` re-exports both. The positivity limiter's back-off step is the third kind-dependent literal, internal to the
limiter.

## Every other constant

Each family's page lists its constants with their values: [WENO-JS](./weno-js#constants), [TENO](./teno#the-kernels),
[Non-polynomial TENO](./nonpolynomial#the-composite-kernels), [Sensors](./sensors), [Limiters](./limiters).

## Coefficient tables

Every coefficient of every scheme is a named constant, written by a generator script:

| Module | Generator | Content |
|---|---|---|
| `seth_weno_js_tables` | `scripts/gen_weno_js_tables.py` | WENO-JS's (ADAM's) candidate, weight and indicator coefficients (`WENO5_JS_A`, `WENO5_JS_P`, `WENO5_JS_D`, ...) |
| `seth_z_teno_tables` | `scripts/gen_z_teno_tables.py` | the descriptors of WENO-Z, TENO and the linear schemes: candidate cells (`*_CLO`, `*_CW`), face coefficients (`*_P`), indicator factors (`*_L`, `*_DQ`, `*_LG`, `*_DG`), linear weights (`*_D`), $\tau$ coefficients (`*_TC`, `*_TG`), linear rows (`LIN<N>_C`) |
| `seth_nonpoly_tables` | `scripts/gen_nonpoly_tables.py` | the rows of stencil targeting and of the RBF correction, DR's steepnesses |

`seth` re-exports the `WENO<N>_JS_D_SCALE` constants of `seth_weno_js_tables`; the other tables are used through
their own modules. What is asserted about them: [C3](/verification/contracts/tables).
