# Interpolation for WENO-Z and TENO: implementation plan

Status: plan (2026-11-15). The owner decided to add the feature (decision 2 of the verification night). Suggested slot:
M5, after the M4 items in `kernel-api-design.md` section 7, or before them if a client needs it first.

## 1. What is missing

SETH reconstructs face values from cell averages with every scheme. It interpolates face values from point values
only with WENO-JS (`seth_weno{1,3,5,7,9}_js_interp[_wr]`, ADAM's `initialize_interpolation_S*` tables, Jiang, Shu &
Zhang 2013). WENO-Z and TENO have no `_interp` kernels. Clients that need them:

- the alternative finite-difference WENO (Shu–Osher's second form; Jiang, Shu & Zhang 2013; Balsara, Bhoriya & Shu)
  interpolates point values to the faces and then applies any Riemann solver;
- WCNS-type schemes (Deng & Zhang; WCNS-Z; WCNS-T of Hiejima 2022) interpolate primitive or characteristic variables
  from the nodes to the faces.

## 2. Design: what changes and what does not

A Z/TENO kernel is a descriptor (candidate cells, face coefficients `P`, indicator factors `L`/`DQ` and global
`LG`/`DG`, linear weights `D`, tau coefficients `TC`/`TG`) instantiated into the shared body
`seth_z_teno_body.INC`. **An interpolation kernel differs only in `P` and `D`.** The body, the mirroring
(`st(m) = vm(1-m)`: the cell centres `m-1/2` are symmetric about the face just as the cells are) and the indicators
stay the same. The work is mostly generator, tests and docs; the Fortran is about 120 lines of instantiation.

### 2.1 Face coefficients and linear weights (derived, exact)

Point value of cell `m` at its centre `x = m - 1/2`, face at `x = 0`. Candidate `k` is the Lagrange interpolant of its
points evaluated at the face; `D` are the weights for which the combination is exact to the highest degree (the
generator's `max_order_weights` on point-value moments). Computed exactly (`notes/interpolation-plan-evidence.txt`):

| Scheme | Interpolation linear weights | Order | Cross-check |
|---|---|---|---|
| WENO S = 2..5 (Z3/5/7/9) | 1/4, 3/4; 1/16, 5/8, 5/16; 1/64 ... 7/64; 1/256 ... 9/256 | 3, 5, 7, 9 | equal to ADAM's `initialize_interpolation_S*` weights **and** face coefficients, all S |
| TENO5 | 5/8, 5/16, 1/16 | 5 | equal to Hiejima 2022 (WCNS-T, K = 5) and Chu, Kurganov & Na 2022 (llm-wiki weno-teno `weighted-compact-nonlinear-scheme.md`) |
| TENO6 | 15/32, 5/16, 1/32, 3/16 | 6 | **derived only**: no published TENO6 interpolation weights found (WCNS-T uses odd stencils). All positive, so no negative-weight splitting is needed |

### 2.2 Smoothness indicators: reuse the reconstruction ones (recommended)

There are two candidate choices.

- **(a) Reuse the reconstruction indicators.** These are the Jiang–Shu integrals of the averages' reconstruction
  polynomial: the same quadratic forms, applied to point values. This is ADAM's choice for WENO-JS interpolation,
  and that of Chu, Kurganov & Na 2022 and Shu 2009, eq. 2.9.
- **(b) Use the Jiang–Shu integrals of the point-value interpolant's derivatives** (Zhou, Ding, Huang & Luo 2023).

For width 2 and 3 the two choices coincide. For width 4 and above they differ: the interpolant and the reconstruction
polynomial differ by `h^2/24 q''`, and from width 4 on `q''` is not constant. The exact tau orders were measured on
point values of a generic smooth function, with generic rational Taylor coefficients and no exponential:

| tau | (a) on point values | (b) on point values | (a) on averages (shipped) |
|---|---|---|---|
| Z3, Z5, Z9, TENO5 | 3, 5, 7, 5 | 3, 5, 7, 5 | 3, 5, 7, 5 |
| Z7 (Don–Borges) | **7** | 5 | 7 |
| TENO6 | **6** | 4 | 6 |

The tau combinations were built for the reconstruction indicators' expansions. With (b), Z7 and TENO6 lose two orders
each; recovering them would need new tau coefficients. Zhou et al. 2023 do exactly that for their six-point candidate,
with tau6 = |beta3 - (beta0 + 6 beta1 + beta2)/8|. **(a) keeps every tau order and keeps the body and the indicator
tables unchanged.** Hiejima's Deng–Zhang difference indicators are a third choice, not needed.

### 2.3 API

The new kernels follow the reconstruction signatures exactly, with the same constants (`SETH_TENO*_CT`, `SETH_TENO_EPS`):

```fortran
call seth_wenoz<3|5|7|9>_interp(vp, vm, eps, p, rp, rm)
call seth_teno<5|6>_interp     (vp, vm, eps, ct, rp, rm)
```

- No `_wr` variants: the Z/TENO reconstruction kernels have none.
- **TENO6-opt and NENO6: no `_interp`.** Their linear weights are tuned for the reconstruction's dispersion, and on
  interpolation candidates they have no published basis. An interpolation-tuned set would be new research, not a port.
- Additions only: no existing kernel, table value or golden value changes.

## 3. Steps and gates

Run as one night (night protocol), one snapshot per step.

1. **Generator** (`gen_z_teno_tables.py`).
   - Add a point-value basis, emit `<NAME>_INTERP_P` and `<NAME>_INTERP_D` for Z3/5/7/9, TENO5 and TENO6.
   - New cross-checks:
     - the WENO interpolation face coefficients and weights equal ADAM's `initialize_interpolation_S*`;
     - TENO5 = (5/8, 5/16, 1/16);
     - every interpolation weight is positive and reaches the design order;
     - the tau orders on point values.
   - **Hardening found while planning:** check every tau order (reconstruction included) on a generic function
     instead of `exp`.
     - For exp, cell averages are a constant multiple of point values, so exp cannot distinguish averages from
       point values at all.
     - It also cancels terms of the form f'f'''' - f''f''', which can overstate an order by one.
     - Re-measured on a generic function, every shipped order holds (Z3 3, Z5 5, Z7 7, Z9 7, TENO5 5, TENO6 6), so
       nothing shipped is wrong.
   - *Gate:* the generator passes; the regenerated tables differ from the committed ones only by added arrays
     (`git diff` shows no `-` line in the existing arrays); `fobis build` passes for the gnu, nvf and intel tests modes.
2. **Kernels.**
   - Add six kernels to `seth_wenoz.F90` and `seth_teno.F90`, instantiating the same body with the `_INTERP_`
     constants, and re-export them from `seth.F90`.
   - *Gate:* all three compilers build; the `tests-nvf-oac` build compiles them as `acc routine seq`; the existing
     15 tests stay green (no regression: additions only).
3. **Tests.** The registry `seth_test_schemes` gets the six schemes (an `interp` attribute; the families stay Z and
   TENO, so the family-driven expectations apply), and every registry-driven test picks them up. Then:
   - **order**: the design order at the face on point values of a generic function, both states;
   - **properties**:
     - exact mirroring;
     - no overshoot across a jump;
     - TENO returns its linear interpolation on smooth data to round-off;
     - a 1e30 jump stays finite;
   - **critical points**: pointwise interpolation error at first-order critical points, asserted against the
     literature where it exists (WCNS-Z expansions, Chen & Deng 2023) and otherwise printed and logged as measured,
     not published;
   - **whole-domain convergence**: the Linf interpolation error over every node of Henrick's field. The
     flux-difference operator test does not apply: differencing interpolated point values is not the alternative
     scheme, which is why `FAM_JS_INTERP` is already excluded;
   - **golden**: a new `seth_golden_v2.txt` covering every scheme. v1 is never edited. Gate: the v2 entries of the
     v1 schemes are bitwise equal to v1 (same compiler and flags), then the test's default path moves to v2;
   - **device**: regenerate `seth_test_device` to cover the six kernels on the GPU;
   - **tables**: re-verify the `_INTERP_` tables at run time;
   - **mutation**: add three mutants (a TENO6 interpolation weight, a Z5 interpolation face coefficient, an
     interpolation kernel instantiated with the reconstruction `P`). All must be killed.
   - *Gate:* the 15 tests (plus any new program) are green on the six CPU modes and fpm; the mutation run kills every
     mutant.
4. **Docs.** Update `api.md`, `schemes.md` (the interpolation section, the derived TENO6 weights flagged as derived,
   the indicator choice and its evidence), `verification.md` and `changelog.md`.
   - *Gate:* `fobis rule --ex makedoc` builds; `fobis scaffold status --strict` passes.
5. **Device run.** Run `tests-nvf-oac` on WSL GPU1, after checking `nvidia-smi` shows it free.
   - *Gate:* the device test is green.

Estimate: one night. The risk sits in step 3's critical-point expectations, where the literature is thin, not in the
kernels.

## 4. Decisions (owner, 2026-11-15: 1 reuse the reconstruction indicators, 2 scope agreed, 3 new golden v2, 4 version deferred)

1. **Indicators:** (a) reuse the reconstruction ones (recommended, section 2.2), or (b) interpolant-derivative
   indicators with new tau coefficients (research, not a port).
2. **Scope:** Z3/5/7/9 and TENO5/6 interpolation; TENO6-opt and NENO6 excluded (recommended).
3. **Golden:** a v2 file superseding v1, with v1 kept as history and v1 ⊂ v2 verified bitwise (recommended); or
   one golden file per feature, read together.
4. **Version:** additions only; the bump is the owner's call (M1, M2 and the test suite are already unreleased on
   top of v0.0.1).
