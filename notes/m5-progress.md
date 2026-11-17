# M5 progress: the non-polynomial reconstructions (issue #6)

One line per gate met, as issue #6 asks: step, what passed, where the evidence is. Commits are the owner's (night
protocol: each step is a snapshot in `.night/<night>/`, committed in the morning); the commit column is filled when
they land.

| Step | Gate | Evidence | Commit |
|---|---|---|---|
| 0.1 | primary reads: every P2 of the plan resolved to P1 or corrected | `notes/m5-sources.md` (42 claims checked; five corrections to the plan) | |
| 0.2 | design note | `notes/kernel-api-design.md` §11 | |
| 0.3 | harness `NG = 6`, registry `-5:6`; every result unchanged | the suite (37/37) and `run_vv.sh --no-figures` on v0.2.0 and on the widened tree, data directories identical (`diff -r`) | |
| 0.4 | generator, aborting on mismatch | `scripts/gen_nonpoly_tables.py`: 36 rationals identified, 6 rows = Eq. 25, 6 truncation errors = Eq. 29, the pairs, Eq. 31, Eq. 16 to 1e-60; a mutated coefficient of Eq. 25 aborts it (nothing written); `--check` up to date | |
| 0.5 | verification group | `docs/verification/nonpolynomial/` (index, eight cards); `npm run docs:build` green | |
| 1.1 | correctness (specs 1–4) | `seth_test_nonpoly`: masks = the cutoff of `seth_teno5` on 10^6 stencils; Tables 1–3 digit for digit (R4P with eps = 1e-12, see the finding); mirror; never 0. gfortran (debug, -O3), nvfortran -fast, ifx, gfortran R4P | |
| 1.2 | correctness | `seth_test_nonpoly`: both detectors equal the printed conditions on every mask combination (8^6, 8^7); the levels of Tables 1–3 | |
| 2.1, 2.2 | correctness (specs 1–5 of 2.2; 2.1 with Takagi's guard) | `seth_test_nonpoly` specs 7–11 on 10^6 triples, R8P and R4P (an R4P-only failure of the first gate was a test bug: the bound's tolerance) | |
| 2.2 | spec 7 (G7: Sod, Lax at N = 100) | `seth_vv_nonpoly`: TENO5-DR L1(rho) 3.27e-3 / 8.27e-3 against TENO5's 5.47e-3 / 1.61e-2; 3 / 2 cells across the contact against 4 | |
| 2.3 | N1, N2 | `seth_test_nonpoly` spec 12–13: bitwise equal to the two-pass stored-mask reference on nine fields; inert on smooth ones | |
| 2.3 | G5, G7, G8, G9 | the registry: properties (mirror, dimensional scaling), every shock case (report tables), device cases 52–53, golden v8 / r4_v3 (old lines identical); ten builds green | |
| 2.4 | E1 (PERF1, PERF2) | `notes/m5-bench-data/40-nonpoly-*`, `docs/guide/performance.md`: DR 2.72 (GPU) / 2.87 (gfortran) / 3.05 (nvfortran) times TENO5; THINC 3.46 / 3.74 / 3.48 (**over 3**); stored masks 1.65–1.94; mask 0.81 / 0.87 / 0.53 (budget 0.70). Stored-mask variant implemented (harness, bench) and compared on Sod, Lax, Shu–Osher, blast (`seth_vv_nonpoly`); the choice left to the owner | |

## Findings

- **Constant data does not give the mask 7** (step 1.1). The indicators of a constant that is not exact in binary are
  round-off, far above eps = 1e-40, and the cutoff drops candidates by the accident of rounding: gfortran gives a
  mask other than 7 on 65 % of 1000 constants in [-10, 10]; nvfortran -fast gives the mask 5 (zeta = 1) on every
  one. In R4P (eps = 1e-30) even data equal to 1 gives the mask 5, and Tables 1–3 are reproduced only with eps >=
  1e-14. States are unaffected; detector activations in uniform regions are not.
- **The DR gain on Sod is a splitting effect** (step 2.3, the registry): TENO5-DR's L1(rho) is 0.51 of TENO5's with
  flux splitting and 0.94 on the Roe path; on Lax it survives on both.
- **The basis of a neighbour's mask matters on systems** (E1): per face against stored masks, Sod 1.39e-3 / 1.24e-3,
  Lax 5.39e-3 / 5.83e-3, the blast waves' density differing by up to 1.9 locally; neither better everywhere.
- **THINC misses the cost budget of 3** on every device (3.46–3.74 times TENO5): three downwind masks per state.
