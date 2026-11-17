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

## Findings

- **Constant data does not give the mask 7** (step 1.1). The indicators of a constant that is not exact in binary are
  round-off, far above eps = 1e-40, and the cutoff drops candidates by the accident of rounding: gfortran gives a
  mask other than 7 on 65 % of 1000 constants in [-10, 10]; nvfortran -fast gives the mask 5 (zeta = 1) on every
  one. In R4P (eps = 1e-30) even data equal to 1 gives the mask 5, and Tables 1–3 are reproduced only with eps >=
  1e-14. States are unaffected; detector activations in uniform regions are not.
