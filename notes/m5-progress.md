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
| 1.3 | the published static tests | `seth_test_nonpoly` spec 19: (a), (b) no activation; (c), (d) every jump flagged by DR and nothing else, N = 100 and 400; (e) reported (it matches the paper's Fig. 5: the jumps at 1.17 and 1.56 missed) — the frozen spec for (e) left to the owner | |
| 2.1, 2.2 | correctness (specs 1–5 of 2.2; 2.1 with Takagi's guard) | `seth_test_nonpoly` specs 7–11 on 10^6 triples, R8P and R4P (an R4P-only failure of the first gate was a test bug: the bound's tolerance) | |
| 2.2 | spec 7 (G7: Sod, Lax at N = 100) | `seth_vv_nonpoly`: TENO5-DR L1(rho) 3.27e-3 / 8.27e-3 against TENO5's 5.47e-3 / 1.61e-2; 3 / 2 cells across the contact against 4 | |
| 2.3 | N1, N2 | `seth_test_nonpoly` spec 12–13: bitwise equal to the two-pass stored-mask reference on nine fields; inert on smooth ones | |
| 2.3 | G5, G7, G8, G9 | the registry: properties (mirror, dimensional scaling), every shock case (report tables), device cases 52–53, golden v8 / r4_v3 (old lines identical); ten builds green | |
| 2.4 | E1 (PERF1, PERF2) | `notes/m5-bench-data/40-nonpoly-*`, `docs/guide/performance.md`: DR 2.72 (GPU) / 2.87 (gfortran) / 3.05 (nvfortran) times TENO5; THINC 3.46 / 3.74 / 3.48 (**over 3**); stored masks 1.65–1.94; mask 0.81 / 0.87 / 0.53 (budget 0.70). Stored-mask variant implemented (harness, bench) and compared on Sod, Lax, Shu–Osher, blast (`seth_vv_nonpoly`); the choice left to the owner | |

| 2.5 | E2, E3, E4 (N4, N6) | `seth_vv_nonpoly` section 3–4: conservation 1.8e-14 with 4312 / 7303 overrides; two blocks with 6 ghost cells bitwise, with the sixth unfilled different; E2 contact 4–5 cells (DR) against 12–15 (TENO5) after 10 periods; E3 gradients >= 2 cells wide bitwise TENO5's | |
| 3.3 | E5 | `seth_vv_nonpoly` section 6: Shu–Osher scaled by 1/8, TENO5 bitwise equal, TENO5-RBF 0.47 in density | |
| 3.3 | E6 (100 periods) | `seth_vv_nonpoly` section 7: inside the range (s = 0.075) the growth rate equals the prediction (6.6e-5 / 6.8e-5); outside, up to 4 times it; the amplitude nearly triples at s = 0.3, CFL 0.05 | |
| 3.4 | the complete scheme, N5 | `seth_teno5_rbf_dr`; device case 56; bench 39; the ablation of seven schemes on Sod, Lax, Shu–Osher, the blast waves (report mode) | |
| 3.1, 3.2, 3.4 | PERF1, PERF2 | `notes/m5-bench-data/41-nonpoly-*`: ST 0.88 / 0.86 / 0.97 × TENO5 (GPU / gfortran / nvfortran; budget 1.05); RBF **1.27** / 0.91 / 1.04 (1.25); RBF+DR 2.97 / 2.77 / **3.14** (3.0); bench_compare ok | |
| 3.1 | correctness (spec 14, rows of spec 16) | `seth_test_nonpoly`: order 5.008 / 4.993; equal to `seth_linear5` bitwise on smooth data; 4.000 beside a jump where TENO5 gives 3.011; overshoot 2e-16; the rows exact | |
| 3.2 | correctness (specs 1–4, 6 partly) | `seth_test_nonpoly` specs 15–16: s = 0 is ST bitwise; order 5.0 with alpha fixed; Eq. 29's error constant within 0.03 %; the kernel equal to a reference (sensor, row) on 10^5 stencils | |
| 3.2, 3.3 | G6, E6 (the band, D7) | spec 18: band edges 0.1640 / 0.3266 against the predicted 0.1643 / 0.3286, peaks within 4 %; validated range s <= 0.123 (TENO6's 2e-4) | |
| 3.3 | E7 | spec 17: the mask (0,1,1) only on noise (0.05 %), never on smooth data, a jump or a kink | |
| 3.1–3.4 | G9 (ST, RBF, RBF+DR) | `seth_test_golden`: golden v9 / r4_v4 add TENO5-ST, TENO5-RBF, TENO5-ST-DR and TENO5-RBF-DR (rbf_dx = 0.01) on the 64 stencils; v8's and r4_v3's 2691 lines their first lines, character for character | |

| 1.1–3.4 | G10 (mutants) | `notes/mutation-run-2026-11-17-nonpoly.txt`: the 38 mutants of issue #6 (mask bits, every term of Eqs. 24, 32 and 20, the branches of zeta, the steepness levels, the guards with their sign and rmu, THINC's beta/2, both mirrors, the early exit, the row by mask, the s^2 rows' sign and table, the sensor's c_r, nodes and s2) all killed | |

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
- **The RBF correction grows long waves** (step 3.2, spec 18): the band theta^2 < 30 s^2 and its peak (200/3) s^6
  are confirmed at small s; the validated range is s <= 0.12, and the paper's Shu-Osher run (s = 0.15 to 0.30) lies
  outside it.
- **TENO5 does not cut at a kink** (E7): on 10^7 kink stencils every mask is (1,1,1).
- **The overrides do not steepen smooth gradients** (E3): a gradient two cells wide or wider is reconstructed
  bitwise as by TENO5; a first version of E3 used a single tanh step on a periodic domain, whose wrap-around jump
  made the composites look better on smooth data: corrected before any page used it.
- **Deviation from a frozen spec (step 3.2, spec 1)**: "s1 = s2 = 0 bitwise equal to seth_teno5_st" holds with
  gfortran and ifx; nvfortran -fast reassociates the row sum, 61 % of 10^5 stencils bitwise and 2.6 ulps at worst, so
  the test asserts 4 eps (found by the step-09 gate, 2026-11-17).
- **What each piece buys** (N5): the DR override the gain on the shock tubes (Sod halved), the RBF term the gain on
  Shu–Osher (halved) at s = 0.15, outside its validated range; stencil targeting alone less than 1 %.
- **The RBF scheme depends on the length unit** (E5): 0.47 in density on Shu–Osher between two units; TENO5 bitwise
  equal (the rescaling must be a power of 2: by 10 the rounding of dx/10 flips a TENO cutoff and moves the profile by
  1e-4).
- **The RBF growth is visible in time only with a small time step** (E6): at CFL 0.4 SSP-RK3's dissipation hides
  part of it; at CFL 0.05 the measured rate matches Eq. 29's lambda term inside the validated range.
