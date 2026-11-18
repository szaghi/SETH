# M5 progress: the non-polynomial reconstructions (issue #6)

One line per gate met, as issue #6 asks: step, what passed, where the evidence is. Commits are the owner's (night
protocol: each step is a snapshot in `.night/<night>/`, committed in the morning); the commit column is filled when
they land.

| Step | Gate | Evidence | Commit |
|---|---|---|---|
| 0.1 | primary reads: every P2 of the plan resolved to P1 or corrected | `notes/m5-sources.md` (42 claims checked; five corrections to the plan) | 6438a61 |
| 0.2 | design note | `notes/kernel-api-design.md` §11 | f9569ae |
| 0.3 | harness `NG = 6`, registry `-5:6`; every result unchanged | the suite (37/37) and `run_vv.sh --no-figures` on v0.2.0 and on the widened tree, data directories identical (`diff -r`) | b9a36ec |
| 0.4 | generator, aborting on mismatch | `scripts/gen_nonpoly_tables.py`: 36 rationals identified, 6 rows = Eq. 25, 6 truncation errors = Eq. 29, the pairs, Eq. 31, Eq. 16 to 1e-60; a mutated coefficient of Eq. 25 aborts it (nothing written); `--check` up to date | 0d2a3ba |
| 0.5 | verification group | `docs/verification/nonpolynomial/` (index, eight cards); `npm run docs:build` green | fb4e630 |
| 1.1 | correctness (specs 1–4) | `seth_test_nonpoly`: masks = the cutoff of `seth_teno5` on 10^6 stencils; Tables 1–3 digit for digit (R4P with eps = 1e-12, see the finding); mirror; never 0. gfortran (debug, -O3), nvfortran -fast, ifx, gfortran R4P | 9f7c24c |
| 1.2 | correctness | `seth_test_nonpoly`: both detectors equal the printed conditions on every mask combination (8^6, 8^7); the levels of Tables 1–3 | 9f7c24c |
| 1.3 | the published static tests | `seth_test_nonpoly` spec 19: (a), (b) no activation; (c), (d) every jump flagged by DR and nothing else, N = 100 and 400; (e) reported (it matches the paper's Fig. 5: the jumps at 1.17 and 1.56 missed) — the frozen spec for (e) left to the owner | fa76aa0 |
| 2.1, 2.2 | correctness (specs 1–5 of 2.2; 2.1 with Takagi's guard) | `seth_test_nonpoly` specs 7–11 on 10^6 triples, R8P and R4P (an R4P-only failure of the first gate was a test bug: the bound's tolerance) | 9f7c24c |
| 2.2 | spec 7 (G7: Sod, Lax at N = 100) | `seth_vv_nonpoly`: TENO5-DR L1(rho) 3.27e-3 / 8.27e-3 against TENO5's 5.47e-3 / 1.61e-2; 3 / 2 cells across the contact against 4 | 5e8ebf1 |
| 2.3 | N1, N2 | `seth_test_nonpoly` spec 12–13: bitwise equal to the two-pass stored-mask reference on nine fields; inert on smooth ones | 9f7c24c, 2d9f75f |
| 2.3 | G5, G7, G8, G9 | the registry: properties (mirror, dimensional scaling), every shock case (report tables), device cases 52–53, golden v8 / r4_v3 (old lines identical); ten builds green | 9f7c24c, 2d9f75f |
| 2.4 | E1 (PERF1, PERF2) | `notes/m5-bench-data/40-nonpoly-*`, `docs/guide/performance.md`: DR 2.72 (GPU) / 2.87 (gfortran) / 3.05 (nvfortran) times TENO5; THINC 3.46 / 3.74 / 3.48 (**over 3**); stored masks 1.65–1.94; mask 0.81 / 0.87 / 0.53 (budget 0.70). Stored-mask variant implemented (harness, bench) and compared on Sod, Lax, Shu–Osher, blast (`seth_vv_nonpoly`); the choice left to the owner | 5e8ebf1 |

| 2.5 | E2, E3, E4 (N4, N6) | `seth_vv_nonpoly` section 3–4: conservation 1.8e-14 with 4312 / 7303 overrides; two blocks with 6 ghost cells bitwise, with the sixth unfilled different; E2 contact 4–5 cells (DR) against 12–15 (TENO5) after 10 periods; E3 gradients >= 2 cells wide bitwise TENO5's | fa76aa0 |
| 3.3 | E5 | `seth_vv_nonpoly` section 6: Shu–Osher scaled by 1/8, TENO5 bitwise equal, TENO5-RBF 0.47 in density | fa76aa0 |
| 3.3 | E6 (100 periods) | `seth_vv_nonpoly` section 7: inside the range (s = 0.075) the growth rate equals the prediction (6.6e-5 / 6.8e-5); outside, up to 4 times it; the amplitude nearly triples at s = 0.3, CFL 0.05 | fa76aa0 |
| 3.4 | the complete scheme, N5 | `seth_teno5_rbf_dr`; device case 56; bench 39; the ablation of seven schemes on Sod, Lax, Shu–Osher, the blast waves (report mode) | fa76aa0 |
| 3.1, 3.2, 3.4 | PERF1, PERF2 | `notes/m5-bench-data/41-nonpoly-*`: ST 0.88 / 0.86 / 0.97 × TENO5 (GPU / gfortran / nvfortran; budget 1.05); RBF **1.27** / 0.91 / 1.04 (1.25); RBF+DR 2.97 / 2.77 / **3.14** (3.0); bench_compare ok | f4f1a93 |
| 3.1 | correctness (spec 14, rows of spec 16) | `seth_test_nonpoly`: order 5.008 / 4.993; equal to `seth_linear5` bitwise on smooth data; 4.000 beside a jump where TENO5 gives 3.011; overshoot 2e-16; the rows exact | fa76aa0 |
| 3.2 | correctness (specs 1–4, 6 partly) | `seth_test_nonpoly` specs 15–16: s = 0 is ST bitwise; order 5.0 with alpha fixed; Eq. 29's error constant within 0.03 %; the kernel equal to a reference (sensor, row) on 10^5 stencils | fa76aa0 |
| 3.2, 3.3 | G6, E6 (the band, D7) | spec 18: band edges 0.1640 / 0.3266 against the predicted 0.1643 / 0.3286, peaks within 4 %; validated range s <= 0.123 (TENO6's 2e-4) | fa76aa0 |
| 3.3 | E7 | spec 17: the mask (0,1,1) only on noise (0.05 %), never on smooth data, a jump or a kink | fa76aa0 |
| 3.1–3.4 | G9 (ST, RBF, RBF+DR) | `seth_test_golden`: golden v9 / r4_v4 add TENO5-ST, TENO5-RBF, TENO5-ST-DR and TENO5-RBF-DR (rbf_dx = 0.01) on the 64 stencils; v8's and r4_v3's 2691 lines their first lines, character for character | c7c3553 |

| 1.1–3.4 | G10 (mutants) | `notes/mutation-run-2026-11-17-nonpoly.txt`: the 38 mutants of issue #6 (mask bits, every term of Eqs. 24, 32 and 20, the branches of zeta, the steepness levels, the guards with their sign and rmu, THINC's beta/2, both mirrors, the early exit, the row by mask, the s^2 rows' sign and table, the sensor's c_r, nodes and s2) all killed | 0c430f6 |
| 4.3 | the 2-D harness (D8) | `seth_vv_2d` in the suite (N = 48, TENO5, TENO5-DR, RBF+DR): L R = I, positivity, conservation, the contact symmetric (configuration 3's symmetry reported: nvfortran breaks it at 2.6e-8) on ten builds | 97e4f5b |
| 4.3 | the 2-D fields (NP9) | report run N = 400, seven schemes, 21 runs positive, symmetric bit for bit (gfortran); `docs/verification/nonpolynomial/two-d.md`: the contact 2-3 cells (DR, THINC) against 6; ripples behind configuration 3's shock, 382-608 local extrema against TENO5's 22 | d9e334f |


### Completion night 2026-11-18 (`.night/2026-11-18-m6-complete`)

A gate-by-gate audit of the issue against the code (two read-only audits, every gap re-checked) found the gates below
unmet or met in another form; each line closes one. Commits: the night's steps (`night commit`).

| Step | Gate | Evidence | Commit |
|---|---|---|---|
| 2.1 | conservation by quadrature (frozen) | `seth_test_nonpoly` spec 9: ten-point Gauss-Legendre on each quarter of the cell, 3.5e-15 (over the whole cell ten points reach 4e-13: tanh's poles); saturated R4P faces counted and skipped (2 of 10^6) | b6088eb |
| 2.2 | specs 4, 5 (frozen) | `seth_test_nonpoly`: 1 + tanh(A) B > 0 on every monotone triple; the offset to 1e-12 |mu|, |mu| in [1, 10] (worst 8.7e-15) | b6088eb |
| 2.2 | G2 | `seth_test_tables`: SETH_DR_T1, SETH_DR_T2 within one ulp of tanh(0.91), tanh(1.23) (0 ulp) | b6088eb |
| 2.2 | spec 6 (Table 4) | `seth_vv_nonpoly` section 8, `nonpoly-pulse.dat`: bitwise TENO5 and zero overrides at N = 80, 160; at N >= 320 overrides on the 1e-15 tails (85 at 320), the solution within 7.8e-15 of TENO5's: asserted to NP_NOISE_FLOOR (deviation, left to the owner); L_inf within 0.3-0.9 % of the paper's from N = 160 | b6088eb |
| 2.3, 3.4 | N2 on Euler | `seth_vv_euler_order`: DR, THINC, RBF+DR against their base schemes on the density wave, N = 20..160: not bitwise (1.8e-15..1.8e-14, the constant acoustic fields' round-off), asserted to NP_NOISE_FLOOR (deviation, left to the owner) | b6088eb |
| 3.2 | spec 3 (all six rows) | `seth_test_nonpoly` spec 16: Richardson on the masks 7, 3, 5, 1, 2, 4 within 0.007-0.04 % of Eq. 29 | b6088eb |
| 3.1, 3.2 | G2 (rows) | `seth_test_tables`: every row zero outside its stencil and exact below its width; the pairs equal to q0/2 + q1/2 and 3 q0/4 + q2/4; D4's row; every s^2 row summing to zero | b6088eb |
| 3.4 | N1, N2, N4, N6 | `seth_test_nonpoly` spec 12-13 (RBF+DR against its two-pass reference, bitwise; inert against TENO5-RBF); `seth_vv_nonpoly` section 3 (conservation 1.8e-14, 6 ghost cells bitwise, 5 not) | b6088eb |
| 3.1, 3.2, 3.4 | G4 | `seth_test_critical_points` over all 38 schemes: the six non-polynomial ones at TENO5's 5.00 / 2.00; TENO5's literature bound asserted | ac61d7d |
| 3.1, 3.2, 3.4 | G5 | `seth_test_properties` over all 38 schemes (ST+DR, RBF+DR descaled); a 1e30 step at every boundary for every TENO scheme | ac61d7d |
| 3.1, 3.2 | G6 | `seth_vv_adr`: ST = linear5 for theta <= pi/4 (5.6e-16), never anti-dissipative; RBF's band asserted at s = 0.03, 0.06 (edges 0.157, 0.314; peaks 4.6e-8, 3.5e-6), at 0.15, 0.3 reported (peak 22 and 4 times the linear prediction) | ac61d7d |
| 2.2, 3.1, 3.2, 3.4 | G7 | the six schemes in `seth_vv_sod`, `_lax`, `_test123`, `_shu_osher`, `_blast`, `_burgers`, `_advection`: positive, within TENO5's recorded bounds (closest: DR on the blast, 96 %); ST's pairs cost no robustness behind the Burgers shock | ac61d7d, 213d389 |
| 1.2, 2.1, 2.2 | G8 | `seth_test_device` cases 57-59: the detectors, zeta, both overrides; bitwise on GPU 1; both detectors fire there | 88289e4 |
| 4.1 | PERF4 | `side_level`'s run-time extent made constant (bitwise); `notes/gpu-kernel-audit.md` row; `-Minfo=accel`: no device allocation (`notes/gpu-audit-data/minfo-nonpoly-2026-11-18.txt`) | 90deaa0 |
| 1.1-3.4 | PERF1-PERF3 | `notes/m5-bench-data/42-*`, `42b-*` (GPU, four quiet runs, ratios within 1 %), `docs/guide/performance.md`: mask one bias 0.41 / 0.44 / 0.28 (0.70); detectors 0.05 / 0.11 / 0.12 (**0.10**); overrides forced on every face DR 0.19 / 0.05 / 0.32 (**0.15**), THINC 0.48 / 0.32 / 0.73 (**0.30**); ST 0.88 / 0.87 / 0.98; RBF **1.27** / 0.91 / 1.06; DR 2.74 / 2.85 / **3.05**; THINC **3.49 / 3.72 / 3.49**; RBF+DR 2.97 / 2.76 / **3.14**; stored masks 1.6-2.1 (GPU / gfortran / nvfortran; misses left to the owner) | 5f3d13d |
| 2.4, 3.4 | PERF7 | `42-fields-*`: per-face DR / THINC from 2.7 / 3.4 (smooth) to 5.3 / 5.8 (a jump every 8 cells) on the GPU, 2.9 / 3.8 to 3.8 / 4.5 on one core; stored masks 1.7-2.0 throughout | 5f3d13d |
| 4.1 | PERF3 | `bench_compare.py`: ok on the gfortran core; on the GPU every old variant, ADAM's included, 4.1-4.6 % slower than the 2026-11-15 baseline in series 42 and 42b alike, and in the control of HEAD 538a471 the benchmark of the unchanged base commit 538a471, run in the same hour, is 4.4-4.7 % slower too (42c-head-gpu-run-*): the device, not the code | 5f3d13d |
| 1.1-3.4 | G10 | three mutants (the composite's mirror, RBF+DR's base states, a pair-row coefficient) and the re-anchored early exit: 4 of 4 killed (`notes/mutation-run-2026-11-18-nonpoly-completion.txt`); 72 mutants in all | 449c0a6 |
| 3.4 | N5 | `nonpoly-ablation.dat`: seven schemes, one at a time, on Sod, Lax, 123, Shu-Osher, the blast waves (L1, overshoot, contact cells, detector activations, limited faces) and, from `seth_vv_burgers` and `seth_vv_advection`, on Burgers and the composite profile | 213d389 |
| 0-4 | DOC1-DOC13 | `schemes.md` (the detector conditions, both override formulas, the sensor, D4, costs, round-off); the compiled reference example; the stencil-targeting figure; the card gallery; ten figures and five animations; the user findings; the stale counts corrected (40 / 37 programs, golden v9, 59 device cases, 72 mutants) | 213d389 |

### Owner decisions D9-D14 and step 2.6 (2026-11-18, `.night/2026-11-18-stored-masks`)

The owner took the decisions of the completion audit on 2026-11-18 (issue #6, section "Decisions found by the
completion audit"): D9, D10 (the round-off form of the bitwise specs), D11 (step 1.3 (e): counts reported), D13 (the
budget misses accepted as recorded) and D14 (step 1.2's signature `(-3:3)`) as recommended; D12: the self-contained
kernels stay the default and the stored-mask variant ships as a kernel pair as well (step 2.6, the lines below).

| Step | Gate | Evidence | Commit |
|---|---|---|---|
| 2.6 | N1, N2 | `seth_test_nonpoly` specs 12-13: `seth_teno5_face_masks` equals the test-code reference's stored masks bitwise and the three `*_stored` kernels equal the self-contained ones bitwise on the nine periodic fields; spec 20: on 10^5 stencils with independent left and right data (about 15000 DR detections per side), on constant data and a jump of BIG_JUMP, bitwise and finite | dada8fe |
| 2.6 | G8 | `seth_test_device` cases 60-63: the four kernels on the GPU, the masks of 61-63 computed there; bitwise | dada8fe |
| 2.6 | G9, G4-G6 | no new golden lines: on scalar data the pair is bitwise the self-contained kernels (N1), which golden v9 and the sweeps already pin | dada8fe |
| 2.6 | G7 (systems), E1 | `seth_vv_euler`'s stored-mask path through the library pair (ST+DR and RBF+DR dispatched to `seth_teno5_rbf_dr_stored`); `seth_vv_nonpoly` E1 with RBF+DR: both variants positive on Sod, Lax, Shu-Osher, the blast waves; DR and THINC rows and every other report file identical to the hand-built two-pass to 16 digits | 07bcad5 |
| 2.6 | PERF2, PERF3 | bench variants 35, 36, 46 through the library pair; series 43 (`notes/m5-bench-data/43-stored-*`, 4 GPU, 3 gfortran, 3 nvfortran runs, checks ok): DR / THINC / RBF+DR stored 1.84 / 1.86 / 2.11 x TENO5 on the GPU, 1.98 / 1.97 / 1.86 with gfortran (the hand-built two-pass of series 42: 1.82 / 1.83 / 2.09, 1.95 / 1.95 / 1.86), 1.95 / 1.95 / 2.02 with nvfortran against 1.64 / 1.64 / 1.70: 12 ns per face in `seth_teno5_face_masks` (interleaved A/B, finding below); bench_compare ok on the GPU (old kernels +2.7..3.1 %) and the gfortran core | b22c6d8 |
| 2.6 | PERF4 | `notes/gpu-kernel-audit.md` row; `-Minfo=accel`: the four routines `acc routine seq`, no device allocation (`notes/gpu-audit-data/minfo-nonpoly-stored-2026-11-18.txt`) | b22c6d8 |
| 2.6 | G10 | four mutants (pass 1's mirror, THINC's detector and RBF+DR's base states in pass 2, the right state's masks): 4 of 4 killed, each by `seth_test_nonpoly` specs 12 and 20 (`notes/mutation-run-2026-11-18-stored-masks.txt`); 76 mutants in all | da19177 |
| 2.6 | DOC1-DOC13 | the reference page's section and compiled example (`reference_nonpoly`, the pair on a 16-cell field, 17 of 17 faces bitwise), `schemes.md`, `notes/kernel-api-design.md` 11.4, the E1 table, the user finding (DOC8), the decisions D9-D13 written into the pages, the device and mutation counts | 1624b5b |
| 2.6 | G11 | green on all ten builds (R8P 40/40 x5, R4P 37/37 x5; candidate tree of steps 01-05) | 1624b5b |

## Findings

- **The stored-mask pair's first pass is 12 ns per face slower than the hand-built one with nvfortran on one core**
  (step 2.6, 2026-11-18). An interleaved A/B on one pinned core (the bench of 0bccbb2 against the pair's) and a
  pass-by-pass micro-benchmark: pass 2 costs the same (41-42 ns per face), pass 1 `seth_teno5_face_masks` 33 ns
  against 20 for the same two `seth_teno5_mask` calls at the call site. The same per-face sequence written at the
  call site with the mirror as five scalar assignments runs at 20; inside the library it does not (tried, reverted).
  The routine is not inlined and not vectorised (`-Minfo=all`). Cause not isolated; gfortran and the GPU show no
  difference. The pair still costs 1.95 times TENO5 there, against 3 to 3.5 for the self-contained kernels.

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
- **Two frozen specs do not hold bit for bit, at round-off** (2026-11-18): the Gaussian pulse of Table 4 (spec 6:
  zero overrides, bitwise TENO5) fails from N = 320, the detector firing on the 1e-15 tails; N2 on the Euler density
  wave fails at every N, on the round-off of the constant acoustic fields. Both now assert agreement to
  NP_NOISE_FLOOR (1e-12 in R8P) and print the bitwise result. Accepting the deviation, or wanting the kernels to
  ignore flags below a magnitude, is the owner's decision.
- **Outside its validated range the nonlinear RBF scheme is far more anti-dissipative than its row** (G6, 2026-11-18):
  the sensor raises s above s1 on a resolved wave's crests; at s = 0.15 the peak of Im theta~ is 1.6e-2, 22 times the
  linear prediction.
- **On the 123 problem the overrides raise the limiter's work** (N5, 2026-11-18): 120-154 limited faces against
  TENO5's 76, the overshoot up from 3e-3 to 5-6e-3; on the blast waves they lower it from ~500 to 70-110.
