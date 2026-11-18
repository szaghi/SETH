# Mutation testing

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| C6 | Contracts | C | CPU (gfortran) | no | `scripts/mutation_test.py` |

## Case and scope

A test suite is only as good as the defects it catches. The mutation run plants one defect at a time in a scratch
copy of the tree, rebuilds (with `--force-compile`: FoBiS does not rebuild the dependents of a changed cpp
`#include`, [FoBiS#212](https://github.com/szaghi/FoBiS/issues/212)), runs the suite, and records which tests fail. A
defect no test catches is a hole in the suite.

## Setup

| Mutant | What it changes |
|---|---|
| `js5-face-coefficient`, `js5-linear-weight`, `js7-indicator`, `js3-interp-weight` | one coefficient of a generated WENO-JS table |
| `js-exponent-s3`, `js-mirror`, `js-eps-ignored`, `js-weights-ratio` | the WENO-JS body: the weight exponent, the right-state mirroring, the use of eps, the weights ratio |
| `descaling` | the scale-invariant descaler |
| `teno6-linear-weight`, `teno6-tau`, `teno6-indicator-pivot`, `wenoz5-tau` | one coefficient of a WENO-Z/TENO descriptor |
| `z-exponent-p2`, `z-teno-mirror`, `teno-cutoff-10x`, `teno-power-4`, `teno-global-indicator-off` | the shared Z/TENO body |
| `teno6-interp-weight`, `wenoz5-interp-coefficient`, `teno5-interp-uses-recon-p` | the interpolations: a weight, a face coefficient, a kernel wired to the reconstruction tables |
| `positivity-pressure-quadratic`, `positivity-density-bound`, `positivity-pressure-floor`, `positivity-right-state` | the positivity limiter: the pressure quadratic, the density bound, the floor, a one-sided state |

```bash
python3 scripts/mutation_test.py --mode tests-gnu
```

## Results

2026-11-15, the final run with all 26 programs: **25 of 25 mutants killed** (`notes/mutation-run-2026-11-15-device-arrays.txt`, after the limiter lost its automatic arrays; before it
`notes/mutation-run-2026-11-15-positivity.txt`;
the earlier run of the 21 kernel mutants: `notes/mutation-run-2026-11-15-final.txt`).
2026-11-16, the centred WENO: **6 of 6 new mutants killed** (`notes/mutation-run-2026-11-16-weno-c.txt`): the three
central rows, the blend's `max`, its descaler, its first pair. Two candidates were left out as equivalent mutants:
the mirrored state's cells in reverse order, and the central row read unmirrored (the row is symmetric, and the
blend does not see the swap).
2026-11-17, the non-polynomial kernels of issue #6: **38 of 38 new mutants killed**
(`notes/mutation-run-2026-11-17-nonpoly.txt`): the mask's bit weights and stencil, every term of the detectors and
overrides (Liang et al. 2022, Eqs. 24 and 32; Takagi et al. 2022, Eq. 20), each branch of $\zeta$, the steepness
levels, both guards with their sign and `rmu`, THINC's $\beta/2$ and both mirrors, the early exit, stencil
targeting's row by mask, the RBF rows' sign and table, the sensor's $c_r$, nodes and $s_2$.
2026-11-18, the completion of issue #6: **3 of 3 new mutants killed**, and the early-exit mutant, re-anchored on the
constant loop extent of `side_level`, killed again (`notes/mutation-run-2026-11-18-nonpoly-completion.txt`): the
right state's mirrored data in `seth_teno5_dr`, TENO5's states in place of the RBF states inside `seth_teno5_rbf_dr`,
a coefficient of the pair row $(1,0,1)$ of `NP_RBF_P`. The same day, step 2.6, the stored-mask pair: **4 of 4
killed** (`notes/mutation-run-2026-11-18-stored-masks.txt`), each by `seth_test_nonpoly` (the pair against the
self-contained kernels, specs 12 and 20): the mirror of `seth_teno5_face_masks`, the DR detector in place of THINC's
in `seth_teno5_thinc_stored`, TENO5's states in place of the RBF states in `seth_teno5_rbf_dr_stored`, the left
state's masks read for the right state in `seth_teno5_dr_stored`. In all, 76 mutants, each caught by at least one
test.
Each line of the record names the tests that caught the mutant. The interpolation mutants are caught by the golden,
interpolation, tables, TENO, WENO-Z and properties tests, and now also by the nonlinear cases (Sod, Lax, 123,
Shu–Osher, blast, Euler order): a kernel wired to the wrong tables changes every downstream result. The four positivity-limiter
mutants are all caught by `seth_test_positivity`; the pressure-quadratic one also by `seth_vv_test123`.

## Findings

- The first runs had survivors: two TENO mutants (a $\tau$ coefficient, an indicator pivot) that change only cut
  decisions no stencil sat near. `seth_test_tables` was added and kills them.
- A mutant that edits a comment is an equivalent mutant (it once survived that way): the anchors must be in code.
- Without `--force-compile` a mutated `.INC` body stayed compiled into later builds and the mutants after it were
  misattributed; reported upstream as FoBiS#212.
