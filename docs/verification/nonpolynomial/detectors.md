# Detectors

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| NP2 | Non-polynomial | P, C | CPU (gfortran, nvfortran, ifx), device | the static tests at $N = 400$ | `seth_test_nonpoly`, `seth_test_device` |

## Case and scope

The two detectors read the flags ([NP1](./flags)) of the cells around a face through $\zeta$ (Liang, Xiao, Shyy & Fu
2022, Eq. 17; Takagi et al. 2022, Eq. 19): 2 where the downwind and the upwind candidates are both cut, 1 where only
the downwind one is, $-1$ where only the upwind one is, 0 otherwise. A discontinuity inside the cell $i$ shows as a
cell upwind with $\zeta \in \{1, 2\}$, the cell itself with $\zeta = 0$ (or 1), a cell downwind with
$\zeta \in \{-1, 2\}$:

| Detector | Kernel | Cells read | Conditions |
|---|---|---|---|
| THINC | `seth_thinc_detect` | $i-3..i+3$ | Takagi et al. 2022, Eq. 20 (Liang et al. 2022, Eq. 18) |
| DR, level 1 | `seth_dr_detect` | $i-3..i+2$ | Liang et al. 2022, Eq. 24 (two conditions) |
| DR, level 2 | `seth_dr_detect` | $i-2..i+2$ | Eq. 32 as well (each of its conditions implies one of Eq. 24's) |

This case verifies the detectors (steps 1.2 and 1.3 of [issue #6](https://github.com/szaghi/SETH/issues/6)).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Exhaustive | every combination of the masks each detector reads, $8^6$ (DR) and $8^7$ (THINC), both states, against a reference written in the test from the printed conditions as set membership | [`seth_test_nonpoly.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_test_nonpoly.F90) |
| Tables 1–3 | the flags of the three profiles (and two mirrored in value): DR level 2, 1, 2 at cell 5 and 0 elsewhere; THINC flags profile 1 only | |
| Static tests | Liang et al. 2022, §3.1, point values at the cell centres of $[0, 2]$, constant continuation beyond the ends, the left state of each face: (a) $\sin 10\pi x$ and (b) $\sin 24\pi x$ at $N = 120$, no activation (the paper's statement); (c) Jiang and Shu's waves and (d) the piecewise function, $N = 100$ and 400: every jump flagged by DR within 1.5 cells, every DR activation within 1.5 cells of a jump of the function or of its derivative; (e) $N = 200$, 400: reported only | |

## Program

```bash
./exe/seth_test_nonpoly
scripts/run_vv.sh --only seth_test_nonpoly
```

## Results

The DR detector equals the printed Eqs. 24 and 32 on all $8^6 = 262144$ combinations of the masks it reads, the
THINC detector Eq. 20 on all $8^7$; on Tables 1–3 they give the levels Liang et al. derive. Where each fires on the
static tests:

<!--@include: ../_tables/nonpoly-static.md-->

![Liang et al.'s static tests (c), (d), (e) at N = 400, where the DR and THINC detectors fire](/verification/nonpoly-static.png)

On the device (`seth_test_device`, case 57) the two detectors and $\zeta$, on the masks of 4096 stencils computed
there, equal the host's bit for bit, and the test asserts that each detector fires on some of them (DR on 264, THINC
on 162).

## Findings

- **The DR detector finds every jump of tests (c) and (d)**, at both resolutions, and nothing else; it does not flag the
  kinks of test (c) (the triangle's peak, which the paper's figure shows flagged), and TENO5's own cutoff does not cut
  at a kink ([NP5](./targeting), E7). The THINC detector misses the jumps of (c) and, at $N = 400$, of (d): the
  defect Liang et al. built the DR detector to correct.
- **Test (e) behaves as the paper's Fig. 5**: at $N = 400$ the DR detector fires at the jumps up to $x = 1.53$ and misses
  the two smallest, at 1.17 and 1.56, as the paper's does; it also fires in the under-resolved sine near the right
  end. The frozen specification of step 1.3, which required every jump flagged and nothing else on (e), contradicts
  the paper's own figure: (e) is reported, not asserted, until the owner decides.
- On flat data the flags, and so the detectors, follow rounding noise ([NP1](./flags)): activation counts in uniform
  regions are not meaningful.
