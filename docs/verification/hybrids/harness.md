# Hybrids in a solver

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| H3 | Hybrids | P | CPU (gfortran, nvfortran, ifx) | — | `seth_vv_hybrid` |

## Case and scope

The pieces of [H1](./sensors) and [H2](./central) assembled in the 1-D Euler harness (characteristic Lax–Friedrichs
flux splitting, SSP-RK3). Li, Chen & Ren's sensor on the density flags faces, dilated by two faces each side; flagged
faces take the nonlinear scheme, the others

- the **cost hybrid**: the linear background of the scheme (`seth_linear*`) on component-wise split fluxes, without
  the characteristic projection and without weights;
- the **KEEP-PE hybrid**: the sixth-order central flux `seth_keep_pe6` of the primitive state.

The harness's hybrids are test code, not library API: the switch, the buffer and the variable fed to the sensor are a
solver's policy (roadmap issue #4, steps 3.4 and 4.3).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Smooth | the density wave of [A4](../accuracy/euler-order), $N = 20 \dots 160$; the sensor must flag nothing | [`seth_vv_hybrid.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_hybrid.F90), [`seth_vv_euler.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_euler.F90) |
| J1 junction | the same wave, KEEP-PE on the left half of the faces and TENO5 on the right, whatever the sensor | |
| J3 conservation | 64 periodic cells, KEEP-PE or TENO5 at random faces, 50 steps: the totals of mass, momentum, energy | |
| Shocks | Sod and Lax (200 cells, exact solution), Shu–Osher (400 cells, 16× WENO5-JS reference, the left end a fixed supersonic inflow); TENO5, WENO-Z5, TENO6, each pure and in both hybrids | |
| Asserted | no flagged face on the smooth wave; the hybrids at the order of their smooth branch; J1 L1 order $\ge 4$; J3 to $10^{-13}$; the cost hybrid within 5% of the pure scheme on Sod and Lax for TENO5 and WENO-Z5; the KEEP-PE hybrid positive | |

## Results

<!--@include: ../_tables/hybrid.md-->

Smooth wave: no face flagged; the cost hybrid of TENO5 at order 5.01, the KEEP-PE hybrid at 6.00. J1: L1 order 4.99,
$L_\infty$ order 4.98 (the junction costs nothing beyond the fifth order of the TENO5 side). J3: mass, momentum and
energy conserved to $2\cdot10^{-15}$.

## Findings

- **The cost hybrid works for the fifth-order schemes**: with about 11% of faces flagged on Sod and Lax and 5% on
  Shu–Osher, TENO5 and WENO-Z5 stay within 1.5% of their pure $L_1$ (WENO-Z5's Shu–Osher error even drops 4%).
- **The plan's 5% holds for TENO6 except on Lax**, where it is 5.4% off; on Shu–Osher it is 4.4% off, the entropy
  waves behind the shock (about 8 cells per wavelength at $N = 400$) going to TENO6's nearly central smooth branch.
  The case prints these, it does not assert them.
- **The KEEP-PE hybrid is as good as the pure scheme on Sod and Lax, and on Shu–Osher** 0.91 (TENO5), 0.87
  (WENO-Z5) and 1.07 (TENO6) times its error: the non-dissipative flux on the smooth faces lowers the fifth-order
  schemes' error and costs TENO6 7%.
- **A central flux needs a well-posed inflow.** An earlier version of this page reported 1.81 (TENO6's cost hybrid)
  and 1.4 to 2.6 (the KEEP-PE hybrid) times the pure error on Shu–Osher, and blamed the entropy waves. It was the
  boundary: with a transmissive (zero-gradient) left end at Shu–Osher's supersonic inflow, the central flux, which
  reads cells downstream, perturbs the first cell, and the ghost cells copy it into the inflow state, which drifts
  from round-off to 0.1–0.15 in density by $t = 1.8$ (the upwind pure schemes damp it). The harness now holds the
  inflow state fixed (`BC_INFLOW`); the drift, and its spread between compilers and precisions, are gone.
- The junction of a central and an upwind flux, and random switching, cost neither order nor conservation.
