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
| Shocks | Sod and Lax (200 cells, exact solution), Shu–Osher (400 cells, 16× WENO5-JS reference); TENO5, WENO-Z5, TENO6, each pure and in both hybrids | |
| Asserted | no flagged face on the smooth wave; the hybrids at the order of their smooth branch; J1 L1 order $\ge 4$; J3 to $10^{-13}$; the cost hybrid within 5% of the pure scheme on Sod and Lax for TENO5 and WENO-Z5; the KEEP-PE hybrid positive | |

## Results

<!--@include: ../_tables/hybrid.md-->

Smooth wave: no face flagged; the cost hybrid of TENO5 at order 5.01, the KEEP-PE hybrid at 6.00. J1: L1 order 4.99,
$L_\infty$ order 4.98 (the junction costs nothing beyond the fifth order of the TENO5 side). J3: mass, momentum and
energy conserved to $2\cdot10^{-15}$.

## Findings

- **The cost hybrid works for the fifth-order schemes**: with about 11% of faces flagged on Sod and Lax and 5% on
  Shu–Osher, TENO5 and WENO-Z5 stay within 1.5% of their pure $L_1$ (WENO-Z5's Shu–Osher error even drops 4%).
- **The plan's 5% does not hold for TENO6**, and the case says so instead of asserting it: 5.4% off on Lax, 1.81
  times the pure error on Shu–Osher. The sensor leaves the entropy waves behind the shock (about 8 cells per
  wavelength at $N = 400$) to the smooth branch, which for TENO6 is nearly the central scheme. A sensor that sees
  discontinuities, not resolution, is not enough when the smooth branch has little dissipation.
- **The KEEP-PE hybrid is worse still on Shu–Osher** (1.4 to 2.6 times the pure scheme's error) and as good as the
  pure scheme on Sod and Lax: a non-dissipative flux on under-resolved waves, which SSP-RK3 does not damp. This is the
  regime the roadmap's Phase 5 is about: the sensor, not the flux, decides.
- The junction of a central and an upwind flux, and random switching, cost neither order nor conservation.
