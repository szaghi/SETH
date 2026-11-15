# Verification

SETH's kernels are verified by 25 programs that are, at the same time, its test suite: they run on gfortran 16.2,
nvfortran 26.9 and ifx 2026.2, debug and optimised, and on the GPU, and each asserts its results against something it
shares no code with. The full documentation, case by case and with figures, is the
**[Verification](/verification/)** section; this page is the short version.

| Group | What is verified | Cases |
|---|---|---|
| [Accuracy](/verification/accuracy/) | design order of every reconstruction and interpolation at a point, at critical points (against Henrick, Borges, Fu, Zhang), over whole domains, and through a complete Euler solver | A1–A4 |
| [Shock capturing](/verification/shocks/) | Sod, Lax and 123 against the exact Riemann solution; Shu–Osher and Woodward–Colella against 16× finer references; Burgers against its exact solution with a moving shock; long-time advection of the Jiang–Shu profile; positivity with and without a limiter | S1–S7 |
| [Spectral](/verification/spectral/) | the modified wavenumber of every linear background; the approximate dispersion relation of every nonlinear scheme | F1–F2 |
| [Contracts](/verification/contracts/) | WENO-JS bitwise equal to ADAM; invariances and robustness; the coefficient tables re-derived and re-verified; golden values; the device; mutation testing | C1–C6 |

```bash
fobis build --mode tests-gnu && bash scripts/run_tests.sh   # the suite (the verification cases on their quick set)
scripts/run_vv.sh                                           # every scheme on the report grids: data, tables, figures
python3 scripts/mutation_test.py --mode tests-gnu           # does every planted defect make a test fail?
```

Findings that matter for users of the kernels, each documented on its page:

- The higher-order interpolations need a positivity safeguard near vacuum, and every interpolation on strong blast
  waves; on the Woodward–Colella blast so do WENO9-JS and the WENO-Z (but Z3) and TENO reconstructions
  ([S3](/verification/shocks/test123), [S5](/verification/shocks/blast)).
- WENO-Z9 is not robust at shocks: 28% overshoot behind the Burgers shock, oscillations on the composite profile,
  loss of positivity in the first steps of the blast wave ([S6](/verification/shocks/burgers)).
- With Lax–Friedrichs splitting, use $\alpha$ a little above $\max|f'|$. With $\alpha = \max|f'|$ exactly, the split
  flux has a second-order critical point where $|f'|$ peaks, and TENO falls to second order there
  ([S6](/verification/shocks/burgers)).
- TENO6 is anti-dissipative in a narrow band near $\vartheta \approx 0.47\pi$ and leaves grid-scale packets behind a
  shock undamped ([F2](/verification/spectral/adr), [S6](/verification/shocks/burgers)).
