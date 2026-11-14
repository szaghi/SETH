# About

SETH (**S**moothness-weighted **E**NO, **T**ENO and **H**ybrids) reconstructs and interpolates with the weighted
essentially non-oscillatory family of schemes: WENO, TENO and the hybrids built on them.

It is extracted from the WENO machinery of ADAM and succeeds WenOOF, whose run-time polymorphism in the inner loops
cost too much on HPC hardware. In the Genesis story Seth is the son of Adam who carries the line on after Abel: here
WenOOF is Abel.

## Status

WENO-JS 1/3/5/7/9 (reconstruction and interpolation, extracted from ADAM and bitwise equal to it), WENO-Z3/5/7/9 and
TENO5/6 with the optimised TENO6 variants: one kernel per scheme, coefficients as compile-time constants, `pure` and
callable on the device. See [Schemes](./schemes), [API](./api), and how it is checked in [Verification](./verification).
