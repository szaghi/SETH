# About

SETH (**S**moothness-weighted **E**NO, **T**ENO and **H**ybrids) reconstructs and interpolates with the weighted
essentially non-oscillatory family of schemes: WENO, TENO and the hybrids built on them.

It is extracted from the WENO machinery of ADAM and succeeds WenOOF, whose run-time polymorphism in the inner loops
cost too much on HPC hardware. In the Genesis story Seth is the son of Adam who carries the line on after Abel: here
WenOOF is Abel.

## Status

The repository is being seeded: the build, the tests and the CI are in place, and the reconstruction machinery is
being moved in from ADAM.
