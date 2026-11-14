# SETH

**S**moothness-weighted **E**NO, **T**ENO and **H**ybrids: high-order non-oscillatory reconstruction and
interpolation for HPC, in modern Fortran.

[![CI](https://github.com/szaghi/SETH/actions/workflows/ci.yml/badge.svg)](https://github.com/szaghi/SETH/actions/workflows/ci.yml)
[![Matrix](https://github.com/szaghi/SETH/actions/workflows/matrix.yml/badge.svg)](https://github.com/szaghi/SETH/actions/workflows/matrix.yml)

SETH is extracted from the WENO machinery of ADAM and succeeds WenOOF: concrete kernels, the scheme chosen once
outside the hot loop, no dynamic dispatch where the cells are swept. Tested with gfortran, nvfortran and ifx;
OpenACC builds with nvfortran.

```bash
fobis fetch && fobis build --mode tests-gnu && bash scripts/run_tests.sh
```

Documentation: https://szaghi.github.io/SETH. License: GPL-3.0.
