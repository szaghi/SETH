---
layout: home

hero:
  name: SETH
  text: Smoothness-weighted ENO, TENO and Hybrids
  tagline: High-order non-oscillatory reconstruction and interpolation for HPC, in modern Fortran
  actions:
    - theme: brand
      text: Guide
      link: /guide/
    - theme: alt
      text: API Reference
      link: /api/
    - theme: alt
      text: View on GitHub
      link: https://github.com/szaghi/SETH

features:
  - title: WENO, TENO and hybrids
    details: One library for the smoothness-weighted families of shock-capturing schemes, open to the newer approaches of the same class.
  - title: Built for HPC
    details: Concrete kernels, scheme chosen once outside the hot loop, no dynamic dispatch where the cells are swept.
  - title: CPU and GPU
    details: Tested with gfortran, nvfortran and ifx; OpenACC builds with nvfortran.
---
