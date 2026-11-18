---
layout: home

hero:
  name: SETH
  text: Smoothness-weighted ENO, TENO and Hybrids
  tagline: High-order non-oscillatory reconstruction and interpolation for HPC, in modern Fortran. WENO and TENO kernels, one pure function each, on every CPU and GPU.
  image:
    light: /logo.svg
    dark: /logo-dark.svg
    alt: SETH, the Riemann fan opened like a folding fan on a crackle-glazed tile
  actions:
    - theme: brand
      text: Quick start
      link: /guide/quickstart
    - theme: alt
      text: Tutorial
      link: /tutorial/
    - theme: alt
      text: Reference
      link: /reference/
    - theme: alt
      text: Verification
      link: /verification/
    - theme: alt
      text: View on GitHub
      link: https://github.com/szaghi/SETH

features:
  - icon:
      light: /icons/schemes.svg
      dark: /icons/schemes-dark.svg
      width: 44
      height: 44
    title: Every scheme, one kernel each
    details: WENO-JS 1/3/5/7/9, WENO-Z 3/5/7/9, TENO5, TENO5-opt, TENO5-AZ, TENO5-A, TENO6, TENO6-opt, NENO6, TENO8 and TENO8-A, for reconstruction from cell averages and interpolation from point values; the centred WENO4/6/8-C; the non-polynomial TENO (DR and THINC overrides, stencil targeting, Gaussian RBF); and the pieces of a hybrid (linear schemes, sensors, central fluxes).
    link: /reference/
    linkText: Feature map
  - icon:
      light: /icons/hpc.svg
      dark: /icons/hpc-dark.svg
      width: 44
      height: 44
    title: Called from inside your loops
    details: Every kernel is pure, acc routine seq and omp declare target; the coefficients are compile-time constants, so nothing is set up or copied to the device. On an RTX 4070, WENO-5 reconstructs a face in 0.97 ns and WENO-7 in 1.76 ns, in double precision.
    link: /guide/performance
    linkText: Performance
  - icon:
      light: /icons/weights.svg
      dark: /icons/weights-dark.svg
      width: 44
      height: 44
    title: Exact coefficients
    details: The WENO-Z and TENO coefficients are derived in exact rational arithmetic, cross-checked against the papers and re-verified as compiled.
    link: /verification/contracts/tables
    linkText: Tables
  - icon:
      light: /icons/smoothness.svg
      dark: /icons/smoothness-dark.svg
      width: 44
      height: 44
    title: Positivity near vacuum
    details: seth_positivity_euler, the flux limiter of Hu, Adams & Shu (2013) for the Euler equations in 1-D to 3-D. The blending factor in closed form, the largest that keeps density and pressure above their floors, on the host and on the GPU.
    link: /reference/limiters#positivity-limiter
    linkText: Positivity limiter
  - icon:
      light: /icons/verification.svg
      dark: /icons/verification-dark.svg
      width: 44
      height: 44
    title: Verification you can read
    details: Forty programs, also the test suite. Orders at a point, at critical points and over whole domains, invariances, golden values, Sod, Lax, 123, Shu–Osher and the blast wave, Burgers against its exact solution, spectral properties, the non-polynomial schemes in 1-D and 2-D; figures and animations on every case page.
    link: /verification/
    linkText: Verification
  - icon:
      light: /icons/cpu-gpu.svg
      dark: /icons/cpu-gpu-dark.svg
      width: 44
      height: 44
    title: Three compilers, CPU and GPU
    details: gfortran, nvfortran and ifx, debug and optimised, plus the OpenACC build on the GPU; built with FoBiS or fpm. Mutation runs catch every planted defect, 69 of 69 in three runs.
    link: /guide/installation
    linkText: Installation
---

<script setup>
import { withBase } from 'vitepress'
import { VPTeamMembers } from 'vitepress/theme'

const authors = [
  {
    avatar: 'https://github.com/szaghi.png',
    name: 'Stefano Zaghi',
    title: 'Chief Yak Shaver & Serial Library Resurrector',
    desc: 'Wrote WenOOF, watched it die of run-time polymorphism in the inner loops, and named its successor after the son Adam had once Abel was gone.',
    links: [{ icon: 'github', link: 'https://github.com/szaghi' }],
  },
  {
    avatar: 'https://github.com/giacrossi.png',
    name: 'Giacomo Rossi',
    title: 'High Priest of the GPU & Exorcist of Data Races',
    desc: 'Co-author of the GPU WENO machinery SETH grew from.',
    links: [{ icon: 'github', link: 'https://github.com/giacrossi' }],
  },
  {
    avatar: withBase('/avatars/claude.svg'),
    name: 'Claude',
    org: 'Anthropic',
    orgLink: 'https://claude.ai',
    title: 'Omniscient Code Oracle & Tireless Rubber Duck',
    desc: "AI pair programmer, responsible for deriving coefficients in exact rational arithmetic and for running mutants all night so humans don't have to.",
  },
]
</script>

<div class="seth-home">

## Quick start

<p class="seth-lead">Build the library and run the verification programs, or pull SETH into your own project.</p>

::: code-group

```bash [Build and test]
git clone https://github.com/szaghi/SETH && cd SETH
fobis fetch                                                   # PENF, into src/third_party/
fobis build --mode static-gnu                                 # lib/gnu/static/libseth.a
fobis build --mode tests-gnu && bash scripts/run_tests.sh     # the 40 test programs
```

```ini [FoBiS dependency]
# fobos of your project; then: fobis fetch && fobis build
[dependencies]
deps_dir = src/third_party
PENF     = https://github.com/szaghi/PENF :: branch=master   # SETH's own dependency: fobis fetch is not transitive
SETH     = https://github.com/szaghi/SETH :: branch=main
```

```toml [fpm dependency]
# fpm.toml of your project
[dependencies]
SETH = { git = "https://github.com/szaghi/SETH", branch = "main" }
```

:::

<p class="seth-lead">The two states at the face i+1/2, from the values around it (cells by offset from the face: 0 left, 1 right):</p>

```fortran
use seth
real(R_P) :: vp(-2:3), vm(-2:3), rp, rm

call seth_weno5_js(vp, vm, SETH_WENO_JS_EPS, 1._R_P, rp, rm)   ! WENO-JS 5: rp left, rm right
call seth_teno5(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)  ! TENO5
```

<p class="seth-lead">The <a href="./tutorial/">tutorial</a> goes from this face to the Euler equations on a GPU;
the <a href="./reference/">reference</a> lists every kernel, its arguments and its constants.</p>

## Authors

</div>

<VPTeamMembers size="medium" :members="authors" />

<div class="seth-home">

## License

<p class="seth-lead">SETH is available under your choice of license.</p>

| License | SPDX | File |
|---|---|---|
| GNU General Public License v3 | `GPL-3.0-or-later` | [LICENSE.gpl3.md](https://github.com/szaghi/SETH/blob/main/LICENSE.gpl3.md) |
| BSD 2-Clause | `BSD-2-Clause` | [LICENSE.bsd-2.md](https://github.com/szaghi/SETH/blob/main/LICENSE.bsd-2.md) |
| BSD 3-Clause | `BSD-3-Clause` | [LICENSE.bsd-3.md](https://github.com/szaghi/SETH/blob/main/LICENSE.bsd-3.md) |
| MIT | `MIT` | [LICENSE.mit.md](https://github.com/szaghi/SETH/blob/main/LICENSE.mit.md) |

<p class="seth-lead">© <a href="https://github.com/szaghi">Stefano Zaghi</a></p>

</div>
