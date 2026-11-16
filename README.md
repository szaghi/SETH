<div align="center">

<picture>
  <source media="(prefers-color-scheme: dark)"  srcset="docs/public/logo-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="docs/public/logo.svg">
  <img src="docs/public/logo.svg" alt="SETH, the Riemann fan opened like a folding fan on a crackle-glazed tile" width="320">
</picture>

#### *the WENO and TENO kernels of ADAM, one pure function each, on every CPU and GPU*

### **S**moothness-weighted **E**NO, **T**ENO and **H**ybrids

> High-order non-oscillatory reconstruction and interpolation for HPC, in modern Fortran. Extracted from the WENO
> machinery of ADAM and successor of WenOOF: one `pure`, device-callable kernel per scheme, coefficients as
> compile-time constants, no dynamic dispatch where the cells are swept.

[![Version](https://img.shields.io/github/v/release/szaghi/SETH?label=version)](https://github.com/szaghi/SETH/releases)
[![CI](https://github.com/szaghi/SETH/actions/workflows/ci.yml/badge.svg)](https://github.com/szaghi/SETH/actions/workflows/ci.yml)
[![Matrix](https://github.com/szaghi/SETH/actions/workflows/matrix.yml/badge.svg)](https://github.com/szaghi/SETH/actions/workflows/matrix.yml)
[![Fortran](https://img.shields.io/badge/Fortran-modern-734f96?logo=fortran&logoColor=white)](https://fortran-lang.org/)
[![License: GPL v3](https://img.shields.io/badge/license-GPL%20v3-blue)](LICENSE.gpl3.md)
[![License: BSD-2](https://img.shields.io/badge/license-BSD--2--Clause-blue)](LICENSE.bsd-2.md)
[![License: BSD-3](https://img.shields.io/badge/license-BSD--3--Clause-blue)](LICENSE.bsd-3.md)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue)](LICENSE.mit.md)

<div>
<table>
<tr>
<td><b>🧮 Every scheme, one kernel each</b><br><sub>WENO-JS 1/3/5/7/9, WENO-Z 3/5/7/9, TENO5, TENO5-opt, TENO5-AZ, TENO5-A, TENO6, TENO6-opt, NENO6, TENO8 and TENO8-A, for reconstruction from cell averages and interpolation from point values, optionally returning the weights ratios. <a href="https://szaghi.github.io/SETH/guide/schemes">Schemes</a> · <a href="https://szaghi.github.io/SETH/guide/api">API</a></sub></td>
<td><b>⚡ Called from inside your loops</b><br><sub>Every kernel is <code>pure</code>, <code>!$acc routine seq</code> and <code>!$omp declare target</code>; the coefficients are compile-time constants, so nothing is set up or copied to the device. On an RTX 4070 the WENO-5/7 kernels run 1.26–1.31× faster than ADAM's. <a href="https://szaghi.github.io/SETH/guide/performance">Performance</a></sub></td>
</tr>
<tr>
<td><b>🎯 Bitwise ADAM, exact coefficients</b><br><sub>WENO-JS is bitwise equal to ADAM's kernels. The WENO-Z and TENO coefficients are derived in exact rational arithmetic, cross-checked against the papers and re-verified as compiled. <a href="https://szaghi.github.io/SETH/verification/contracts/adam-equivalence">Equivalence with ADAM</a> · <a href="https://szaghi.github.io/SETH/verification/contracts/tables">Tables</a></sub></td>
<td><b>🛡️ Positivity near vacuum</b><br><sub><code>seth_positivity_euler</code>, the flux limiter of Hu, Adams &amp; Shu (2013) for the Euler equations in 1-D to 3-D: the blending factor in closed form, the largest that keeps density and pressure above their floors, on the host and on the GPU. <a href="https://szaghi.github.io/SETH/guide/api#positivity-limiter">Positivity limiter</a></sub></td>
</tr>
<tr>
<td><b>🔬 Verification you can read</b><br><sub>Twenty-seven programs, also the test suite: orders at a point, at critical points and over whole domains, invariances, golden values, Sod, Lax, 123, Shu–Osher and the blast wave in a 1-D Euler harness, Burgers against its exact solution, spectral properties; figures and animations on every case page. <a href="https://szaghi.github.io/SETH/verification/">Verification</a></sub></td>
<td><b>🧰 Three compilers, CPU and GPU</b><br><sub>gfortran, nvfortran and ifx, debug and optimised, plus the OpenACC build on the GPU; built with FoBiS or fpm. A mutation run catches 25 of 25 planted defects. <a href="https://szaghi.github.io/SETH/guide/installation">Installation</a> · <a href="https://szaghi.github.io/SETH/verification/contracts/mutation">Mutation testing</a></sub></td>
</tr>
</table>
</div>

**[Full documentation](https://szaghi.github.io/SETH/)**: [guide](https://szaghi.github.io/SETH/guide/) · [verification](https://szaghi.github.io/SETH/verification/) · [API reference](https://szaghi.github.io/SETH/api/)

</div>

---

## Quick start

```bash
git clone https://github.com/szaghi/SETH && cd SETH
fobis fetch                                                   # PENF, into src/third_party/
fobis build --mode static-gnu                                 # lib/gnu/static/libseth.a
fobis build --mode tests-gnu && bash scripts/run_tests.sh     # the 27 verification programs
```

as a FoBiS dependency, in the `fobos` of your project (then `fobis fetch && fobis build`):

```ini
[dependencies]
deps_dir = src/third_party
PENF     = https://github.com/szaghi/PENF :: branch=master   # SETH's own dependency: fobis fetch is not transitive
SETH     = https://github.com/szaghi/SETH :: branch=main
```

or as an fpm dependency, in your `fpm.toml`:

```toml
[dependencies]
SETH = { git = "https://github.com/szaghi/SETH", branch = "main" }
```

The two states at the face i+1/2, from the values around it (cells by offset from the face: 0 left, 1 right):

```fortran
use seth
real(R8P) :: vp(-2:3), vm(-2:3), rp, rm

call seth_weno5_js(vp, vm, SETH_WENO_JS_EPS, 1._R8P, rp, rm)   ! WENO-JS 5: rp left, rm right
call seth_teno5(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)  ! TENO5
```

The [quick start](https://szaghi.github.io/SETH/guide/quickstart) lists every kernel and its parameters;
[conventions](https://szaghi.github.io/SETH/guide/conventions) explains the stencil layout and the kinds.

## Authors

**[Stefano Zaghi](https://github.com/szaghi)** · stefano.zaghi@gmail.com
> *Chief Yak Shaver & Serial Library Resurrector* — wrote WenOOF, watched it die of run-time polymorphism in the inner loops, and named its successor after the son Adam had once Abel was gone

**[Giacomo Rossi](https://github.com/giacrossi)**
> *High Priest of the GPU & Exorcist of Data Races* — ADAM co-author, from whose WENO machinery SETH was carved out, rib by rib

**[Claude](https://claude.ai)** (Anthropic)
> *Omniscient Code Oracle & Tireless Rubber Duck* — AI pair programmer, responsible for deriving coefficients in exact rational arithmetic and for running mutants all night so humans don't have to

## License

SETH is available under your choice of license:

| License | SPDX | File |
|---|---|---|
| GNU General Public License v3 | `GPL-3.0-or-later` | [LICENSE.gpl3.md](LICENSE.gpl3.md) |
| BSD 2-Clause | `BSD-2-Clause` | [LICENSE.bsd-2.md](LICENSE.bsd-2.md) |
| BSD 3-Clause | `BSD-3-Clause` | [LICENSE.bsd-3.md](LICENSE.bsd-3.md) |
| MIT | `MIT` | [LICENSE.mit.md](LICENSE.mit.md) |

© [Stefano Zaghi](https://github.com/szaghi)
