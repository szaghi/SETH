---
title: Tooling
---

# Tooling

The build modes of `fobos` and the scripts of the repository, each with what it is for.

## Build modes

`fobis build --mode <mode>`; `fobis build --lmodes` lists them.

| Modes | What they build |
|---|---|
| `static-gnu`, `static-nvf`, `static-intel` | the static library, `lib/<compiler>/static/libseth.a`, gfortran, nvfortran, ifx, optimised |
| `shared-gnu`, `shared-nvf`, `shared-intel` | the shared library, `libseth.so` |
| `static-gnu-r4`, `static-nvf-r4`, `static-intel-r4`, `shared-gnu-r4`, `shared-nvf-r4`, `shared-intel-r4` | the same in single precision, `R_P = R4P` ([Precision](./precision)) |
| `example-t11-nvf-oac` | tutorial chapter 11 on the GPU (`exe-gpu/t11_gpu`): nvfortran with OpenACC on `$GPU_CC`, SELAH and FUNDAL in their device build (`-DDEV_OAC -DCOMPILER_NVF`) |
| `static-selah-gnu` | [SELAH](https://szaghi.github.io/SELAH/)'s library from `src/third_party/SELAH`, gfortran: the time integrator of the docs examples, linked by `scripts/docs_examples.sh` after `libseth.a` |
| `tests-gnu`, `tests-nvf`, `tests-intel` | the test programs, optimised, into `exe/` |
| `tests-gnu-debug`, `tests-nvf-debug`, `tests-intel-debug` | the test programs with every runtime check |
| `tests-gnu-coverage` | the test programs instrumented for coverage (`fobis rule --ex makecoverage-analysis`) |
| `tests-nvf-oac` | the test programs with OpenACC, on the GPU of compute capability `$GPU_CC` (`--varset cc80`, `cc89`, `cc90`) |
| `tests-gnu-r4`, `tests-gnu-debug-r4`, `tests-nvf-r4`, `tests-nvf-oac-r4`, `tests-intel-r4` | the test programs in single precision, into `exe-r4/` |
| `bench-nvf-oac`, `bench-nvf`, `bench-gnu` | the kernel benchmark against ADAM's kernels, on the GPU or one CPU core ([Performance](./performance)) |
| `bench-arrays-nvf-oac`, `bench-arrays-nvf`, `bench-arrays-gnu` | the benchmark of the device-kernel rules: automatic arrays and loop extents |
| `bench-kinds-nvf-oac`, `bench-kinds-gnu` | the kernels' benchmark in `R_P` (double precision), the reference of the single-precision comparison |
| `bench-nvf-oac-r4`, `bench-gnu-r4` | the same in single precision |

## Scripts

| Script | What it does |
|---|---|
| `scripts/run_tests.sh` | runs the test programs of `exe/` (`--exe-dir exe-r4`: the single-precision build) |
| `scripts/run_vv.sh` | runs the verification cases in report mode and regenerates the documentation's data, tables, figures and animations |
| `scripts/docs_examples.sh` | builds and runs the documentation's example programs, regenerating the code and output the pages include |
| `scripts/check_docs.py` | checks that the documentation covers the code and names nothing the code does not have |
| `scripts/check_papers.py` | checks SETH's tables and constants against the papers they come from, in exact arithmetic |
| `scripts/check_literals.py` | finds the kind-sensitive real literals outside `seth_kind_parameters` |
| `scripts/gen_weno_js_tables.py` | generates `seth_weno_js_tables` from ADAM's WENO initialisers |
| `scripts/gen_z_teno_tables.py` | generates `seth_z_teno_tables`, the descriptors of WENO-Z, TENO and the linear schemes |
| `scripts/gen_nonpoly_tables.py` | generates `seth_nonpoly_tables`, the constants of the non-polynomial schemes |
| `scripts/mutation_test.py` | mutation testing: does every planted defect make some test fail? ([C6](/verification/contracts/mutation)) |
| `scripts/bench_compare.py` | compares runs of the benchmark with a recorded baseline |
| `scripts/compute-coverage.sh` | the coverage percentage of a coverage build |
| `scripts/install.sh` | downloads and/or builds the project from GitHub |
| `scripts/release.sh` | bumps the version, updates the changelog, commits and tags |
| `scripts/ci/matrix.sh` | the project's pipeline on one toolchain, for the CI matrix |
| `docs/scripts/make_doc_figures.py` | draws the verification figures and writes the result tables from the report data |
| `docs/scripts/make_doc_animations.py` | renders the verification animations |
| `docs/scripts/make_convention_figures.py` | draws the figures and the animation of [Conventions](./conventions) |
| `docs/scripts/make_performance_figures.py` | draws the figures of [Performance](./performance) |
| `docs/scripts/figure_theme.py` | the light and dark styles of the figures |
| `docs/scripts/seth_palette.py` | the colours of the figures and animations |
