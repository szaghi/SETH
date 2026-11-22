# Installation

## FoBiS

```bash
git clone https://github.com/szaghi/SETH && cd SETH
fobis fetch                      # PENF, and for the tests FUNDAL and SELAH, into src/third_party/
fobis build --mode static-gnu    # lib/gnu/static/libseth.a (shared-gnu: libseth.so)
```

`fobis build --lmodes` lists the modes: static and shared libraries and tests for gfortran (`-gnu`), nvfortran
(`-nvf`) and ifx (`-intel`), and `tests-nvf-oac`, the OpenACC build (`--varset cc80|cc89|cc90` picks the compute
capability; cc89 by default; `iac`: cc80 with CUDA 12.9, the CNR-IAC A30 nodes). `tests-amd` and `tests-amd-omp` build
the tests with amdflang, on the host and with OpenMP offload to an AMD GPU (`--varset gfx90a|gfx942|gfx950`; gfx942 by
default); they have not been compiled by the maintainers yet.

The tests have two more dependencies, which the library modes never compile: [SELAH](https://github.com/szaghi/SELAH),
the source of the time integrators a verification case needs (SETH implements space operators only), and
[FUNDAL](https://github.com/szaghi/FUNDAL), which SELAH uses. Every test mode adds FUNDAL's include directory
(SELAH's kernels include `fundal.H`) and `-DFUNDAL_NO_MPI`.

### As a FoBiS dependency

Declare SETH in the `fobos` of your project, then `fobis fetch && fobis build`:

```ini
[dependencies]
deps_dir = src/third_party
PENF     = https://github.com/szaghi/PENF :: branch=master
SETH     = https://github.com/szaghi/SETH :: branch=main
```

PENF is SETH's own dependency and must be declared too: `fobis fetch` does not fetch the dependencies of a
dependency. Both are cloned as sources and scanned with yours; only the modules your targets use are compiled.

## fpm

```toml
[dependencies]
SETH = { git = "https://github.com/szaghi/SETH", branch = "main" }
```

## Tests

```bash
fobis build --mode tests-gnu && bash scripts/run_tests.sh   # the test programs, gfortran
fobis rule --ex makecoverage-analysis                       # coverage into docs/public/coverage.json
```

What the programs check: [Test programs](/verification/unit-tests).

## Documentation examples

The example programs of the documentation are in `docs/examples/src/`.
`bash scripts/docs_examples.sh` builds the library, compiles each program against it, runs it, and regenerates the
code and output the pages include; CI fails when they drift from the committed ones.
