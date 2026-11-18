# Installation

## FoBiS

```bash
git clone https://github.com/szaghi/SETH && cd SETH
fobis fetch                      # PENF, into src/third_party/
fobis build --mode static-gnu    # lib/gnu/static/libseth.a (shared-gnu: libseth.so)
```

`fobis build --lmodes` lists the modes: static and shared libraries and tests for gfortran (`-gnu`), nvfortran
(`-nvf`) and ifx (`-intel`), and `tests-nvf-oac`, the OpenACC build (`--varset cc80|cc89|cc90` picks the compute
capability; cc89 by default).

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
