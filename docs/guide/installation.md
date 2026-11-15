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
