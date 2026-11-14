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

## fpm

```toml
[dependencies]
SETH = { git = "https://github.com/szaghi/SETH", branch = "main" }
```
