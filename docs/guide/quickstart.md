# Quick Start

```fortran
use seth, only : I4P, R8P   ! the working kinds, re-exported from PENF
```

The reconstruction and interpolation API is being extracted from ADAM; this page follows it.

## Tests

```bash
fobis fetch
fobis build --mode tests-gnu && bash scripts/run_tests.sh
fobis rule --ex makecoverage-analysis   # what CI runs: coverage into docs/public/coverage.json
```
