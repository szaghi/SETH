#!/usr/bin/env bash
# The project's pipeline on one toolchain, run by .github/workflows/matrix.yml on each toolchain of the fobos
# [scaffold] matrix, with the fortran-toolchains environment contract: FC, CC, CXX, TOOLCHAIN, TOOLCHAIN_OS and
# TOOLCHAIN_FAMILY (gnu, nvf, intel). Created once by fobis scaffold, then the project's: edit it freely. The same
# script runs locally, in `docker run <image>` and in CI.
set -euo pipefail
echo "toolchain: $TOOLCHAIN ($TOOLCHAIN_OS, $TOOLCHAIN_FAMILY), FC=$FC"
# the separation of time (SELAH) and space (GIDEON, SETH) operators: no copy in the tests (SELAH issue #2)
python3 scripts/check_separation.py --self-test
python3 scripts/check_separation.py
fobis fetch
fobis build --mode "tests-$TOOLCHAIN_FAMILY"
bash scripts/run_tests.sh
