#!/usr/bin/env bash
# Run SETH's verification cases in report mode and regenerate the documentation's data, tables, figures and animations.
#
# Usage: scripts/run_vv.sh [--mode MODE] [--data DIR] [--no-build] [--no-figures] [--only NAME ...]
#   --mode MODE    fobis mode that builds the cases (default: tests-gnu, gfortran -O2)
#   --data DIR     where the report data goes (default: exe/vv-data; build products, not committed)
#   --no-build     use the executables already in exe/
#   --no-figures   run the cases only
#   --only NAME    run only these programs (seth_vv_sod, seth_test_spectral, ...); repeatable
#
# Every verification program asserts its bounds in both modes; with SETH_VV_DIR set it also runs every scheme on the
# report grids and writes its data there. The figures (light and dark twins, docs/scripts/figure_theme.py) and the
# animations go to docs/public/verification/, the result tables to docs/verification/_tables/ (included by the case
# pages). Needs the Python of exe/.vv-venv (numpy, matplotlib) or SETH_VV_PYTHON, and ffmpeg for the animations.
set -euo pipefail
cd "$(dirname "$0")/.."
MODE=tests-gnu
DATA=exe/vv-data
BUILD=1
FIGURES=1
ONLY=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) MODE="$2"; shift 2 ;;
    --data) DATA="$2"; shift 2 ;;
    --no-build) BUILD=0; shift ;;
    --no-figures) FIGURES=0; shift ;;
    --only) ONLY+=("$2"); shift 2 ;;
    *) printf 'unknown argument: %s\n' "$1" >&2; exit 2 ;;
  esac
done
PY="${SETH_VV_PYTHON:-exe/.vv-venv/bin/python}"
[[ $BUILD == 1 ]] && fobis build --mode "$MODE"
PROGRAMS=(seth_vv_sod seth_vv_lax seth_vv_test123 seth_vv_shu_osher seth_vv_blast seth_vv_burgers seth_vv_euler_order
          seth_vv_advection seth_vv_adr seth_test_spectral seth_test_convergence seth_test_interpolation
          seth_test_weno_js_order seth_test_wenoz seth_test_teno seth_test_critical_points seth_test_sensors seth_vv_hybrid
          seth_vv_roe)
[[ ${#ONLY[@]} -gt 0 ]] && PROGRAMS=("${ONLY[@]}")
mkdir -p "$DATA"
DATA_ABS="$(cd "$DATA" && pwd)"
status=0
for p in "${PROGRAMS[@]}"; do
  printf '%-28s ' "$p"
  start=$(date +%s)
  if SETH_VV_DIR="$DATA_ABS" "./exe/$p" > "$DATA_ABS/$p.log" 2>&1; then
    printf 'ok    %5ss\n' "$(( $(date +%s) - start ))"
  else
    printf 'FAIL  %5ss  (see %s)\n' "$(( $(date +%s) - start ))" "$DATA/$p.log"; status=1
  fi
done
if [[ $FIGURES == 1 ]]; then
  "$PY" docs/scripts/make_doc_figures.py --data "$DATA_ABS" --out docs/public/verification --tables docs/verification/_tables
  "$PY" docs/scripts/make_doc_animations.py --data "$DATA_ABS" --out docs/public/verification/casts
fi
exit $status
