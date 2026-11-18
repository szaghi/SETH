#!/usr/bin/env bash
# Build and run the documentation examples, regenerating everything the pages include from them (FUNDAL's scheme).
#
#   docs/examples/src/*.F90    the example programs (hand-written), with marker comments:
#                                !run ID COMMAND        a run shown in the pages (!run -s: with its exit status)
#                                !region NAME ... !endregion NAME   a part of the program included on its own
#   docs/examples/snippets/    generated: <program>.F90 without the markers, and <program>-<region>.F90 (dedented)
#   docs/examples/output/      generated: <ID>.txt, "$ COMMAND" then its standard output (and [exit status N] with -s)
#   docs/examples/device/      captured by hand, never touched here: the outputs of the GPU chapter on a device
#                              (its command is in the file's first line)
#
# The library is built by FoBiS (mode static-gnu: gfortran, host, the archive carries PENF) and every example is
# compiled against it exactly as a user program would be: `gfortran -I<mods> prog.F90 libseth.a` (with -fcheck=bounds:
# an example that reads out of bounds fails the build). OpenACC and OpenMP
# directives in an example are comments here: the programs run on the host. The examples must print only values that
# do not depend on the machine (explicit formats, a few significant digits; no timings, no device names). Each run
# happens in a scratch directory with a minimal environment and no standard input. The docs-examples check fails when
# the committed snippets or outputs differ from the regenerated ones.
#
# Usage: bash scripts/docs_examples.sh      (FC=/path/to/gfortran: another gfortran; it must be the one fobis uses)
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ex=$root/docs/examples
build=$root/build/docs-examples # git-ignored
lib=$root/lib/gnu/static        # fobis build --mode static-gnu
run_dir=$build/run
home=/home/user                 # how the run directory is shown
fc=${FC:-gfortran}

[ -d "$root/src/third_party/PENF/src" ] || (cd "$root" && fobis fetch) # a fresh clone: PENF
(cd "$root" && fobis build --mode static-gnu) > "$root/build-docs-examples.log" 2>&1 || {
  cat "$root/build-docs-examples.log"; echo "docs_examples: library build failed" >&2; exit 1; }
rm -f "$root/build-docs-examples.log"
rm -rf -- "$build"
mkdir -p "$build/bin" "$build/mod" "$run_dir" "$ex/snippets" "$ex/output"
rm -f -- "$ex"/snippets/*.F90 "$ex"/output/*.txt

markers='^ *!(run|region|endregion) '
dedent() { awk '{l[NR]=$0; if ($0 ~ /[^ ]/) {match($0, /^ */); if (m == "" || RLENGTH < m) m = RLENGTH}}
                END {for (i = 1; i <= NR; i++) print substr(l[i], m + 1)}' "$1"; }
for src in "$ex"/src/*.F90; do
  name=$(basename "$src" .F90)
  grep -Ev "$markers" "$src" > "$ex/snippets/$name.F90" || true
  for region in $(sed -n 's/^ *!region \([A-Za-z0-9_-]*\).*/\1/p' "$src"); do
    awk -v r="$region" '$1 == "!endregion" && $2 == r {on = 0}
                        on && $0 !~ /^ *!(run|region|endregion) / {print}
                        $1 == "!region" && $2 == r {on = 1}' "$src" > "$build/region.F90"
    dedent "$build/region.F90" > "$ex/snippets/$name-$region.F90"
  done
done

for src in "$ex"/src/*.F90; do
  "$fc" -cpp -O2 -fcheck=bounds,do,mem,pointer -I"$lib/mod" -J"$build/mod" -o "$build/bin/$(basename "$src" .F90)" "$src" "$lib/libseth.a"
done

run() { # run [-s] ID COMMAND
  local show=0 status=0
  if [ "$1" = -s ]; then show=1; shift; fi
  local id=$1; shift
  local cmd="$*"
  {
    printf '$ %s\n' "$cmd"
    (cd "$run_dir" && env -i HOME="$run_dir" PATH="$build/bin:/usr/bin:/bin" SHELL=/bin/bash LC_ALL=C \
                     LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}" OMP_NUM_THREADS=1 GFORTRAN_UNBUFFERED_PRECONNECTED=y \
                     GFORTRAN_ERROR_BACKTRACE=0 bash -c "$cmd" < /dev/null 2>&1) || status=$?
    if [ $show = 1 ]; then printf '[exit status %d]\n' "$status"; fi
  } | sed -E -e "s|$run_dir|$home|g" > "$ex/output/$id.txt"
}
for src in "$ex"/src/*.F90; do
  while IFS= read -r line; do
    line=${line#*!run }
    if [ "${line%% *}" = -s ]; then line=${line#-s }; run -s "${line%% *}" "${line#* }"
    else run "${line%% *}" "${line#* }"; fi
  done < <(grep -E '^ *!run ' "$src" || true)
done
echo "docs_examples: $(ls "$ex"/src/*.F90 | wc -l) programs, $(ls "$ex"/output/*.txt | wc -l) runs"
