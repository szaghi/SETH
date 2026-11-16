#!/usr/bin/env bash
# Build and run the centred-WENO packet experiment (notes/weno-c-burgers-packet.md) and check that it reproduces
# seth_vv_burgers. Run from the repository root, with gfortran on the PATH: bash notes/experiments/weno-c-packet/run.sh [outdir]
set -euo pipefail
here=notes/experiments/weno-c-packet
out=$(realpath -m "${1:-${TMPDIR:-/tmp}/weno-c-packet}")
mkdir -p "$out"
echo "=== build tests-gnu (the library and the test modules)"
fobis build --mode tests-gnu > "$out/build.log" 2>&1 || { tail -5 "$out/build.log"; exit 1; }
# the objects of every module (not of the programs) of the tests-gnu build
progs=$(grep -l '^program' src/tests/*.F90 | xargs -n1 basename | sed 's/\.F90$/.o/')
objs=()
for o in exe/obj/*.o; do
  grep -qx "$(basename "$o")" <<< "$progs" || objs+=("$o")
done
gfortran -O3 -funroll-loops -std=f2018 -Iexe/mod -J "$out" "$here/weno_c_packet.F90" "${objs[@]}" -o "$out/weno_c_packet"
echo "=== the harness in report mode (seth_vv_burgers, every scheme)"
SETH_VV_DIR="$out/harness" ./exe/seth_vv_burgers > "$out/harness.log" 2>&1 || { tail -5 "$out/harness.log"; exit 1; }
echo "=== the experiment"
(cd "$out" && ./weno_c_packet > experiment.txt 2>&1) || { tail -5 "$out/experiment.txt"; exit 1; }
echo "=== adam against the harness: L1 before, L1 after, L1 away, overshoot within 1e-8 relative + 1e-14"
# Not bitwise: the two programs are compiled apart and differ in the last bits, which the shock amplifies; 1e-14 is
# the round-off floor of an L1 (the pre-shock L1 of WENO8-C at 640 cells is 9e-13).
python3 - "$out" <<'EOF'
import sys
out = sys.argv[1]
def rows(path):
    r = {}
    for line in open(path):
        if line.startswith('#'):
            continue
        w = line.split()
        r[(w[0], int(w[1]))] = [float(x) for x in w[2:]]
    return r
pre = rows(f'{out}/harness/burgers-smooth.dat')
post = rows(f'{out}/harness/burgers-shock.dat')
adam = rows(f'{out}/packet-adam.dat')
worst = 0.0
worst_away = 0.0
for (name, n), (e_pre, e_all, e_away, osh) in sorted(adam.items()):
    ref = [pre[(name, n)][0], post[(name, n)][0], post[(name, n)][1], post[(name, n)][2]]
    for a, b in zip((e_pre, e_all, e_away, osh), ref):
        worst = max(worst, abs(a - b) / (1e-8 * abs(b) + 1e-14))
    worst_away = max(worst_away, abs(e_away - ref[2]) / ref[2])
print(f'   {len(adam)} runs; largest difference over the bound {worst:.2f}; L1 away: largest relative difference '
      f'{worst_away:.1e}')
sys.exit(0 if worst <= 1.0 else 1)
EOF
echo "run.sh: ok"
