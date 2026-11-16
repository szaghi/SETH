#!/usr/bin/env python3
"""Compare runs of seth_bench_weno with a recorded baseline: does any kernel run slower than before?

A baseline is a directory of raw outputs of seth_bench_weno, several runs of the same build on the same machine
(notes/bench-baseline/<machine>-<build>/run-*.txt). For every variant the baseline gives the median of its runs' ns
per face and their spread (max - min). A new set of runs passes when, for every variant present in both, the median
of the new runs is at most

    base_median * (1 + tolerance) + base_spread        (tolerance 3 % by default)

Variants only in the new runs (kernels added since the baseline) are listed, not judged; variants missing from the
new runs fail. Exit status 0 when nothing regressed, 1 otherwise.

Usage:
    bench_compare.py BASELINE_DIR NEW_RUN.txt [NEW_RUN.txt ...] [--tolerance 0.03]
    bench_compare.py --summary RUN.txt [RUN.txt ...]     (print the per-variant median and spread only)
"""

from __future__ import annotations

import argparse
import re
import statistics
import sys
from pathlib import Path

ROW = re.compile(r'^(?P<name>.{31})\s*(?P<med>[0-9.]+)\s+(?P<min>[0-9.]+)\s+(?P<ns>[0-9.]+)(\s+[-0-9.]+)?\s*$')


def parse(path: Path) -> dict[str, float]:
    """ns per face of every variant of one run (the table after the 'variant' header)."""
    out: dict[str, float] = {}
    table = False
    for line in path.read_text().splitlines():
        if line.lstrip().startswith('variant') and 'ns/face' in line:
            table = True
            continue
        if table:
            m = ROW.match(line)
            if not m:
                if out:
                    break
                continue
            out[m['name'].strip()] = float(m['ns'])
    if not out:
        sys.exit(f'{path}: no benchmark table found')
    return out


def summarise(paths: list[Path]) -> dict[str, tuple[float, float, int]]:
    """Per variant: median of the runs' ns per face, spread (max - min), number of runs."""
    runs = [parse(p) for p in paths]
    names = list(dict.fromkeys(n for r in runs for n in r))
    out = {}
    for n in names:
        vals = [r[n] for r in runs if n in r]
        out[n] = (statistics.median(vals), max(vals) - min(vals), len(vals))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('baseline', nargs='?', type=Path, help='directory of baseline runs')
    ap.add_argument('runs', nargs='*', type=Path, help='new runs')
    ap.add_argument('--tolerance', type=float, default=0.03, help='relative tolerance (default 0.03)')
    ap.add_argument('--summary', nargs='+', type=Path, help='only print the summary of these runs')
    args = ap.parse_args()
    if args.summary:
        for n, (med, spread, k) in summarise(args.summary).items():
            print(f'{n:31s} {med:10.3f} ns/face  spread {spread:8.3f}  runs {k}')
        return 0
    if args.baseline is None or not args.runs:
        ap.error('a baseline directory and at least one new run are required')
    base_files = sorted(args.baseline.glob('run-*.txt'))
    if not base_files:
        sys.exit(f'{args.baseline}: no run-*.txt')
    base = summarise(base_files)
    new = summarise(args.runs)
    failed = False
    print(f'{"variant":31s} {"base ns":>10s} {"spread":>8s} {"new ns":>10s} {"ratio":>7s}  verdict')
    for n, (bmed, bspread, _k) in base.items():
        if n not in new:
            print(f'{n:31s} {bmed:10.3f} {bspread:8.3f} {"-":>10s} {"-":>7s}  MISSING')
            failed = True
            continue
        nmed = new[n][0]
        limit = bmed * (1 + args.tolerance) + bspread
        ok = nmed <= limit
        failed |= not ok
        print(f'{n:31s} {bmed:10.3f} {bspread:8.3f} {nmed:10.3f} {nmed / bmed:7.3f}  {"ok" if ok else "SLOWER"}')
    for n, (nmed, _s, _k) in new.items():
        if n not in base:
            print(f'{n:31s} {"-":>10s} {"-":>8s} {nmed:10.3f} {"-":>7s}  new')
    print('bench_compare: ' + ('REGRESSION' if failed else 'ok'))
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
