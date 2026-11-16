#!/usr/bin/env python3
"""The R4P speed-up of SETH's kernels: S = t(R8P) / t(R4P) per kernel of seth_bench_kinds, and its median.

DIR holds the raw outputs of interleaved runs, DIR/r8p/run-*.txt and DIR/r4p/run-*.txt. Per kernel: the median over
the runs of each kind (ns per face) and their spread (max - min); S is the ratio of the two medians. Every run must
have passed its checks.

Usage: speedup.py DIR
"""

from __future__ import annotations

import re
import statistics
import sys
from pathlib import Path

ROW = re.compile(r'^(?P<name>.{31})\s*(?P<med>[0-9.]+)\s+(?P<min>[0-9.]+)\s+(?P<ns>[0-9.]+)\s*$')


def load(directory: Path) -> list[dict[str, float]]:
    """ns per face per kernel, one dictionary per run."""
    runs = []
    for path in sorted(directory.glob('run-*.txt')):
        text = path.read_text()
        if 'checks ok' not in text:
            sys.exit(f'{path}: the checks did not pass')
        rows = (ROW.match(line) for line in text.splitlines())
        runs.append({m['name'].strip(): float(m['ns']) for m in rows if m})
    return runs


def main() -> int:
    """Print the per-kernel table and the median S."""
    root = Path(sys.argv[1])
    r8, r4 = load(root / 'r8p'), load(root / 'r4p')
    print(f'{len(r8)} R8P runs, {len(r4)} R4P runs; ns per face: median over the runs (spread max - min)')
    print(f'{"kernel":31s} {"R8P":>16s} {"R4P":>16s} {"S":>6s}')
    speedups = []
    for kernel in r8[0]:
        a = [run[kernel] for run in r8]
        b = [run[kernel] for run in r4]
        ma, mb = statistics.median(a), statistics.median(b)
        speedups.append(ma / mb)
        print(f'{kernel:31s} {ma:8.3f} ({max(a) - min(a):5.3f}) {mb:8.3f} ({max(b) - min(b):5.3f}) {ma / mb:6.2f}')
    print(f'median S over {len(speedups)} kernels: {statistics.median(speedups):.2f} '
          f'(min {min(speedups):.2f}, max {max(speedups):.2f})')
    return 0


if __name__ == '__main__':
    sys.exit(main())
