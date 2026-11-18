#!/usr/bin/env python3
"""Draw the figures of the Performance page of SETH's guide (docs/guide/performance.md) from the benchmark outputs.

Why this script exists: the numbers of the page are measurements, kept as the raw outputs of the two benchmark
programs; the figures are drawn from those files, not typed in, so that a rerun of a benchmark only needs a rerun of
this script.

Reads
    notes/m3-bench-data/bench-gpu1-run{1,2,3,4-sumsquares}.txt, bench-cpu-{gnu,nvf}.txt
                                        src/bench/seth_bench_weno.F90: SETH against ADAM (RTX 4070, one CPU core)
    notes/gpu-audit-data/gpu1-after-fix-with-auto.txt, cpu-{gfortran,nvfortran}.txt
                                        src/bench/seth_bench_device_arrays.F90: automatic arrays, loop extents

Writes, into --out (default docs/public/guide), each with its silver-ink twin <name>-dark.png (docs/scripts/
figure_theme.py, the harness of the ADAM and FLUME docs):

    perf-kernels.png        ns per face, SETH against ADAM's per-width and generic kernels, GPU and one CPU core
    perf-device-rules.png   ns per face (log scale) of the same kernel following or breaking the two device rules
    perf-nonpoly.png        the non-polynomial kernels of issue #6: ns per face over seth_teno5's, GPU and both CPU
                            compilers, against their budgets; and the composites' cost against the spacing of the
                            jumps of the field (notes/m5-bench-data/42-nonpoly-*, 42-fields-*)

GPU: the median over the runs of the median of each run, with the spread of the runs as an error bar (WENO-JS and
ADAM: runs 1 to 4; WENO-Z5 and TENO6: run 4, after their indicators became sums of squares). CPU: WENO-5 and WENO-7
only, since the CPU runs predate the sum-of-squares indicators of WENO-Z and TENO.

Usage: make_performance_figures.py [--out DIR] [--theme {light,dark,both}]
Needs numpy and matplotlib.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib as mpl
import numpy as np

mpl.use('Agg')
mpl.rcParams['figure.facecolor'] = '#f3ede1'
mpl.rcParams['savefig.facecolor'] = '#f3ede1'
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_theme import THEMES, theme, themes  # noqa: E402
import seth_palette as P  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
M3 = ROOT / 'notes/m3-bench-data'
M5 = ROOT / 'notes/m5-bench-data'
AUDIT = ROOT / 'notes/gpu-audit-data'
DPI = 130
# seth_palette: SETH in the identity's blue, ADAM in neutral greys; follows the rules blue, half-follows light blue,
# automatic arrays umber
SETH, WIDTH, GENERIC, HOST = P.BLUE, '0.45', '0.72', '0.88'
SLOW, MILD, GOOD = P.UMBER, '#8db6e0', P.BLUE
ROW = re.compile(r'^(?P<name>\S.*?)\s+(?P<med>[\d.]+)\s+(?P<min>[\d.]+)\s+(?P<ns>[\d.]+)\s+(?P<ratio>[\d.]+)\s*$')


def table(path: Path) -> dict[str, float]:
    """The ns-per-face column of a benchmark's result table, by variant name (whitespace collapsed)."""
    out: dict[str, float] = {}
    for line in path.read_text().splitlines():
        m = ROW.match(line)
        if m and not line.lstrip().startswith(('variant', 'faces')):
            out[' '.join(m['name'].split())] = float(m['ns'])
    return out


def spread(values: list[float]) -> tuple[float, float, float]:
    """Median and the distances to the extremes (an asymmetric error bar)."""
    v = np.array(values)
    med = float(np.median(v))
    return med, med - v.min(), v.max() - med


# ------------------------------------------------------------------------------------------------- figures
def fig_kernels(out: Path) -> None:
    """SETH against ADAM: WENO-5 and WENO-7, GPU and one CPU core; SETH's WENO-Z5 and TENO6 for reference."""
    runs = [table(M3 / f) for f in ('bench-gpu1-run1.txt', 'bench-gpu1-run2.txt', 'bench-gpu1-run3.txt',
                                     'bench-gpu1-run4-sumsquares.txt')]
    gnu, nvf = table(M3 / 'bench-cpu-gnu.txt'), table(M3 / 'bench-cpu-nvf.txt')
    groups = [('WENO-5', 'SETH weno5_js', 'ADAM width s3', 'ADAM generic S=3', 'ADAM host S=3 (FLUME CPU)'),
              ('WENO-7', 'SETH weno7_js', 'ADAM width s4', 'ADAM generic S=4', 'ADAM host S=4 (FLUME CPU)')]
    bars = (('SETH', SETH), ("ADAM per-width kernel", WIDTH), ('ADAM generic kernel', GENERIC),
            ("ADAM host routine (FLUME's CPU path)", HOST))
    fig, axs = plt.subplots(1, 3, figsize=(10.4, 3.9), gridspec_kw={'width_ratios': [1.35, 1.0, 1.0]})
    # GPU
    ax = axs[0]
    for g, (lab, *names) in enumerate(groups):
        for b, name in enumerate(names[:3]):
            med, lo, hi = spread([r[name] for r in runs])
            x = g + (b - 1) * 0.26
            ax.bar(x, med, 0.24, color=bars[b][1], edgecolor='k', lw=0.4, label=bars[b][0] if g == 0 else None)
            ax.errorbar(x, med, yerr=[[lo], [hi]], fmt='none', ecolor='k', elinewidth=0.7, capsize=2)
            ax.text(x, med + hi + 0.04, f'{med:.2f}', ha='center', va='bottom', fontsize=7)
    for g, name in ((2, 'SETH wenoz5 (p=1)'), (3, 'SETH teno6 (C_T=1e-7)')):
        v = runs[3][name]
        ax.bar(g, v, 0.24, color=SETH, edgecolor='k', lw=0.4, hatch='//',
               label='SETH, no ADAM counterpart' if g == 2 else None)
        ax.text(g, v + 0.04, f'{v:.2f}', ha='center', va='bottom', fontsize=7)
    ax.set_xticks(range(4), ['WENO-5', 'WENO-7', 'WENO-Z5\n(p = 1)', 'TENO6'])
    ax.set_ylabel('ns per face')
    ax.set_title('RTX 4070 (GPU), nvfortran -acc', fontsize=10)
    ax.set_ylim(0, 2.9)
    # CPU, two compilers
    for ax, data, comp in ((axs[1], gnu, 'gfortran -O3'), (axs[2], nvf, 'nvfortran -fast')):
        for g, (lab, *names) in enumerate(groups):
            for b, name in enumerate(names):
                x = g + (b - 1.5) * 0.2
                v = data[name]
                ax.bar(x, v, 0.19, color=bars[b][1], edgecolor='k', lw=0.4,
                       label=bars[b][0] if (g == 0 and ax is axs[1]) else None)
                ax.text(x, v + 1.0, f'{v:.0f}', ha='center', va='bottom', fontsize=7)
        # no WENO-Z5/TENO6 here: the CPU runs predate their sum-of-squares indicators
        ax.set_xticks(range(2), ['WENO-5', 'WENO-7'])
        ax.set_title(f'one i9-13900KF core, {comp}', fontsize=10)
        ax.set_ylim(0, 98)
    for ax in axs:
        ax.grid(axis='y', alpha=0.25)
        ax.set_axisbelow(True)
    handles, labels = axs[1].get_legend_handles_labels()
    h0, l0 = axs[0].get_legend_handles_labels()
    handles.append(h0[-1])
    labels.append(l0[-1])
    fig.legend(handles, labels, loc='lower center', ncol=5, fontsize=7.5, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(out / 'perf-kernels.png', dpi=DPI)
    plt.close(fig)


NP_KERNELS = [  # bench variant, label, budget against seth_teno5 (issue #6, PERF1; None: no budget)
    ('SETH teno5_mask (one bias)', 'mask, one bias', 0.70),
    ('SETH dr_detect (stored masks)', 'DR detector', 0.10),
    ('SETH thinc_detect (stored masks)', 'THINC detector', 0.10),
    ('SETH dr_override (level 1)', 'DR override', 0.15),
    ('SETH thinc_override (level 1)', 'THINC override', 0.30),
    ('SETH teno5_st', 'TENO5-ST', 1.05),
    ('SETH teno5_rbf (s1 = s2 = 0.09)', 'TENO5-RBF', 1.25),
    ('SETH teno5_dr', 'TENO5-DR', 3.0),
    ('SETH teno5_thinc', 'TENO5-THINC', 3.0),
    ('SETH teno5_rbf_dr', 'TENO5-RBF+DR', 3.0),
    ('SETH teno5_dr, stored masks', 'DR, stored masks', None),
    ('SETH teno5_thinc, stored masks', 'THINC, stored masks', None),
    ('SETH teno5_rbf_dr, stored masks', 'RBF+DR, stored masks', None),
]
TENO5 = 'SETH teno5 (C_T=1e-5)'


def lookup(run: dict[str, float], name: str) -> float:
    """A variant's ns per face; the benchmark truncates its names to 31 characters (30 after collapsing spaces)."""
    hits = [v for k, v in run.items() if k == name or (len(k) >= 28 and name.startswith(k))]
    return hits[0] if hits else float('nan')


def np_ratios(files: list[Path]) -> dict[str, float]:
    """Per variant, the median over the runs of (ns per face / seth_teno5's ns per face of the same run)."""
    runs = [table(f) for f in files]
    return {name: float(np.median([lookup(r, name) / r[TENO5] for r in runs])) for name, _, _ in NP_KERNELS}


def fig_nonpoly(out: Path) -> None:
    """Issue #6: the cost of every non-polynomial kernel against its budget, and against the density of jumps."""
    dev = {'GPU (RTX 4070)': sorted(M5.glob('42b-nonpoly-gpu-run-*.txt')),
           'gfortran core': sorted(M5.glob('42-nonpoly-cpu-gnu-run-*.txt')),
           'nvfortran core': sorted(M5.glob('42-nonpoly-cpu-nvf-run-*.txt'))}
    if not all(dev.values()):
        print('perf-nonpoly: skipped, no series 42')
        return
    rat = {d: np_ratios(f) for d, f in dev.items()}
    cols = (SETH, MILD, SLOW)
    fig, axs = plt.subplots(1, 2, figsize=(12.0, 4.2), gridspec_kw={'width_ratios': [2.4, 1.0]})
    ax = axs[0]
    x = np.arange(len(NP_KERNELS))
    for j, (d, r) in enumerate(rat.items()):
        ax.bar(x + (j - 1) * 0.27, [r[n] for n, _, _ in NP_KERNELS], 0.26, color=cols[j], edgecolor='k', lw=0.3,
               label=d)
    for i, (_, _, b) in enumerate(NP_KERNELS):
        if b is not None:
            ax.plot([i - 0.42, i + 0.42], [b, b], '-', color='k', lw=1.4)
    ax.plot([], [], '-', color='k', lw=1.4, label='budget')
    ax.set_xticks(x, [lab for _, lab, _ in NP_KERNELS], rotation=40, ha='right', fontsize=7.5)
    ax.set_yscale('log')
    ax.set_ylabel('ns per face / seth_teno5\'s')
    ax.set_title('the pieces and the kernels, a jump every 2048 cells', fontsize=10)
    ax.legend(fontsize=7, ncol=4, loc='upper left')
    # PERF7: against the spacing of the jumps
    ax = axs[1]
    fields = {'GPU (RTX 4070)': 'gpu', 'gfortran core': 'cpu-gnu'}
    spacings = [0, 2048, 64, 8]
    for j, (d, tag) in enumerate(fields.items()):
        for name, lab, ls in (('SETH teno5_dr', 'DR', '-'), ('SETH teno5_thinc', 'THINC', '--'),
                              ('SETH teno5_dr, stored masks', 'DR, stored', ':')):
            ys = []
            for js in spacings:
                f = M5 / (f'42-fields-{tag}-js{js}.txt' if js != 2048 else f'42-nonpoly-{tag}-run-1.txt')
                t = table(f) if f.exists() else {}
                ys.append(lookup(t, name) / t[TENO5] if TENO5 in t else np.nan)
            ax.plot(range(len(spacings)), ys, ls, marker='o', ms=3, color=cols[j], lw=1.1, label=f'{lab}, {d}')
    ax.set_xticks(range(len(spacings)), ['smooth', '2048', '64', '8'])
    ax.set_xlabel('cells between jumps')
    ax.set_title('the composites against the density of jumps', fontsize=10)
    ax.legend(fontsize=6.5)
    for a in axs:
        a.grid(axis='y', alpha=0.25, which='both')
        a.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(out / 'perf-nonpoly.png', dpi=DPI)
    plt.close(fig)


def fig_device_rules(out: Path) -> None:
    """The same kernel following or breaking the two device rules: GPU and one CPU core, log scale."""
    gpu = table(AUDIT / 'gpu1-after-fix-with-auto.txt')
    gnu, nvf = table(AUDIT / 'cpu-gfortran.txt'), table(AUDIT / 'cpu-nvfortran.txt')
    rows = [  # (label, variant name in the outputs, colour)
        ('WENO5-JS, library: constant arrays, constant loops', 'weno5_js library (const arrays, const loops)', GOOD),
        ('WENO5-JS, S_MAX arrays, run-time loops', 'weno5_js S_MAX arrays, run-time loops', MILD),
        ('WENO5-JS, automatic arrays, constant loops', 'weno5_js automatic arrays, const loops', SLOW),
        ('WENO5-JS, automatic arrays, run-time loops', 'weno5_js automatic arrays, run-time loops', SLOW),
        (None, None, None),
        ('limiter 1-D, library: NV_MAX arrays, loops to nv', 'positivity nv=3 library (NV_MAX, loops to nv)', GOOD),
        ('limiter 1-D, nv constant', 'positivity nv=3 const arrays, const loops', GOOD),
        ('limiter 1-D, first version: automatic arrays', 'positivity nv=3 automatic arrays (ea4959b)', SLOW),
        (None, None, None),
        ('limiter 3-D, library: NV_MAX arrays, loops to nv', 'positivity nv=5 library (NV_MAX, loops to nv)', GOOD),
        ('limiter 3-D, nv constant', 'positivity nv=5 const arrays, const loops', GOOD),
        ('limiter 3-D, first version: automatic arrays', 'positivity nv=5 automatic arrays (ea4959b)', SLOW),
    ]
    fig, axs = plt.subplots(1, 3, figsize=(10.4, 4.6), sharey=True)
    y = np.arange(len(rows))[::-1]
    for ax, data, title in ((axs[0], gpu, 'RTX 4070 (GPU)'), (axs[1], gnu, 'one CPU core, gfortran'),
                            (axs[2], nvf, 'one CPU core, nvfortran')):
        for yy, (lab, name, col) in zip(y, rows):
            if name is None:
                continue
            v = data[name]
            ax.barh(yy, v, 0.7, color=col, edgecolor='k', lw=0.4)
            ax.text(v * 1.12, yy, f'{v:.3g}', va='center', ha='left', fontsize=7)
        ax.set_xscale('log')
        ax.set_title(title, fontsize=10)
        ax.set_xlabel('ns per face (log)')
        ax.grid(axis='x', alpha=0.25, which='both')
        ax.set_axisbelow(True)
    axs[0].set_xlim(0.5, 9000)
    axs[1].set_xlim(10, 600)
    axs[2].set_xlim(10, 600)
    axs[0].set_yticks([yy for yy, r in zip(y, rows) if r[0]], [r[0] for r in rows if r[0]], fontsize=8)
    fig.tight_layout()
    fig.savefig(out / 'perf-device-rules.png', dpi=DPI)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--out', type=Path, default=ROOT / 'docs/public/guide')
    ap.add_argument('--theme', choices=THEMES, default='both')
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    global SETH, SLOW, MILD, GOOD  # noqa: PLW0603 - the drawing functions read the module colours
    for th in themes(args.theme):
        dark = th == 'dark'
        SETH = GOOD = P.BLUE_D if dark else P.BLUE
        SLOW = P.UMBER_D if dark else P.UMBER
        MILD = '#9cc6f0' if dark else '#8db6e0'
        with theme(th):
            fig_kernels(args.out)
            fig_device_rules(args.out)
            fig_nonpoly(args.out)
        print(f'figures ({th}): {args.out}/perf-{{kernels,device-rules,nonpoly}}.png')
    return 0


if __name__ == '__main__':
    sys.exit(main())
