#!/usr/bin/env python3
"""Draw the figures and write the result tables of SETH's verification documentation from the report data.

Why this script exists: the verification pages show real results of the verification programs (src/tests/seth_vv_*
and the report-capable seth_test_* programs), not sketches. scripts/run_vv.sh runs them with SETH_VV_DIR set; they
write their data there, and this script turns it into

* figures, docs/public/verification/<name>.png, each with its silver-ink twin <name>-dark.png drawn by the same code
  inside figure_theme.theme("dark") (the harness of the ADAM and FLUME documentation, copied in docs/scripts/);
* result tables, docs/verification/_tables/<name>.md, included by the case pages (so the numbers on the pages are
  always those of the last report run).

A figure or table whose data is missing is skipped with a message.

Usage: make_doc_figures.py --data DIR [--out DIR] [--tables DIR] [--theme {light,dark,both}] [--only NAME ...]
Needs numpy and matplotlib.
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import numpy as np

mpl.use('Agg')
# the washi ground of the ADAM docs identity: the figures read as part of the page, as ADAM's and FLUME's do
mpl.rcParams['figure.facecolor'] = '#f3ede1'
mpl.rcParams['savefig.facecolor'] = '#f3ede1'
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_theme import THEMES, theme, themes  # noqa: E402
from seth_palette import SCHEME_STYLE, SCHEME_STYLE_DARK  # noqa: E402

DPI = 110
# one colour and marker per scheme family, the same on every figure (seth_palette: the Neon Moegi identity)
STYLE = dict(SCHEME_STYLE)  # swapped for SCHEME_STYLE_DARK while the dark twins are drawn (main)
LABEL = {'weno1_js': 'WENO1-JS', 'weno3_js': 'WENO3-JS', 'weno5_js': 'WENO5-JS', 'weno7_js': 'WENO7-JS', 'weno9_js': 'WENO9-JS',
         'wenoz3': 'WENO-Z3', 'wenoz5': 'WENO-Z5', 'wenoz7': 'WENO-Z7', 'wenoz9': 'WENO-Z9', 'teno5': 'TENO5',
         'teno6': 'TENO6', 'teno6_opt': 'TENO6-opt', 'teno6_neno6': 'TENO6-NENO6', 'teno5_opt': 'TENO5-opt', 'teno5_az': 'TENO5-AZ', 'teno8': 'TENO8', 'teno5_a': 'TENO5-A', 'teno8_a': 'TENO8-A', 'teno5_dr': 'TENO5-DR', 'teno5_thinc': 'TENO5-THINC', 'teno5_st': 'TENO5-ST', 'teno5_rbf': 'TENO5-RBF', 'teno5_st_dr': 'TENO5-ST+DR', 'teno5_rbf_dr': 'TENO5-RBF+DR',
         'weno4_c': 'WENO4-C', 'weno6_c': 'WENO6-C', 'weno8_c': 'WENO8-C'}


def base(s: str) -> str:
    return s[:-7] if s.endswith('_interp') else s


def style(s: str) -> tuple[str, str]:
    return STYLE.get(base(s), ('0.4', '.'))


def label(s: str) -> str:
    return LABEL.get(base(s), s) + (' interp.' if s.endswith('_interp') else '')


def fmt(x: float) -> str:
    """A number for a table: 4 significant digits, scientific below 1e-3."""
    if not math.isfinite(x) or x < 0:
        return '—'
    return f'{x:.3e}' if (x != 0 and abs(x) < 1e-3) or abs(x) >= 1e4 else f'{x:.4g}'


class Figures:
    def __init__(self, data: Path, out: Path, tables: Path) -> None:
        self.data, self.out, self.tables = data, out, tables
        self.out.mkdir(parents=True, exist_ok=True)
        self.tables.mkdir(parents=True, exist_ok=True)
        self.dark = False

    # ---------------------------------------------------------------------------------------------------------- io
    def has(self, *names: str) -> bool:
        missing = [n for n in names if not (self.data / n).exists()]
        if missing:
            print(f'  skipped: missing {", ".join(missing)}')
        return not missing

    def load(self, name: str) -> np.ndarray:
        return np.loadtxt(self.data / name, ndmin=2)

    def rows(self, name: str) -> list[list[str]]:
        return [ln.split() for ln in (self.data / name).read_text().splitlines() if ln and not ln.startswith('#')]

    def save(self, fig: plt.Figure, name: str) -> None:
        fig.savefig(self.out / f'{name}.png', dpi=DPI, bbox_inches='tight')
        plt.close(fig)

    def table(self, name: str, header: list[str], body: list[list[str]], note: str = '') -> None:
        if self.dark:
            return  # tables once, in the light pass
        lines = ['| ' + ' | '.join(header) + ' |', '|' + '---|' * len(header)]
        lines += ['| ' + ' | '.join(r) + ' |' for r in body]
        if note:
            lines += ['', note]
        (self.tables / f'{name}.md').write_text('\n'.join(lines) + '\n')

    # ----------------------------------------------------------------------------------------------- shock tubes
    def riemann(self, case: str, title: str, show: list[str], zoom: tuple[float, float] | None) -> None:
        files = [f'{case}-exact.dat'] + [f'{case}-{s}.dat' for s in show]
        if not self.has(*files):
            return
        ex = self.load(f'{case}-exact.dat')
        fig, axs = plt.subplots(1, 3, figsize=(12, 3.6))
        for ax, j, nm in zip(axs, (1, 2, 3), ('density', 'velocity', 'pressure')):
            ax.plot(ex[:, 0], ex[:, j], 'k-', lw=1.0, label='exact')
            for s in show:
                d = self.load(f'{case}-{s}.dat')
                c, m = style(s)
                interp = s.endswith('_interp')  # interpolation path: filled, smaller markers
                ax.plot(d[:, 0], d[:, j], m, ms=2.2 if interp else 3, mfc=c if interp else 'none', mec=c, mew=0.8,
                        label=label(s))
            ax.set_xlabel('$x$')
            ax.set_title(nm, fontsize=10)
            ax.grid(alpha=0.3)
        if zoom:
            ins = axs[0].inset_axes([0.55, 0.5, 0.42, 0.42])
            ins.plot(ex[:, 0], ex[:, 1], 'k-', lw=1.0)
            for s in show:
                d = self.load(f'{case}-{s}.dat')
                c, m = style(s)
                ins.plot(d[:, 0], d[:, 1], m, ms=3, mfc='none', mec=c, mew=0.8)
            sel = (ex[:, 0] >= zoom[0]) & (ex[:, 0] <= zoom[1])
            lo, hi = ex[sel, 1].min(), ex[sel, 1].max()
            ins.set_xlim(*zoom)
            ins.set_ylim(lo - 0.08 * (hi - lo), hi + 0.08 * (hi - lo))
            ins.tick_params(labelsize=7)
            axs[0].indicate_inset_zoom(ins, edgecolor='0.4')
        axs[2].legend(fontsize=8, loc='best')
        fig.suptitle(title, fontsize=11)
        fig.tight_layout()
        self.save(fig, case)
        self.riemann_table(case)

    def riemann_table(self, case: str) -> None:
        if not self.has(f'{case}-l1.dat', f'{case}-summary.dat'):
            return
        l1 = defaultdict(dict)
        for s, n, e in self.rows(f'{case}-l1.dat'):
            l1[s][int(n)] = float(e)
        summ = {r[0]: r[1:] for r in self.rows(f'{case}-summary.dat')}
        body = []
        for s in l1:
            ns = sorted(l1[s])
            o = (math.log2(l1[s][ns[-2]] / l1[s][ns[-1]])
                 if len(ns) > 1 and l1[s][ns[-1]] > 0 and l1[s][ns[-2]] > 0 else float('nan'))
            pos, nlim, _ = summ.get(s, ['?', '?', '0'])
            body.append([label(s), fmt(l1[s].get(200, -1)), fmt(l1[s].get(1600, -1)),
                         f'{o:.2f}' if math.isfinite(o) else '—', 'yes' if pos == 'T' else '**no**', nlim])
        self.table(case, ['scheme', '$L_1(\\rho)$, $N=200$', '$L_1(\\rho)$, $N=1600$', 'order 800→1600',
                          'positive without limiter', 'limited faces ($N=200$)'], body,
                   '$L_1$ with the positivity limiter on (it is inactive where the last column is 0).')

    def riemann_convergence(self) -> None:
        cases = [c for c in ('sod', 'lax', '123') if (self.data / f'{c}-l1.dat').exists()]
        if not cases:
            print('  skipped: no <case>-l1.dat')
            return
        fig, axs = plt.subplots(1, len(cases), figsize=(4.2 * len(cases), 3.6), squeeze=False)
        for ax, c in zip(axs[0], cases):
            l1 = defaultdict(dict)
            for s, n, e in self.rows(f'{c}-l1.dat'):
                if float(e) > 0:
                    l1[s][int(n)] = float(e)
            for s, d in l1.items():
                col, m = style(s)
                ns = sorted(d)
                ax.loglog(ns, [d[n] for n in ns], m + ('--' if s.endswith('_interp') else '-'), color=col, ms=3,
                          lw=0.8, label=label(s))
            n0 = np.array([100, 1600])
            ref = min(min(d.values()) for d in l1.values())
            ax.loglog(n0, ref * 3 * (n0 / 1600.0) ** -1, 'k:', lw=0.8, label='slope 1')
            ax.set_title({'sod': 'Sod', 'lax': 'Lax', '123': '123 problem'}[c], fontsize=10)
            ax.set_xlabel('$N$')
            ax.grid(alpha=0.3, which='both')
        axs[0][0].set_ylabel('$L_1(\\rho)$')
        axs[0][-1].legend(fontsize=6, ncol=2, loc='lower left')
        fig.tight_layout()
        self.save(fig, 'riemann-convergence')

    # ------------------------------------------------------------------------------------- reference problems
    def reference_case(self, case: str, title: str, show: list[str], zoom: tuple[float, float], var: int = 1) -> None:
        files = [f'{case}-reference.dat'] + [f'{case}-{s}.dat' for s in show]
        if not self.has(*files):
            return
        ref = self.load(f'{case}-reference.dat')
        fig, axs = plt.subplots(1, 2, figsize=(12, 3.8), gridspec_kw={'width_ratios': [1.3, 1]})
        for k, ax in enumerate(axs):
            ax.plot(ref[:, 0], ref[:, var], 'k-', lw=0.9, label='reference')
            for s in show:
                d = self.load(f'{case}-{s}.dat')
                c, m = style(s)
                ls = '--' if s.endswith('_interp') else '-'
                ax.plot(d[:, 0], d[:, var], m + ls, color=c, ms=2.5 if k else 0, lw=0.8, mfc='none', label=label(s))
            if k:
                sel = (ref[:, 0] >= zoom[0]) & (ref[:, 0] <= zoom[1])
                lo, hi = ref[sel, var].min(), ref[sel, var].max()
                ax.set_xlim(*zoom)
                ax.set_ylim(lo - 0.08 * (hi - lo), hi + 0.08 * (hi - lo))
                ax.set_title('zoom', fontsize=10)
            else:
                ax.set_title(title, fontsize=10)
            ax.set_xlabel('$x$')
            ax.grid(alpha=0.3)
        axs[0].set_ylabel('density')
        axs[1].legend(fontsize=8)
        fig.tight_layout()
        self.save(fig, case)

    def reference_table(self, case: str) -> None:
        if not self.has(f'{case}-summary.dat'):
            return
        body = []
        for r in self.rows(f'{case}-summary.dat'):
            if len(r) == 4:
                body.append([label(r[0]), fmt(float(r[3])), 'yes' if r[1] == 'T' else '**no**', r[2]])
            else:
                body.append([label(r[0]), fmt(float(r[1])), '—', '—'])
        self.table(case, ['scheme', '$L_1(\\rho)$ against the reference', 'positive without limiter',
                          'limited faces'], body)

    def positivity(self) -> None:
        cases = [c for c in ('blast', '123') if (self.data / f'{c}-summary.dat').exists()]
        if not cases:
            print('  skipped: no summaries')
            return
        fig, axs = plt.subplots(1, len(cases), figsize=(6 * len(cases), 3.8), squeeze=False)
        for ax, c in zip(axs[0], cases):
            rows = self.rows(f'{c}-summary.dat')
            names = [r[0] for r in rows]
            nl = [max(int(r[2]), 0.8) for r in rows]
            cols = [style(n)[0] for n in names]
            bars = ax.bar(range(len(rows)), nl, color=cols, edgecolor='k', linewidth=0.4)
            for b, r in zip(bars, rows):
                if r[1] != 'T':
                    b.set_hatch('///')
            ax.set_yscale('log')
            ax.set_xticks(range(len(rows)))
            ax.set_xticklabels([label(n) for n in names], rotation=60, ha='right', fontsize=7)
            ax.set_ylabel('limited faces')
            ax.set_title({'blast': 'Woodward–Colella blast, $N=400$', '123': '123 problem, $N=200$'}[c], fontsize=10)
            ax.grid(alpha=0.3, axis='y')
        fig.text(0.01, 0.01, 'hatched: loses positivity without the limiter', fontsize=8)
        fig.tight_layout()
        self.save(fig, 'positivity')

    # -------------------------------------------------------------------------------------------------- Burgers
    def burgers(self) -> None:
        if not self.has('burgers-smooth.dat', 'burgers-shock.dat'):
            return
        sm = defaultdict(dict)
        for s, n, e1, ei in self.rows('burgers-smooth.dat'):
            sm[s][int(n)] = float(e1)
        sh = defaultdict(dict)
        for s, n, e, es, osh in self.rows('burgers-shock.dat'):
            sh[s][int(n)] = (float(e), float(es), float(osh))
        fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
        for s, d in sm.items():
            c, m = style(s)
            ns = sorted(d)
            axs[0].loglog(ns, [d[n] for n in ns], m + '-', color=c, ms=3, lw=0.8, label=label(s))
            axs[1].loglog(ns, [sh[s][n][1] for n in ns], m + '-', color=c, ms=3, lw=0.8)
        axs[0].set_title('before the shock, $t = 0.5/\\pi$: $L_1$', fontsize=10)
        axs[1].set_title('after, $t = 1.5/\\pi$: $L_1$ farther than 0.1 from the shock', fontsize=10)
        for ax in axs[:2]:
            ax.set_xlabel('$N$')
            ax.grid(alpha=0.3, which='both')
        axs[0].legend(fontsize=6, ncol=2)
        for s in ('weno5_js', 'teno5', 'teno6'):
            if (self.data / f'burgers-post-{s}.dat').exists():
                d = self.load(f'burgers-post-{s}.dat')
                c, m = style(s)
                axs[2].semilogy(d[:, 0], np.maximum(abs(d[:, 1] - d[:, 2]), 1e-17), m + '-', color=c, ms=2.5, lw=0.7,
                                label=label(s))
        axs[2].set_title('pointwise error after the shock, $N = 80$', fontsize=10)
        axs[2].set_xlabel('$x$')
        axs[2].legend(fontsize=7)
        axs[2].grid(alpha=0.3)
        fig.tight_layout()
        self.save(fig, 'burgers')
        body = []
        for s, d in sm.items():
            ns = sorted(d)
            pairs = [(a, b) for a, b in zip(ns, ns[1:]) if d[b] > 1e-11]  # as the program: above round-off
            o = math.log2(d[pairs[-1][0]] / d[pairs[-1][1]]) if pairs else float('nan')
            body.append([label(s), fmt(d[ns[-1]]), f'{o:.2f}', fmt(sh[s][ns[-1]][0]), fmt(sh[s][ns[-1]][1]),
                         fmt(max(v[2] for v in sh[s].values()))])
        nmax = max(max(d) for d in sm.values())
        self.table('burgers', ['scheme', f'$L_1$ before, $N={nmax}$', 'order (finest pair above 1e-11)',
                               f'$L_1$ after, $N={nmax}$',
                               '$L_1$ after, away from the shock', 'max overshoot'], body)

    # --------------------------------------------------------------------------------------------- Euler order
    def euler_order(self) -> None:
        if not self.has('euler-order.dat'):
            return
        d = defaultdict(dict)
        for s, n, e1, ei in self.rows('euler-order.dat'):
            d[s][int(n)] = float(e1)
        fig, axs = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
        for s, e in d.items():
            ax = axs[1 if s.endswith('_interp') else 0]
            c, m = style(s)
            ns = sorted(e)
            ax.loglog(ns, [e[n] for n in ns], m + '-', color=c, ms=3, lw=0.8, label=label(s))
        for ax, t in zip(axs, ('splitting path (reconstruction)', 'interpolation path (HLLC + flux corrections)')):
            ax.set_title(t, fontsize=10)
            ax.set_xlabel('$N$')
            ax.grid(alpha=0.3, which='both')
            ax.legend(fontsize=6, ncol=2)
        axs[0].set_ylabel('$L_1(\\rho)$')
        fig.tight_layout()
        self.save(fig, 'euler-order')
        body = []
        for s, e in d.items():
            ns = sorted(e)
            pairs = [(a, b) for a, b in zip(ns, ns[1:]) if e[b] > 1e-10]
            o = math.log2(e[pairs[-1][0]] / e[pairs[-1][1]]) if pairs else float('nan')
            body.append([label(s), fmt(e[ns[-1]]), f'{o:.2f}'])
        self.table('euler-order', ['scheme', f'$L_1(\\rho)$, $N={max(ns)}$', 'observed order (finest pair above 1e-10)'],
                   body)

    # ------------------------------------------------------------------------------------------------ advection
    def advection(self) -> None:
        show = ['weno5_js', 'wenoz5', 'teno5', 'teno6']
        if not self.has(*[f'advection-{s}.dat' for s in show], 'advection.dat'):
            return
        fig, axs = plt.subplots(1, 3, figsize=(13, 3.8), gridspec_kw={'width_ratios': [2, 1, 1]})
        ex = self.load('advection-weno5_js.dat')
        for k, (ax, lim) in enumerate(zip(axs, (None, (-0.45, -0.15), (0.35, 0.65)))):
            ax.plot(ex[:, 0], ex[:, 2], 'k-', lw=0.9, label='exact')
            for s in show:
                d = self.load(f'advection-{s}.dat')
                c, m = style(s)
                ax.plot(d[:, 0], d[:, 1], m + '-', color=c, ms=2 if k else 0, lw=0.8, mfc='none', label=label(s))
            if lim:
                ax.set_xlim(*lim)
                ax.set_ylim(-0.05, 1.12)
            ax.set_xlabel('$x$')
            ax.grid(alpha=0.3)
        axs[0].set_title('composite profile after four periods ($t = 8$, $N = 200$)', fontsize=10)
        axs[1].set_title('square wave', fontsize=10)
        axs[2].set_title('half ellipse', fontsize=10)
        axs[0].legend(fontsize=7)
        fig.tight_layout()
        self.save(fig, 'advection')
        body = [[label(r[0]), fmt(float(r[1])), fmt(float(r[2]) - 1.0), fmt(-float(r[3]))]
                for r in self.rows('advection.dat')]
        self.table('advection', ['scheme', '$L_1$', 'max $u$ − 1', '−min $u$'], body)

    # ------------------------------------------------------------------------------------------------------ ADR
    def adr(self) -> None:
        names = [p.name[4:-4] for p in sorted(self.data.glob('adr-*.dat'))]
        if not names:
            print('  skipped: no adr-*.dat')
            return
        fig, axs = plt.subplots(1, 2, figsize=(12, 4))
        th = np.linspace(0, np.pi, 200)
        axs[0].plot(th / np.pi, th, 'k-', lw=0.8, label='exact')
        axs[1].axhline(0, color='k', lw=0.8)
        for s in names:
            d = self.load(f'adr-{s}.dat')
            c, m = style(s)
            axs[0].plot(d[:, 0] / np.pi, d[:, 1], '-', color=c, lw=0.9, label=label(s))
            axs[1].plot(d[:, 0] / np.pi, d[:, 2], '-', color=c, lw=0.9)
            if s in ('teno5', 'teno6'):
                axs[0].plot(d[:, 0] / np.pi, d[:, 3], ':', color=c, lw=0.8)
                axs[1].plot(d[:, 0] / np.pi, d[:, 4], ':', color=c, lw=0.8)
        axs[0].set_title('dispersion, Re $\\tilde\\vartheta$ (dotted: linear background of TENO5/6)', fontsize=10)
        axs[1].set_title('dissipation, Im $\\tilde\\vartheta$', fontsize=10)
        for ax in axs:
            ax.set_xlabel('$\\vartheta/\\pi$')
            ax.grid(alpha=0.3)
        axs[0].legend(fontsize=6, ncol=2)
        fig.tight_layout()
        self.save(fig, 'adr')

    def hybrid(self) -> None:
        if not self.has('hybrid-summary.dat'):
            return
        rows = self.rows('hybrid-summary.dat')
        pure = {(r[0], r[1]): float(r[3]) for r in rows if r[2] == 'pure'}
        body = []
        for r in rows:
            ratio = float(r[3]) / pure[(r[0], r[1])]
            body.append([r[0], label(r[1]), r[2], f'{float(r[3]):.4g}', f'{ratio:.3f}', f'{100 * float(r[4]):.1f} %',
                         'yes' if r[5] == 'T' else '**no**'])
        self.table('hybrid', ['case', 'scheme', 'mode', '$L_1(\\rho)$', 'against pure', 'faces flagged', 'positive'],
                   body)

    def sensors(self) -> None:
        if not self.has('sensors-k.dat', 'sensors-b.dat'):
            return
        names = ['TENO-based (Fu)', 'Ren et al.', 'Li, Chen & Ren', 'Wan & Xia']
        colours = [style('teno6')[0], style('weno5_js')[0], style('wenoz5')[0], style('teno5')[0]]
        k = self.load('sensors-k.dat')
        b = self.load('sensors-b.dat')
        fig, axs = plt.subplots(1, 2, figsize=(12, 4))
        for i, (n, c) in enumerate(zip(names, colours)):
            axs[0].plot(k[:, 0], k[:, i + 1], '-', color=c, lw=1.0, label=n)
            axs[1].plot(b[:, 0], 1.0 - b[:, i + 2] + 1.25 * i, '-', color=c, lw=0.9, label=n)
        for x, txt in ((0.64, '0.64'), (2.06, '2.06')):
            axs[0].axvline(x, color='0.5', lw=0.6, ls=':')
            axs[0].text(x, 0.05, txt, fontsize=7, color='0.4')
        axs[0].set_title('sine of wavenumber $k$, 192 points: fraction of faces judged smooth', fontsize=10)
        axs[0].set_xlabel('$k$')
        axs[0].set_ylim(-0.02, 1.05)
        ax2 = axs[1].twinx()
        ax2.plot(b[:, 0], b[:, 1], '-', color='0.55', lw=0.7)
        ax2.set_ylabel('$f$', color='0.45')
        axs[1].set_title('Li, Chen & Ren test (b): $1 - s$ per sensor (offset), and the data (grey)', fontsize=10)
        axs[1].set_xlabel('$x$')
        axs[1].set_yticks([])
        axs[1].legend(fontsize=6, loc='upper left')
        for ax in axs:
            ax.grid(alpha=0.3)
        axs[0].legend(fontsize=7)
        fig.tight_layout()
        self.save(fig, 'sensors')

    # ------------------------------------------------------------------------------ non-polynomial (issue #6)
    FLAG_PROFILES = {'profile1': 'profile 1 (Table 1)', 'profile2': 'profile 2 (Table 2)',
                     'profile2m': 'profile 2, mirrored in value', 'profile3': 'profile 3 (Table 3)',
                     'profile3m': 'profile 3, mirrored in value'}

    def nonpoly_flags(self) -> None:
        """The flags of Liang et al. 2022, Tables 1-3, against seth_teno5_mask (seth_test_nonpoly, spec 2)."""
        if self.has('nonpoly-flags.dat'):
            rows = self.rows('nonpoly-flags.dat')
            prof: dict[str, list[list[str]]] = defaultdict(list)
            for r in rows:
                prof[r[0]].append(r)
            body = []
            for key, title in self.FLAG_PROFILES.items():
                rs = prof[key]
                f = ' '.join(fmt(float(r[2])) for r in rs)
                got = [''.join(r[3 + k] for r in rs[2:9]) for k in range(3)]
                ref = [''.join(r[6 + k] for r in rs[2:9]) for k in range(3)]
                body.append([title, f'`{f}`'] + [f'`{g}`' for g in got] + ['yes' if got == ref else '**no**'])
            self.table('nonpoly-flags', ['profile (cells 0..10)', 'values', r'$\delta_0$ (cells 2..8)',
                                         r'$\delta_1$', r'$\delta_2$', 'as printed'], body)
            keys = ['profile1', 'profile2', 'profile3']
            fig, axs = plt.subplots(2, 3, figsize=(12, 3.8), gridspec_kw={'height_ratios': [2.2, 1.2]}, sharex=True)
            for j, key in enumerate(keys):
                rs = prof[key]
                c = np.array([int(r[1]) for r in rs])
                f = np.array([float(r[2]) for r in rs])
                ax = axs[0, j]
                ax.step(c, f, where='mid', color=style('teno5')[0], lw=1.2)
                ax.plot(c, f, 'o', color=style('teno5')[0], ms=3)
                ax.set_title(self.FLAG_PROFILES[key], fontsize=10)
                ax.set_ylim(-0.15, 1.15)
                ax.grid(alpha=0.3)
                sx = axs[1, j]
                for k, nm in enumerate(('$\\delta_0$', '$\\delta_1$', '$\\delta_2$')):
                    d = np.array([int(r[3 + k]) for r in rs])
                    ref = np.array([int(r[6 + k]) for r in rs])
                    y = 2 - k
                    for ci, di, ri in zip(c, d, ref):
                        fc = style('teno6')[0] if di else 'none'
                        ec = style('teno6')[0] if di == ri else 'red'
                        sx.add_patch(mpl.patches.Rectangle((ci - 0.4, y - 0.35), 0.8, 0.7, facecolor=fc,
                                                           edgecolor=ec, lw=0.8))
                    sx.text(-0.9, y, nm, fontsize=8, ha='right', va='center')
                sx.set_ylim(-0.6, 2.6)
                sx.set_yticks([])
                sx.set_xlim(-1.6, 10.6)
                sx.set_xlabel('cell')
            fig.suptitle('TENO5 flags per cell ($C_T = 10^{-5}$): filled = candidate kept; asserted on cells 2..8',
                         fontsize=10)
            fig.tight_layout()
            self.save(fig, 'nonpoly-flags')
        if self.has('nonpoly-masks.dat'):
            rows = self.rows('nonpoly-masks.dat')
            body = [[r[0]] + r[1:] for r in rows]
            self.table('nonpoly-masks', ['stencils'] + [f'{m} ({m & 1}{m >> 1 & 1}{m >> 2 & 1})' for m in range(8)],
                       body, note='Mask $m = \\delta_0 + 2\\delta_1 + 4\\delta_2$, in parentheses '
                                  '$(\\delta_0\\delta_1\\delta_2)$. *random*: the $10^6$ stencils of spec 1 '
                                  '(smooth cubics, independent values, jumps, kinks); *constant*: 1000 constants drawn '
                                  'in $[-10, 10]$.')

    def nonpoly_census(self) -> None:
        """seth_test_nonpoly, spec 17 (E7): the masks of TENO5's cutoff by data."""
        if self.has('nonpoly-census.dat'):
            body = []
            for r in self.rows('nonpoly-census.dat'):
                n = int(r[1])
                body.append([r[0], f'{n:.0e}'.replace('e+0', 'e')] +
                            [f'{int(c) / n:.2%}' if int(c) else '0' for c in r[2:10]])
            self.table('nonpoly-census', ['data', 'stencils'] + [f'{m} ({m & 1}{m >> 1 & 1}{m >> 2 & 1})'
                                                                 for m in range(8)], body,
                       note='Mask $m = \\delta_0 + 2\\delta_1 + 4\\delta_2$, in parentheses '
                            '$(\\delta_0\\delta_1\\delta_2)$; the fraction of the stencils.')

    def nonpoly_st(self) -> None:
        """seth_test_nonpoly, spec 14: stencil targeting's orders."""
        if self.has('nonpoly-st.dat'):
            v = {r[0]: r[1:] for r in self.rows('nonpoly-st.dat')}
            body = [['smooth data (the generic $g$), left / right state', f"{float(v['orders-smooth'][0]):.3f} / "
                     f"{float(v['orders-smooth'][1]):.3f}", 'equal to `seth_linear5` within '
                     f"{fmt(float(v['st-minus-linear5'][0]))}"],
                    ['beside a jump, mask (1,1,0)', f"{float(v['orders-jump'][0]):.3f}",
                     f"TENO5: {float(v['orders-jump'][1]):.3f}"],
                    ['across a unit jump', '—', f"overshoot {fmt(float(v['overshoot'][0]))}"],
                    ['RBF, $s = 3\\,dx$, left / right state', f"{float(v['orders-smooth'][2]):.3f} / "
                     f"{float(v['orders-smooth'][3]):.3f}", '$\\alpha$ fixed, $dx \\to 0$']]
            self.table('nonpoly-st', ['case', 'observed order', ''], body)

    def nonpoly_adr(self) -> None:
        """seth_test_nonpoly, spec 18 (G6, E6): the long-wave growth of the linear RBF row, against the prediction."""
        if self.has('nonpoly-adr.dat'):
            txt = (self.data / 'nonpoly-adr.dat').read_text().splitlines()
            rows = [ln.split() for ln in txt if ln and not ln.startswith('#')]
            valid = [ln.split()[-1] for ln in txt if ln.startswith('# validated-s')]
            body = [[fmt(float(r[0])), fmt(float(r[0]) ** 0.5), f'{float(r[1]):.4f}', f'{float(r[2]):.4f}',
                     fmt(float(r[3])), fmt(float(r[4]))] for r in rows]
            note = (f'The largest $s$ whose peak stays at or below $2\\cdot10^{{-4}}$ (TENO6\'s): {valid[0]}.'
                    if valid else '')
            self.table('nonpoly-adr', ['$s^2$', '$s$', 'band edge $\\theta$', 'predicted $\\sqrt{30}\\,s$',
                                       'peak of Im $\\tilde\\theta$', 'predicted $200 s^6/3$'], body, note=note)

    def nonpoly_advection(self) -> None:
        """seth_vv_nonpoly, sections 3-4: N4, N6 (the halo), E2 (contact spreading), E3 (smooth gradients)."""
        names = ['teno5', 'teno5_dr', 'teno5_thinc']
        if self.has('nonpoly-halo.dat'):
            body = [[label(r[0]), r[1], fmt(float(r[2])), 'yes' if r[3] == 'T' else '**no**',
                     'yes' if r[4] == 'T' else '**no**'] for r in self.rows('nonpoly-halo.dat')]
            self.table('nonpoly-halo', ['scheme', 'states overridden (block 1)', 'sum of $u$: relative change',
                                        'two blocks, 6 ghost cells: bitwise', 'two blocks, 5 ghost cells: differs'],
                       body)
        if self.has('nonpoly-e2.dat'):
            d = self.load('nonpoly-e2.dat')
            self.table('nonpoly-e2', ['period'] + [label(n) for n in names],
                       [[str(int(r[0]))] + [str(int(x)) for x in r[1:]] for r in d],
                       note='Cells whose value lies between 5 % and 95 % of the unit jump, both jumps of the square '
                            'wave, $N = 100$, CFL 0.4.')
            fig, ax = plt.subplots(figsize=(6.4, 3.4))
            for j, n in enumerate(names):
                c, mk = style(n)
                ax.plot(d[:, 0], d[:, j + 1], '-', color=c, marker=mk, ms=4, lw=1.2, label=label(n))
            ax.set_xlabel('periods')
            ax.set_ylabel('cells across the two jumps')
            ax.set_ylim(0, None)
            ax.grid(alpha=0.3)
            ax.legend(fontsize=8)
            ax.set_title('E2: a square wave advected 10 periods, $N = 100$', fontsize=10)
            fig.tight_layout()
            self.save(fig, 'nonpoly-e2')
        if self.has('nonpoly-e3.dat'):
            d = self.load('nonpoly-e3.dat')
            self.table('nonpoly-e3', ['gradient width (cells)'] + [f'$L_1$ {label(n)}' for n in names],
                       [[f'{r[0]:g}'] + [fmt(x) for x in r[1:]] for r in d],
                       note='A smooth periodic plateau, $0.5(\\tanh((x - 0.25)/(w h)) - \\tanh((x - 0.75)/(w h)))$, '
                            'advected one period, $N = 100$.')

    def nonpoly_static(self) -> None:
        """seth_test_nonpoly, spec 19: Liang et al. 2022's static tests, where each detector fires."""
        if self.has('nonpoly-static.dat'):
            act: dict[tuple[str, str], dict[str, list[str]]] = {}
            for r in self.rows('nonpoly-static.dat'):
                act.setdefault((r[0], r[1]), {})[r[2]] = r[3:]
            body = []
            for (case, n), d in act.items():
                dr, th = d.get('dr', []), d.get('thinc', [])
                body.append([f'({case})', n, str(len(dr)), ' '.join(f'{float(x):.3f}' for x in dr) or '—',
                             str(len(th)), ' '.join(f'{float(x):.3f}' for x in th) or '—'])
            self.table('nonpoly-static', ['test', '$N$', 'DR', 'where ($x$)', 'THINC', 'where ($x$)'], body)

    def nonpoly_ablation(self) -> None:
        """seth_vv_nonpoly, sections 5-6: N5, the ablation, and E5, the length unit."""
        if self.has('nonpoly-ablation.dat'):
            rows = self.rows('nonpoly-ablation.dat')
            cases = list(dict.fromkeys(r[0] for r in rows))
            schemes = list(dict.fromkeys(r[1] for r in rows))
            val = {(r[0], r[1]): r for r in rows}
            body = []
            for sch in schemes:
                line = [label(sch)]
                for c in cases:
                    r = val[(c, sch)]
                    line.append(fmt(float(r[2])) + ('' if r[4] == 'T' else ' **not positive**'))
                line.append(val[('blast', sch)][3] if ('blast', sch) in val else '—')
                body.append(line)
            self.table('nonpoly-ablation', ['scheme'] + [f'$L_1(\\rho)$ {self.CASE_NAMES.get(c, c)}' for c in cases]
                       + ['blast: faces limited'], body,
                       note='Sod, Lax: $N = 200$ against the exact solution; Shu–Osher: $N = 200$ on $[-5, 5]$ against '
                            'WENO5-JS on 3200 cells; blast waves: $N = 400$, the positivity limiter on, against WENO5-JS '
                            'on 6400 cells. The RBF kernels with $s = 3\\,\\Delta x$ (Liang & Fu\'s $\\alpha = 3$): '
                            '0.15 on Shu–Osher, outside the validated range $s \\le 0.12$.')
        if self.has('nonpoly-units.dat'):
            body = [[label(r[0]), fmt(float(r[1])), f'{float(r[2]):.3f}', f'{float(r[3]):.4f}']
                    for r in self.rows('nonpoly-units.dat')]
            self.table('nonpoly-units', ['scheme', 'max density difference', '$s$ on $[-5, 5]$', '$s$ scaled by 1/8'],
                       body)

    def nonpoly_e6(self) -> None:
        """seth_vv_nonpoly, section 7 (E6): the growth of a sine over 100 periods under the RBF scheme."""
        if self.has('nonpoly-e6.dat'):
            body = []
            for r in self.rows('nonpoly-e6.dat'):
                s = float(r[3])
                body.append([r[0], r[1], r[2], f'{s:.3f}' + ('' if s <= 0.12 else ' (outside)'), f'{float(r[4]):.3e}',
                             f'{float(r[5]):.3e}', f'{float(r[6]):.3e}'])  # signed: a negative rate is a decay
            self.table('nonpoly-e6', ['$\\lambda$', '$N$', 'CFL', '$s$', 'growth rate, RBF', 'growth rate, TENO5',
                                      'predicted $\\lambda$ term'], body,
                       note='$\\sin 2\\pi x$ advected 100 periods, SSP-RK3; the growth rate is $\\log(A(t)/A(0))/t$; '
                            'the prediction is the $\\lambda$ term of Eq. 29, $\\Delta x^5 k^4 \\lambda^2/2$.')

    NP2D = ['teno5', 'teno5_st', 'teno5_rbf', 'teno5_thinc', 'teno5_dr', 'teno5_st_dr', 'teno5_rbf_dr']
    NP2D_CASES = {'contact': 'oblique contact, one period', 'kh': 'Kelvin–Helmholtz, $t = 0.5$',
                  'config3': '2-D Riemann problem, configuration 3, $t = 0.3$'}
    NP2D_WINDOW = {'config3': (0.15, 0.65)}  # the region of the small-scale structures (the lobe on the diagonal)

    def nonpoly_2d(self) -> None:
        """seth_vv_2d in report mode: the density fields at N = 400, one figure per case, a panel per scheme."""
        if self.has('nonpoly2d-summary.dat'):
            short = {'contact': 'oblique contact', 'kh': 'Kelvin–Helmholtz', 'config3': 'configuration 3'}
            body = [[short.get(r[0], r[0]), label(r[1]), r[2], 'yes' if r[3] == 'T' else '**no**',
                     fmt(float(r[4])) if r[0] != 'config3' else '— (transmissive)',
                     fmt(float(r[5])) if r[0] != 'kh' else '—']
                    for r in self.rows('nonpoly2d-summary.dat')]
            self.table('nonpoly-2d', ['case', 'scheme', '$N$', 'positive', 'conservation', 'symmetry'], body)
        self.nonpoly_2d_fields()
        for case, title in self.NP2D_CASES.items():
            have = [sc for sc in self.NP2D if (self.data / f'nonpoly2d-{case}-{sc}.dat').exists()]
            if not have:
                continue
            fields = {sc: self.load(f'nonpoly2d-{case}-{sc}.dat') for sc in have}
            lo = min(float(d.min()) for d in fields.values())
            hi = max(float(d.max()) for d in fields.values())
            ncol = min(4, len(have))
            nrow = -(-len(have) // ncol)
            fig, axs = plt.subplots(nrow, ncol, figsize=(3.2 * ncol + 0.8, 3.3 * nrow + 0.4), squeeze=False)
            for ax in axs.flat:
                ax.set_axis_off()
            for ax, (sc, d) in zip(axs.flat, fields.items()):
                im = ax.imshow(d, origin='lower', extent=(0, 1, 0, 1), vmin=lo, vmax=hi,
                               cmap='viridis' if not self.dark else 'magma', interpolation='nearest')
                ax.set_title(label(sc), fontsize=9)
                if case in self.NP2D_WINDOW:
                    ax.set_xlim(*self.NP2D_WINDOW[case])
                    ax.set_ylim(*self.NP2D_WINDOW[case])
            # one colour scale for every scheme: the panels compare by eye
            fig.colorbar(im, ax=axs.ravel().tolist(), shrink=0.8, label='density')
            n = next(iter(fields.values())).shape[0]
            win = self.NP2D_WINDOW.get(case)
            where = f', window $[{win[0]}, {win[1]}]^2$' if win else ''
            fig.suptitle(f'{title}: density, $N = {n}${where}', fontsize=10)
            self.save(fig, f'nonpoly-2d-{case}')

    def nonpoly_2d_fields(self) -> None:
        """Measures of the N = 400 density fields: the contact's width and bounds, configuration 3's ripples."""
        body = []
        for sc in self.NP2D:
            con, cf3 = (self.data / f'nonpoly2d-{c}-{sc}.dat' for c in ('contact', 'config3'))
            if not (con.exists() and cf3.exists()):
                continue
            d = self.load(con.name)
            # cells strictly between 5 % and 95 % of the jump 1 -> 1.5, per row (each row crosses two jumps)
            width = ((d > 1.025) & (d < 1.475)).sum(axis=1).mean() / 2
            over = max(float(d.max()) - 1.5, 0.0) + max(1.0 - float(d.min()), 0.0)
            q = self.load(cf3.name)
            n = q.shape[0]
            q = q[int(0.55 * n):int(0.80 * n), int(0.55 * n):int(0.80 * n)]  # behind the shock, smooth
            ext = 0
            for ax in (0, 1):  # local extrema along x and y: sign changes of slopes above the data's 5 digits
                g = np.diff(q, axis=ax)
                sg = np.sign(np.where(np.abs(g) > 1e-4, g, 0.0))
                ext += int((np.abs(np.diff(sg, axis=ax)) == 2).sum())
            k = 7  # ripple: rms of the field minus its 7 x 7 running mean
            pad = np.pad(q, k // 2, mode='edge')
            mean = np.lib.stride_tricks.sliding_window_view(pad, (k, k)).mean(axis=(-1, -2))
            rip = float(np.sqrt(((q - mean) ** 2).mean()))
            body.append([label(sc), f'{width:.1f}', fmt(over) if over > 0 else '0', str(ext), f'{rip:.1e}'])
        if body:
            self.table('nonpoly-2d-fields', ['scheme', 'contact: cells per jump', 'contact: over/undershoot',
                                             'config. 3: local extrema', 'config. 3: ripple (rms)'], body,
                       note='Contact after one period: the cells strictly between 5 % and 95 % of the jump, per row '
                            'and jump; the over/undershoot to the 5 digits of the stored fields. Configuration 3, '
                            'the box $[0.55, 0.80]^2$ behind the shock, where the flow is smooth: the local extrema '
                            'along $x$ and $y$ (slopes above $10^{-4}$) and the rms of the density minus its '
                            '$7\\times 7$ running mean.')

    CASE_NAMES = {'sod': 'Sod', 'lax': 'Lax', 'shu-osher': 'Shu–Osher', 'blast': 'blast waves'}

    def nonpoly_solver(self) -> None:
        """seth_vv_nonpoly: the contact at N = 100 (step 2.2, spec 7) and experiment E1 (the two mask variants)."""
        if self.has('nonpoly-contact.dat'):
            body = [[self.CASE_NAMES.get(r[0], r[0]), label(r[1]), fmt(float(r[2])), r[3], 'yes' if r[4] == 'T' else '**no**']
                    for r in self.rows('nonpoly-contact.dat')]
            self.table('nonpoly-contact', ['case', 'scheme', '$L_1(\\rho)$', 'cells across the contact', 'positive'],
                       body)
        if self.has('nonpoly-e1.dat'):
            body = []
            for r in self.rows('nonpoly-e1.dat'):
                a, b = float(r[2]), float(r[3])
                body.append([self.CASE_NAMES.get(r[0], r[0]), label(r[1]), fmt(a) if a >= 0 else '—', fmt(b) if b >= 0 else '—',
                             fmt(float(r[4])), 'yes' if r[5] == 'T' and r[6] == 'T' else '**no**'])
            self.table('nonpoly-e1', ['case', 'scheme', '$L_1(\\rho)$ per face', '$L_1(\\rho)$ stored masks',
                                      'max $\\lvert\\Delta\\rho\\rvert$ between them', 'positive'], body)

    SENSOR_NAMES = {'teno': 'TENO-based (Fu)', 'ren': 'Ren et al.', 'lcr': 'Li, Chen & Ren', 'tau': 'Wan & Xia',
                    'ducros': 'Ducros', 'ducros-c': 'Ducros, compression'}

    def sensor_static(self) -> None:
        """Tables of the critical points and of the scaled data (seth_test_sensors, sections 4-5)."""
        if self.has('sensors-cp.dat'):
            rows = self.rows('sensors-cp.dat')
            ms = sorted({int(r[1]) for r in rows})
            body = []
            for ncp in ('1', '2'):
                for j, key in enumerate(('teno', 'ren', 'lcr', 'tau')):
                    fl = {int(r[1]): r[2 + j] == '1' for r in rows if r[0] == ncp}
                    body.append([ncp, self.SENSOR_NAMES[key]] + ['**flagged**' if fl[m] else '·' for m in ms])
            self.table('sensors-cp', ['$n_{cp}$', 'sensor'] + [str(m) for m in ms], body)
        if self.has('sensors-scaling.dat'):
            rows = self.rows('sensors-scaling.dat')
            body = []
            for r in rows:
                if float(r[0]) == 1.0 and float(r[1]) == 0.0:
                    continue
                body.append([f'{float(r[0]):.0e}', f'{float(r[1]):g}', r[2]] + [f'{100 * float(x):.1f} %' for x in r[3:7]])
            self.table('sensors-scaling', ['$a$', '$b$', 'normalisation', 'TENO-based', 'Ren', 'Li, Chen & Ren',
                                           'Wan & Xia'], body)

    def sensor_flows(self) -> None:
        """x-t flag fields of the four sensors on Shu-Osher and the blast (pure TENO5), and the hybrid table."""
        cases = [('shu-osher', 'Shu–Osher'), ('blast', 'Woodward–Colella blast')]
        if self.has(*[f'sensors-{c}-flags.dat' for c, _ in cases]):
            keys = ['teno', 'ren', 'lcr', 'tau']
            fig, axs = plt.subplots(4, 2, figsize=(12, 8.4), sharex='col')
            ink = '0.85' if self.dark else '0.15'
            cmap = mpl.colors.ListedColormap(['none', ink])
            for col, (case, title) in enumerate(cases):
                lines = [ln for ln in (self.data / f'sensors-{case}-flags.dat').read_text().splitlines()
                         if ln and not ln.startswith('#')]
                x = np.array([float(t) for t in lines[0].split()])
                rec = [ln.split() for ln in lines[1:]]
                for row, key in enumerate(keys):
                    sel = [r for r in rec if int(r[1]) == row + 1]
                    t = np.array([float(r[0]) for r in sel])
                    z = np.array([[int(v) for v in r[2:]] for r in sel])
                    ax = axs[row, col]
                    ax.pcolormesh(x, t, z, cmap=cmap, shading='nearest', rasterized=True)
                    ax.set_ylabel(f'{self.SENSOR_NAMES[key]}\n$t$', fontsize=8)
                    if row == 0:
                        ax.set_title(f'{title}, pure TENO5, 400 cells: faces flagged ($s < 1$)', fontsize=10)
                axs[3, col].set_xlabel('$x$')
            fig.tight_layout()
            self.save(fig, 'sensor-flows')
        if self.has('sensors-flows.dat'):
            body = []
            for r in self.rows('sensors-flows.dat'):
                if r[2] == 'pure':
                    continue
                body.append([r[0], label(r[1]), self.SENSOR_NAMES[r[2]], f'{float(r[3]):.4g}', f'{float(r[4]):.3f}',
                             f'{100 * float(r[5]):.1f} %', 'yes' if r[6] == 'T' else '**no**'])
            self.table('sensor-flows', ['case', 'scheme', 'sensor', '$L_1(\\rho)$', 'against pure', 'faces flagged',
                                        'positive'], body)
        if self.has('sensor-fields.dat'):
            rows = self.rows('sensor-fields.dat')
            keys = ['teno', 'ren', 'lcr', 'tau', 'ducros', 'ducros-c']
            par = {'taylor-green': '$N$ = {:g}', 'shock': 'Mach {:g}', 'shock-vortex': '$\\omega_0/\\lvert\\theta\\rvert$ = {:g}',
                   'expansion': '', 'contact': ''}
            groups: dict[tuple[str, str], dict[str, float]] = {}
            for r in rows:
                groups.setdefault((r[0], r[1]), {})[r[2]] = float(r[3])
            body = []
            for (case, p), v in groups.items():
                body.append([case, par[case].format(float(p))] + [f'{v[k]:.3f}' for k in keys])
            self.table('sensor-fields', ['field', 'parameter'] + [self.SENSOR_NAMES[k] for k in keys], body,
                       note='Taylor–Green, expansion, contact: the fraction of all interior cells flagged; shock: of '
                            'the cells next to it; shock in a vortex: of the shock cells within $r_0$ of its centre.')

    # ---------------------------------------------------------------------------------------- unit-test reports
    def spectral(self) -> None:
        names = [p.name[9:-4] for p in sorted(self.data.glob('spectral-*.dat'))]
        if not names:
            print('  skipped: no spectral-*.dat')
            return
        fig, axs = plt.subplots(1, 2, figsize=(12, 4))
        th = np.linspace(0, np.pi, 200)
        axs[0].plot(th / np.pi, th, 'k-', lw=0.8, label='exact')
        axs[1].axhline(0, color='k', lw=0.8)
        for s in names:
            d = self.load(f'spectral-{s}.dat')
            c, m = style(s)
            axs[0].plot(d[:, 0] / np.pi, d[:, 1], '-', color=c, lw=0.9, label=label(s))
            axs[1].plot(d[:, 0] / np.pi, d[:, 2], '-', color=c, lw=0.9)
        axs[0].set_title('linear backgrounds: dispersion', fontsize=10)
        axs[1].set_title('linear backgrounds: dissipation', fontsize=10)
        for ax in axs:
            ax.set_xlabel('$\\xi/\\pi$')
            ax.grid(alpha=0.3)
        axs[0].legend(fontsize=6, ncol=2)
        fig.tight_layout()
        self.save(fig, 'spectral')

    # ------------------------------------------------------------------------------------------ accuracy tables
    def order_tables(self) -> None:
        """A1: the observed orders of the three pointwise order tests."""
        body = []
        if (self.data / 'orders-js.dat').exists():
            o = defaultdict(dict)
            for S, interp, side, order, design in self.rows('orders-js.dat'):
                o[(int(S), interp == 'T')][side] = float(order)
            for (S, it), d in sorted(o.items()):
                body.append([f'WENO{2 * S - 1}-JS' + (' interp.' if it else ''), f'{d["+"]:.2f}', f'{d["-"]:.2f}',
                             str(2 * S - 1)])
        if (self.data / 'orders-wenoz.dat').exists():
            o = defaultdict(dict)
            for kind, order_, p, side, val in self.rows('orders-wenoz.dat'):
                o[(kind, int(order_), int(p))][side] = float(val)
            for (kind, od, p), d in sorted(o.items(), key=lambda t: (t[0][0] != 'recon', t[0][1], t[0][2])):
                body.append([f'WENO-Z{od}' + (' interp.' if kind == 'interp' else '') + f', $p = {p}$',
                             f'{d["+"]:.2f}', f'{d["-"]:.2f}', str(od)])
        if (self.data / 'orders-teno.dat').exists():
            o = defaultdict(dict)
            des = {}
            for name, side, val, design in self.rows('orders-teno.dat'):
                o[name][side] = float(val)
                des[name] = design
            for name, d in o.items():
                body.append([name.replace('-interp', ' interp.'), f'{d["+"]:.2f}', f'{d["-"]:.2f}',
                             des[name].rstrip('0').rstrip('.')])
        if body:
            self.table('order-of-accuracy', ['kernel', 'order, left state', 'order, right state', 'design'], body,
                       'Observed on the finest pair of $h$ above round-off; asserted $\\ge$ design $- 0.5$.')
        else:
            print('  skipped: no orders-*.dat')

    def critical_tables(self) -> None:
        """A2: orders at critical points, reconstructions and interpolations."""
        if (self.data / 'critical-points.dat').exists():
            body = [[label(r[0]) + (f', $p = {r[1]}$' if base(r[0]).startswith('wenoz') else ''),
                     f'{float(r[2]):.2f}', f'{float(r[3]):.2f}', f'{float(r[4]):.2f}', f'{float(r[5]):.2f}']
                    for r in self.rows('critical-points.dat')]
            self.table('critical-points', ['scheme', 'first order, $a > 0$', 'first order, $a < 0$',
                                           'second order, $a > 0$', 'second order, $a < 0$'], body)
        if (self.data / 'interpolation-orders.dat').exists():
            body = [[label(r[0]) + (f', $p = {r[1]}$' if base(r[0]).startswith('wenoz') else '')] +
                    [f'{float(x):.2f}' for x in r[2:8]] for r in self.rows('interpolation-orders.dat')]
            self.table('interpolation-orders', ['interpolation', 'domain, left', 'domain, right', 'first-order cp, left',
                                                'first-order cp, right', 'second-order cp, left',
                                                'second-order cp, right'], body)
        if (self.data / 'convergence-orders.dat').exists():
            body = [[label(r[0]) + (f', $p = {r[1]}$' if base(r[0]).startswith('wenoz') else ''), f'{float(r[2]):.2f}',
                     f'{float(r[3]):.2f}' if float(r[3]) > 0 else '—'] for r in self.rows('convergence-orders.dat')]
            self.table('convergence-orders', ['scheme', 'operator order', 'RK3 advection order'], body)

    def weno_c(self) -> None:
        """A5: the centred WENO's orders, its step overshoot against its upwind part, its equivalence with ADAM."""
        if (self.data / 'weno-c-orders.dat').exists():
            body = [[label(r[0]), 'left' if r[1] == '+' else 'right', f'{float(r[2]):.2f}', fmt(float(r[3]))]
                    for r in self.rows('weno-c-orders.dat')]
            self.table('weno-c-orders', ['scheme', 'state', 'observed order', 'error, finest grid'], body)
        if (self.data / 'weno-c-step.dat').exists():
            body = [[label(r[0]), fmt(float(r[1])), fmt(float(r[2]))] for r in self.rows('weno-c-step.dat')]
            self.table('weno-c-step', ['scheme', 'worst overshoot, unit step', 'its upwind WENO-JS'], body)
        if (self.data / 'weno-c-adam.dat').exists():
            body = [[label(f'weno{2 * int(r[0])}_c'), fmt(float(r[1])), fmt(float(r[2])), f'{int(r[3])} of {int(r[4])}',
                     fmt(float(r[5]))] for r in self.rows('weno-c-adam.dat')]
            self.table('weno-c-adam', ['scheme', 'against ADAM with D1 fixed, both states',
                                       'against ADAM as it is, right state', 'left state different',
                                       'largest left difference'], body,
                       'Differences relative to the stencil\'s magnitude; random stencils at magnitudes 1e-6..1e6.')

    def convergence(self, name: str, file: str, title: str, ylabel: str) -> None:
        if not self.has(file):
            return
        d = defaultdict(dict)
        for r in self.rows(file):
            s, p, n, e = r[0], r[1], int(r[2]), float(r[3])
            key = s + (f' p={p}' if base(s).startswith('wenoz') else '')
            d[key][n] = e
        fig, ax = plt.subplots(figsize=(7, 4.4))
        for k, e in d.items():
            s = k.split()[0]
            c, m = style(s)
            ns = sorted(e)
            ax.loglog(ns, [max(e[n], 1e-17) for n in ns], m + (':' if k.endswith('p=1') else '-'), color=c, ms=3,
                      lw=0.8, label=label(s) + (' ' + k.split()[1] if ' ' in k else ''))
        ax.set_title(title, fontsize=10)
        ax.set_xlabel('$N$')
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3, which='both')
        ax.legend(fontsize=6, ncol=2)
        fig.tight_layout()
        self.save(fig, name)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--data', type=Path, required=True)
    ap.add_argument('--out', type=Path, default=Path('docs/public/verification'))
    ap.add_argument('--tables', type=Path, default=Path('docs/verification/_tables'))
    ap.add_argument('--theme', choices=THEMES, default='both')
    ap.add_argument('--only', nargs='*')
    args = ap.parse_args()
    f = Figures(args.data, args.out, args.tables)
    jobs = {
        'sod': lambda: f.riemann('sod', 'Sod shock tube, $N = 200$, $t = 0.2$', ['weno5_js', 'wenoz5', 'teno6',
                                 'teno6_interp'], (0.62, 0.78)),
        'lax': lambda: f.riemann('lax', 'Lax shock tube, $N = 200$, $t = 0.13$', ['weno5_js', 'wenoz5', 'teno6',
                                 'teno6_interp'], (0.6, 0.85)),
        '123': lambda: f.riemann('123', '123 problem, $N = 200$, $t = 0.15$ (positivity limiter on)',
                                 ['weno5_js', 'teno6', 'teno6_interp'], (0.35, 0.65)),
        'riemann-convergence': f.riemann_convergence,
        'shu-osher': lambda: (f.reference_case('shu-osher', 'Shu–Osher, $N = 400$, $t = 1.8$',
                                               ['weno5_js', 'wenoz5', 'teno6', 'teno6_interp'], (0.4, 2.6)),
                              f.reference_table('shu-osher')),
        'blast': lambda: (f.reference_case('blast', 'Woodward–Colella blast, $N = 400$, $t = 0.038$',
                                           ['weno5_js', 'wenoz5', 'teno6', 'teno6_interp'], (0.6, 0.86)),
                          f.reference_table('blast')),
        'positivity': f.positivity,
        'burgers': f.burgers,
        'euler-order': f.euler_order,
        'advection': f.advection,
        'adr': f.adr,
        'sensors': lambda: (f.sensors(), f.sensor_static()),
        'sensor-flows': f.sensor_flows,
        'hybrid': f.hybrid,
        'spectral': f.spectral,
        'order-tables': f.order_tables,
        'critical-tables': f.critical_tables,
        'weno-c': f.weno_c,
        'nonpoly-flags': f.nonpoly_flags,
        'nonpoly-solver': f.nonpoly_solver,
        'nonpoly-census': f.nonpoly_census,
        'nonpoly-st': f.nonpoly_st,
        'nonpoly-adr': f.nonpoly_adr,
        'nonpoly-advection': f.nonpoly_advection,
        'nonpoly-static': f.nonpoly_static,
        'nonpoly-ablation': f.nonpoly_ablation,
        'nonpoly-e6': f.nonpoly_e6,
        'nonpoly-2d': f.nonpoly_2d,
        'convergence': lambda: f.convergence('convergence', 'convergence.dat',
                                             'whole-domain operator error, Henrick field (critical points included)',
                                             '$L_\\infty(D - u\')$'),
        'interpolation': lambda: f.convergence('interpolation', 'interpolation.dat',
                                               'whole-domain interpolation error, Henrick field', '$L_\\infty$'),
    }
    for th in themes(args.theme):
        f.dark = th == 'dark'
        STYLE.clear()
        STYLE.update(SCHEME_STYLE_DARK if f.dark else SCHEME_STYLE)
        with theme(th):
            for name, job in jobs.items():
                if args.only and name not in args.only:
                    continue
                print(f'{name} ({th})')
                job()
    return 0


if __name__ == '__main__':
    sys.exit(main())
