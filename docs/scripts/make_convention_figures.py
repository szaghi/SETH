#!/usr/bin/env python3
"""Draw the figures and the animation of the Conventions page of SETH's guide (docs/guide/conventions.md).

Why this script exists: the stencil convention (cells by offset from the face, the two biased arrays, the mirror, the
mapping to ADAM's rows) is the one thing a caller must get right, and a picture says it faster than a table. The
figures are drawn from SETH's own data: the stencil map reads the bounds and the candidate stencils of every scheme
from src/lib/seth_z_teno_tables.F90 (and S for WENO-JS), and the animation runs WENO5-JS with SETH's coefficients
(src/lib/seth_weno_js_tables.F90: the Jiang-Shu weights with ADAM's exponent S = 3 and eps = 1e-6). Its states equal
those of seth_weno5_js to round-off (checked when the figures were made: worst difference 3.3e-16 over the sweep).

Writes, into --out (default docs/public/guide):

    conv-face.png         one face, the two arrays vp and vm of WENO5, the cells each reads, its candidates
    conv-schemes.png      every scheme: the bounds of vp and vm and the cells each reads; TENO's candidates
    conv-mirror.png       vm(m) = vp(1 - m) gives rm = rp
    conv-adam-rows.png    SETH's arrays and ADAM's rows v(1:2, 1-S:S-1)
    conv-nonpoly.png      the -5:6 footprint of the composite non-polynomial kernels: the masks each state reads
    casts/conv-sweep-{light,dark}.{webm,mp4,png}   WENO5 sweeping a field with a jump: the candidates, the weights
                                                    and the two states at every face

each static figure with its silver-ink twin <name>-dark.png (docs/scripts/figure_theme.py, the ADAM/FLUME harness).

Usage: make_convention_figures.py [--out DIR] [--theme {light,dark,both}] [--no-cast] [--fps N] [--width PX]
Needs numpy, matplotlib and (for the animation) ffmpeg on PATH.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib as mpl
import numpy as np

mpl.use('Agg')
mpl.rcParams['figure.facecolor'] = '#f3ede1'
mpl.rcParams['savefig.facecolor'] = '#f3ede1'
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, Rectangle  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_theme import THEMES, theme, themes  # noqa: E402
from matplotlib.colors import to_rgba  # noqa: E402
import seth_palette as P  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
DPI = 130
# neutral greys ('k', '0.42'), which figure_theme maps onto the silver-ink ladder in the dark twin; translucent fills,
# which read on either ground; the two data colours from seth_palette
INK, GREY = 'k', '0.42'
LEFT, RIGHT = P.LEFT, P.RIGHT  # vp / left state blue, vm / right state moegi (seth_palette)
LEFT_SOFT, RIGHT_SOFT = to_rgba(LEFT, 0.2), to_rgba(RIGHT, 0.2)

# ------------------------------------------------------------------------------------------------- SETH's data
# WENO5-JS (src/lib/seth_weno_js_tables.F90, row 2 = the left-biased reconstruction): candidate k reads cells k-2..k
P5 = np.array([[1 / 3, -7 / 6, 11 / 6], [-1 / 6, 5 / 6, 1 / 3], [1 / 3, 5 / 6, -1 / 6]])
D5 = np.array([0.1, 0.6, 0.3])
EPS5, WEXP5 = 1.0e-6, 3


def weno5_left(u: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    """WENO5-JS left state from u = values at the cells -2..2 of a face: (state, candidate values, weights)."""
    q = np.array([P5[k] @ u[k:k + 3] for k in range(3)])
    b = np.array([
        13 / 12 * (u[0] - 2 * u[1] + u[2]) ** 2 + 1 / 4 * (u[0] - 4 * u[1] + 3 * u[2]) ** 2,
        13 / 12 * (u[1] - 2 * u[2] + u[3]) ** 2 + 1 / 4 * (u[1] - u[3]) ** 2,
        13 / 12 * (u[2] - 2 * u[3] + u[4]) ** 2 + 1 / 4 * (3 * u[2] - 4 * u[3] + u[4]) ** 2,
    ])
    a = D5 / (EPS5 + b) ** WEXP5
    w = a / a.sum()
    return float(w @ q), q, w


def parse_tables() -> dict[str, dict[str, object]]:
    """Bounds and candidate stencils of the WENO-Z and TENO kernels, from src/lib/seth_z_teno_tables.F90."""
    src = (ROOT / 'src/lib/seth_z_teno_tables.F90').read_text()
    out: dict[str, dict[str, object]] = {}
    for name, key, val in re.findall(r'parameter :: (\w+?)_(LO|UMAX|NC) = (-?\d+)_I4P', src):
        out.setdefault(name, {})[key] = int(val)
    for name, key, val in re.findall(r'parameter :: (\w+?)_(CLO|CW)\(0:\d+\) = \[([^\]]*)\]', src):
        out.setdefault(name, {})[key] = [int(v.replace('_I4P', '')) for v in val.split(',')]
    return out


# ----------------------------------------------------------------------------------------------- drawing aids
def cell(ax, m: int, y: float, h: float, fill: str | None, edge: str = GREY, hatch: str | None = None,
         lw: float = 0.8) -> None:
    ax.add_patch(Rectangle((m - 0.5, y), 1.0, h, facecolor=fill or 'none', edgecolor=edge, lw=lw, hatch=hatch))


def bracket(ax, lo: int, hi: int, y: float, color: str, text: str, below: bool = True, fs: float = 8) -> None:
    """A bracket under (or over) the cells lo..hi, with a label."""
    s = -1 if below else 1
    ax.plot([lo - 0.42, lo - 0.42, hi + 0.42, hi + 0.42], [y - s * 0.08, y, y, y - s * 0.08], color=color, lw=1.1)
    ax.text((lo + hi) / 2, y + s * 0.07, text, ha='center', va='top' if below else 'bottom', fontsize=fs, color=color)


def face_line(ax, y0: float, y1: float, label: str | None = 'face $i+1/2$') -> None:
    ax.plot([0.5, 0.5], [y0, y1], '--', color=INK, lw=1.0)
    if label:
        ax.text(0.5, y1 + 0.05, label, ha='center', va='bottom', fontsize=9, color=INK)


def arrow(ax, x0: float, y0: float, x1: float, y1: float, color: str, rad: float = 0.0, lw: float = 1.3) -> None:
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='-|>', mutation_scale=11, color=color, lw=lw,
                                 connectionstyle=f'arc3,rad={rad}'))


def bare(ax) -> None:
    ax.set_axis_off()
    ax.set_aspect('equal')


# ------------------------------------------------------------------------------------------------------ figures
def fig_face(out: Path) -> None:
    """One face of WENO5: the mesh, vp and vm with the cells they read, the candidates and their linear weights."""
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    bare(ax)
    ym, yp, yv = 4.4, 3.0, 0.9  # rows: mesh, vp, vm
    for m in range(-2, 4):
        cell(ax, m, ym, 0.45, None, edge=INK)
        ax.text(m, ym + 0.22, f'{m}', ha='center', va='center', fontsize=10, color=INK)
    ax.text(-2.75, ym + 0.22, 'offset $m$', ha='right', va='center', fontsize=9, color=GREY)
    face_line(ax, -0.05, ym + 0.6)
    for y, arr, col, soft, read, lab in ((yp, 'vp', LEFT, LEFT_SOFT, range(-2, 3), 'vp(-2:3)\nleft-biased'),
                                         (yv, 'vm', RIGHT, RIGHT_SOFT, range(-1, 4), 'vm(-2:3)\nright-biased')):
        for m in range(-2, 4):
            r = m in read
            cell(ax, m, y, 0.45, soft if r else None, edge=col if r else GREY)
            ax.text(m, y + 0.22, f'{arr}({m})' if r else 'not read', ha='center', va='center', fontsize=7.5,
                    color=INK if r else GREY, style='normal' if r else 'italic')
        ax.text(-2.75, y + 0.22, lab, ha='right', va='center', fontsize=9, color=col)
    for k, (lo, d) in enumerate(zip((-2, -1, 0), D5)):
        yy = yp - 0.15 - 0.3 * k
        bracket(ax, lo, lo + 2, yy, LEFT, '')
        ax.text(3.65, yy, f'candidate {k}: cells {lo}..{lo + 2}, $d_{k}$ = {d:g}', ha='left', va='center',
                fontsize=8, color=LEFT)
    for k, (lo, d) in enumerate(zip((1, 0, -1), D5)):
        yy = yv - 0.15 - 0.3 * k
        bracket(ax, lo, lo + 2, yy, RIGHT, '')
        ax.text(3.65, yy, f'mirror of candidate {k}: cells {lo}..{lo + 2}, $d_{k}$ = {d:g}', ha='left',
                va='center', fontsize=8, color=RIGHT)
    arrow(ax, -0.6, yp + 0.75, 0.46, yp + 0.75, LEFT)
    ax.text(-0.65, yp + 0.75, 'rp: the left state', ha='right', va='center', fontsize=9, color=LEFT)
    arrow(ax, 1.6, yv + 0.75, 0.54, yv + 0.75, RIGHT)
    ax.text(1.65, yv + 0.75, 'rm: the right state', ha='left', va='center', fontsize=9, color=RIGHT)
    ax.set_xlim(-4.6, 8.4)
    ax.set_ylim(-0.1, ym + 1.0)
    fig.savefig(out / 'conv-face.png', dpi=DPI, bbox_inches='tight', pad_inches=0.08)
    plt.close(fig)


def fig_schemes(out: Path) -> None:
    """Every scheme: the array bounds, the cells vp and vm read, and the candidates of TENO."""
    t = parse_tables()
    rows: list[tuple[str, int, int, tuple[int, int], tuple[int, int], list[tuple[int, int]] | None]] = []
    for s in range(1, 6):
        rows.append((f'WENO{2 * s - 1}-JS', 1 - s, s, (1 - s, s - 1), (2 - s, s), None))
    for name, lab in (('WENOZ3', 'WENO-Z3'), ('WENOZ5', 'WENO-Z5'), ('WENOZ7', 'WENO-Z7'), ('WENOZ9', 'WENO-Z9'),
                      ('TENO5', 'TENO5'), ('TENO6', 'TENO6')):
        d = t[name]
        lo, umax = int(d['LO']), int(d['UMAX'])
        cand = [(c, c + w - 1) for c, w in zip(d['CLO'], d['CW'])] if name.startswith('TENO') else None
        rows.append((lab, lo, 1 - lo, (lo, umax), (1 - umax, 1 - lo), cand))
    n = len(rows)
    heights = [1.4 if r[5] else 0.95 for r in rows]  # TENO rows carry their candidates above the cells
    fig, ax = plt.subplots(figsize=(8.6, 0.65 * sum(heights) + 1.2))
    bare(ax)
    top = sum(heights)
    for m in range(-4, 6):
        ax.text(m, top + 0.35, f'{m}', ha='center', va='center', fontsize=9, color=INK)
    face_line(ax, -0.3, top + 0.55, 'face')
    for r, (lab, lo, hi, (pl, ph), (ml, mh), cand) in enumerate(rows):
        y = top - sum(heights[:r + 1]) + 0.2
        ax.text(-4.9, y + 0.33, lab, ha='right', va='center', fontsize=9, color=INK)
        ax.text(6.0, y + 0.33, f'({lo}:{hi})', ha='left', va='center', fontsize=8, color=GREY, family='monospace')
        for m in range(lo, hi + 1):
            cell(ax, m, y + 0.33, 0.3, LEFT_SOFT if pl <= m <= ph else None, edge=LEFT if pl <= m <= ph else GREY)
            cell(ax, m, y, 0.3, RIGHT_SOFT if ml <= m <= mh else None, edge=RIGHT if ml <= m <= mh else GREY)
        if cand:
            for k, (c0, c1) in enumerate(cand):
                yy = y + 0.74 + 0.11 * k
                ax.plot([c0 - 0.4, c1 + 0.4], [yy, yy], '-', color=LEFT, lw=1.6, solid_capstyle='butt')
                ax.text(c0 - 0.5, yy, f'{k}', ha='right', va='center', fontsize=6.5, color=LEFT)
    ax.text(-4.9, top + 0.35, 'offset', ha='right', va='center', fontsize=9, color=GREY)
    ax.add_patch(Rectangle((1.6, -0.75), 0.5, 0.25, facecolor=LEFT_SOFT, edgecolor=LEFT, lw=0.8))
    ax.text(2.2, -0.62, 'read by vp (left state)', ha='left', va='center', fontsize=8, color=INK)
    ax.add_patch(Rectangle((1.6, -1.1), 0.5, 0.25, facecolor=RIGHT_SOFT, edgecolor=RIGHT, lw=0.8))
    ax.text(2.2, -0.97, 'read by vm (right state)', ha='left', va='center', fontsize=8, color=INK)
    ax.plot([-4.3, -3.6], [-0.62, -0.62], '-', color=LEFT, lw=1.0)
    ax.text(-3.5, -0.62, 'TENO candidates of vp', ha='left', va='center', fontsize=8, color=INK)
    ax.add_patch(Rectangle((-4.3, -1.1), 0.5, 0.25, facecolor='none', edgecolor=GREY, lw=0.8))
    ax.text(-3.7, -0.97, 'in the bounds, not read', ha='left', va='center', fontsize=8, color=INK)
    ax.set_xlim(-7.4, 7.6)
    ax.set_ylim(-1.3, top + 0.8)
    fig.savefig(out / 'conv-schemes.png', dpi=DPI, bbox_inches='tight', pad_inches=0.08)
    plt.close(fig)


def fig_mirror(out: Path) -> None:
    """vm(m) = vp(1 - m): the right state of the mirrored data equals the left state of the original."""
    m = np.arange(-2, 4)
    vp = np.array([1.0, 1.0, 0.95, 0.25, 0.2, 0.2])
    vm = vp[::-1]  # vm(m) = vp(1 - m): m = -2..3 reads vp(3..-2)
    rp, _, _ = weno5_left(vp[:5])
    rm, _, _ = weno5_left(vm[::-1][:5])  # the right state is the left state of the mirrored right-biased values
    fig, axs = plt.subplots(1, 2, figsize=(8.6, 3.0), sharey=True)
    for ax, v, col, soft, name, r, side in ((axs[0], vp, LEFT, LEFT_SOFT, 'vp', rp, 'rp'),
                                            (axs[1], vm, RIGHT, RIGHT_SOFT, 'vm', rm, 'rm')):
        ax.bar(m, v, width=0.92, color=soft, edgecolor=col, lw=0.8)
        ax.axvline(0.5, ls='--', color=INK, lw=1.0)
        ax.plot([0.5], [r], 'o', color=col, ms=7, zorder=5)
        ax.annotate(f'{side} = {r:.6f}', (0.5, r), xytext=(18 if side == 'rp' else -18, 22), textcoords='offset points',
                    ha='left' if side == 'rp' else 'right', fontsize=9, color=col,
                    arrowprops={'arrowstyle': '-', 'color': col, 'lw': 0.8})
        ax.set_xticks(m)
        ax.set_xlabel('offset $m$')
        ax.set_title(f'{name}(m)', fontsize=10, color=col)
        ax.set_ylim(0, 1.25)
        ax.grid(axis='y', alpha=0.25)
    axs[0].set_ylabel('value')
    fig.text(0.5, -0.04, 'vm(m) = vp(1 − m): the same data seen from the other side; WENO5-JS gives rm = rp',
             ha='center', fontsize=9, color=INK)
    fig.tight_layout()
    fig.savefig(out / 'conv-mirror.png', dpi=DPI, bbox_inches='tight', pad_inches=0.08)
    plt.close(fig)


def fig_nonpoly(out: Path) -> None:
    """The footprint of the composite non-polynomial kernels (issue #6): the masks a face reads, both states."""
    fig, ax = plt.subplots(figsize=(9.6, 7.4))
    bare(ax)
    ym, yp, yv = 7.6, 6.6, 2.7  # rows: mesh, vp, vm
    for m in range(-5, 7):
        cell(ax, m, ym, 0.45, None, edge=INK)
        ax.text(m, ym + 0.22, f'{m}', ha='center', va='center', fontsize=9, color=INK)
    ax.text(-5.75, ym + 0.22, 'offset $m$', ha='right', va='center', fontsize=9, color=GREY)
    face_line(ax, 0.0, ym + 0.6)
    for y, arr, col, soft, dr, th, lab in ((yp, 'vp', LEFT, LEFT_SOFT, range(-5, 5), [5], 'vp(-5:6)\nleft state'),
                                           (yv, 'vm', RIGHT, RIGHT_SOFT, range(-3, 7), [-4], 'vm(-5:6)\nright state')):
        for m in range(-5, 7):
            if m in dr:
                cell(ax, m, y, 0.45, soft, edge=col)
            elif m in th:
                cell(ax, m, y, 0.45, None, edge=col, hatch='///')
            else:
                cell(ax, m, y, 0.45, None, edge=GREY)
        ax.text(-5.75, y + 0.22, lab, ha='right', va='center', fontsize=9, color=col)
    # the masks the left state reads: cells -3..2 (DR), -3..3 (THINC), each from its five values c-2..c+2
    for j, c in enumerate(range(-3, 4)):
        yy = yp - 0.25 - 0.42 * j
        bracket(ax, c - 2, c + 2, yy, LEFT, '')
        ax.plot([c], [yy + 0.13], 'o', color=LEFT, ms=4)
        ax.text(6.8, yy, f'mask of cell {c}' + (' (THINC only)' if c == 3 else ''), ha='left', va='center',
                fontsize=8, color=LEFT)
    # the right state: the same construction on w(m) = vm(1 - m), so the mirrored cells 1 - c
    for j, c in enumerate(range(-3, 4)):
        cm = 1 - c
        yy = yv - 0.25 - 0.42 * j
        bracket(ax, cm - 2, cm + 2, yy, RIGHT, '')
        ax.plot([cm], [yy + 0.13], 'o', color=RIGHT, ms=4)
        ax.text(6.8, yy, f'mask of cell {cm} (mirrored {c})' + (' (THINC only)' if c == 3 else ''), ha='left',
                va='center', fontsize=8, color=RIGHT)
    ax.text(0.5, -0.35, 'shaded: read by TENO5-DR; hatched: read by TENO5-THINC only. Both states together: twelve cells, '
            '-5..6, so 6 ghost cells', ha='center', va='top', fontsize=8.5, color=INK)
    ax.set_xlim(-7.4, 11.2)
    ax.set_ylim(-0.7, ym + 1.0)
    fig.savefig(out / 'conv-nonpoly.png', dpi=DPI, bbox_inches='tight', pad_inches=0.08)
    plt.close(fig)


def fig_adam_rows(out: Path) -> None:
    """SETH's vp(1-S:S), vm(1-S:S) against ADAM's v(1:2, 1-S:S-1), S = 3."""
    fig, ax = plt.subplots(figsize=(8.6, 4.2))
    bare(ax)
    y_vp, y_r2, y_r1, y_vm = 3.6, 2.3, 1.55, 0.25  # SETH vp, ADAM row 2, ADAM row 1, SETH vm
    for m in range(-2, 4):
        cell(ax, m, y_vp, 0.45, LEFT_SOFT if m <= 2 else None, edge=LEFT if m <= 2 else GREY)
        ax.text(m, y_vp + 0.22, f'vp({m})', ha='center', va='center', fontsize=7.5, color=INK if m <= 2 else GREY)
        cell(ax, m, y_vm, 0.45, RIGHT_SOFT if m >= -1 else None, edge=RIGHT if m >= -1 else GREY)
        ax.text(m, y_vm + 0.22, f'vm({m})', ha='center', va='center', fontsize=7.5, color=INK if m >= -1 else GREY)
    for m in range(-2, 3):
        cell(ax, m, y_r2, 0.45, LEFT_SOFT, edge=LEFT)
        ax.text(m, y_r2 + 0.22, f'v(2,{m})', ha='center', va='center', fontsize=7.5, color=INK)
        cell(ax, m, y_r1, 0.45, RIGHT_SOFT, edge=RIGHT)
        ax.text(m, y_r1 + 0.22, f'v(1,{m})', ha='center', va='center', fontsize=7.5, color=INK)
        arrow(ax, m, y_vp - 0.02, m, y_r2 + 0.5, LEFT, lw=0.9)                 # vp(m) -> v(2,m)
        arrow(ax, m + 1, y_vm + 0.47, m + 0.08, y_r1 - 0.04, RIGHT, lw=0.9)    # vm(m+1) -> v(1,m)
    ax.text(-2.75, y_vp + 0.22, 'SETH  vp(-2:3)', ha='right', va='center', fontsize=9, color=LEFT)
    ax.text(-2.75, y_r2 + 0.22, 'ADAM  row 2', ha='right', va='center', fontsize=9, color=LEFT)
    ax.text(-2.75, y_r1 + 0.22, 'ADAM  row 1', ha='right', va='center', fontsize=9, color=RIGHT)
    ax.text(-2.75, y_vm + 0.22, 'SETH  vm(-2:3)', ha='right', va='center', fontsize=9, color=RIGHT)
    ax.text(2.75, y_r2 + 0.22, 'left state, vr(2)', ha='left', va='center', fontsize=8.5, color=LEFT)
    ax.text(2.75, y_r1 + 0.22, 'right state, vr(1): vm shifted one cell left', ha='left', va='center',
            fontsize=8.5, color=RIGHT)
    face_line(ax, 0.1, 4.2)
    ax.text(0.5, -0.3, r'$v(2,m) = v_p(m)$,   $v(1,m) = v_m(m+1)$,   $m = 1-S, \dots, S-1$', ha='center',
            va='center', fontsize=10, color=INK)
    ax.set_xlim(-5.2, 7.4)
    ax.set_ylim(-0.6, 4.6)
    fig.savefig(out / 'conv-adam-rows.png', dpi=DPI, bbox_inches='tight', pad_inches=0.08)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------- animation
CAST_THEMES = (  # seth_palette: blue for vp and the left state, moegi for vm and the right state
    ('light', {'bg': P.WASHI, 'bg_alt': P.WASHI_ALT, 'text': P.INK_LIGHT, 'mesh': P.MESH_LIGHT, 'left': P.LEFT,
               'right': P.RIGHT, 'left_soft': '#c7d6e8', 'right_soft': '#d3dfbf', 'avg': P.MESH_LIGHT}),
    ('dark', {'bg': P.SILVER, 'bg_alt': P.SILVER_ALT, 'text': P.INK_DARK, 'mesh': P.MESH_DARK, 'left': P.LEFT_D,
              'right': P.RIGHT_D, 'left_soft': '#17304a', 'right_soft': '#22340d', 'avg': '#a9b1bb'}),
)
NCELL = 32


def sweep_field() -> tuple[np.ndarray, np.ndarray]:
    """Cell averages on [0, 1], NCELL cells: a smooth wave with a jump at x = 0.53 (exact averages)."""
    xe = np.linspace(0.0, 1.0, NCELL + 1)
    h = xe[1] - xe[0]

    def prim(x: np.ndarray) -> np.ndarray:  # integral of 0.45 sin(2 pi x) + 0.5 + 0.6 H(x - 0.53)
        return -0.45 / (2 * np.pi) * np.cos(2 * np.pi * x) + 0.5 * x + 0.6 * np.maximum(x - 0.53, 0.0)

    return 0.5 * (xe[1:] + xe[:-1]), (prim(xe[1:]) - prim(xe[:-1])) / h


def quad_through_averages(ubar: np.ndarray, x0: float, h: float) -> tuple[np.ndarray, np.ndarray]:
    """The quadratic whose averages over the three cells [x0 + j h, x0 + (j+1) h] are ubar: its curve on them."""
    xs = np.linspace(x0, x0 + 3 * h, 60)
    # average of 1, x, x^2 over each cell -> 3x3 system
    a = np.array([[1.0, x0 + (j + 0.5) * h, ((x0 + (j + 1) * h) ** 3 - (x0 + j * h) ** 3) / (3 * h)] for j in range(3)])
    c = np.linalg.solve(a, ubar)
    return xs, c[0] + c[1] * xs + c[2] * xs ** 2


def render_sweep(out: Path, fps: int, width: int) -> None:
    xc, ub = sweep_field()
    h = xc[1] - xc[0]
    faces = list(range(2, NCELL - 3))  # face i+1/2 between cells i and i+1 (0-based): stencils inside the domain
    states = {}
    for i in faces:
        rp, qp, wp = weno5_left(ub[i - 2:i + 3])
        rm, qm, wm = weno5_left(ub[i + 3:i - 2 if i >= 2 else None:-1][:5])  # mirrored: cells i+3..i-1
        states[i] = (rp, qp, wp, rm, qm, wm)
    hold = 4  # frames per face
    for name, c in CAST_THEMES:
        tmp = Path(tempfile.mkdtemp(prefix=f'seth-conv-{name}-'))
        try:
            rc = {'figure.facecolor': c['bg'], 'savefig.facecolor': c['bg'], 'axes.facecolor': c['bg_alt'],
                  'axes.edgecolor': c['mesh'], 'axes.labelcolor': c['text'], 'text.color': c['text'],
                  'xtick.color': c['mesh'], 'ytick.color': c['mesh'], 'xtick.labelcolor': c['text'],
                  'ytick.labelcolor': c['text'], 'legend.facecolor': c['bg'], 'legend.edgecolor': c['mesh']}
            with mpl.rc_context(rc):
                f = 0
                for n, i in enumerate(faces):
                    for _ in range(hold if n < len(faces) - 1 else 4 * hold):
                        draw_sweep_frame(c, xc, ub, h, faces[:n + 1], states, i, tmp / f'f{f:04d}.png')
                        f += 1
            out.mkdir(parents=True, exist_ok=True)
            base = out / f'conv-sweep-{name}'
            scale = f'scale={width}:-2:flags=lanczos'
            subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', str(fps), '-i', str(tmp / 'f%04d.png'),
                            '-vf', scale, '-c:v', 'libvpx-vp9', '-crf', '36', '-b:v', '0', '-pix_fmt', 'yuv420p',
                            f'{base}.webm'], check=True)
            subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', str(fps), '-i', str(tmp / 'f%04d.png'),
                            '-vf', scale, '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-movflags', '+faststart',
                            f'{base}.mp4'], check=True)
            poster = faces.index(17) * hold  # the face next to the jump
            subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-i', str(tmp / f'f{poster:04d}.png'), '-vf', scale,
                            f'{base}.png'], check=True)
            print(f'  {base}.{{webm,mp4,png}}')
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


def draw_sweep_frame(c: dict[str, str], xc: np.ndarray, ub: np.ndarray, h: float, done: list[int],
                     states: dict, i: int, path: Path) -> None:
    rp, qp, wp, rm, qm, wm = states[i]
    xf = xc[i] + 0.5 * h
    fig = plt.figure(figsize=(7.2, 4.6), dpi=100)
    ax = fig.add_axes([0.1, 0.36, 0.88, 0.6])
    # the two stencils
    ax.axvspan(xc[i - 2] - 0.5 * h, xc[i + 2] + 0.5 * h, ymin=0.52, ymax=1.0, color=c['left_soft'], lw=0, alpha=0.55)
    ax.axvspan(xc[i - 1] - 0.5 * h, xc[i + 3] + 0.5 * h, ymin=0.0, ymax=0.48, color=c['right_soft'], lw=0, alpha=0.55)
    # cell averages
    for j in range(xc.size):
        ax.plot([xc[j] - 0.5 * h, xc[j] + 0.5 * h], [ub[j], ub[j]], '-', color=c['avg'], lw=1.4,
                solid_capstyle='butt')
    # candidate polynomials of the left state (cells i-2+k .. i+k)
    for k in range(3):
        xs, ys = quad_through_averages(ub[i - 2 + k:i + 1 + k], xc[i - 2 + k] - 0.5 * h, h)
        ax.plot(xs, ys, '-', color=c['left'], lw=0.5 + 1.5 * wp[k], alpha=0.3 + 0.6 * wp[k])
        ax.plot([xf], [qp[k]], 'o', mfc='none', mec=c['left'], ms=4, mew=0.8)
    for k in range(3):  # right state: candidates on cells i+1-k .. i+3-k
        lo = i + 1 - k
        xs, ys = quad_through_averages(ub[lo:lo + 3], xc[lo] - 0.5 * h, h)
        ax.plot(xs, ys, '-', color=c['right'], lw=0.5 + 1.5 * wm[k], alpha=0.3 + 0.6 * wm[k])
        ax.plot([xf], [qm[k]], 's', mfc='none', mec=c['right'], ms=4, mew=0.8)
    # the states of the faces swept so far, and of this one
    ax.plot([xc[j] + 0.5 * h for j in done], [states[j][0] for j in done], '<', color=c['left'], ms=3.5, lw=0)
    ax.plot([xc[j] + 0.5 * h for j in done], [states[j][3] for j in done], '>', color=c['right'], ms=3.5, lw=0)
    ax.axvline(xf, ls='--', color=c['text'], lw=0.9)
    ax.plot([xf - 0.006], [rp], '<', color=c['left'], ms=9, mec=c['bg'], mew=0.6, zorder=6)
    ax.plot([xf + 0.006], [rm], '>', color=c['right'], ms=9, mec=c['bg'], mew=0.6, zorder=6)
    ax.set_xlim(0, 1)
    ax.set_ylim(0.1, 1.6)
    ax.set_xticks([])
    ax.set_ylabel('averages, states')
    ax.text(0.01, 1.56, 'vp: cells $i-2..i+2$, left state rp ◀', color=c['left'], fontsize=8, va='top')
    ax.text(0.01, 1.47, 'vm: cells $i-1..i+3$, right state rm ▶', color=c['right'], fontsize=8, va='top')
    # weights
    for col, w, x0, lab, cells in ((c['left'], wp, 0.08, 'weights of rp', ('i−2..i', 'i−1..i+1', 'i..i+2')),
                                   (c['right'], wm, 0.55, 'weights of rm', ('i+1..i+3', 'i..i+2', 'i−1..i+1'))):
        bx = fig.add_axes([x0 + 0.02, 0.07, 0.38, 0.2])
        bx.bar(range(3), w, color=col, width=0.6)
        bx.plot(range(3), D5, '_', color=c['text'], ms=22, mew=1.2)
        bx.set_ylim(0, 1.0)
        bx.set_xticks(range(3), [f'{cl}' for cl in cells], fontsize=7)
        bx.set_yticks([0, 0.5, 1.0])
        bx.tick_params(labelsize=7)
        bx.set_title(f'{lab} (bars) and linear weights (ticks)', fontsize=7.5, color=c['text'], pad=3)
    fig.text(0.1, 0.31, f'face $i+1/2$, i = {i + 1} of {xc.size}', fontsize=8.5, color=c['text'])
    fig.savefig(path)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--out', type=Path, default=ROOT / 'docs/public/guide')
    ap.add_argument('--theme', choices=THEMES, default='both')
    ap.add_argument('--no-cast', action='store_true')
    ap.add_argument('--fps', type=int, default=12)
    ap.add_argument('--width', type=int, default=720)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    global LEFT, RIGHT, LEFT_SOFT, RIGHT_SOFT  # noqa: PLW0603 - the drawing functions read the module colours
    for th in themes(args.theme):
        dark = th == 'dark'
        LEFT, RIGHT = (P.LEFT_D, P.RIGHT_D) if dark else (P.LEFT, P.RIGHT)
        LEFT_SOFT, RIGHT_SOFT = to_rgba(LEFT, 0.22), to_rgba(RIGHT, 0.22)
        with theme(th):
            fig_face(args.out)
            fig_schemes(args.out)
            fig_mirror(args.out)
            fig_adam_rows(args.out)
            fig_nonpoly(args.out)
        print(f'figures ({th}): {args.out}/conv-{{face,schemes,mirror,adam-rows,nonpoly}}.png')
    if not args.no_cast:
        if shutil.which('ffmpeg') is None:
            print('ffmpeg not found: no animation rendered')
            return 1
        print('conv-sweep')
        render_sweep(args.out / 'casts', args.fps, args.width)
    return 0


if __name__ == '__main__':
    sys.exit(main())
