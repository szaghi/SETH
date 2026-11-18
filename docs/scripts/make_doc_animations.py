#!/usr/bin/env python3
"""Render the animations of SETH's verification documentation from the frames of the verification programs.

Modelled on ADAM's docs/scripts/make_gallery.py (the landing-page animations of the ADAM and FLUME docs): each
animation is drawn twice, on the two grounds of the docs identity, washi (light) and silver ink (dark), with no
title in the frame, a small time label under the panel, and no margin of another colour. SETH's cases are 1-D, so a
frame is a set of curves (the exact or reference solution and two or more schemes) on fixed axes. Per animation and
theme it writes, into --out:

    <name>-<theme>.webm   VP9, constant quality (crf 36), 720 px wide
    <name>-<theme>.mp4    H.264, yuv420p, faststart
    <name>-<theme>.png    the poster: the last frame

Animations and the frames they read (written by scripts/run_vv.sh into --data):

    shu-osher   shu-osher-frames-{weno5_js,teno6}.dat     density, the shock running into the entropy waves
    blast       blast-frames-{weno5_js,teno6}.dat         density, the two blast waves and their collision
    burgers     burgers-frames-{weno5_js,teno6}.dat       u and the exact solution while the shock forms
    advection   advection-frames-{weno5_js,wenoz5,teno5,teno6}.dat   the composite profile over four periods

The non-polynomial reconstructions of issue #6 (docs/verification/nonpolynomial/):

    nonpoly-adr-morph   no data: the linear row P + s^2 C of the Gaussian-RBF TENO5 (the five-point row of
                        src/lib/seth_nonpoly_tables.F90), its modified wavenumber as s = lambda dx grows from 0 to 0.3
    nonpoly-square      nonpoly-square-frames-{teno5,teno5_dr,teno5_thinc}.dat: the square wave over 10 periods, a
                        strip per composite under it marking the faces where its state differs from TENO5's
    nonpoly-sod         nonpoly-sod-frames-{teno5,teno5_dr,teno5_thinc}.dat: Sod's density, the contact zoomed
    nonpoly-lambda      nonpoly-shu-osher-frames-teno5_rbf.dat: Shu-Osher's density under TENO5-RBF and, under it,
                        lambda = s/dx of the scheme's sensor (Liang & Fu 2024, Eqs. 30-32) evaluated on the density
    nonpoly-blast       nonpoly-blast-{frames,flags}-teno5_dr.dat: the blast waves under TENO5-DR, the DR detector's
                        level on the density under them

Usage: make_doc_animations.py --data DIR [--out DIR] [--only NAME ...] [--fps N] [--width PX]
Needs numpy, matplotlib and ffmpeg on PATH.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import matplotlib as mpl
import numpy as np

mpl.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import seth_palette as P  # noqa: E402


@dataclass(frozen=True)
class Theme:
    """The colours of one ground of the docs identity, and the line colours (seth_palette)."""

    name: str
    bg: str
    bg_alt: str
    text: str
    mesh: str
    lines: tuple[str, ...]  # exact/reference first, then the schemes


THEMES = (  # the Neon Moegi identity (seth_palette): exact/reference in ink, then blue, moegi, plum, umber
    Theme('light', P.WASHI, P.WASHI_ALT, P.INK_LIGHT, P.MESH_LIGHT, (P.INK_LIGHT, *P.LINES_LIGHT)),
    Theme('dark', P.SILVER, P.SILVER_ALT, P.INK_DARK, P.MESH_DARK, (P.INK_DARK, *P.LINES_DARK)),
)
LABEL = {'weno5_js': 'WENO5-JS', 'wenoz5': 'WENO-Z5', 'teno5': 'TENO5', 'teno6': 'TENO6', 'teno5_dr': 'TENO5-DR',
         'teno5_thinc': 'TENO5-THINC'}
# the colour follows the scheme, as on the figures (index into Theme.lines after the exact solution)
SLOT = {'WENO5-JS': 0, 'WENO-Z5': 1, 'TENO5': 2, 'TENO6': 3, 'TENO5-DR': 0, 'TENO5-THINC': 1}  # the composites: WENO's slots


def read_frames(path: Path, euler: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]:
    """x, t, the field per frame (density for Euler, u for scalar), the exact field per frame (scalar) or None."""
    rows = [ln for ln in path.read_text().splitlines() if ln and not ln.startswith('#')]
    x = np.array(rows[0].split(), float)
    n = x.size
    data = np.array([r.split() for r in rows[1:]], float)
    t = data[:, 0]
    field = data[:, 1:1 + n]
    exact = None if euler else data[:, 1 + n:1 + 2 * n]
    return x, t, field, exact


def render(name: str, x: np.ndarray, t: np.ndarray, curves: list[tuple[str, np.ndarray]], exact: np.ndarray | None,
           ylabel: str, xlim: tuple[float, float], ylim: tuple[float, float], out: Path, fps: int, width: int,
           tlabel: str) -> None:
    """One animation in both themes: frames to PNG, then ffmpeg to webm and mp4, and the poster."""
    for th in THEMES:
        tmp = Path(tempfile.mkdtemp(prefix=f'seth-{name}-'))
        try:
            with mpl.rc_context({'figure.facecolor': th.bg, 'savefig.facecolor': th.bg, 'axes.facecolor': th.bg_alt,
                                 'axes.edgecolor': th.mesh, 'axes.labelcolor': th.text, 'text.color': th.text,
                                 'xtick.color': th.mesh, 'ytick.color': th.mesh, 'xtick.labelcolor': th.text,
                                 'ytick.labelcolor': th.text, 'grid.color': th.mesh, 'legend.facecolor': th.bg,
                                 'legend.edgecolor': th.mesh}):
                for f in range(t.size):
                    fig = plt.figure(figsize=(7.2, 3.2), dpi=100)
                    ax = fig.add_axes([0.08, 0.2, 0.9, 0.76])
                    if exact is not None:
                        ax.plot(x, exact[f], '-', color=th.lines[0], lw=1.0, label='exact')
                    for lab, field in curves:
                        ax.plot(x, field[f], '-', color=th.lines[1 + SLOT[lab]], lw=1.1, label=lab)
                    ax.set_xlim(*xlim)
                    ax.set_ylim(*ylim)
                    ax.set_ylabel(ylabel)
                    ax.grid(alpha=0.25)
                    ax.legend(loc='upper right', fontsize=8, frameon=True)
                    fig.text(0.08, 0.04, f'{tlabel} = {t[f]:.4g}', fontsize=9, color=th.text)
                    fig.savefig(tmp / f'f{f:04d}.png')
                    plt.close(fig)
            out.mkdir(parents=True, exist_ok=True)
            base = out / f'{name}-{th.name}'
            scale = f'scale={width}:-2:flags=lanczos'
            subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', str(fps), '-i', str(tmp / 'f%04d.png'),
                            '-vf', scale, '-c:v', 'libvpx-vp9', '-crf', '36', '-b:v', '0', '-pix_fmt', 'yuv420p',
                            f'{base}.webm'], check=True)
            subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', str(fps), '-i', str(tmp / 'f%04d.png'),
                            '-vf', scale, '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-movflags', '+faststart',
                            f'{base}.mp4'], check=True)
            subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-i', str(tmp / f'f{t.size - 1:04d}.png'), '-vf',
                            scale, f'{base}.png'], check=True)
            print(f'  {base}.{{webm,mp4,png}}')
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


def encode(tmp: Path, out: Path, name: str, th: Theme, fps: int, width: int, last: int) -> None:
    """The frames tmp/f%04d.png to webm, mp4 and the poster (the last frame)."""
    out.mkdir(parents=True, exist_ok=True)
    base = out / f'{name}-{th.name}'
    scale = f'scale={width}:-2:flags=lanczos'
    subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', str(fps), '-i', str(tmp / 'f%04d.png'),
                    '-vf', scale, '-c:v', 'libvpx-vp9', '-crf', '36', '-b:v', '0', '-pix_fmt', 'yuv420p',
                    f'{base}.webm'], check=True)
    subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', str(fps), '-i', str(tmp / 'f%04d.png'),
                    '-vf', scale, '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-movflags', '+faststart',
                    f'{base}.mp4'], check=True)
    subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-i', str(tmp / f'f{last:04d}.png'), '-vf', scale,
                    f'{base}.png'], check=True)
    print(f'  {base}.{{webm,mp4,png}}')


def animate(name: str, nframes: int, draw, out: Path, fps: int, width: int, size: tuple[float, float]) -> None:
    """One animation in both themes, each frame drawn by draw(fig, frame, theme)."""
    for th in THEMES:
        tmp = Path(tempfile.mkdtemp(prefix=f'seth-{name}-'))
        try:
            with mpl.rc_context({'figure.facecolor': th.bg, 'savefig.facecolor': th.bg, 'axes.facecolor': th.bg_alt,
                                 'axes.edgecolor': th.mesh, 'axes.labelcolor': th.text, 'text.color': th.text,
                                 'xtick.color': th.mesh, 'ytick.color': th.mesh, 'xtick.labelcolor': th.text,
                                 'ytick.labelcolor': th.text, 'grid.color': th.mesh, 'legend.facecolor': th.bg,
                                 'legend.edgecolor': th.mesh}):
                for f in range(nframes):
                    fig = plt.figure(figsize=size, dpi=100)
                    draw(fig, f, th)
                    fig.savefig(tmp / f'f{f:04d}.png')
                    plt.close(fig)
            encode(tmp, out, name, th, fps, width, nframes - 1)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


def np_rows() -> tuple[np.ndarray, np.ndarray]:
    """The five-point rows (mask 7) of NP_RBF_P and NP_RBF_C, read from src/lib/seth_nonpoly_tables.F90."""
    src = (Path(__file__).resolve().parents[2] / 'src/lib/seth_nonpoly_tables.F90').read_text()

    def table(name: str) -> np.ndarray:
        body = re.search(name + r'\(-2:2,0:7\) = reshape\(\[(.*?)\], \[5,8\]\)', src, re.S).group(1)
        body = re.sub(r'!.*', '', body).replace('&', ' ').replace('_R_P', '')
        vals = [float(Fraction(*(int(float(x)) for x in v.split('/'))) if '/' in v else float(v))
                for v in (t.strip() for t in body.split(',')) if v]
        return np.array(vals).reshape(8, 5)[7]  # column-major: column 7 = mask 7, cells -2..2

    return table('NP_RBF_P'), table('NP_RBF_C')


def adr_morph(out: Path, fps: int, width: int) -> None:
    """The linear Gaussian-RBF row's modified wavenumber, s = 0 .. 0.3 (Liang & Fu 2024, Eqs. 25, 29)."""
    p, c = np_rows()
    th_ = np.linspace(1e-3, np.pi, 400)
    m = np.arange(-2, 3)
    s_all = np.linspace(0.0, 0.3, 61)

    def xit(row: np.ndarray) -> np.ndarray:
        return -1j * (row[None, :] * np.exp(1j * np.outer(th_, m))).sum(1) * (1 - np.exp(-1j * th_))

    lin5 = xit(p)

    def draw(fig, f: int, th: Theme) -> None:
        s = s_all[f]
        z = xit(p + s * s * c)
        a1 = fig.add_axes([0.09, 0.2, 0.38, 0.74])
        a2 = fig.add_axes([0.61, 0.2, 0.37, 0.74])
        a1.plot(th_, th_, '-', color=th.lines[0], lw=0.8, label='exact')
        a1.plot(th_, lin5.real, '--', color=th.lines[3], lw=1.0, label='linear5 ($s = 0$)')
        a1.plot(th_, z.real, '-', color=th.lines[2], lw=1.3, label='RBF row')
        a1.set_xlim(0, np.pi)
        a1.set_ylim(0, 2.2)
        a1.set_xlabel('$\\theta$')
        a1.set_ylabel('Re $\\tilde\\theta$')
        a1.legend(fontsize=7, loc='upper left')
        a2.axhline(0, color=th.lines[0], lw=0.6)
        a2.plot(th_, lin5.imag, '--', color=th.lines[3], lw=1.0)
        a2.plot(th_, z.imag, '-', color=th.lines[2], lw=1.3)
        if s > 0:
            a2.axvline(np.sqrt(30) * s, color=th.mesh, lw=0.8, ls=':')
        a2.set_xlim(0, 1.8)
        a2.set_ylim(-0.04, 0.03)
        a2.set_xlabel('$\\theta$')
        a2.set_ylabel('Im $\\tilde\\theta$ (> 0: growth)')
        for a in (a1, a2):
            a.grid(alpha=0.25)
        tag = 'inside the validated range' if s <= 0.12 else 'outside the validated range ($s > 0.12$)'
        fig.text(0.07, 0.04, f'$s = \\lambda \\Delta x$ = {s:.3f}: {tag}; dotted: the predicted band edge '
                 f'$\\sqrt{{30}}\\,s$', fontsize=8.5, color=th.text)

    print('nonpoly-adr-morph')
    animate('nonpoly-adr-morph', s_all.size, draw, out, fps, width, (7.2, 3.2))


def read_square(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """x, t, u, u_exact, flags per frame of a nonpoly-square-frames file."""
    rows = [ln for ln in path.read_text().splitlines() if ln and not ln.startswith('#')]
    x = np.array(rows[0].split(), float)
    n = x.size
    d = np.array([r.split() for r in rows[1:]], float)
    return x, d[:, 0], d[:, 1:1 + n], d[:, 1 + n:1 + 2 * n], d[:, 1 + 2 * n:1 + 3 * n]


def nonpoly_square(data: Path, out: Path, fps: int, width: int) -> None:
    names = ['teno5', 'teno5_dr', 'teno5_thinc']
    files = [data / f'nonpoly-square-frames-{n}.dat' for n in names]
    if not all(p.exists() for p in files):
        print('nonpoly-square: skipped, missing frames')
        return
    sets = [read_square(p) for p in files]
    x, t = sets[0][0], sets[0][1]
    lab = {'teno5': 'TENO5', 'teno5_dr': 'TENO5-DR', 'teno5_thinc': 'TENO5-THINC'}
    slot = {'teno5': 2, 'teno5_dr': 0, 'teno5_thinc': 1}

    def draw(fig, f: int, th: Theme) -> None:
        ax = fig.add_axes([0.14, 0.40, 0.84, 0.56])
        ax.plot(x, sets[0][3][f], '-', color=th.lines[0], lw=1.0, label='exact')
        for n, st in zip(names, sets):
            ax.plot(x, st[2][f], '-', color=th.lines[1 + slot[n]], lw=1.1, label=lab[n])
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.15, 1.2)
        ax.set_ylabel('u')
        ax.grid(alpha=0.25)
        ax.legend(loc='upper right', fontsize=7, ncol=4)
        sx = fig.add_axes([0.14, 0.16, 0.84, 0.18], sharex=ax)
        for j, (n, st) in enumerate(zip(names[1:], sets[1:])):
            fl = st[4][f] > 0
            sx.plot(x[fl] + 0.5 * (x[1] - x[0]), np.full(fl.sum(), 1 - j), '|', color=th.lines[1 + slot[n]], ms=9)
            sx.text(-0.01, 1 - j, lab[n], ha='right', va='center', fontsize=7, color=th.text)
        sx.set_ylim(-0.6, 1.6)
        sx.set_yticks([])
        sx.set_xlabel('x')
        fig.text(0.02, 0.02, f'periods = {t[f]:.1f}', fontsize=9, color=th.text)

    print('nonpoly-square')
    animate('nonpoly-square', t.size, draw, out, fps, width, (7.2, 3.6))


def nonpoly_sod(data: Path, out: Path, fps: int, width: int) -> None:
    names = ['teno5', 'teno5_dr', 'teno5_thinc']
    files = [data / f'nonpoly-sod-frames-{n}.dat' for n in names]
    if not all(p.exists() for p in files):
        print('nonpoly-sod: skipped, missing frames')
        return
    curves = []
    for n, p in zip(names, files):
        x, t, field, _ = read_frames(p, True)
        curves.append((LABEL.get(n, n), field))
    render('nonpoly-sod', x, t, curves, None, 'density', (0.0, 1.0), (0.05, 1.08), out, fps, width, 't')


def sensor_lambda(rho: np.ndarray, a1: float = 3.0, a2: float = 3.0, cr: float = 0.4) -> np.ndarray:
    """lambda = s/dx of Liang & Fu 2024's sensor (Eqs. 30-32) at every face i+1/2 of rho (constant continuation)."""
    e = 1.0e-40
    v = np.concatenate([[rho[0]] * 2, rho, [rho[-1]] * 2])
    d = np.diff(v)  # d[j] = v[j+1] - v[j]
    eta = (np.abs(2 * d[:-1] * d[1:]) + e) / (d[:-1] ** 2 + d[1:] ** 2 + e)  # eta at the nodes 1..n+2 of v
    etamin = np.minimum(np.minimum(eta[:-2], eta[1:-1]), eta[2:])  # nodes i-1..i+1 of cell i
    m = 1 - np.minimum(1.0, etamin / cr)
    g = (1 - m) ** 4 * (1 + 4 * m)
    return a1 + a2 * (1 - g)


def nonpoly_lambda(data: Path, out: Path, fps: int, width: int) -> None:
    p = data / 'nonpoly-shu-osher-frames-teno5_rbf.dat'
    if not p.exists():
        print('nonpoly-lambda: skipped, missing frames')
        return
    x, t, rho, _ = read_frames(p, True)

    def draw(fig, f: int, th: Theme) -> None:
        ax = fig.add_axes([0.10, 0.44, 0.88, 0.53])
        ax.plot(x, rho[f], '-', color=th.lines[3], lw=1.1, label='TENO5-RBF, N = 200')
        ax.set_xlim(-5, 5)
        ax.set_ylim(0.6, 5.0)
        ax.set_ylabel('density')
        ax.grid(alpha=0.25)
        ax.legend(loc='upper right', fontsize=7)
        lx = fig.add_axes([0.10, 0.16, 0.88, 0.22], sharex=ax)
        lam = sensor_lambda(rho[f])
        lx.plot(x, lam[:x.size], '-', color=th.lines[2], lw=1.0)
        lx.set_ylim(2.8, 6.2)
        lx.set_yticks([3, 6])
        lx.set_ylabel('$\\lambda$')
        lx.grid(alpha=0.25)
        lx.set_xlabel('x')
        fig.text(0.02, 0.02, f't = {t[f]:.3f}', fontsize=9, color=th.text)

    print('nonpoly-lambda')
    animate('nonpoly-lambda', t.size, draw, out, fps, width, (7.2, 3.6))


def nonpoly_blast(data: Path, out: Path, fps: int, width: int) -> None:
    pf, pl = data / 'nonpoly-blast-frames-teno5_dr.dat', data / 'nonpoly-blast-flags-teno5_dr.dat'
    if not (pf.exists() and pl.exists()):
        print('nonpoly-blast: skipped, missing frames')
        return
    x, t, rho, _ = read_frames(pf, True)
    fl = np.array([ln.split() for ln in pl.read_text().splitlines() if ln and not ln.startswith('#')], float)[:, 1:]

    def draw(fig, f: int, th: Theme) -> None:
        ax = fig.add_axes([0.12, 0.40, 0.86, 0.56])
        ax.plot(x, rho[f], '-', color=th.lines[1], lw=1.1, label='TENO5-DR, N = 400 (limiter on)')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 7)
        ax.set_ylabel('density')
        ax.grid(alpha=0.25)
        ax.legend(loc='upper left', fontsize=7)
        sx = fig.add_axes([0.12, 0.16, 0.86, 0.18], sharex=ax)
        xf = x + 0.5 * (x[1] - x[0])
        for lev, col in ((1, th.lines[2]), (2, th.lines[1])):
            sel = fl[f] == lev
            sx.plot(xf[sel], np.full(sel.sum(), 2 - lev), '|', color=col, ms=9)
            sx.text(-0.01, 2 - lev, f'level {lev}', ha='right', va='center', fontsize=7, color=th.text)
        sx.set_ylim(-0.6, 1.6)
        sx.set_yticks([])
        sx.set_xlabel('x')
        fig.text(0.02, 0.02, f't = {t[f]:.4f}', fontsize=9, color=th.text)

    print('nonpoly-blast')
    animate('nonpoly-blast', t.size, draw, out, fps, width, (7.2, 3.6))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--data', type=Path, required=True)
    ap.add_argument('--out', type=Path, default=Path('docs/public/verification/casts'))
    ap.add_argument('--only', nargs='*')
    ap.add_argument('--fps', type=int, default=15)
    ap.add_argument('--width', type=int, default=720)
    args = ap.parse_args()
    if shutil.which('ffmpeg') is None:
        print('ffmpeg not found: no animation rendered')
        return 1
    jobs = {
        'shu-osher': (True, ['weno5_js', 'teno6'], 'density', (-5.0, 5.0), (0.6, 5.0), 't'),
        'blast': (True, ['weno5_js', 'teno6'], 'density', (0.0, 1.0), (0.0, 7.0), 't'),
        'burgers': (False, ['weno5_js', 'teno6'], 'u', (-1.0, 1.0), (-0.65, 1.65), 't'),
        'advection': (False, ['weno5_js', 'wenoz5', 'teno5', 'teno6'], 'u', (-1.0, 1.0), (-0.15, 1.2), 't'),
    }
    if not args.only or 'nonpoly-adr-morph' in args.only:
        adr_morph(args.out, args.fps, args.width)
    for name, fn in (('nonpoly-square', nonpoly_square), ('nonpoly-sod', nonpoly_sod), ('nonpoly-lambda', nonpoly_lambda),
                     ('nonpoly-blast', nonpoly_blast)):
        if not args.only or name in args.only:
            fn(args.data, args.out, args.fps, args.width)
    for name, (euler, schemes, ylabel, xlim, ylim, tlabel) in jobs.items():
        if args.only and name not in args.only:
            continue
        files = [args.data / f'{name}-frames-{s}.dat' for s in schemes]
        if not all(p.exists() for p in files):
            print(f'{name}: skipped, missing frames')
            continue
        print(name)
        curves, exact, x, t = [], None, None, None
        for s, p in zip(schemes, files):
            x, t, field, ex = read_frames(p, euler)
            curves.append((LABEL.get(s, s), field))
            if ex is not None:
                exact = ex
        render(name, x, t, curves, exact, ylabel, xlim, ylim, args.out, args.fps, args.width, tlabel)
    return 0


if __name__ == '__main__':
    sys.exit(main())
