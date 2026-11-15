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

Usage: make_doc_animations.py --data DIR [--out DIR] [--only NAME ...] [--fps N] [--width PX]
Needs numpy, matplotlib and ffmpeg on PATH.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
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
LABEL = {'weno5_js': 'WENO5-JS', 'wenoz5': 'WENO-Z5', 'teno5': 'TENO5', 'teno6': 'TENO6'}
# the colour follows the scheme, as on the figures (index into Theme.lines after the exact solution)
SLOT = {'WENO5-JS': 0, 'WENO-Z5': 1, 'TENO5': 2, 'TENO6': 3}


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
