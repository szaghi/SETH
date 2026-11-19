#!/usr/bin/env python3
"""Check that a library's tests do not re-implement another library's operators (SELAH issue #2).

Why this script exists: the time operators live in SELAH only, the space operators in GIDEON and SETH only. A test
that writes its own Runge-Kutta step (in GIDEON or SETH) or its own stencil (in SELAH) is a second source of truth that
drifts from the first: issue #2 removed eleven such copies. This check, the same file in the three repositories, stops
one from coming back unnoticed. It reads the fobos [project] name to know which rule applies:

- SELAH (time operators): no stencil in src/tests/ and docs/examples/src/, i.e. no assignment to an array element
  that combines neighbour values a(i+-m), a(j+-m), a(k+-m) with an arithmetic operator, outside the procedures named
  fill_ghosts (a ghost fill copies values, it is the test problem's boundary condition) and outside calls (a call
  gathers a window for GIDEON); a scalar on the left (a norm, a total variation) is a diagnostic, not a residual;
- GIDEON and SETH (space operators): no Runge-Kutta combination in src/tests/ and docs/examples/src/: the SSP-RK3
  stages (0.75 u + 0.25 ..., u/3 + 2/3 ...) or the classic RK4 sum (k1 + 2 k2 + 2 k3 + k4).

The patterns are regular expressions, not a Fortran parser: a false positive is waived by name in WAIVERS below, with
the reason. `--self-test` plants a copy of each kind in a scratch file and checks that the rule catches it.

Usage: python3 scripts/check_separation.py [--root DIR] [--self-test]     (exit status 1 on a finding)
"""

from __future__ import annotations

import argparse
import re
import sys
import tempfile
from pathlib import Path

# a neighbour index i+1, j-2, k+m ... inside an array reference (an integer or a named offset: i - 0.5_R_P, a cell
# centre, is not one)
NEIGHBOUR = re.compile(r'\(\s*[^()]*\b[ijk]\s*[-+]\s*(?:\d+(?![.\d_])|[A-Za-z]\w*)[^()]*\)')
# an arithmetic operator applied to an array element: a(...) op  or  op a(...)
ARITH = re.compile(r'\)\s*[-+*/]\s*[\w.(]|[-+*/]\s*\w+\s*\(')
RK3 = re.compile(r'0\.75_R_P\s*\*.*0\.25_R_P|/\s*3\._R_P\s*\+\s*2\._R_P\s*/\s*3\._R_P|\(2\._R_P\s*/\s*3\._R_P\)\s*\*')
RK4 = re.compile(r'\bk1\b\s*\+\s*2\w*\s*\*\s*k2\b')

# file name -> (pattern that a waived line matches, reason)
WAIVERS: dict[str, list[tuple[str, str]]] = {
    'selah_adam_reference.F90': [(r'bk\(i,k\+1\)', "ADAM's tableau conversion, copied verbatim: k is a stage "
                                                    'index, not a neighbour')],
    't03_advection.F90': [('.', "SETH tutorial chapter 3 writes its SSP-RK3 by hand to show what a time step is (the "
                                 'later chapters take it from SELAH)')],
    't11_gpu.F90': [('.', "SETH tutorial chapter 11 runs on device memory: its step waits for SELAH's device build")],
    'seth_test_tables.F90': [(r'NP_RBF_P|pair\s*=', 'a row of the non-polynomial tables is a 3/4, 1/4 average of '
                                                       'two stencils: a space-operator coefficient, not a time step')],
}


def project(root: Path) -> str:
    m = re.search(r'(?m)^\s*name\s*=\s*(\S+)', (root / 'fobos').read_text())
    if not m:
        sys.exit('check_separation: no [project] name in fobos')
    return m.group(1)


def statements(path: Path):
    """(line number, statement, enclosing procedure) of a Fortran source, continuation lines joined, comments cut."""
    proc = ''
    buf, start = '', 0
    for n, raw in enumerate(path.read_text(errors='replace').splitlines(), 1):
        code = raw.split('!')[0].rstrip()
        if not buf:
            start = n
        if code.endswith('&'):
            buf += code[:-1] + ' '
            continue
        stmt = (buf + code).strip()
        buf = ''
        m = re.match(r'(?i)(?:\w+\s+)*(subroutine|function)\s+(\w+)', stmt)
        if m and not stmt.lower().startswith('end'):
            proc = m.group(2).lower()
        for part in stmt.split(';'):
            if part.strip():
                yield start, part.strip(), proc


def waived(path: Path, stmt: str) -> bool:
    return any(re.search(p, stmt) for p, _ in WAIVERS.get(path.name, []))


def stencils(path: Path) -> list[str]:
    out = []
    for n, s, proc in statements(path):
        if proc == 'fill_ghosts' or re.match(r'(?i)call\b', s) or '=' not in s or s.lower().startswith('do '):
            continue
        lhs, rhs = s.split('=', 1)
        if '(' not in lhs:  # a scalar (a norm, a sum, a coordinate): a diagnostic, not a residual array
            continue
        if NEIGHBOUR.search(rhs) and ARITH.search(rhs) and not waived(path, s):
            out.append(f'{path}:{n}: stencil: {s[:110]}')
    return out


def rk_steps(path: Path) -> list[str]:
    return [f'{path}:{n}: Runge-Kutta step: {s[:110]}' for n, s, _ in statements(path)
            if (RK3.search(s) or RK4.search(s)) and not waived(path, s)]


def check(root: Path) -> list[str]:
    name = project(root)
    if name == 'SELAH':
        files = sorted((root / 'src/tests').glob('*.[fF]90')) + sorted((root / 'docs/examples/src').glob('*.[fF]90'))
        return [f for p in files for f in stencils(p)]
    if name in ('GIDEON', 'SETH'):
        files = sorted((root / 'src/tests').glob('*.[fF]90')) + sorted((root / 'docs/examples/src').glob('*.[fF]90'))
        return [f for p in files for f in rk_steps(p)]
    sys.exit(f'check_separation: project {name} has no rule (SELAH, GIDEON, SETH)')


def self_test() -> int:
    """Plant a stencil and a Runge-Kutta step in scratch files: each rule must catch its copy."""
    planted = {
        'stencil.F90': ('      r(1,i,1,1,b) = -(u(1,i+1,1,1,b) - u(1,i-1,1,1,b)) / (2._R_P*dx)\n', stencils),
        'gather.F90': ('      w = [u(1,i-1,1,1,b), u(1,i,1,1,b)]\n      call gideon_fv_face_c2(w, fr)\n', stencils),
        'centres.F90': ('      q(1,i,1,1,1) = sin(6.28_R_P * (i - 0.5_R_P) * dx)\n', stencils),
        'tv.F90': ('      v = v + abs(q(1,i+1,1,1,1) - q(1,i,1,1,1))\n', stencils),
        'rk3.F90': ('      u2 = 0.75_R_P*u + 0.25_R_P*(u1 + dt*r)\n', rk_steps),
        'rk4.F90': ('      u = u + DT*(k1 + 2*k2 + 2*k3 + k4)/6\n', rk_steps),
    }
    expect = {'stencil.F90': 1, 'gather.F90': 0, 'centres.F90': 0, 'tv.F90': 0, 'rk3.F90': 1, 'rk4.F90': 1}
    bad = 0
    with tempfile.TemporaryDirectory() as d:
        for name, (src, rule) in planted.items():
            p = Path(d) / name
            p.write_text('program p\n' + src + 'endprogram p\n')
            got = len(rule(p))
            ok = got == expect[name]
            bad += not ok
            print(f'self-test {name:12s} findings {got} (expected {expect[name]}): {"ok" if ok else "FAIL"}')
    return 1 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--root', default=str(Path(__file__).resolve().parent.parent))
    ap.add_argument('--self-test', action='store_true')
    args = ap.parse_args()
    if args.self_test:
        return self_test()
    root = Path(args.root)
    findings = check(root)
    for f in findings:
        print(f)
    print(f'check_separation ({project(root)}): {len(findings)} finding(s)')
    return 1 if findings else 0


if __name__ == '__main__':
    sys.exit(main())
