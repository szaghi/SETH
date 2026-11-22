#!/usr/bin/env python3
"""Generate src/lib/seth_weno_js_tables.F90 from ADAM's WENO initialisers.

The coefficient expressions are copied verbatim from ADAM's `adam_weno_object.F90` (initialize_upwind_S<N> and
initialize_interpolation_S<N>), so every constant folds to the same double as ADAM's run-time tables. Before writing,
each table is checked in exact rational arithmetic:

- the optimal weights of each interface sum to 1;
- every candidate polynomial reproduces constants;
- the optimal combination is exact for polynomials of degree 2S-2 at its face (cell averages for the reconstruction,
  point values for the interpolation);
- the smoothness indicators vanish on constant data and equal 1 on data of unit slope (S >= 2).

Usage: gen_weno_js_tables.py <path to adam_weno_object.F90> <output .F90> <adam commit>

ADAM deleted these initialisers in its commit ba37c4cb (ADAM issue 139, ADAM calls SETH's kernels): take the file of
an earlier ADAM commit, e.g. `git -C adam show 302f093e:src/lib/common/adam_weno_object.F90`, the one SETH's tables
were generated from.
"""

from __future__ import annotations

import re
import sys
import itertools
from fractions import Fraction
from pathlib import Path

ASSIGN = re.compile(r'\b([apd])\(\s*(\d+)\s*,\s*(\d+)\s*(?:,\s*(\d+)\s*)?,\s*S\s*\)\s*=\s*([^;!]+)')
NUMBER = re.compile(r'(\d+\.\d*|\d*\.\d+|\d+)')


def rational(expr: str) -> Fraction:
    """Exact value of an ADAM coefficient expression such as ``-7._R8P/6._R8P`` (emitted as R_P)."""
    text = expr.replace('_R8P', '').replace('_R_P', '').replace(' ', '')
    text = NUMBER.sub(lambda m: f'Fraction("{m.group(1)}")', text)
    return eval(text, {'Fraction': Fraction})  # noqa: S307 - input is our own source tree, digits and /-+ only


def routine(source: str, name: str) -> str:
    match = re.search(rf'subroutine {name}\(self\)(.*?)endsubroutine {name}\b', source, re.S)
    if match is None:
        sys.exit(f'routine {name} not found')
    return match.group(1)


def parse(body: str) -> dict[str, dict[tuple[int, ...], str]]:
    tables: dict[str, dict[tuple[int, ...], str]] = {'a': {}, 'p': {}, 'd': {}}
    for line in body.splitlines():
        code = line.split('!')[0]
        for m in ASSIGN.finditer(code):
            idx = tuple(int(g) for g in m.groups()[1:4] if g is not None)
            tables[m.group(1)][idx] = m.group(5).strip()
    return tables


def cell_average(m: int, n: int) -> Fraction:
    """Average of x**n over the cell [m-1/2, m+1/2]."""
    lo, hi = Fraction(2 * m - 1, 2), Fraction(2 * m + 1, 2)
    return (hi ** (n + 1) - lo ** (n + 1)) / (n + 1)


def check(S: int, t: dict, interpolation: bool, label: str) -> Fraction:
    """Check the tables; return the common scale of the smoothness indicators (1 for the normalised ones)."""
    a, p, d = t['a'], t['p'], t['d']
    need_a = {(f, k) for f in (1, 2) for k in range(S)}
    need_p = {(f, c, k) for f in (1, 2) for c in range(S) for k in range(S)}
    if set(a) != need_a or set(p) != need_p:
        sys.exit(f'{label}: incomplete a/p tables')
    if not interpolation and set(d) != {(i, j, k) for i in range(S) for j in range(S) for k in range(S)}:
        sys.exit(f'{label}: incomplete d table')
    for f in (1, 2):
        if sum(rational(a[f, k]) for k in range(S)) != 1:
            sys.exit(f'{label}: weights of interface {f} do not sum to 1')
        for k in range(S):
            if sum(rational(p[f, c, k]) for c in range(S)) != 1:
                sys.exit(f'{label}: candidate {k} of interface {f} does not reproduce constants')
        face = Fraction(-1, 2) if f == 1 else Fraction(1, 2)  # row 1: face i-1/2, row 2: face i+1/2 (unshifted)
        for n in range(2 * S - 1):
            value = Fraction(0)
            for k in range(S):
                for c in range(S):
                    m = k - c
                    sample = Fraction(m) ** n if interpolation else cell_average(m, n)
                    value += rational(a[f, k]) * rational(p[f, c, k]) * sample
            if value != face ** n:
                sys.exit(f'{label}: interface {f} not exact for degree {n}')
    scale = Fraction(1)
    if not interpolation and S >= 2:
        slopes = set()
        for k in range(S):
            const = sum(rational(d[i, j, k]) for i in range(S) for j in range(S))
            slope = sum(rational(d[i, j, k]) * (k - i) * (k - j) for i in range(S) for j in range(S))
            if const != 0 or slope <= 0:
                sys.exit(f'{label}: smoothness indicator {k} wrong (constant {const}, unit slope {slope})')
            slopes.add(slope)
        if len(slopes) != 1:
            sys.exit(f'{label}: candidates scaled differently: {sorted(slopes)}')
        scale = slopes.pop()
    return scale


def fortran_array(name: str, shape: list[tuple[int, int]], values: list[str], doc: str) -> str:
    bounds = ','.join(f'{lo}:{hi}' for lo, hi in shape)
    extents = ','.join(str(hi - lo + 1) for lo, hi in shape)
    rows = []
    for i in range(0, len(values), 4):
        rows.append('      ' + ', '.join(f"{v.replace('_R8P', '_R_P'):>18}" for v in values[i:i + 4]))
    body = ', &\n'.join(rows)
    return (f'real(R_P), parameter :: {name}({bounds}) = reshape([ & !< {doc}\n'
            f'{body}], [{extents}])\n')


def ordered(table: dict, ranges: list[range]) -> list[str]:
    """Values in Fortran array-element order (first index fastest)."""
    out: list[str] = []
    for combo in itertools.product(*reversed(ranges)):
        out.append(table[tuple(reversed(combo))])
    return out


def main() -> None:
    src_path, out_path, commit = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
    source = src_path.read_text()
    blocks = []
    for S in range(1, 6):
        order = 2 * S - 1
        if S == 1:  # initialize_upwind_S1 assigns the whole arrays to 1
            rec = {'a': {(f, 0): '1._R_P' for f in (1, 2)},
                   'p': {(f, 0, 0): '1._R_P' for f in (1, 2)},
                   'd': {(0, 0, 0): '1._R_P'}}
        else:
            rec = parse(routine(source, f'initialize_upwind_S{S}'))
        itp = parse(routine(source, f'initialize_interpolation_S{S}'))
        scale = check(S, rec, interpolation=False, label=f'reconstruction S={S}')
        check(S, itp, interpolation=True, label=f'interpolation S={S}')
        r2, rs = range(1, 3), range(S)
        name = f'WENO{order}_JS'
        blocks.append(f'! WENO-{order} (S = {S})\n')
        blocks.append(f'real(R_P), parameter :: {name}_D_SCALE = {scale}._R_P !< IS of unit-slope data (1: normalised)\n')
        blocks.append(fortran_array(f'{name}_A', [(1, 2), (0, S - 1)], ordered(rec['a'], [r2, rs]),
                                    'optimal weights (interface, candidate)'))
        blocks.append(fortran_array(f'{name}_P', [(1, 2), (0, S - 1), (0, S - 1)],
                                    ordered(rec['p'], [r2, rs, rs]),
                                    'reconstruction coefficients (interface, cell, candidate)'))
        blocks.append(fortran_array(f'{name}_D', [(0, S - 1), (0, S - 1), (0, S - 1)],
                                    ordered(rec['d'], [rs, rs, rs]), 'smoothness-indicator coefficients'))
        blocks.append(fortran_array(f'{name}_INTERP_A', [(1, 2), (0, S - 1)], ordered(itp['a'], [r2, rs]),
                                    'interpolation optimal weights'))
        blocks.append(fortran_array(f'{name}_INTERP_P', [(1, 2), (0, S - 1), (0, S - 1)],
                                    ordered(itp['p'], [r2, rs, rs]), 'interpolation coefficients'))
    header = f'''!< SETH, WENO-JS coefficient tables (generated: do not edit).
module seth_weno_js_tables
!< SETH, WENO-JS coefficient tables (generated: do not edit).
!<
!< Generated by scripts/gen_weno_js_tables.py from ADAM's adam_weno_object.F90 (commit {commit}); the expressions are
!< ADAM's, verbatim. Every table was checked in exact rational arithmetic by the generator (weights sum, constant
!< reproduction, exactness of the optimal combination up to degree 2S-2, smoothness indicators on constant and
!< unit-slope data). Layout as ADAM's: interface 1 is the face i-1/2 of the centre cell, interface 2 the face i+1/2;
!< P(interface, cell c, candidate k) multiplies the value of cell k-c. `WENO<N>_JS_D_SCALE` is the value of the
!< smoothness indicators on unit-slope data: 1 up to WENO-5; WENO-7 and WENO-9 keep Balsara & Shu's (2000) integer
!< forms, so the same regulariser eps acts relative to IS/scale (eps/scale in normalised terms).
use penf, only : R_P
implicit none
private
public :: R_P
'''
    names = []
    for S in range(1, 6):
        o = 2 * S - 1
        names += [f'WENO{o}_JS_D_SCALE', f'WENO{o}_JS_A', f'WENO{o}_JS_P', f'WENO{o}_JS_D', f'WENO{o}_JS_INTERP_A', f'WENO{o}_JS_INTERP_P']
    publics = ''.join(f'public :: {n}\n' for n in names)
    out_path.write_text(header + publics + '\n' + '\n'.join(blocks) + 'endmodule seth_weno_js_tables\n')
    print(f'{out_path}: 5 schemes, all tables checked')


if __name__ == '__main__':
    main()
