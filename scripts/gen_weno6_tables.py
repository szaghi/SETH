#!/usr/bin/env python3
"""Generate src/lib/seth_point_weno6_tables.F90: the weights of seth_weno6_point (issue #7, E2).

seth_weno6_point interpolates 6 equispaced point values u(1:6), node k at x = k, at the target x = t4/4, t4 = 3..25
(a node, t4 = 0 mod 4, returns its value and reads no table). Its candidates are the 4-point sub-stencils s = nodes
s..s+3 (s = 1, 2, 3) that contain the target's interval m = clamp(floor(t4/4), 1, 5), i.e. the nodes m and m+1; none
extrapolates but at the two outermost targets, t4 = 3 and 25, where the one candidate is the one-sided cubic of the
end interval. Every weight is derived here from that definition, in exact rational arithmetic:

- the cubic weights W(1:4, s, t4): the Lagrange weights of the cubic through the nodes s..s+3 at x = t4/4;
- the linear weights D(s, t4): the unique weights that make the convex combination of the candidates' cubics equal to
  the Lagrange interpolant of the candidates' union at x (an overdetermined linear system, solved and then checked on
  every node of the union); every one is asserted positive, and zero marks a sub-stencil that is not a candidate.

Then the checks: each D column sums to 1 and each W row to 1; the combination reproduces every polynomial of degree up
to the size of the union minus 1 (6 with three candidates, 5 with two, 4 with one) at x; and the merged table equals
ADAM's four tables (adam_seam_interpolation_library.F90 at ADAM commit 2ea71524, frozen below): the midpoints
t4 = 4 code + 6 are ADAM's restriction (SEAM_W_CUBIC_MID, SEAM_D_WENO6), the odd t4 = 4 p + 2 sub - 3 its fill
(SEAM_W_CUBIC_Q, SEAM_D_FILL_WENO6). Every entry is written as ADAM writes it, an integer numerator over ADAM's
denominator of that table (16, 128 or 320), so that each constant folds to the same double as ADAM's parameter:
seth_test_weno6_point_adam compares the compiled tables and the kernel's results with ADAM's, bit for bit.

Usage: gen_weno6_tables.py [--out src/lib/seth_point_weno6_tables.F90] [--check]
"""

from __future__ import annotations

import argparse
import sys
from fractions import Fraction as Fr
from pathlib import Path

HERE = Path(__file__).resolve().parent
T4 = [t for t in range(3, 26) if t % 4]  # the targets that read the tables (a node returns its value)
NODES = range(1, 7)
FAILS: list[str] = []

# ADAM's tables at ADAM commit 2ea71524 (adam_seam_interpolation_library.F90, lines 179-232), numerators as printed
ADAM_W_CUBIC_MID = ([5, 15, -5, 1], [-1, 9, 9, -1], [1, -5, 15, 5])  # / 16, interval 1..3 of a sub-stencil
ADAM_D_WENO6 = ([16, 0, 0], [10, 6, 0], [3, 10, 3], [0, 6, 10], [0, 0, 16])  # / 16, code 0..4
ADAM_W_CUBIC_Q = ([195, -117, 65, -15], [77, 77, -33, 7], [15, 135, -27, 5], [-7, 105, 35, -5],
                  [-5, 35, 105, -7], [5, -27, 135, 15], [7, -33, 77, 77], [-15, 65, -117, 195])  # / 128, r = 1..8
ADAM_D_FILL_WENO6 = ([320, 0, 0], [320, 0, 0], [320, 0, 0], [220, 100, 0], [180, 140, 0], [77, 198, 45],
                     [45, 198, 77], [0, 140, 180], [0, 100, 220], [0, 0, 320], [0, 0, 320], [0, 0, 320])  # / 320
DEN_D = {0: 320, 2: 16}  # ADAM's denominator by the target's kind: t4 odd (the fill), t4 = 2 mod 4 (the restriction)
DEN_W = {0: 128, 2: 16}


def check(ok: bool, msg: str) -> None:
    if not ok:
        print(f'  FAIL  {msg}')
        FAILS.append(msg)


def lagrange(nodes: list[int], x: Fr) -> list[Fr]:
    """The Lagrange weights of the interpolant through `nodes` at x."""
    out = []
    for j in nodes:
        w = Fr(1)
        for i in nodes:
            if i != j:
                w *= (x - i) / Fr(j - i)
        out.append(w)
    return out


def solve(a: list[list[Fr]], b: list[Fr]) -> list[Fr]:
    """Gauss-Jordan on a square system, exact."""
    n = len(a)
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for c in range(n):
        piv = next(r for r in range(c, n) if m[r][c] != 0)
        m[c], m[piv] = m[piv], m[c]
        for r in range(n):
            if r != c and m[r][c] != 0:
                f = m[r][c] / m[c][c]
                m[r] = [x - f * y for x, y in zip(m[r], m[c])]
    return [m[i][n] / m[i][i] for i in range(n)]


def interval(t4: int) -> int:
    return min(max(t4 // 4, 1), 5)


def derive(t4: int) -> tuple[list[Fr], list[list[Fr]]]:
    """D(1:3) and W(1:4, 1:3) at the target t4, from the definition."""
    x = Fr(t4, 4)
    m = interval(t4)
    cands = [s for s in (1, 2, 3) if s <= m and s + 3 >= m + 1]
    check(bool(cands), f't4 = {t4}: no candidate')
    extrapolates = [s for s in cands if not (s <= x <= s + 3)]
    check(extrapolates == ([] if t4 not in (3, 25) else cands), f't4 = {t4}: a candidate extrapolates')
    w = [[Fr(0)] * 4 for _ in range(3)]
    for s in cands:
        w[s - 1] = lagrange(list(range(s, s + 4)), x)
    union = list(range(cands[0], cands[-1] + 4))
    ell = lagrange(union, x)
    # sum_s D_s W_s(k) = ell(k) on every node k of the union: normal equations, then the residual on every node
    rows = [[w[s - 1][k - s] if s <= k <= s + 3 else Fr(0) for s in cands] for k in union]
    ata = [[sum(r[i] * r[j] for r in rows) for j in range(len(cands))] for i in range(len(cands))]
    atb = [sum(r[i] * e for r, e in zip(rows, ell)) for i in range(len(cands))]
    dc = solve(ata, atb)
    check(all(sum(r[i] * dc[i] for i in range(len(cands))) == e for r, e in zip(rows, ell)),
          f't4 = {t4}: no linear weights reproduce the Lagrange interpolant of the union')
    check(all(v > 0 for v in dc), f't4 = {t4}: a linear weight is not positive')
    d = [Fr(0)] * 3
    for s, v in zip(cands, dc):
        d[s - 1] = v
    return d, w


def properties(t4: int, d: list[Fr], w: list[list[Fr]]) -> None:
    x = Fr(t4, 4)
    check(sum(d) == 1, f't4 = {t4}: the linear weights do not sum to 1')
    cands = [s for s in (1, 2, 3) if d[s - 1] > 0]
    for s in cands:
        check(sum(w[s - 1]) == 1, f't4 = {t4}: the cubic of sub-stencil {s} does not reproduce constants')
    degree = cands[-1] - cands[0] + 3  # the union's size minus 1
    for n in range(degree + 1):
        value = sum(d[s - 1] * sum(w[s - 1][j] * Fr(s + j) ** n for j in range(4)) for s in cands)
        check(value == x ** n, f't4 = {t4}: x^{n} not reproduced')


def adam(t4: int) -> tuple[list[Fr], list[list[Fr]]]:
    """ADAM's D and W at the target t4, as seam_restriction_weno6_1d and seam_fill_weno6_1d select them."""
    w = [[Fr(0)] * 4 for _ in range(3)]
    if t4 % 4 == 2:
        code = (t4 - 6) // 4
        d = [Fr(n, 16) for n in ADAM_D_WENO6[code]]
        for s in (1, 2, 3):
            if d[s - 1] > 0:
                w[s - 1] = [Fr(n, 16) for n in ADAM_W_CUBIC_MID[code - s + 1]]
    else:
        j = (t4 - 1) // 2
        d = [Fr(n, 320) for n in ADAM_D_FILL_WENO6[j - 1]]
        for s in (1, 2, 3):
            if d[s - 1] > 0:
                w[s - 1] = [Fr(n, 128) for n in ADAM_W_CUBIC_Q[(t4 - 4 * s + 3) // 2 - 1]]
    return d, w


def literal(x: Fr, den: int, what: str) -> str:
    """x as ADAM writes it: an integer numerator over ADAM's denominator of the table (0 as 0)."""
    if x == 0:
        return '0._R_P'
    n = x * den
    check(n.denominator == 1, f'{what}: {x} is not an integer over {den}')
    return f'{n.numerator}._R_P/{den}._R_P'


def describe(t4: int, d: list[Fr]) -> str:
    cands = [str(s) for s in (1, 2, 3) if d[s - 1] > 0]
    kind = 'midpoint' if t4 % 4 == 2 else 'quarter'
    edge = ', one-sided' if t4 in (3, 25) else ''
    return f't4 = {t4:2d}, x = {t4}/4 ({kind}{edge}): sub-stencil{"s" if len(cands) > 1 else ""} {", ".join(cands)}'


def emit(tables: dict[int, tuple[list[Fr], list[list[Fr]]]]) -> str:
    head = '''!< SETH, weights of the WENO-6 point interpolation (generated: do not edit).
module seth_point_weno6_tables
!< SETH, weights of the WENO-6 point interpolation (generated: do not edit).
!<
!< Generated by scripts/gen_weno6_tables.py, which derives every weight from its definition in exact rational
!< arithmetic and checks it (issue #7, E2). Six point values, node k at x = k (k = 1..6), the target at x = t4/4,
!< t4 = 3..25; sub-stencil s = nodes s..s+3 (s = 1..3). WENO6_POINT_W(1:4, s, t4): the Lagrange weights of the cubic
!< of sub-stencil s at the target. WENO6_POINT_D(s, t4): the linear weight of sub-stencil s, positive for the
!< sub-stencils that contain the target's interval (the nodes m, m+1, m = clamp(floor(t4/4), 1, 5)), zero for the
!< others; the positive ones make the combination of the candidates' cubics the Lagrange interpolant of their union
!< (order 6, 5 or 4 with three, two or one candidates). The columns of a node (t4 = 0 mod 4) are zero and never read.
!< Each entry is written as ADAM writes it (adam_seam_interpolation_library.F90 at ADAM commit 2ea71524), an integer
!< over 16 at the midpoints (t4 = 2 mod 4) and over 128 (W) or 320 (D) at the quarters (t4 odd): the same doubles.
use penf, only : R_P
implicit none
private
public :: WENO6_POINT_D, WENO6_POINT_W

'''
    dl, wl = [], []
    for t4 in range(3, 26):
        d, w = tables.get(t4, ([Fr(0)] * 3, [[Fr(0)] * 4 for _ in range(3)]))
        note = describe(t4, d) if t4 % 4 else f't4 = {t4:2d}, node {t4 // 4}: not read'
        kind = 2 if t4 % 4 == 2 else 0
        dd, dw = DEN_D[kind], DEN_W[kind]
        dl.append((', '.join(literal(v, dd, f'D t4 = {t4}') for v in d), note))
        for s in (1, 2, 3):
            role = 'not read' if t4 % 4 == 0 or d[s - 1] == 0 else 'candidate'
            wl.append((', '.join(literal(v, dw, f'W t4 = {t4}') for v in w[s - 1]),
                       f't4 = {t4:2d}, sub-stencil {s}: {role}'))

    def array(name: str, bounds: str, extents: str, lines: list[tuple[str, str]], doc: str) -> str:
        out = f'real(R_P), parameter :: {name}({bounds}) = reshape([ & !< {doc}\n'
        width = max(len(v) for v, _ in lines)
        for i, (vals, note) in enumerate(lines):
            sep = ', &' if i < len(lines) - 1 else '], &'
            out += f'   {vals:<{width}}{sep} ! {note}\n'
        return out + f'   [{extents}])\n'

    body = array('WENO6_POINT_D', '1:3,3:25', '3,23', dl, 'linear weights (sub-stencil, target t4)')
    body += array('WENO6_POINT_W', '1:4,1:3,3:25', '4,3,23', wl,
                  'cubic weights (node, sub-stencil, target t4)')
    return head + body + 'endmodule seth_point_weno6_tables\n'


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--out', type=Path, default=HERE.parent / 'src/lib/seth_point_weno6_tables.F90')
    ap.add_argument('--check', action='store_true', help='compare with the file on disk instead of writing it')
    args = ap.parse_args()
    tables = {}
    for t4 in T4:
        d, w = derive(t4)
        properties(t4, d, w)
        da, wa = adam(t4)
        check(d == da, f't4 = {t4}: the linear weights differ from ADAM\'s ({d} != {da})')
        check(w == wa, f't4 = {t4}: the cubic weights differ from ADAM\'s')
        tables[t4] = (d, w)
    print(f'  {len(T4)} targets: the weights derived, the properties and ADAM\'s four tables checked exactly')
    text = emit(tables)
    if FAILS:
        print(f'{len(FAILS)} check(s) failed: nothing written')
        return 1
    if args.check:
        same = args.out.read_text() == text
        print(f'{args.out}: ' + ('up to date' if same else 'DIFFERS from the generator'))
        return 0 if same else 1
    args.out.write_text(text)
    print(f'written {args.out}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
