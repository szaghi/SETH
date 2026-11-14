#!/usr/bin/env python3
"""Generate src/lib/seth_z_teno_tables.F90: the descriptors of SETH's WENO-Z and TENO schemes.

Everything is derived in exact rational arithmetic from the definitions listed in notes/m2-sources.md:

- each candidate's face coefficients: the polynomial of degree w-1 matching the cell averages of its w cells,
  evaluated at the face;
- each Jiang-Shu smoothness indicator (Li, Fu & Adams, eq. 10): sum_{j=1}^{w-1} int_{cell 0} (d^j p/dx^j)^2 dx, as a
  quadratic form in the cell averages (unit cells, so the dx powers are absorbed in undivided values);
- the max-order linear weights: the combination of the candidates exact to the highest degree.

Cells are indexed by their offset from the face (cell 0 = [-1, 0] left of it, cell 1 = [0, 1]), the convention of
SETH's kernels. Cross-checks run before anything is written (the script exits non-zero on any failure):

1. every candidate is exact for polynomials up to degree w-1;
2. the TENO candidates' face coefficients equal Li, Fu & Adams (2020) Table 1;
3. the indicators of the upwind 3-, 4- and 5-point stencils equal ADAM's WENO-5/7/9 tables divided by 1, 240, 5040;
4. the max-order linear weights are exact to degree K-1 and are the WENO weights where K = 2S-1;
5. the order of every tau on smooth data, measured exactly on a Taylor polynomial of exp.

Usage: gen_z_teno_tables.py <adam_weno_object.F90> <output .F90>
"""

from __future__ import annotations

import importlib.util
import math
import sys
from fractions import Fraction as Fr
from pathlib import Path

sys.dont_write_bytecode = True  # importing the sibling generator must not leave scripts/__pycache__ behind
HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location('gen_weno_js_tables', HERE / 'gen_weno_js_tables.py')
_js = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_js)

# ---------------------------------------------------------------------------------------------------------- polynomials
Poly = list  # coefficients, lowest degree first


def padd(a: Poly, b: Poly) -> Poly:
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def pmul(a: Poly, b: Poly) -> Poly:
    out = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def pder(a: Poly) -> Poly:
    return [i * a[i] for i in range(1, len(a))] or [Fr(0)]


def pint(a: Poly, lo: Fr, hi: Fr) -> Fr:
    return sum(c * (hi ** (i + 1) - lo ** (i + 1)) / (i + 1) for i, c in enumerate(a))


def peval(a: Poly, x: Fr) -> Fr:
    return sum(c * x ** i for i, c in enumerate(a))


def cell_avg_monomial(m: int, n: int) -> Fr:
    """Average of x**n over cell m = [m-1, m]."""
    return (Fr(m) ** (n + 1) - Fr(m - 1) ** (n + 1)) / (n + 1)


def solve(a: list[list[Fr]], b: list[Fr]) -> list[Fr]:
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


def rank(a: list[list[Fr]]) -> int:
    m = [row[:] for row in a]
    r = 0
    for c in range(len(m[0]) if m else 0):
        piv = next((i for i in range(r, len(m)) if m[i][c] != 0), None)
        if piv is None:
            continue
        m[r], m[piv] = m[piv], m[r]
        for i in range(len(m)):
            if i != r and m[i][c] != 0:
                f = m[i][c] / m[r][c]
                m[i] = [x - f * y for x, y in zip(m[i], m[r])]
        r += 1
    return r


def basis(cells: list[int]) -> list[Poly]:
    """Lagrange-like basis: basis[j] is the polynomial whose cell averages are 1 on cells[j], 0 on the others."""
    w = len(cells)
    a = [[cell_avg_monomial(m, n) for n in range(w)] for m in cells]
    out = []
    for j in range(w):
        rhs = [Fr(1) if i == j else Fr(0) for i in range(w)]
        out.append(solve(a, rhs))
    return out


def face_coefficients(cells: list[int]) -> list[Fr]:
    return [peval(b, Fr(0)) for b in basis(cells)]


def indicator(cells: list[int]) -> list[list[Fr]]:
    """Symmetric matrix B of the Jiang-Shu indicator on cell 0: beta = sum_ab B[a][b] u[a] u[b]."""
    bs = basis(cells)
    w = len(cells)
    mat = [[Fr(0)] * w for _ in range(w)]
    for j in range(1, w):
        ders = []
        for b in bs:
            d = b
            for _ in range(j):
                d = pder(d)
            ders.append(d)
        for a in range(w):
            for c in range(w):
                mat[a][c] += pint(pmul(ders[a], ders[c]), Fr(-1), Fr(0))
    return mat


def quad(mat: list[list[Fr]], u: list[Fr]) -> Fr:
    return sum(mat[a][c] * u[a] * u[c] for a in range(len(u)) for c in range(len(u)))


# ------------------------------------------------------------------------------------------------------------- schemes
def weno_cells(S: int) -> list[list[int]]:
    """WENO candidate k: cells k-S+1..k (beta_0 on the most upwind stencil)."""
    return [list(range(k - S + 1, k + 1)) for k in range(S)]


TENO_CELLS = [[-1, 0, 1], [0, 1, 2], [-2, -1, 0], [0, 1, 2, 3]]  # Li, Fu & Adams 2020, Table 1 (i -> 0)
TABLE1 = [[Fr(-1, 6), Fr(5, 6), Fr(2, 6)], [Fr(2, 6), Fr(5, 6), Fr(-1, 6)], [Fr(2, 6), Fr(-7, 6), Fr(11, 6)],
          [Fr(3, 12), Fr(13, 12), Fr(-5, 12), Fr(1, 12)]]
# TENO6-opt of Fu, Hu & Adams 2016, Table 5 (llm-wiki weno-teno targeted-eno.md, high-order-capturing teno.md)
TENO6_OPT_D = ['0.462_R8P', '0.3_R8P', '0.054_R8P', '0.184_R8P']
# the optimised weights of the NENO6 paper (arXiv:2207.08500, section 2.2), a later and different optimisation
TENO6_NENO6_D = ['0.4294317718960898_R8P', '0.1727270875843552_R8P', '0.0855682281039113_R8P', '0.312272912415645_R8P']

SCHEMES = {
    # name: (candidate cells, tau coefficients on the candidates, global stencil or None, global coefficient,
    #        weights 'max' or quoted strings, family, tau order required, description)
    'WENOZ3': (weno_cells(2), [1, -1], None, 0, 'max', 'z', 3, 'WENO-Z3: tau_3 = |beta_0 - beta_1|'),
    'WENOZ5': (weno_cells(3), [1, 0, -1], None, 0, 'max', 'z', 5, 'WENO-Z5: tau_5 = |beta_0 - beta_2|'),
    # tau_7 orders measured exactly: Don & Borges 7, Castro et al. 6 (Rathan & Naga Raju state 8 and 7)
    'WENOZ7': (weno_cells(4), [1, 3, -3, -1], None, 0, 'max', 'z', 7,
               'WENO-Z7: tau_7 = |beta_0 + 3 beta_1 - 3 beta_2 - beta_3| (Don & Borges)'),
    # Castro et al.'s rule tau = |beta_0 - beta_{r-1}| for odd r (Li et al. 2022b, eq. 15), here r = 5: O(dx^7)
    'WENOZ9': (weno_cells(5), [1, 0, 0, 0, -1], None, 0, 'max', 'z', 7, 'WENO-Z9: tau_9 = |beta_0 - beta_4|'),
    'TENO5': (TENO_CELLS[:3], [0, 1, -1], None, 0, 'max', 'teno', 5, 'TENO5: tau_5 = |beta_1 - beta_2|'),
    'TENO6': (TENO_CELLS, [Fr(-4, 6), Fr(-1, 6), Fr(-1, 6), 0], list(range(-2, 4)), 1, 'max', 'teno', 6,
              'TENO6: tau_6 = beta_6 - (beta_1 + beta_2 + 4 beta_0)/6, max-order linear weights'),
    'TENO6_OPT': (TENO_CELLS, [Fr(-4, 6), Fr(-1, 6), Fr(-1, 6), 0], list(range(-2, 4)), 1, TENO6_OPT_D, 'teno', 6,
                  'TENO6-opt: TENO6 with the optimised linear weights of Fu, Hu & Adams 2016 (5th order)'),
    'TENO6_NENO6': (TENO_CELLS, [Fr(-4, 6), Fr(-1, 6), Fr(-1, 6), 0], list(range(-2, 4)), 1, TENO6_NENO6_D, 'teno', 6,
                    'TENO6 with the optimised linear weights of the NENO6 paper (4th order)'),
}


def max_order_weights(cands: list[list[int]], label: str) -> tuple[list[Fr], int]:
    """Linear weights exact to the highest degree; returns them and that degree + 1 (the order)."""
    nc = len(cands)
    coefs = [face_coefficients(c) for c in cands]

    def moment(k: int, n: int) -> Fr:
        return sum(coefs[k][j] * cell_avg_monomial(m, n) for j, m in enumerate(cands[k]))

    # exactness for degrees 0, 1, 2, ...: the low degrees hold for every candidate (rows dependent on degree 0), so take
    # the first nc independent rows in degree order; the weights are then unique, and the order is the first degree
    # they miss
    rows, rhs, n = [], [], 0
    while len(rows) < nc:
        row = [moment(k, n) for k in range(nc)]
        trial = rows + [row]
        if rank(trial) == len(trial):
            rows.append(row)
            rhs.append(Fr(1) if n == 0 else Fr(0))
        elif n > 0 and any(row):  # dependent row: it must already hold, else no consistent max-order weights
            pass
        n += 1
        if n > 4 * nc + 8:
            sys.exit(f'{label}: cannot find {nc} independent exactness rows')
    d = solve(rows, rhs)
    order = 0
    while sum(d[k] * moment(k, order) for k in range(nc)) == (1 if order == 0 else 0):
        order += 1
    if any(x <= 0 for x in d):
        sys.exit(f'{label}: a max-order linear weight is not positive: {d}')
    return d, order


def taylor_exp_averages(h: Fr, cells: range, degree: int = 16) -> dict[int, Fr]:
    """Exact cell averages of the degree-`degree` Taylor polynomial of exp on cells [(m-1)h, mh]."""
    coef = [Fr(1, math.factorial(n)) for n in range(degree + 1)]
    prim = [Fr(0)] + [c / (n + 1) for n, c in enumerate(coef)]
    return {m: (peval(prim, m * h) - peval(prim, (m - 1) * h)) / h for m in cells}


def check_fail(msg: str) -> None:
    sys.exit('CHECK FAILED: ' + msg)


def run_checks(adam_source: str) -> None:
    # 1-2: candidates and Table 1
    for cells in [c for s in (2, 3, 4, 5) for c in weno_cells(s)] + TENO_CELLS + [list(range(-2, 4))]:
        coefs = face_coefficients(cells)
        for n in range(len(cells)):
            if sum(c * cell_avg_monomial(m, n) for c, m in zip(coefs, cells)) != (1 if n == 0 else 0):
                check_fail(f'candidate {cells} not exact to degree {n}')
    for k, cells in enumerate(TENO_CELLS):
        if face_coefficients(cells) != TABLE1[k]:
            check_fail(f'TENO candidate {k} differs from Table 1: {face_coefficients(cells)}')
    # 3: indicators vs ADAM's tables (row 2, stencil s1 = k: cells k-S+1..k, IS = sum d(s3,s2,k) v(k-s3) v(k-s2))
    for S, scale in ((3, 1), (4, 240), (5, 5040)):
        d = _js.parse(_js.routine(adam_source, f'initialize_upwind_S{S}'))['d']
        for k, cells in enumerate(weno_cells(S)):
            mine = indicator(cells)
            adam = [[Fr(0)] * S for _ in range(S)]
            for s3 in range(S):
                for s2 in range(S):
                    a, c = cells.index(k - s3), cells.index(k - s2)
                    val = _js.rational(d[s3, s2, k])
                    adam[a][c] += val / 2
                    adam[c][a] += val / 2
            if any(adam[a][c] != scale * mine[a][c] for a in range(S) for c in range(S)):
                check_fail(f'indicator of S={S} candidate {k} differs from ADAM / {scale}')
    # 4: max-order weights
    for S, expect in ((2, [Fr(1, 3), Fr(2, 3)]), (3, [Fr(1, 10), Fr(6, 10), Fr(3, 10)])):
        d, order = max_order_weights(weno_cells(S), f'WENO S={S}')
        if d != expect or order != 2 * S - 1:
            check_fail(f'WENO S={S} linear weights {d} order {order}')
    d, order = max_order_weights(TENO_CELLS[:3], 'TENO5')
    if d != [Fr(6, 10), Fr(3, 10), Fr(1, 10)] or order != 5:
        check_fail(f'TENO5 linear weights {d} order {order}')
    d, order = max_order_weights(TENO_CELLS, 'TENO6')
    if d != [Fr(9, 20), Fr(3, 10), Fr(1, 20), Fr(1, 5)] or order != 6:  # = Fu et al. 2016 Table 5: 0.45, 0.3, 0.05, 0.2
        check_fail(f'TENO6 linear weights {d} reach order {order}, not Fu et al.\'s 6th-order ones')
    # the derived WENO linear weights equal ADAM's (row 2: the left-biased state at i+1/2, candidate k = cells k-S+1..k)
    for S in range(2, 6):
        a = _js.parse(_js.routine(adam_source, f'initialize_upwind_S{S}'))['a']
        d, order = max_order_weights(weno_cells(S), f'WENO S={S}')
        if d != [_js.rational(a[2, k]) for k in range(S)] or order != 2 * S - 1:
            check_fail(f'WENO S={S}: derived linear weights {d} differ from ADAM\'s')
    # the order of the quoted optimised TENO6 weights: Fu's lose one order, NENO6's two
    for label, wq, want in (('TENO6-opt', TENO6_OPT_D, 5), ('TENO6-NENO6', TENO6_NENO6_D, 4)):
        dq = [Fr(w.replace('_R8P', '')) for w in wq]
        coefs = [face_coefficients(c) for c in TENO_CELLS]
        got = 0
        for n in range(8):
            r = sum(dq[k] * sum(coefs[k][j] * cell_avg_monomial(m, n) for j, m in enumerate(TENO_CELLS[k]))
                    for k in range(4)) - (1 if n == 0 else 0)
            if abs(float(r)) > 1e-12:
                break
            got += 1
        print(f'  {label}: linear-scheme order {got} (expected {want})')
        if got != want:
            check_fail(f'{label}: linear-scheme order {got}, expected {want}')
    # 5: tau orders
    for name, (cands, tc, glob, tg, _w, _f, need, _doc) in SCHEMES.items():
        mats = [indicator(c) for c in cands]
        gmat = indicator(glob) if glob else None
        taus = []
        for e in range(4, 8):
            h = Fr(1, 2 ** e)
            avg = taylor_exp_averages(h, range(-4, 6))
            tau = sum(Fr(c) * quad(m, [avg[x] for x in cl]) for c, m, cl in zip(tc, mats, cands))
            if gmat:
                tau += tg * quad(gmat, [avg[x] for x in glob])
            taus.append(abs(tau))
        order = math.log2(taus[-2] / taus[-1])
        print(f'  {name}: tau order {order:.2f} (required {need})')
        if order < need - 0.1:
            check_fail(f'{name}: tau of order {order:.2f} < {need}')


# ---------------------------------------------------------------------------------------------------------- emission
def fstr(x) -> str:
    if isinstance(x, str):
        return x
    x = Fr(x)
    if abs(x.numerator) >= 2 ** 53 or x.denominator >= 2 ** 53:
        check_fail(f'{x} is not a ratio of exactly representable integers')
    return f'{x.numerator}._R8P' if x.denominator == 1 else f'{x.numerator}._R8P/{x.denominator}._R8P'


def farray(name: str, bounds: str, extents: str, values: list[str], doc: str, kind: str = 'real(R8P)') -> str:
    rows = [', '.join(f'{v:>26}' for v in values[i:i + 4]) for i in range(0, len(values), 4)]
    body = ', &\n      '.join(rows)
    return f'{kind}, parameter :: {name}({bounds}) = reshape([ & !< {doc}\n      {body}], [{extents}])\n'


def upper(mat: list[list[Fr]]) -> list[list[Fr]]:
    """The quadratic form as an upper-triangular matrix: U[a][a] = B[a][a], U[a][c] = 2 B[a][c] for a < c."""
    w = len(mat)
    return [[mat[a][c] if a == c else (2 * mat[a][c] if a < c else Fr(0)) for c in range(w)] for a in range(w)]


def emit(name: str) -> tuple[str, list[str]]:
    cands, tc, glob, tg, wts, family, _need, doc = SCHEMES[name]
    nc = len(cands)
    maxw = max(len(c) for c in cands)
    lo = min(min(c) for c in cands)
    umax = max(max(c) for c in cands)
    d = max_order_weights(cands, name)[0] if wts == 'max' else wts
    out = [f'! {doc}\n']
    names = []

    def scalar(n: str, kind: str, val: str, doc_: str) -> None:
        out.append(f'{kind}, parameter :: {name}_{n} = {val} !< {doc_}\n')
        names.append(f'{name}_{n}')

    scalar('NC', 'integer(I4P)', f'{nc}_I4P', 'candidate stencils')
    scalar('LO', 'integer(I4P)', f'{lo}_I4P', 'first cell of the stencil arrays (bounds LO:1-LO)')
    scalar('UMAX', 'integer(I4P)', f'{umax}_I4P', 'last cell any candidate reads (left bias)')
    scalar('MAXW', 'integer(I4P)', f'{maxw}_I4P', 'widest candidate')
    scalar('TG', 'real(R8P)', fstr(tg), 'coefficient of the global indicator in tau (0: none)')
    out.append(f'integer(I4P), parameter :: {name}_CLO(0:{nc - 1}) = [{", ".join(f"{min(c)}_I4P" for c in cands)}]'
               ' !< first cell of each candidate\n')
    out.append(f'integer(I4P), parameter :: {name}_CW(0:{nc - 1}) = [{", ".join(f"{len(c)}_I4P" for c in cands)}]'
               ' !< width of each candidate\n')
    names += [f'{name}_CLO', f'{name}_CW']
    pvals = []
    uvals = []
    for k in range(nc):  # Fortran order: first index fastest
        coefs = face_coefficients(cands[k]) + [Fr(0)] * (maxw - len(cands[k]))
        pvals += [fstr(x) for x in coefs]
    for k in range(nc):
        u = upper(indicator(cands[k]))
        w = len(cands[k])
        full = [[u[a][c] if a < w and c < w else Fr(0) for c in range(maxw)] for a in range(maxw)]
        uvals += [fstr(full[a][c]) for c in range(maxw) for a in range(maxw)]
    out.append(farray(f'{name}_P', f'0:{maxw - 1},0:{nc - 1}', f'{maxw},{nc}', pvals,
                      'face coefficients P(j,k) of cell CLO(k)+j'))
    out.append(farray(f'{name}_U', f'0:{maxw - 1},0:{maxw - 1},0:{nc - 1}', f'{maxw},{maxw},{nc}', uvals,
                      'indicators, upper-triangular U(a,b,k) on cells CLO(k)+a, CLO(k)+b'))
    out.append(farray(f'{name}_D', f'0:{nc - 1}', f'{nc}', [fstr(x) for x in d], 'linear weights'))
    out.append(farray(f'{name}_TC', f'0:{nc - 1}', f'{nc}', [fstr(x) for x in tc], 'tau coefficients on the beta_k'))
    names += [f'{name}_P', f'{name}_U', f'{name}_D', f'{name}_TC']
    g = glob if glob else [lo]
    gm = upper(indicator(glob)) if glob else [[Fr(0)]]
    gw = len(g)
    scalar('GLO', 'integer(I4P)', f'{g[0]}_I4P', 'first cell of the global stencil')
    scalar('GW', 'integer(I4P)', f'{gw}_I4P', 'width of the global stencil')
    out.append(farray(f'{name}_UG', f'0:{gw - 1},0:{gw - 1}', f'{gw},{gw}',
                      [fstr(gm[a][c]) for c in range(gw) for a in range(gw)], 'global indicator, upper-triangular'))
    names.append(f'{name}_UG')
    return ''.join(out), names


def main() -> None:
    adam = Path(sys.argv[1]).read_text()
    out_path = Path(sys.argv[2])
    print('checks:')
    run_checks(adam)
    blocks, names = [], []
    for n in SCHEMES:
        b, ns = emit(n)
        blocks.append(b)
        names += ns
    header = '''!< SETH, WENO-Z and TENO scheme descriptors (generated: do not edit).
module seth_z_teno_tables
!< SETH, WENO-Z and TENO scheme descriptors (generated: do not edit).
!<
!< Generated by scripts/gen_z_teno_tables.py, every coefficient derived in exact rational arithmetic from the
!< definitions in notes/m2-sources.md and cross-checked (candidates exact, TENO Table 1, ADAM's WENO-5/7/9 indicators,
!< max-order weights, tau orders). Cells by offset from the face i+1/2 (cell 0 left of it); left-biased stencils.
use penf, only : I4P, R8P
implicit none
private
'''
    publics = ''.join(f'public :: {n}\n' for n in names)
    out_path.write_text(header + publics + '\n' + '\n'.join(blocks) + 'endmodule seth_z_teno_tables\n')
    print(f'{out_path}: {len(SCHEMES)} schemes written')


if __name__ == '__main__':
    main()
