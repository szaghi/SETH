#!/usr/bin/env python3
"""Check SETH's tables and constants against the papers they come from, in exact arithmetic.

The numbers below are transcribed from the papers (read in the paper; the PDFs' text) and checked twice: against an
exact derivation (the algebra of scripts/gen_z_teno_tables.py: Jiang-Shu indicators from their definition, the
candidate reconstructions, the maximum-order weights) and against what SETH compiles (its source files, parsed).

- Balsara & Shu (2000), J. Comput. Phys. 160:405, section III.a-b (r = 4, 5): the smoothness indicators IS_k^r, the
  candidate coefficients a_{k,i}^r (Tables I-II), the optimal weights C_k^r and the optimal combination u^L_{j+1/2};
  SETH's WENO7_JS_* and WENO9_JS_* (src/lib/seth_weno_js_tables.F90).
- Huang, Liang & Fu (2023), Commun. Comput. Phys. (arXiv 2303.10020), Eqs. (18)-(20), the standard TENO5-A of Fu,
  Hu & Adams (AIAA J. 57(3):1143, 2019): C_r, xi, alpha_1, alpha_2 (src/lib/seth_teno.F90).
- Fu (2023), Arch. Comput. Methods Eng. 30:2493, section 3.4, Eqs. (21)-(23), the TENO-based detector of Fu (Commun.
  Comput. Phys. 26:973, 2019): the power q and the threshold C_T (src/lib/seth_sensors.F90).

Usage (from the repository root): python3 scripts/check_papers.py [--tables FILE] [--teno FILE] [--sensors FILE]
Exit status 0 when every check passes. No ADAM input.
"""

from __future__ import annotations

import argparse
import importlib.util
import re
import sys
from fractions import Fraction as Fr
from pathlib import Path

sys.dont_write_bytecode = True  # importing the sibling generator must not leave scripts/__pycache__ behind
HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location('gen_z_teno_tables', HERE / 'gen_z_teno_tables.py')
_zt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_zt)

FAILED: list[str] = []


def fail(msg: str) -> None:
    print(f'  FAIL: {msg}')
    FAILED.append(msg)


def F(s: str) -> Fr:
    return Fr(s)


# ---------------------------------------------------------------------------------------------------------------------
# Balsara & Shu 2000, section III. Cells by offset from j (the cell left of the face j+1/2): -3 = u_{j-3}.
# Indicators as printed: {(p, q): coefficient of u_p u_q}, p <= q.
BS_IS = {
    4: [
        {(-3, -3): 547, (-3, -2): -3882, (-3, -1): 4642, (-3, 0): -1854, (-2, -2): 7043, (-2, -1): -17246,
         (-2, 0): 7042, (-1, -1): 11003, (-1, 0): -9402, (0, 0): 2107},
        {(-2, -2): 267, (-2, -1): -1642, (-2, 0): 1602, (-2, 1): -494, (-1, -1): 2843, (-1, 0): -5966,
         (-1, 1): 1922, (0, 0): 3443, (0, 1): -2522, (1, 1): 547},
        {(-1, -1): 547, (-1, 0): -2522, (-1, 1): 1922, (-1, 2): -494, (0, 0): 3443, (0, 1): -5966, (0, 2): 1602,
         (1, 1): 2843, (1, 2): -1642, (2, 2): 267},
        {(0, 0): 2107, (0, 1): -9402, (0, 2): 7042, (0, 3): -1854, (1, 1): 11003, (1, 2): -17246, (1, 3): 4642,
         (2, 2): 7043, (2, 3): -3882, (3, 3): 547},
    ],
    5: [
        {(-4, -4): 22658, (-4, -3): -208501, (-4, -2): 364863, (-4, -1): -288007, (-4, 0): 86329,
         (-3, -3): 482963, (-3, -2): -1704396, (-3, -1): 1358458, (-3, 0): -411487,
         (-2, -2): 1521393, (-2, -1): -2462076, (-2, 0): 758823, (-1, -1): 1020563, (-1, 0): -649501,
         (0, 0): 107918},
        {(-3, -3): 6908, (-3, -2): -60871, (-3, -1): 99213, (-3, 0): -70237, (-3, 1): 18079,
         (-2, -2): 138563, (-2, -1): -464976, (-2, 0): 337018, (-2, 1): -88297,
         (-1, -1): 406293, (-1, 0): -611976, (-1, 1): 165153, (0, 0): 242723, (0, 1): -140251, (1, 1): 22658},
        {(-2, -2): 6908, (-2, -1): -51001, (-2, 0): 67923, (-2, 1): -38947, (-2, 2): 8209,
         (-1, -1): 104963, (-1, 0): -299076, (-1, 1): 179098, (-1, 2): -38947,
         (0, 0): 231153, (0, 1): -299076, (0, 2): 67923, (1, 1): 104963, (1, 2): -51001, (2, 2): 6908},
        {(-1, -1): 22658, (-1, 0): -140251, (-1, 1): 165153, (-1, 2): -88297, (-1, 3): 18079,
         (0, 0): 242723, (0, 1): -611976, (0, 2): 337018, (0, 3): -70237,
         (1, 1): 406293, (1, 2): -464976, (1, 3): 99213, (2, 2): 138563, (2, 3): -60871, (3, 3): 6908},
        {(0, 0): 107918, (0, 1): -649501, (0, 2): 758823, (0, 3): -411487, (0, 4): 86329,
         (1, 1): 1020563, (1, 2): -2462076, (1, 3): 1358458, (1, 4): -288007,
         (2, 2): 1521393, (2, 3): -1704396, (2, 4): 364863, (3, 3): 482963, (3, 4): -208501, (4, 4): 22658},
    ],
}
# Tables I-II: a_{k,i}, i = 1..r, on the cells k-r+1 .. k (left to right)
BS_A = {
    4: [['-1/4', '13/12', '-23/12', '25/12'], ['1/12', '-5/12', '13/12', '1/4'],
        ['-1/12', '7/12', '7/12', '-1/12'], ['1/4', '13/12', '-5/12', '1/12']],
    5: [['1/5', '-21/20', '137/60', '-163/60', '137/60'], ['-1/20', '17/60', '-43/60', '77/60', '1/5'],
        ['1/30', '-13/60', '47/60', '9/20', '-1/20'], ['-1/20', '9/20', '47/60', '-13/60', '1/30'],
        ['1/5', '77/60', '-43/60', '17/60', '-1/20']],
}
BS_C = {4: ['1/35', '12/35', '18/35', '4/35'], 5: ['1/126', '10/63', '10/21', '20/63', '5/126']}
# u^L_{j+1/2} with the optimal weights, cells -r+1 .. r-1
BS_UL = {4: ['-1/140', '5/84', '-101/420', '319/420', '107/210', '-19/210', '1/105'],
         5: ['1/630', '-41/2520', '199/2520', '-641/2520', '1879/2520', '275/504', '-61/504', '11/504', '-1/504']}
BS_SCALE = {4: 240, 5: 5040}  # the printed integers are the indicator times this (SETH's WENO<2r-1>_JS_D_SCALE)

# Huang, Liang & Fu 2023, Eqs. (19)-(20): the standard TENO5-A
HLF_TENO5_A = {'SETH_TENO5_A_CR': F('0.24'), 'SETH_TENO5_A_XI': F('1e-3'), 'SETH_TENO5_A_A1': F('10'),
               'SETH_TENO5_A_A2': F('5')}
# Fu 2023, Eqs. (21)-(23): gamma_k = (beta_k + eps)^-q, q = 6; smooth iff every chi_k > C_T = 5e-4
FU_Q = 6
FU_CT = F('5e-4')


# ---------------------------------------------------------------------------------------------------------------------
def fortran_arrays(text: str) -> dict[str, tuple[list[tuple[int, int]], list[Fr]]]:
    """Every `real(R_P), parameter :: NAME(lo:hi, ...) = reshape([...], [...])` of a SETH source, exact values."""
    out = {}
    pat = re.compile(r'real\(R_P\), parameter :: (\w+)\(([^)]*)\) = reshape\(\[(.*?)\], \[[^\]]*\]\)', re.S)
    for m in pat.finditer(text):
        dims = []
        for d in m.group(2).split(','):
            lo, hi = d.split(':')
            dims.append((int(lo), int(hi)))
        body = re.sub(r'!<[^\n]*', '', m.group(3)).replace('&', ' ')
        vals = [_zt._js.rational(v) for v in body.split(',') if v.strip()]
        out[m.group(1)] = (dims, vals)
    return out


def element(arr: tuple[list[tuple[int, int]], list[Fr]], *idx: int) -> Fr:
    """Column-major element of a parsed Fortran array."""
    dims, vals = arr
    off, stride = 0, 1
    for (lo, hi), i in zip(dims, idx):
        off += (i - lo) * stride
        stride *= hi - lo + 1
    return vals[off]


def scalar(text: str, name: str) -> Fr:
    m = re.search(rf'real\(R_P\), parameter :: {name}\s*=\s*([^!\n]+)', text)
    if m is None:
        fail(f'{name} not found')
        return Fr(0)
    text = m.group(1).strip()
    if '/' in text:
        return _zt._js.rational(text)
    return Fr(text.replace('_R_P', ''))  # exact decimal, exponents included (1.e-3)


def as_matrix(form: dict[tuple[int, int], int], cells: list[int]) -> list[list[Fr]]:
    """Symmetric matrix of a printed quadratic form sum c_pq u_p u_q on the given cells."""
    n = len(cells)
    mat = [[Fr(0)] * n for _ in range(n)]
    for (p, q), c in form.items():
        a, b = cells.index(p), cells.index(q)
        if a == b:
            mat[a][a] += c
        else:
            mat[a][b] += Fr(c, 2)
            mat[b][a] += Fr(c, 2)
    return mat


def check_balsara_shu(tables: str) -> None:
    arrays = fortran_arrays(tables)
    for r in (4, 5):
        name = f'WENO{2 * r - 1}_JS'
        cands = _zt.weno_cells(r)  # candidate k: cells k-r+1 .. k
        derived_c, order = _zt.max_order_weights(cands, f'WENO r={r}')
        if derived_c != [F(c) for c in BS_C[r]] or order != 2 * r - 1:
            fail(f'r={r}: Balsara & Shu C_k {BS_C[r]} differ from the derived maximum-order weights {derived_c}')
        row = [Fr(0)] * (2 * r - 1)
        for k, cells in enumerate(cands):
            # the indicator: printed = scale x derived from the Jiang-Shu definition = SETH's D
            printed = as_matrix(BS_IS[r][k], cells)
            derived = _zt.indicator(cells)
            if any(printed[a][b] != BS_SCALE[r] * derived[a][b] for a in range(r) for b in range(r)):
                fail(f'r={r} k={k}: Balsara & Shu IS differs from {BS_SCALE[r]} x the Jiang-Shu indicator')
            seth = [[Fr(0)] * r for _ in range(r)]
            for s3 in range(r):
                for s2 in range(r):
                    a, b = cells.index(k - s3), cells.index(k - s2)
                    v = element(arrays[f'{name}_D'], s3, s2, k)
                    seth[a][b] += v / 2
                    seth[b][a] += v / 2
            if seth != printed:
                fail(f'r={r} k={k}: SETH {name}_D differs from Balsara & Shu IS_{k}')
            # the candidate: printed = derived face coefficients = SETH's P(2, c, k) on cell k-c
            a_bs = [F(x) for x in BS_A[r][k]]
            if a_bs != _zt.face_coefficients(cells):
                fail(f'r={r} k={k}: Table a_k differs from the derived candidate reconstruction')
            a_seth = [element(arrays[f'{name}_P'], 2, k - cell, k) for cell in cells]
            if a_seth != a_bs:
                fail(f'r={r} k={k}: SETH {name}_P differs from Balsara & Shu a_k')
            if element(arrays[f'{name}_A'], 2, k) != F(BS_C[r][k]):
                fail(f'r={r} k={k}: SETH {name}_A differs from Balsara & Shu C_k')
            for cell, a in zip(cells, a_bs):
                row[cell + r - 1] += F(BS_C[r][k]) * a
        if row != [F(x) for x in BS_UL[r]]:
            fail(f'r={r}: Balsara & Shu u^L differs from sum C_k a_k')
        if scalar(tables, f'{name}_D_SCALE') != BS_SCALE[r]:
            fail(f'{name}_D_SCALE is not {BS_SCALE[r]}')
        print(f'  Balsara & Shu r={r}: indicators, Table {"I" if r == 4 else "II"}, C_k, u^L checked against the '
              f'derivation and SETH {name}_*')


def check_teno5_a(teno: str) -> None:
    for name, want in HLF_TENO5_A.items():
        got = scalar(teno, name)
        if got != want:
            fail(f'{name} = {got}, Huang, Liang & Fu 2023 Eqs. (19)-(20) give {want}')
    print('  Huang, Liang & Fu 2023: TENO5-A C_r, xi, alpha_1, alpha_2 checked')


def check_sensor_teno(sensors: str) -> None:
    if scalar(sensors, 'SETH_SENSOR_TENO_CT') != FU_CT:
        fail(f'SETH_SENSOR_TENO_CT differs from Fu 2023 Eq. (23) ({FU_CT})')
    body = re.search(r'subroutine seth_sensor_teno\((.*?)endsubroutine seth_sensor_teno', sensors, re.S)
    power = body and re.search(r'g\(k\) = (\(g\(k\)\*g\(k\)\)(?:\*\(g\(k\)\*g\(k\)\))*)', body.group(1))
    q = power.group(1).count('g(k)') if power else 0
    if q != FU_Q:
        fail(f'seth_sensor_teno raises gamma to the power {q}, Fu 2023 Eq. (21) has q = {FU_Q}')
    print('  Fu 2023: TENO detector q and C_T checked')


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--tables', type=Path, default=Path('src/lib/seth_weno_js_tables.F90'))
    ap.add_argument('--teno', type=Path, default=Path('src/lib/seth_teno.F90'))
    ap.add_argument('--sensors', type=Path, default=Path('src/lib/seth_sensors.F90'))
    args = ap.parse_args()
    check_balsara_shu(args.tables.read_text())
    check_teno5_a(args.teno.read_text())
    check_sensor_teno(args.sensors.read_text())
    if FAILED:
        print(f'check_papers: {len(FAILED)} check(s) failed')
        return 1
    print('check_papers: ok')
    return 0


if __name__ == '__main__':
    sys.exit(main())
