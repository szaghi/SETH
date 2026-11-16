#!/usr/bin/env python3
"""Generate src/lib/seth_z_teno_tables.F90: the descriptors of SETH's WENO-Z and TENO schemes.

Everything is derived in exact rational arithmetic from the definitions listed in notes/m2-sources.md:

- each candidate's face coefficients: the polynomial of degree w-1 matching the cell averages of its w cells,
  evaluated at the face;
- each Jiang-Shu smoothness indicator (Li, Fu & Adams, eq. 10): sum_{j=1}^{w-1} int_{cell 0} (d^j p/dx^j)^2 dx, as a
  quadratic form in the cell averages (unit cells, so the dx powers are absorbed in undivided values);
- the max-order linear weights: the combination of the candidates exact to the highest degree;
- for interpolation (notes/interpolation-plan.md): the same with point values at the cell centres; the indicators and
  tau are the reconstruction ones.

Cells are indexed by their offset from the face (cell 0 = [-1, 0] left of it, cell 1 = [0, 1]), the convention of
SETH's kernels. Cross-checks run before anything is written (the script exits non-zero on any failure):

1. every candidate is exact for polynomials up to degree w-1;
2. the TENO candidates' face coefficients equal Li, Fu & Adams (2020) Table 1;
3. the indicators of the upwind 3-, 4- and 5-point stencils equal ADAM's WENO-5/7/9 tables divided by 1, 240, 5040;
4. the max-order linear weights are exact to degree K-1 and are the WENO weights where K = 2S-1;
5. the order of every tau on smooth data, measured exactly on a generic polynomial (not exp: its averages are a
   constant multiple of its point values), on cell averages and, for the interpolating schemes, on point values;
6. the interpolation tables: the WENO ones equal ADAM's initialize_interpolation_S*, TENO5's weights equal Hiejima
   2022 and Chu, Kurganov & Na 2022, every interpolation weight is positive and reaches the design order.

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


def pv_basis(cells: list[int]) -> list[Poly]:
    """Lagrange basis of point values: basis[j] is 1 at the centre of cells[j] (x = m - 1/2), 0 at the others."""
    w = len(cells)
    a = [[Fr(2 * m - 1, 2) ** n for n in range(w)] for m in cells]
    return [solve(a, [Fr(1) if i == j else Fr(0) for i in range(w)]) for j in range(w)]


def pv_face_coefficients(cells: list[int]) -> list[Fr]:
    """Interpolation: the face value (x = 0) of the polynomial through the point values of the cells."""
    return [peval(b, Fr(0)) for b in pv_basis(cells)]


def sample(m: int, n: int, interp: bool) -> Fr:
    """The datum of cell m for the monomial x^n: its point value at the centre (interpolation) or its average."""
    return Fr(2 * m - 1, 2) ** n if interp else cell_avg_monomial(m, n)


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


def ldl(mat: list[list[Fr]]) -> tuple[list[list[Fr]], list[Fr]]:
    """Exact LDL^T of a symmetric positive semi-definite matrix: mat = L diag(D) L^T, L unit lower-triangular, D >= 0.

    The smoothness indicators are written as sum_i D_i (sum_a L[a][i] u_a)^2: a sum of non-negative terms, so beta >= 0
    in floating point by construction (the upper-triangular form cancels and can round to small negative values,
    which with eps = 1e-40 sends tau/(beta + eps) to -inf: NaN on constant data). A zero pivot (the indicators vanish
    on constants, so one pivot is zero) leaves its column of L zero.
    """
    w = len(mat)
    lo = [[Fr(0)] * w for _ in range(w)]
    d = [Fr(0)] * w
    for i in range(w):
        lo[i][i] = Fr(1)
        d[i] = mat[i][i] - sum(lo[i][j] ** 2 * d[j] for j in range(i))
        if d[i] < 0:
            sys.exit('LDL: negative pivot, the indicator is not positive semi-definite')
        for k in range(i + 1, w):
            num = mat[k][i] - sum(lo[k][j] * lo[i][j] * d[j] for j in range(i))
            if d[i] == 0:
                if num != 0:
                    sys.exit('LDL: zero pivot with a non-zero column')
                lo[k][i] = Fr(0)
            else:
                lo[k][i] = num / d[i]
    for a in range(w):  # exact check: L D L^T == mat
        for c in range(w):
            if sum(lo[a][i] * d[i] * lo[c][i] for i in range(w)) != mat[a][c]:
                sys.exit('LDL: reconstruction differs from the indicator matrix')
    return lo, d


def quad(mat: list[list[Fr]], u: list[Fr]) -> Fr:
    return sum(mat[a][c] * u[a] * u[c] for a in range(len(u)) for c in range(len(u)))


# ------------------------------------------------------------------------------------------------------------- schemes
def weno_cells(S: int) -> list[list[int]]:
    """WENO candidate k: cells k-S+1..k (beta_0 on the most upwind stencil)."""
    return [list(range(k - S + 1, k + 1)) for k in range(S)]


TENO_CELLS = [[-1, 0, 1], [0, 1, 2], [-2, -1, 0], [0, 1, 2, 3]]  # Li, Fu & Adams 2020, Table 1 (i -> 0)
TABLE1 = [[Fr(-1, 6), Fr(5, 6), Fr(2, 6)], [Fr(2, 6), Fr(5, 6), Fr(-1, 6)], [Fr(2, 6), Fr(-7, 6), Fr(11, 6)],
          [Fr(3, 12), Fr(13, 12), Fr(-5, 12), Fr(1, 12)],
          # k = 4, 5 of the eight-point TENO: Fu, Hu & Adams 2016, Table 1
          [Fr(-3, 12), Fr(13, 12), Fr(-23, 12), Fr(25, 12)],
          [Fr(12, 60), Fr(77, 60), Fr(-43, 60), Fr(17, 60), Fr(-3, 60)]]
TENO8_CELLS = TENO_CELLS + [[-3, -2, -1, 0], [0, 1, 2, 3, 4]]  # Fu, Hu & Adams 2016, Table 1, k = 0..5
# TENO6-opt of Fu, Hu & Adams 2016, Table 5 (llm-wiki weno-teno targeted-eno.md, high-order-capturing teno.md)
TENO6_OPT_D = ['0.462_R_P', '0.3_R_P', '0.054_R_P', '0.184_R_P']
# The optimised backgrounds of Fu, Hu & Adams 2016, section 4, Eq. (31): eta g1 + (1 - eta) g2, two sub-schemes of
# order K-1 on the K-point stencil: Eq. (32) for K = 5 (eta = 0.2: TENO5-opt), Eq. (33) for K = 6 (eta = 0.54:
# TENO6-opt). Cells by offset (i -> 0), coefficients as printed.
# The eight-point TENO of Liang & Fu 2024b (Comput. Phys. Commun. 302:109236), section 3.2: its sixth-order background,
# the face coefficients of Eq. (46) and the linear weights printed after it (16 digits, from a numerical optimisation:
# quoted, and checked below against Eq. (46), Fu et al.'s Table 1 candidates and the order)
TENO8_D = ['0.4130855804023061_R_P', '0.2193140179474727_R_P', '0.06459052081827904_R_P', '0.1439236310125986_R_P',
           '0.02746675592227204_R_P', '0.1316194938970745_R_P']
TENO8_EQ46 = [Fr('-0.006866688980568011'), Fr('0.05128582585522106'), Fr('-0.1968478198727312'),
              Fr('0.6552858258552222'), Fr('0.6452858258552218'), Fr('-0.1908478198727314'),
              Fr('0.04928582585522098'), Fr('-0.006580974694853729')]  # f_{i-3} .. f_{i+4}
# tau_8 = |beta_8 - (beta_1 + beta_2 + 4 beta_0)/6| as printed after Eq. (15), over 62270208000: the terms with f_{i+4},
# which only beta_8 contains (f_{i+4} f_{i+4}, f_{i+3} f_{i+4}, f_{i+2} f_{i+4}, f_{i+1} f_{i+4}, f_i f_{i+4},
# f_{i-1} f_{i+4}, f_{i-2} f_{i+4}, f_{i-3} f_{i+4})
TAU8_PRINTED = {(4, 4): 75349098471, (3, 4): -1078504915264, (2, 4): 3263178215782, (1, 4): -5401061230160,
                (0, 4): 5274436892970, (-1, 4): -3038037798592, (-2, 4): 956371298594, (-3, 4): -127080660272}
OPT_SUBSCHEMES = {
    5: ([-2, -1, 0, 1], [Fr(1, 12), Fr(-5, 12), Fr(13, 12), Fr(3, 12)],
        [-1, 0, 1, 2], [Fr(-1, 12), Fr(7, 12), Fr(7, 12), Fr(-1, 12)]),
    6: ([-2, -1, 0, 1, 2], [Fr(2, 60), Fr(-13, 60), Fr(47, 60), Fr(27, 60), Fr(-3, 60)],
        [-1, 0, 1, 2, 3], [Fr(-3, 60), Fr(27, 60), Fr(47, 60), Fr(-13, 60), Fr(2, 60)]),
}
# the optimised weights of the NENO6 paper (arXiv:2207.08500, section 2.2), a later and different optimisation
TENO6_NENO6_D = ['0.4294317718960898_R_P', '0.1727270875843552_R_P', '0.0855682281039113_R_P', '0.312272912415645_R_P']

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
    'TENO5_OPT': (TENO_CELLS[:3], [0, 1, -1], None, 0, 'opt5', 'teno', 5,
                  'TENO5-opt: TENO5 with the optimised linear weights of Fu, Hu & Adams 2016 (eta = 0.2, 4th order)'),
    # Tang 2024, Eq. (19): tau = |beta_0 + beta_2 - 2 beta_1| in WENO labels (beta_0 upwind, beta_1 central, beta_2
    # downwind), i.e. |beta_S2 + beta_S1 - 2 beta_S0| in TENO labels; the cutoff min(C_T, tau) is the kernel's
    'TENO5_AZ': (TENO_CELLS[:3], [-2, 1, 1], None, 0, 'max', 'teno', 4,
                 'TENO5-AZ (Tang 2024): tau = |beta_1 + beta_2 - 2 beta_0| (TENO labels), cutoff min(C_T, tau)'),
    'TENO8': (TENO8_CELLS, [Fr(-4, 6), Fr(-1, 6), Fr(-1, 6), 0, 0, 0], list(range(-3, 5)), 1, TENO8_D, 'teno', 6,
              'TENO8: eight-point TENO with the sixth-order background of Liang & Fu 2024b'),
    'TENO6_OPT': (TENO_CELLS, [Fr(-4, 6), Fr(-1, 6), Fr(-1, 6), 0], list(range(-2, 4)), 1, TENO6_OPT_D, 'teno', 6,
                  'TENO6-opt: TENO6 with the optimised linear weights of Fu, Hu & Adams 2016 (5th order)'),
    'TENO6_NENO6': (TENO_CELLS, [Fr(-4, 6), Fr(-1, 6), Fr(-1, 6), 0], list(range(-2, 4)), 1, TENO6_NENO6_D, 'teno', 6,
                    'TENO6 with the optimised linear weights of the NENO6 paper (4th order)'),
}

# the schemes that also get interpolation tables (owner decision 2026-11-15: not TENO6-opt nor NENO6, whose weights are
# tuned for the reconstruction)
INTERP = ('WENOZ3', 'WENOZ5', 'WENOZ7', 'WENOZ9', 'TENO5', 'TENO6')


def opt_weights(cands: list[list[int]], k: int, eta: Fr) -> list[Fr]:
    """Linear weights of the candidates whose combination is eta g1 + (1 - eta) g2 (Fu, Hu & Adams 2016, Eq. (31)):
    solved on the face coefficients over the union of the cells, and checked on every cell (exit if inconsistent)."""
    c1, g1, c2, g2 = OPT_SUBSCHEMES[k]
    cells = sorted({x for c in cands for x in c} | set(c1) | set(c2))
    target = {m: Fr(0) for m in cells}
    for m, g in zip(c1, g1):
        target[m] += eta * g
    for m, g in zip(c2, g2):
        target[m] += (1 - eta) * g
    a = [[dict(zip(c, face_coefficients(c))).get(m, Fr(0)) for c in cands] for m in cells]
    rows, rhs = [], []
    for m, row in zip(cells, a):
        if rank(rows + [row]) > len(rows):
            rows.append(row)
            rhs.append(target[m])
        if len(rows) == len(cands):
            break
    d = solve(rows, rhs)
    for m, row in zip(cells, a):
        if sum(x * y for x, y in zip(d, row)) != target[m]:
            sys.exit(f'opt weights K={k} eta={eta}: the candidates cannot represent eta g1 + (1 - eta) g2')
    return d


def max_order_weights(cands: list[list[int]], label: str, interp: bool = False) -> tuple[list[Fr], int]:
    """Linear weights exact to the highest degree (of the reconstruction, or of the interpolation with `interp`);
    returns them and that degree + 1 (the order)."""
    nc = len(cands)
    coefs = [(pv_face_coefficients if interp else face_coefficients)(c) for c in cands]

    def moment(k: int, n: int) -> Fr:
        return sum(coefs[k][j] * sample(m, n, interp) for j, m in enumerate(cands[k]))

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


# A generic smooth function for the tau orders: f(x) = sum_n c_n (x + 3/10)^n with rational c_n of no pattern. Not exp:
# the cell averages of exp are a constant multiple of its point values (the two data kinds are indistinguishable), and
# exp satisfies identities such as f' f'''' = f'' f''' that can cancel the leading term of a tau and overstate its order.
GENERIC = [Fr((-1) ** (n * (n + 1) // 2) * (n + 3), (n * n + 7) * math.factorial(n)) for n in range(17)]
GENERIC_SHIFT = Fr(3, 10)


def generic_data(h: Fr, cells: range, interp: bool) -> dict[int, Fr]:
    """Exact data of GENERIC on cells [(m-1)h, mh]: point values at the centres (interp) or cell averages."""
    if interp:
        return {m: sum(c * ((m - Fr(1, 2)) * h + GENERIC_SHIFT) ** n for n, c in enumerate(GENERIC)) for m in cells}

    def prim(x: Fr) -> Fr:
        return sum(c / (n + 1) * (x + GENERIC_SHIFT) ** (n + 1) for n, c in enumerate(GENERIC))

    return {m: (prim(m * h) - prim((m - 1) * h)) / h for m in cells}


def check_fail(msg: str) -> None:
    sys.exit('CHECK FAILED: ' + msg)


def run_checks(adam_source: str) -> None:
    # 1-2: candidates and Table 1
    for cells in [c for s in (2, 3, 4, 5) for c in weno_cells(s)] + TENO_CELLS + [list(range(-2, 4))]:
        coefs = face_coefficients(cells)
        for n in range(len(cells)):
            if sum(c * cell_avg_monomial(m, n) for c, m in zip(coefs, cells)) != (1 if n == 0 else 0):
                check_fail(f'candidate {cells} not exact to degree {n}')
    for k, cells in enumerate(TENO8_CELLS):
        if face_coefficients(cells) != TABLE1[k]:
            check_fail(f'TENO candidate {k} differs from Table 1: {face_coefficients(cells)}')
    # the sub-schemes of Eqs. (32)-(33) are the order K-1 reconstructions on their cells, and the eta combinations give
    # Fu et al.'s optimised weights: Table 3 (eta = 0.2) = Table 5 TENO5-opt, Table 4 (eta = 0.54) = Table 5 TENO6-opt
    for k, (c1, g1, c2, g2) in OPT_SUBSCHEMES.items():
        if face_coefficients(c1) != g1 or face_coefficients(c2) != g2:
            check_fail(f'K={k}: the sub-schemes of Fu et al. Eq. (32)/(33) are not the reconstructions of their cells')
    for k, eta, cands, want in ((5, Fr(1, 5), TENO_CELLS[:3], [Fr(55, 100), Fr(40, 100), Fr(5, 100)]),
                                (6, Fr(54, 100), TENO_CELLS, [Fr(w.replace('_R_P', '')) for w in TENO6_OPT_D])):
        got = opt_weights(cands, k, eta)
        print(f'  K={k}, eta={eta}: optimised weights {[str(x) for x in got]}')
        if got != want:
            check_fail(f'K={k}, eta={eta}: optimised weights {got} differ from Fu et al. Table 5 {want}')
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
    # the eight-point TENO: Fu et al.'s max-order weights (Table 2, K = 8), then Liang & Fu's background (Eq. (46)), its
    # weights (sum, order, face coefficients) and tau_8 (the printed f_{i+4} terms, exact)
    d, order = max_order_weights(TENO8_CELLS, 'TENO8 max-order')
    if d != [Fr(30, 70), Fr(18, 70), Fr(4, 70), Fr(12, 70), Fr(1, 70), Fr(5, 70)] or order != 8:
        check_fail(f'TENO8 max-order weights {d} (order {order}) differ from Fu et al. Table 2, K = 8')
    d8 = [Fr(w.replace('_R_P', '')) for w in TENO8_D]
    if abs(float(sum(d8)) - 1) > 1e-14:  # printed to 16 digits: they sum to 1 + 2.9e-15
        check_fail(f'TENO8 weights sum to {float(sum(d8))!r}')
    face = {m: Fr(0) for m in range(-3, 5)}
    for k, cells in enumerate(TENO8_CELLS):
        for c, m in zip(face_coefficients(cells), cells):
            face[m] += d8[k] * c
    worst = max(abs(float(face[m] - e)) for m, e in zip(range(-3, 5), TENO8_EQ46))
    print(f'  TENO8: face coefficients vs Liang & Fu Eq. (46): max difference {worst:.2e}')
    if worst > 2e-15:
        check_fail(f'TENO8: the weights do not reproduce Eq. (46) ({worst:.2e})')
    gmat = indicator(list(range(-3, 5)))
    for (a, b), num in TAU8_PRINTED.items():
        mine = gmat[a + 3][b + 3] * (1 if a == b else 2)
        if mine != Fr(num, 62270208000):
            check_fail(f'TENO8: beta_8 coefficient of f_{a} f_{b} is {mine}, Liang & Fu print {num}/62270208000')
    print('  TENO8: beta_8 equals the printed tau_8 on the eight f_{i+4} terms (exact)')
    for label, wq, want, cl in (('TENO6-opt', TENO6_OPT_D, 5, TENO_CELLS), ('TENO6-NENO6', TENO6_NENO6_D, 4, TENO_CELLS),
                                ('TENO5-opt', opt_weights(TENO_CELLS[:3], 5, Fr(1, 5)), 4, TENO_CELLS[:3]),
                                ('TENO8 (Liang & Fu)', TENO8_D, 6, TENO8_CELLS)):
        dq = [Fr(w.replace('_R_P', '')) if isinstance(w, str) else w for w in wq]
        coefs = [face_coefficients(c) for c in cl]
        got = 0
        for n in range(8):
            r = sum(dq[k] * sum(coefs[k][j] * cell_avg_monomial(m, n) for j, m in enumerate(cl[k]))
                    for k in range(len(cl))) - (1 if n == 0 else 0)
            if abs(float(r)) > 1e-12:
                break
            got += 1
        print(f'  {label}: linear-scheme order {got} (expected {want})')
        if got != want:
            check_fail(f'{label}: linear-scheme order {got}, expected {want}')
    # 5: tau orders on a generic function, on cell averages (reconstruction) and on point values (interpolation: the
    # same indicators applied to point values, notes/interpolation-plan.md section 2.2)
    for name, (cands, tc, glob, tg, _w, _f, need, _doc) in SCHEMES.items():
        mats = [indicator(c) for c in cands]
        gmat = indicator(glob) if glob else None
        for interp in (False, True) if name in INTERP else (False,):
            taus = []
            for e in range(5, 9):
                dat = generic_data(Fr(1, 2 ** e), range(-4, 6), interp)
                tau = sum(Fr(c) * quad(m, [dat[x] for x in cl]) for c, m, cl in zip(tc, mats, cands))
                if gmat:
                    tau += tg * quad(gmat, [dat[x] for x in glob])
                taus.append(abs(tau))
            order = math.log2(taus[-2] / taus[-1])
            kind = 'point values' if interp else 'cell averages'
            print(f'  {name}: tau order {order:.2f} on {kind} (required {need})')
            if order < need - 0.1:
                check_fail(f'{name}: tau of order {order:.2f} < {need} on {kind}')
    # 6: interpolation tables. The point-value face coefficients and max-order weights of the WENO candidates equal
    # ADAM's interpolation tables (row 2); TENO5's weights are those printed by Hiejima 2022 (WCNS-T, K = 5) and by Chu,
    # Kurganov & Na 2022; every interpolation weight is positive (max_order_weights) and reaches the design order
    # (TENO6's: derived here, no published set)
    for S in range(2, 6):
        t = _js.parse(_js.routine(adam_source, f'initialize_interpolation_S{S}'))
        d, order = max_order_weights(weno_cells(S), f'WENO S={S} interpolation', interp=True)
        if d != [_js.rational(t['a'][2, k]) for k in range(S)] or order != 2 * S - 1:
            check_fail(f'WENO S={S}: interpolation weights {d} (order {order}) differ from ADAM')
        for k, cells in enumerate(weno_cells(S)):
            pv = pv_face_coefficients(cells)
            if any(pv[cells.index(k - s2)] != _js.rational(t['p'][2, s2, k]) for s2 in range(S)):
                check_fail(f'WENO S={S} candidate {k}: interpolation coefficients differ from ADAM')
    d, order = max_order_weights(TENO_CELLS[:3], 'TENO5 interpolation', interp=True)
    if d != [Fr(5, 8), Fr(5, 16), Fr(1, 16)] or order != 5:
        check_fail(f'TENO5 interpolation weights {d} order {order}')
    d, order = max_order_weights(TENO_CELLS, 'TENO6 interpolation', interp=True)
    print(f'  TENO6 interpolation weights {[str(x) for x in d]}, order {order} (derived)')
    if order != 6:
        check_fail(f'TENO6 interpolation weights reach order {order}')


# ---------------------------------------------------------------------------------------------------------- emission
def fstr(x) -> str:
    if isinstance(x, str):
        return x
    x = Fr(x)
    if abs(x.numerator) >= 2 ** 53 or x.denominator >= 2 ** 53:
        # not a ratio of exactly representable integers (the LDL factors of the 6-point indicator): the correctly
        # rounded double, 17 significant digits (round-trips exactly)
        return f'{float(x):.16e}_R_P'
    return f'{x.numerator}._R_P' if x.denominator == 1 else f'{x.numerator}._R_P/{x.denominator}._R_P'


def farray(name: str, bounds: str, extents: str, values: list[str], doc: str, kind: str = 'real(R_P)') -> str:
    rows, row = [], []
    for v in values:  # pack the values by length: free-form lines stay within 132 columns
        if row and len('      ' + ', '.join(row + [v]) + ', &') > 120:
            rows.append(', '.join(row))
            row = []
        row.append(v)
    rows.append(', '.join(row))
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
    if wts == 'max':
        d = max_order_weights(cands, name)[0]
    elif wts == 'opt5':
        d = opt_weights(cands, 5, Fr(1, 5))
    else:
        d = wts
    out = [f'! {doc}\n']
    names = []

    def scalar(n: str, kind: str, val: str, doc_: str) -> None:
        out.append(f'{kind}, parameter :: {name}_{n} = {val} !< {doc_}\n')
        names.append(f'{name}_{n}')

    scalar('NC', 'integer(I4P)', f'{nc}_I4P', 'candidate stencils')
    scalar('LO', 'integer(I4P)', f'{lo}_I4P', 'first cell of the stencil arrays (bounds LO:1-LO)')
    scalar('UMAX', 'integer(I4P)', f'{umax}_I4P', 'last cell any candidate reads (left bias)')
    scalar('MAXW', 'integer(I4P)', f'{maxw}_I4P', 'widest candidate')
    scalar('TG', 'real(R_P)', fstr(tg), 'coefficient of the global indicator in tau (0: none)')
    out.append(f'integer(I4P), parameter :: {name}_CLO(0:{nc - 1}) = [{", ".join(f"{min(c)}_I4P" for c in cands)}]'
               ' !< first cell of each candidate\n')
    out.append(f'integer(I4P), parameter :: {name}_CW(0:{nc - 1}) = [{", ".join(f"{len(c)}_I4P" for c in cands)}]'
               ' !< width of each candidate\n')
    names += [f'{name}_CLO', f'{name}_CW']
    pvals = []
    lvals = []
    dvals = []
    for k in range(nc):  # Fortran order: first index fastest
        coefs = face_coefficients(cands[k]) + [Fr(0)] * (maxw - len(cands[k]))
        pvals += [fstr(x) for x in coefs]
    for k in range(nc):
        lt, dd = ldl(indicator(cands[k]))
        w = len(cands[k])
        full = [[lt[a][i] if a < w and i < w else Fr(0) for i in range(maxw)] for a in range(maxw)]
        lvals += [fstr(full[a][i]) for i in range(maxw) for a in range(maxw)]
        dvals += [fstr(dd[i]) if i < w else fstr(0) for i in range(maxw)]
    out.append(farray(f'{name}_P', f'0:{maxw - 1},0:{nc - 1}', f'{maxw},{nc}', pvals,
                      'face coefficient of cell CLO(k)+j'))
    out.append(farray(f'{name}_L', f'0:{maxw - 1},0:{maxw - 1},0:{nc - 1}', f'{maxw},{maxw},{nc}', lvals,
                      'indicator factors: beta_k = sum_i DQ(i,k) (sum_a L(a,i,k) u_a)^2'))
    out.append(farray(f'{name}_DQ', f'0:{maxw - 1},0:{nc - 1}', f'{maxw},{nc}', dvals,
                      'indicator pivots DQ(i,k) >= 0'))
    out.append(farray(f'{name}_D', f'0:{nc - 1}', f'{nc}', [fstr(x) for x in d], 'linear weights'))
    out.append(farray(f'{name}_TC', f'0:{nc - 1}', f'{nc}', [fstr(x) for x in tc], 'tau coefficients on the beta_k'))
    names += [f'{name}_P', f'{name}_L', f'{name}_DQ', f'{name}_D', f'{name}_TC']
    g = glob if glob else [lo]
    gl, gd = ldl(indicator(glob)) if glob else ([[Fr(0)]], [Fr(0)])
    gw = len(g)
    scalar('GLO', 'integer(I4P)', f'{g[0]}_I4P', 'first cell of the global stencil')
    scalar('GW', 'integer(I4P)', f'{gw}_I4P', 'width of the global stencil')
    out.append(farray(f'{name}_LG', f'0:{gw - 1},0:{gw - 1}', f'{gw},{gw}',
                      [fstr(gl[a][i]) for i in range(gw) for a in range(gw)], 'global indicator, sum of squares: L'))
    out.append(farray(f'{name}_DG', f'0:{gw - 1}', f'{gw}', [fstr(x) for x in gd], 'global indicator pivots'))
    names += [f'{name}_LG', f'{name}_DG']
    if name in INTERP:  # interpolation: other face coefficients and weights, the same indicators and tau
        ivals = []
        for k in range(nc):
            ivals += [fstr(x) for x in pv_face_coefficients(cands[k]) + [Fr(0)] * (maxw - len(cands[k]))]
        di = max_order_weights(cands, name + ' interpolation', interp=True)[0]
        out.append(farray(f'{name}_INTERP_P', f'0:{maxw - 1},0:{nc - 1}', f'{maxw},{nc}', ivals,
                          'interpolation: face coefficient of cell CLO(k)+j'))
        out.append(farray(f'{name}_INTERP_D', f'0:{nc - 1}', f'{nc}', [fstr(x) for x in di],
                          'interpolation linear weights'))
        names += [f'{name}_INTERP_P', f'{name}_INTERP_D']
    return ''.join(out), names


# The linear backgrounds as single face-coefficient rows (cells LO..UMAX of the left-biased stencil): the smooth branch
# of a hybrid. name: (scheme, interp)
LINEARS = {'LIN3': ('WENOZ3', False), 'LIN5': ('WENOZ5', False), 'LIN7': ('WENOZ7', False), 'LIN9': ('WENOZ9', False),
           'LIN3_INTERP': ('WENOZ3', True), 'LIN5_INTERP': ('WENOZ5', True), 'LIN7_INTERP': ('WENOZ7', True),
           'LIN9_INTERP': ('WENOZ9', True), 'LIN6': ('TENO6', False), 'LIN6_INTERP': ('TENO6', True),
           'LIN8': ('TENO8', False)}


def linear_row(scheme: str, interp: bool) -> tuple[int, int, list]:
    """sum_k d_k q_k collapsed to one row c(lo..hi); exact unless the weights are quoted decimals (TENO8: float)."""
    cands, _tc, _g, _tg, wts, _f, _n, _doc = SCHEMES[scheme]
    if wts == 'max':
        d = max_order_weights(cands, scheme, interp=interp)[0]
    else:
        d = [Fr(w.replace('_R_P', '')) for w in wts]
    d = [x / sum(d) for x in d]  # what the TENO kernel computes, d_k / sum d (TENO8's printed weights sum to 1 + 2.9e-15)
    lo, hi = min(min(c) for c in cands), max(max(c) for c in cands)
    row = {m: Fr(0) for m in range(lo, hi + 1)}
    for k, c in enumerate(cands):
        for coef, m in zip((pv_face_coefficients if interp else face_coefficients)(c), c):
            row[m] += d[k] * coef
    vals = [row[m] for m in range(lo, hi + 1)]
    if sum(vals) != 1 and abs(float(sum(vals)) - 1) > 1e-14:
        check_fail(f'{scheme} linear row does not sum to 1')
    return lo, hi, vals


def emit_linear(name: str) -> tuple[str, list[str]]:
    scheme, interp = LINEARS[name]
    lo, hi, vals = linear_row(scheme, interp)
    kind = 'interpolation from point values' if interp else 'reconstruction from cell averages'
    exact = all(isinstance(v, Fr) and v.denominator < 2 ** 53 for v in vals) and SCHEMES[scheme][4] == 'max'
    txt = [f'! {name}: the linear background of {scheme}, {kind}; cells {lo}..{hi}\n',
           f'integer(I4P), parameter :: {name}_LO = {lo}_I4P !< first cell\n',
           f'integer(I4P), parameter :: {name}_HI = {hi}_I4P !< last cell\n',
           farray(f'{name}_C', f'{lo}:{hi}', f'{hi - lo + 1}',
                  [fstr(v) if exact else f'{float(v):.16e}_R_P' for v in vals], 'face coefficient of each cell')]
    return ''.join(txt), [f'{name}_LO', f'{name}_HI', f'{name}_C']


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
    for n in LINEARS:
        b, ns = emit_linear(n)
        blocks.append(b)
        names += ns
    header = '''!< SETH, WENO-Z and TENO scheme descriptors (generated: do not edit).
module seth_z_teno_tables
!< SETH, WENO-Z and TENO scheme descriptors (generated: do not edit).
!<
!< Generated by scripts/gen_z_teno_tables.py, every coefficient derived in exact rational arithmetic from the
!< definitions in notes/m2-sources.md and cross-checked (candidates exact, TENO Table 1, ADAM's WENO-5/7/9 indicators,
!< max-order weights, tau orders, ADAM's interpolation tables). Cells by offset from the face i+1/2 (cell 0 left of
!< it); left-biased stencils. The _INTERP_ arrays interpolate from the point values at the cell centres.
use penf, only : I4P, R_P
implicit none
private
'''
    publics = ''.join(f'public :: {n}\n' for n in names)
    out_path.write_text(header + publics + '\n' + '\n'.join(blocks) + 'endmodule seth_z_teno_tables\n')
    print(f'{out_path}: {len(SCHEMES)} schemes written')


if __name__ == '__main__':
    main()
