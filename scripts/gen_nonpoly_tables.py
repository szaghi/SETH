#!/usr/bin/env python3
"""Generate src/lib/seth_nonpoly_tables.F90: the constants of the non-polynomial reconstructions (issue #6, step 0.4).

Every constant is derived here and checked against the printed formula; the script aborts on any mismatch, and writes
the Fortran module only when every check passes. Sources: notes/m5-sources.md (Takagi et al. 2022; Liang, Xiao, Shyy &
Fu 2022; Liang & Fu 2024).

(a) Liang & Fu 2024, Eq. 25: the six Gaussian-RBF flux rows. R'(x) is the derivative of the Gaussian-RBF interpolant
    of the primitive H through the faces of the stencil (Eqs. 16-24). The rows are computed in 400-digit arithmetic at
    s = lambda dx = 1e-30 and 2e-30 and split into P + C s^2 + O(s^4); P and C are identified as rationals (residual
    below 1e-25) and compared with the printed rows, with the paper's labels R'_2 and R'_3 swapped against the
    stencils of its section 3.1 (notes/m5-sources.md). P is also checked against the polynomial reconstruction on the
    same cells, exactly.
(b) Liang & Fu 2024, Eq. 29: the leading truncation errors of the six rows, exactly (rational Taylor moments).
(c) The pair weights: on i-2..i+1 the polynomial row is 3 q0/4 + q2/4, on i-1..i+2 it is q0/2 + q1/2 (q_k TENO5's
    candidates), exactly.
(d) Liang et al. 2022, Eq. 31 equals m + h tanh(A + atanh B) (the tanh addition formula), checked at 60 digits on
    random monotone triples; the literals tanh(0.91), tanh(1.23).
(e) Takagi et al. 2022, Eq. 16 with A = tanh(alpha beta/2)/tanh(beta/2): the face value of the THINC function of
    Eq. 15 under the conservation constraint, checked at 60 digits against a root-find of d and the integral.
(f) The untruncated Gaussian rows against the truncated ones at s up to 0.3: what the truncation of Eq. 25 drops
    (printed, not asserted: the paper does not say which form its code used).

The rows are emitted by mask m = delta_0 + 2 delta_1 + 4 delta_2 (seth_teno5_mask), cells by offset from the cell i
(-2..2): m = 7 the full stencil, 5 (1,0,1) the cells -2..1, 3 (1,1,0) the cells -1..2, 1, 2, 4 the candidates q0, q1,
q2; m = 6 (0,1,1), a mask no row of the paper covers, takes TENO5's renormalised weights 3 q1/4 + q2/4 (owner decision
D4 of issue #6), with the same combination of the RBF rows of q1 and q2 (SETH's choice); m = 0 never occurs (one
chi_k >= 1/3 > C_T) and repeats m = 7.

Usage: gen_nonpoly_tables.py [--out src/lib/seth_nonpoly_tables.F90] [--check]   (needs mpmath)
"""

from __future__ import annotations

import argparse
import importlib.util
import random
import sys
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location('gen_z_teno_tables', HERE / 'gen_z_teno_tables.py')
_zt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_zt)

FAILS: list[str] = []


def check(ok: bool, msg: str) -> None:
    print(('  ok    ' if ok else '  FAIL  ') + msg)
    if not ok:
        FAILS.append(msg)


# ------------------------------------------------------------------------------------- (a) Gaussian-RBF flux rows
# stencils by SETH's mask (offsets from cell i); the paper's row label for each (Eq. 25)
STENCILS = {7: (-2, 2), 5: (-2, 1), 3: (-1, 2), 1: (-1, 1), 2: (0, 2), 4: (-2, 0)}
PAPER_LABEL = {7: 1, 5: 3, 3: 2, 1: 4, 2: 5, 4: 6}  # R'_2 sits on -1..2, R'_3 on -2..1 (the label swap)
F = Fr
# Eq. 25 as printed: label -> {offset: (polynomial coefficient, s^2 coefficient)}
EQ25 = {
    1: {2: (F(-3, 60), F(-18, 60)), 1: (F(27, 60), F(42, 60)), 0: (F(47, 60), F(-18, 60)),
        -1: (F(-13, 60), F(-18, 60)), -2: (F(2, 60), F(12, 60))},
    2: {2: (F(-1, 12), F(-4, 12)), 1: (F(7, 12), F(4, 12)), 0: (F(7, 12), F(4, 12)), -1: (F(-1, 12), F(-4, 12))},
    3: {1: (F(3, 12), F(9, 12)), 0: (F(13, 12), F(-15, 12)), -1: (F(-5, 12), F(3, 12)), -2: (F(1, 12), F(3, 12))},
    4: {1: (F(2, 6), F(4, 6)), 0: (F(5, 6), F(-2, 6)), -1: (F(-1, 6), F(-2, 6))},
    5: {2: (F(-1, 6), F(-2, 6)), 1: (F(5, 6), F(-2, 6)), 0: (F(2, 6), F(4, 6))},
    6: {0: (F(11, 6), F(-18, 6)), -1: (F(-7, 6), F(18, 6)), -2: (F(2, 6), F(0))},
}


def rbf_row(lo: int, hi: int, s) -> dict[int, mp.mpf]:
    """Weights on f_lo..f_hi of R'(x_{i+1/2}), the derivative of the Gaussian-RBF interpolant of the primitive H
    through the faces x_{lo-1/2}..x_{hi+1/2} (dx = 1, cell i centred at 0, H(x_{l+1/2}) = sum_{j <= l} f_j)."""
    nodes = [mp.mpf(l) + mp.mpf(1) / 2 for l in range(lo - 1, hi + 1)]
    n = len(nodes)
    a = mp.matrix(n, n)
    for p in range(n):
        for q in range(n):
            a[p, q] = mp.exp(-(s ** 2) * (nodes[p] - nodes[q]) ** 2)
    xf = mp.mpf(1) / 2
    dphi = mp.matrix([-2 * s ** 2 * (xf - nodes[q]) * mp.exp(-(s ** 2) * (xf - nodes[q]) ** 2) for q in range(n)])
    y = mp.lu_solve(a.T, dphi)
    return {j: sum(y[p] for p in range(n) if lo - 1 + p >= j) for j in range(lo, hi + 1)}


def rational(x: mp.mpf, what: str) -> Fr:
    r = Fr(mp.nstr(x, 60)).limit_denominator(100000)
    check(abs(x - mp.mpf(r.numerator) / r.denominator) < mp.mpf('1e-25'), f'{what} = {r} (identified as a rational)')
    return r


def rbf_rows() -> dict[int, dict[int, tuple[Fr, Fr]]]:
    mp.mp.dps = 400
    s1, s2 = mp.mpf('1e-30'), mp.mpf('2e-30')
    rows = {}
    print('(a) Gaussian-RBF rows, P + C s^2, against Liang & Fu 2024, Eq. 25')
    for m, (lo, hi) in STENCILS.items():
        w1, w2 = rbf_row(lo, hi, s1), rbf_row(lo, hi, s2)
        row = {}
        for j in range(lo, hi + 1):
            c = (w2[j] - w1[j]) / (s2 ** 2 - s1 ** 2)
            p = w1[j] - c * s1 ** 2
            row[j] = (rational(p, f'mask {m} P[{j:+d}]'), rational(c, f'mask {m} C[{j:+d}]'))
        printed = EQ25[PAPER_LABEL[m]]
        same = all(row[j] == printed.get(j, (Fr(0), Fr(0))) for j in range(lo, hi + 1)) and set(printed) <= set(row)
        check(same, f"mask {m}, cells {lo}..{hi}: equal to the paper's R'_{PAPER_LABEL[m]}")
        poly = _zt.face_coefficients(list(range(lo, hi + 1)))
        check([row[j][0] for j in range(lo, hi + 1)] == poly, f'mask {m}: P is the polynomial reconstruction')
        check(sum(row[j][1] for j in row) == 0, f'mask {m}: the s^2 row sums to zero')
        rows[m] = row
    return rows


# --------------------------------------------------------------------------------- (b) truncation errors, Eq. 29
def truncation(rows: dict[int, dict[int, tuple[Fr, Fr]]]) -> None:
    """R - h(face) = sum_n c_n dx^n [sum_j P_j a_jn] + lambda^2 sum_n c_n dx^(n+2) [sum_j C_j a_jn], with h = sum c_n
    x^n/n! at the face and a_jn = ((j+1)^(n+1) - j^(n+1))/(n+1)! the average over cell j = [j-1/2 .. j+1/2] shifted so
    that the face is at 0: cell j spans [j-1, j] in face coordinates (cell i = j = 0 spans [-1, 0])."""
    print('(b) truncation errors against Liang & Fu 2024, Eq. 29 (coefficients of dx^k h^(n), lambda^0 and lambda^2)')

    def a(j: int, n: int) -> Fr:
        fact = 1
        for k in range(2, n + 2):
            fact *= k
        return Fr(j ** (n + 1) - (j - 1) ** (n + 1), fact)

    def err(row: dict[int, tuple[Fr, Fr]]) -> dict[tuple[int, int, int], Fr]:
        out: dict[tuple[int, int, int], Fr] = {}  # (power of lambda, power of dx, derivative n)
        for n in range(9):
            e0 = sum(p * a(j, n) for j, (p, _) in row.items()) - (1 if n == 0 else 0)
            e2 = sum(c * a(j, n) for j, (_, c) in row.items())
            if e0:
                out[(0, n, n)] = e0
            if e2:
                out[(2, n + 2, n)] = e2
        return out

    # Eq. 29 as printed, the leading terms: (lambda power, dx power, derivative) -> coefficient
    eq29 = {
        7: {(2, 5, 3): F(-30, 60), (0, 5, 5): F(-1, 60), (2, 6, 4): F(12, 60)},
        3: {(2, 4, 2): F(-20, 30), (0, 4, 4): F(-1, 30)},   # the paper's R'_2
        5: {(2, 4, 2): F(120, 120), (0, 4, 4): F(6, 120)},  # the paper's R'_3
        1: {(2, 3, 1): F(360, 360), (0, 3, 3): F(30, 360)},
        2: {(2, 3, 1): F(-360, 360), (0, 3, 3): F(-30, 360)},
        4: {(2, 3, 1): F(-360, 120), (0, 3, 3): F(-30, 120)},
    }
    for m, row in rows.items():
        e = err(row)
        lead = min(k[1] for k in e)
        ok = all(e.get(k, Fr(0)) == v for k, v in eq29[m].items())
        ok = ok and all(k[1] >= lead for k in e) and all(k in eq29[m] for k in e if k[1] == lead)
        terms = ', '.join(f'{v} lambda^{k[0]} dx^{k[1]} h^({k[2]})' for k, v in sorted(e.items(), key=lambda t: t[0][1])
                          if k[1] <= lead + 1)
        check(ok, f'mask {m}: {terms}')


# ------------------------------------------------------------------------------------------- (c) the pair weights
def pairs(rows: dict[int, dict[int, tuple[Fr, Fr]]]) -> None:
    print('(c) the pair weights of stencil targeting at lambda = 0')
    q = {k: dict(zip(range(lo, lo + 3), _zt.face_coefficients(list(range(lo, lo + 3)))))
         for k, lo in ((0, -1), (1, 0), (2, -2))}

    def comb(*wq: tuple[Fr, int]) -> dict[int, Fr]:
        out: dict[int, Fr] = {}
        for w, k in wq:
            for j, c in q[k].items():
                out[j] = out.get(j, Fr(0)) + w * c
        return out

    for m, wq, txt in ((5, ((F(3, 4), 0), (F(1, 4), 2)), '3 q0/4 + q2/4'), (3, ((F(1, 2), 0), (F(1, 2), 1)),
                                                                               'q0/2 + q1/2')):
        want = comb(*wq)
        got = {j: p for j, (p, _) in rows[m].items()}
        check(all(got.get(j, 0) == want.get(j, 0) for j in set(got) | set(want)), f'mask {m}: {txt}')


# ------------------------------------------------------------------------------------------------ (d) DR, Eq. 31
def dr() -> tuple[str, str]:
    print('(d) Liang et al. 2022, Eq. 31 = m + h tanh(A + atanh B); the literals')
    mp.mp.dps = 60
    rng = random.Random(20261117)
    worst = mp.mpf(0)
    for _ in range(2000):
        f = sorted(mp.mpf(rng.uniform(-3, 3)) for _ in range(3))
        if rng.random() < 0.5:
            f = f[::-1]
        fm1, fi, fp1 = f
        if (fp1 - fi) * (fi - fm1) <= mp.mpf('1e-15'):
            continue
        for a in (mp.mpf('0.91'), mp.mpf('1.23')):
            b = (2 * fi - (fp1 + fm1)) / (fp1 - fm1)
            eq31 = (fp1 + fm1) / 2 + (fp1 - fm1) / 2 * (mp.tanh(a) + b) / (1 + mp.tanh(a) * b)
            alt = (fp1 + fm1) / 2 + (fp1 - fm1) / 2 * mp.tanh(a + mp.atanh(b))
            worst = max(worst, abs(eq31 - alt))
    check(worst < mp.mpf('1e-50'), f'Eq. 31 = m + h tanh(A + atanh B), worst {mp.nstr(worst, 3)}')
    t1, t2 = mp.tanh(mp.mpf('0.91')), mp.tanh(mp.mpf('1.23'))
    lit = (f'{float(t1):.17g}', f'{float(t2):.17g}')
    check(abs(float(lit[0]) - float(t1)) == 0 and abs(float(lit[1]) - float(t2)) == 0,
          f'tanh(0.91) = {lit[0]}, tanh(1.23) = {lit[1]} (the nearest doubles)')
    return lit


# --------------------------------------------------------------------------------------------- (e) THINC, Eq. 16
def thinc() -> None:
    print('(e) Takagi et al. 2022, Eq. 16, A = tanh(alpha beta/2)/tanh(beta/2), against the THINC function of Eq. 15')
    mp.mp.dps = 60
    rng = random.Random(20261118)
    beta = mp.mpf('1.6')
    worst = mp.mpf(0)
    for _ in range(200):
        f = sorted(mp.mpf(rng.uniform(-3, 3)) for _ in range(3))
        if rng.random() < 0.5:
            f = f[::-1]
        fm1, fi, fp1 = f
        if (fp1 - fi) * (fi - fm1) <= mp.mpf('1e-15'):
            continue
        fp, fn = (fm1 + fp1) / 2, (fp1 - fm1) / 2
        alpha = (fi - fp) / fn

        def avg(d):  # the cell average of h in the unit coordinate X in [-1/2, 1/2]
            return fp + fn * (mp.log(mp.cosh(beta * (mp.mpf(1) / 2 - d))) - mp.log(mp.cosh(beta * (-mp.mpf(1) / 2 - d)))) / beta
        d = mp.findroot(lambda d: avg(d) - fi, 0)
        exact = fp + fn * mp.tanh(beta * (mp.mpf(1) / 2 - d))
        t = mp.tanh(beta / 2)
        a = mp.tanh(alpha * beta / 2) / t
        eq16 = fp + fn * (a + t) / (1 + a * t)
        b = mp.exp(abs(alpha) * beta)  # the printed form, sgn(alpha)/T (B-1)/(B+1)
        a_printed = mp.sign(alpha) / t * (b - 1) / (b + 1)
        worst = max(worst, abs(eq16 - exact), abs(a - a_printed))
    check(worst < mp.mpf('1e-40'), f'Eq. 16 exact under the conservation constraint, worst {mp.nstr(worst, 3)}')


# --------------------------------------------------------------------------------- (f) what the truncation drops
def untruncated(rows: dict[int, dict[int, tuple[Fr, Fr]]]) -> None:
    print('(f) the untruncated Gaussian rows against P + C s^2 (max over the row, report only)')
    mp.mp.dps = 60
    for s in ('0.03', '0.06', '0.12', '0.15', '0.3'):
        sv = mp.mpf(s)
        out = []
        for m, (lo, hi) in STENCILS.items():
            w = rbf_row(lo, hi, sv)
            d = max(abs(w[j] - (mp.mpf(rows[m][j][0].numerator) / rows[m][j][0].denominator
                                + sv ** 2 * mp.mpf(rows[m][j][1].numerator) / rows[m][j][1].denominator))
                    for j in range(lo, hi + 1))
            out.append(f'm{m} {mp.nstr(d, 2)}')
        print(f'  s = {s:5s}: ' + ', '.join(out))


# ------------------------------------------------------------------------------------------------------- emission
def emit(rows: dict[int, dict[int, tuple[Fr, Fr]]], t: tuple[str, str]) -> str:
    full = dict(rows)
    full[6] = {j: (F(3, 4) * rows[2].get(j, (0, 0))[0] + F(1, 4) * rows[4].get(j, (0, 0))[0],
                   F(3, 4) * rows[2].get(j, (0, 0))[1] + F(1, 4) * rows[4].get(j, (0, 0))[1]) for j in range(-2, 3)}
    full[0] = rows[7]
    pv, cv = [], []
    for m in range(8):
        for j in range(-2, 3):
            p, c = full[m].get(j, (Fr(0), Fr(0)))
            pv.append(_zt.fstr(p))
            cv.append(_zt.fstr(c))
    head = '''!< SETH, constants of the non-polynomial reconstructions (generated: do not edit).
module seth_nonpoly_tables
!< SETH, constants of the non-polynomial reconstructions (generated: do not edit).
!<
!< Generated by scripts/gen_nonpoly_tables.py, which derives every constant and checks it against the printed formula
!< (notes/m5-sources.md): the Gaussian-RBF flux rows of Liang & Fu 2024 (Eq. 25), P + s^2 C, s = lambda dx, by mask
!< m = delta_0 + 2 delta_1 + 4 delta_2 of the cell (seth_teno5_mask), cells by offset from the cell i (-2..2): m = 7
!< the five cells, 5 the cells -2..1, 3 the cells -1..2, 1, 2, 4 TENO5's candidates q0, q1, q2, 6 TENO5's renormalised
!< weights 3 q1/4 + q2/4 (issue #6, D4), 0 never produced (= 7). P is the polynomial reconstruction on the same cells
!< (stencil targeting), every C row sums to zero. The tanh of DR's two steepness constants (Liang, Xiao, Shyy & Fu
!< 2022, Eq. 31) as the nearest doubles.
use penf, only : R_P
implicit none
private
public :: NP_RBF_P, NP_RBF_C, SETH_DR_T1, SETH_DR_T2

'''
    body = _zt.farray('NP_RBF_P', '-2:2,0:7', '5,8', pv, 'polynomial rows by mask (stencil targeting)')
    body += _zt.farray('NP_RBF_C', '-2:2,0:7', '5,8', cv, 'the s^2 rows by mask (Gaussian RBF, Eq. 25)')
    body += f'real(R_P), parameter :: SETH_DR_T1 = {t[0]}_R_P !< tanh(0.91): DR steepness, level 1 (Eq. 31).\n'
    body += f'real(R_P), parameter :: SETH_DR_T2 = {t[1]}_R_P !< tanh(1.23): DR steepness, level 2 (Eqs. 31-32).\n'
    return head + body + 'endmodule seth_nonpoly_tables\n'


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--out', type=Path, default=HERE.parent / 'src/lib/seth_nonpoly_tables.F90')
    ap.add_argument('--check', action='store_true', help='compare with the file on disk instead of writing it')
    args = ap.parse_args()
    rows = rbf_rows()
    truncation(rows)
    pairs(rows)
    t = dr()
    thinc()
    untruncated(rows)
    if FAILS:
        print(f'{len(FAILS)} check(s) failed: nothing written')
        return 1
    text = emit(rows, t)
    if args.check:
        same = args.out.read_text() == text
        print(f'{args.out}: ' + ('up to date' if same else 'DIFFERS from the generator'))
        return 0 if same else 1
    args.out.write_text(text)
    print(f'written {args.out}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
