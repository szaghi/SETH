#!/usr/bin/env python3
"""Evidence for notes/interpolation-plan.md: exact interpolation weights and tau orders on point values.

Usage (repository root): python3 -I notes/interpolation-plan-evidence.py <path to ADAM adam_weno_object.F90>
"""
import sys, math, importlib.util
from fractions import Fraction as Fr
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location('g', 'scripts/gen_z_teno_tables.py'); g = importlib.util.module_from_spec(spec); spec.loader.exec_module(g)
js = g._js
src = open(sys.argv[1]).read()
# point value of cell m sits at its centre x = m - 1/2; face at x = 0
def pv_basis(cells):
    w = len(cells); a = [[Fr(2*m-1, 2)**n for n in range(w)] for m in cells]
    return [g.solve(a, [Fr(int(i == j)) for i in range(w)]) for j in range(w)]
def pv_face(cells): return [g.peval(b, Fr(0)) for b in pv_basis(cells)]
def pv_indicator(cells):
    bs = pv_basis(cells); w = len(cells); mat = [[Fr(0)]*w for _ in range(w)]
    for j in range(1, w):
        ders = []
        for b in bs:
            d = b
            for _ in range(j): d = g.pder(d)
            ders.append(d)
        for a in range(w):
            for c in range(w): mat[a][c] += g.pint(g.pmul(ders[a], ders[c]), Fr(-1), Fr(0))
    return mat
def pv_weights(cands):
    nc = len(cands); co = [pv_face(c) for c in cands]
    mom = lambda k, n: sum(co[k][j]*Fr(2*m-1, 2)**n for j, m in enumerate(cands[k]))
    rows, rhs, n = [], [], 0
    while len(rows) < nc:
        row = [mom(k, n) for k in range(nc)]
        if g.rank(rows+[row]) == len(rows)+1: rows.append(row); rhs.append(Fr(int(n == 0)))
        n += 1
    d = g.solve(rows, rhs); o = 0
    while sum(d[k]*mom(k, o) for k in range(nc)) == (1 if o == 0 else 0): o += 1
    return d, o
for S in range(2, 6):
    d, o = pv_weights(g.weno_cells(S))
    a = js.parse(js.routine(src, f'initialize_interpolation_S{S}'))['a']
    adam = [js.rational(a[2, k]) for k in range(S)]
    print(f'WENO S={S}: interp d {[str(x) for x in d]} order {o}; ADAM interp a(2,:) equal: {d == adam}')
    p = js.parse(js.routine(src, f'initialize_interpolation_S{S}'))['p']
    ok = all(pv_face(c)[c.index(k - s2)] == js.rational(p[2, s2, k]) for k, c in enumerate(g.weno_cells(S)) for s2 in range(S))
    print(f'   face coefficients equal ADAM interp p(2,...): {ok}')
for lab, c in (('TENO5', g.TENO_CELLS[:3]), ('TENO6', g.TENO_CELLS)):
    d, o = pv_weights(c); print(f'{lab}: interp d {[str(x) for x in d]} = {[float(x) for x in d]} order {o}')
# indicators: reconstruction vs interpolant-derivative, per width
for cells in [[-1, 0], [-2, -1, 0], [-3, -2, -1, 0], [-4, -3, -2, -1, 0], [0, 1, 2, 3], list(range(-2, 4))]:
    print(f'indicator {cells}: interp-integral == reconstruction: {pv_indicator(cells) == g.indicator(cells)}')
# tau orders on point values exp(x) at x = (m - 1/2) h
def pvals(h): return {m: sum(Fr(1, math.factorial(n))*((m - Fr(1, 2))*h)**n for n in range(17)) for m in range(-4, 6)}
for name, (cands, tc, glob, tg, *_r) in g.SCHEMES.items():
    if name in ('TENO6_OPT', 'TENO6_NENO6'): continue
    for lab, ind in (('recon-beta', g.indicator), ('interp-beta', pv_indicator)):
        mats = [ind(c) for c in cands]; gm = ind(glob) if glob else None; taus = []
        for e in range(4, 8):
            v = pvals(Fr(1, 2**e))
            t = sum(Fr(c)*g.quad(m, [v[x] for x in cl]) for c, m, cl in zip(tc, mats, cands))
            if gm: t += tg*g.quad(gm, [v[x] for x in glob])
            taus.append(abs(t))
        print(f'  {name:8s} {lab}: tau order on point values {math.log2(taus[-2]/taus[-1]):.2f}')
print('=== generic function (no exp degeneracy): f = sum c_n x^n, c_n = s_n (n+3)/((n*n+7) n!), x0 shift 0.3')
coef = [Fr((-1)**(n*(n+1)//2)*(n+3), (n*n+7)*math.factorial(n)) for n in range(17)]
def fpt(x): return sum(c*(x+Fr(3,10))**n for n, c in enumerate(coef))
def fprim(x):  # primitive of f
    return sum(c*(x+Fr(3,10))**(n+1)/(n+1) for n, c in enumerate(coef))
def gen_pv(h): return {m: fpt((m - Fr(1, 2))*h) for m in range(-4, 6)}
def gen_av(h): return {m: (fprim(m*h) - fprim((m-1)*h))/h for m in range(-4, 6)}
for name, (cands, tc, glob, tg, *_r) in g.SCHEMES.items():
    if name in ('TENO6_OPT', 'TENO6_NENO6'): continue
    for lab, ind, data in (('recon-beta  on averages', g.indicator, gen_av), ('recon-beta  on points  ', g.indicator, gen_pv),
                           ('interp-beta on points  ', pv_indicator, gen_pv)):
        mats = [ind(c) for c in cands]; gm = ind(glob) if glob else None; taus = []
        for e in range(5, 10):
            v = data(Fr(1, 2**e))
            t = sum(Fr(c)*g.quad(m, [v[x] for x in cl]) for c, m, cl in zip(tc, mats, cands))
            if gm: t += tg*g.quad(gm, [v[x] for x in glob])
            taus.append(abs(t))
        print(f'  {name:8s} {lab}: tau orders {[round(math.log2(taus[i]/taus[i+1]), 2) for i in range(len(taus)-1)]}')
