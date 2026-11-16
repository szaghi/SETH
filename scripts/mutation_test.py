#!/usr/bin/env python3
"""Mutation testing of SETH's test suite: does every test-worthy defect make some test fail?

Each mutant changes one thing in a scratch copy of the tree (a coefficient of a generated table, a weight, a tau
coefficient, the TENO cutoff or power, the WENO-JS exponent, the mirroring, the descaler, the weights ratio, the use
of eps, the interpolation tables and their wiring, the positivity limiter), the tests are rebuilt (--force-compile)
and run with scripts/run_tests.sh, and the mutant is *killed* when at least one test fails. A surviving mutant is a
defect the suite cannot see: the report lists them.

Usage (from the repository root, with fobis and the compilers of the mode on PATH, PENF fetched):

    python3 scripts/mutation_test.py [--mode tests-gnu] [--workdir DIR] [--only NAME ...]

The working tree is never modified: the copy goes to --workdir (default: a new temporary directory). Exit status 0
when every mutant is killed, 1 otherwise.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# (name, file, anchor, old, new): replace the first `old` after `anchor` (the anchor itself when old == anchor). Anchor
# in code, never in a comment: a mutated comment is an equivalent mutant (z-teno-mirror once survived that way).
MUTANTS = [
    ('js5-face-coefficient', 'src/lib/seth_weno_js_tables.F90', 'WENO5_JS_P(', '11._R8P/6._R8P', '11.000001_R8P/6._R8P'),
    ('js5-linear-weight', 'src/lib/seth_weno_js_tables.F90', 'WENO5_JS_A(', '0.3_R8P', '0.30001_R8P'),
    ('js7-indicator', 'src/lib/seth_weno_js_tables.F90', 'WENO7_JS_D(', '2107._R8P', '2108._R8P'),
    ('js3-interp-weight', 'src/lib/seth_weno_js_tables.F90', 'WENO3_JS_INTERP_A(', '3._R8P/4._R8P', '0.7501_R8P'),
    ('js-exponent-s3', 'src/lib/seth_weno_js_body.INC', 'case(3)', 'bp = b*(b*b)', 'bp = b*b'),
    ('js-mirror', 'src/lib/seth_weno_js_body.INC', 'v(1,m) = vm(m+1)', 'v(1,m) = vm(m+1)', 'v(1,m) = vm(m)'),
    ('js-eps-ignored', 'src/lib/seth_weno_js_body.INC', 'b = eps + ', 'b = eps + ', 'b = 1.e-6_R8P + '),
    ('js-weights-ratio', 'src/lib/seth_weno_js_body.INC', 'wr(f) = min(', 'wr(f) = min(', 'wr(f) = max('),
    ('descaling-ai-mean', 'src/lib/seth_descaling.F90', 'mu = mu + abs(vp(m) - avp)', 'abs(vp(m) - avp)', 'abs(vp(m))'),
    ('mp-curvature', 'src/lib/seth_mp.F90', 'lc = c + 0.5_R8P*(c - b)', '(4._R8P/3._R8P)*m4m', '(1._R8P/3._R8P)*m4m'),
    ('mp-median', 'src/lib/seth_mp.F90', 'v = max(vmin, min(v, vmax))', 'v = max(vmin, min(v, vmax))', 'v = min(v, vmax)'),
    ('descaling', 'src/lib/seth_descaling.F90', 'mu = mu + abs(vp(m))', 'mu = mu + abs(vp(m))',
     'mu = mu + abs(vp(m+1))'),
    ('teno6-linear-weight', 'src/lib/seth_z_teno_tables.F90', 'TENO6_D(0:3)', '9._R8P/20._R8P', '0.451_R8P'),
    ('teno5-az-tau', 'src/lib/seth_z_teno_tables.F90', 'TENO5_AZ_TC(0:2)', '-2._R8P', '-1._R8P'),
    ('teno5-az-cutoff', 'src/lib/seth_z_teno_body.INC', 'min(ct, (tau*rmu)*rmu)', 'min(ct, (tau*rmu)*rmu)', 'ct'),
    ('teno8-linear-weight', 'src/lib/seth_z_teno_tables.F90', 'TENO8_D(0:5)', '0.4130855804023061_R8P', '0.4131855804023061_R8P'),
    ('teno8-tau', 'src/lib/seth_z_teno_tables.F90', 'TENO8_TC(0:5)', '-2._R8P/3._R8P', '-0.6_R8P'),
    ('teno-a-sensor-min', 'src/lib/seth_z_teno_body.INC', 'etamin = min(etamin, eta)', 'etamin = min(etamin, eta)', 'etamin = max(etamin, eta)'),
    ('teno-a-exponent', 'src/lib/seth_z_teno_body.INC', 'ct = POW10(', 'a1 - a2*(1._R8P - gg)', 'a1 - a2*gg'),
    ('teno5-opt-linear-weight', 'src/lib/seth_z_teno_tables.F90', 'TENO5_OPT_D(0:2)', '11._R8P/20._R8P', '0.5501_R8P'),
    ('teno6-tau', 'src/lib/seth_z_teno_tables.F90', 'TENO6_TC(0:3)', '-2._R8P/3._R8P', '-0.6_R8P'),
    ('teno6-indicator-pivot', 'src/lib/seth_z_teno_tables.F90', 'TENO6_DQ(0:3,0:3)', '4._R8P/3._R8P', '1.4_R8P'),
    ('wenoz5-tau', 'src/lib/seth_z_teno_tables.F90', 'WENOZ5_TC(0:2)', '1._R8P, 0._R8P, -1._R8P',
     '1._R8P, -1._R8P, 0._R8P'),
    ('z-exponent-p2', 'src/lib/seth_z_teno_body.INC', 'case(2)', 'x = x*x', 'x = x*x*x'),
    ('z-teno-mirror', 'src/lib/seth_z_teno_body.INC', 'do m=LO, UMAX\n            st(m) = vm', 'st(m) = vm(1-m)',
     'st(m) = vm(-m)'),
    ('teno-cutoff-10x', 'src/lib/seth_z_teno_body.INC', 'if (g(k)/x >= ct)', 'if (g(k)/x >= ct)',
     'if (g(k)/x >= 10._R8P*ct)'),
    ('teno-power-4', 'src/lib/seth_z_teno_body.INC', 'g(k) = (g(k)*g(k))*(g(k)*g(k))*(g(k)*g(k))',
     'g(k) = (g(k)*g(k))*(g(k)*g(k))*(g(k)*g(k))', 'g(k) = (g(k)*g(k))*(g(k)*g(k))'),
    ('teno-global-indicator-off', 'src/lib/seth_z_teno_body.INC', 'tau = tau + TG*betag', 'tau = tau + TG*betag',
     'tau = tau + 0._R8P*betag'),
    # the interpolations (2026-11-15): a weight, a face coefficient, and a kernel wired to the reconstruction tables
    ('teno6-interp-weight', 'src/lib/seth_z_teno_tables.F90', 'TENO6_INTERP_D(0:3)', '15._R8P/32._R8P', '0.469_R8P'),
    ('wenoz5-interp-coefficient', 'src/lib/seth_z_teno_tables.F90', 'WENOZ5_INTERP_P(0:2,0:2)', '15._R8P/8._R8P',
     '15.0001_R8P/8._R8P'),
    ('teno5-interp-uses-recon-p', 'src/lib/seth_teno.F90', 'seth_teno5_interp(vp', 'CP(0:MAXW-1,0:NC-1) = TENO5_INTERP_P',
     'CP(0:MAXW-1,0:NC-1) = TENO5_P'),
    # the positivity limiter (2026-11-15): the pressure quadratic, the density bound, the floor, a one-sided state
    ('positivity-pressure-quadratic', 'src/lib/seth_positivity.F90', 'qb = s0(1)*d(nv)',
     '- sd - kappa*d(1)', '+ sd - kappa*d(1)'),
    ('positivity-density-bound', 'src/lib/seth_positivity.F90', 'if (s1(1) < eps_rho) t =',
     't = (s0(1) - eps_rho)/(s0(1) - s1(1))', 't = 0.5_R8P*(s0(1) - eps_rho)/(s0(1) - s1(1))'),
    ('positivity-pressure-floor', 'src/lib/seth_positivity.F90', 'ok = (gamma - 1._R8P)*(q(nv)', '>= eps_p',
     '>= 0._R8P'),
    ('positivity-right-state', 'src/lib/seth_positivity.F90', 'b0(k) = qr(k) + lam*f_low(k)',
     'b0(k) = qr(k) + lam*f_low(k)', 'b0(k) = qr(k) - lam*f_low(k)'),
]


def run(cmd: list[str], cwd: Path, timeout: int = 3600) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)


def failing_tests(work: Path, mode: str) -> tuple[bool, list[str]]:
    """Build the tests; return (built, names of the failing test programs)."""
    # --force-compile: FoBiS does not rebuild the dependents of a changed cpp #include (the kernels' .INC bodies), so
    # an incremental build would keep a stale mutant compiled in, or miss a new one
    b = run(['fobis', 'build', '--mode', mode, '--force-compile'], work)
    if b.returncode != 0 or 'Error' in b.stdout and 'error' in b.stdout.lower() and not (work / 'exe').exists():
        return False, []
    t = run(['bash', 'scripts/run_tests.sh'], work)
    # the report lines of run_tests.sh: '  FAIL  <name> (<reason>)' (a test's own 'FAILED: ...' messages are not)
    failed = [line.split()[1].split('/')[-1] for line in t.stdout.splitlines()
              if line.split()[:1] in (['FAIL'], ['XPASS'])]
    if t.returncode != 0 and not failed:
        failed = ['(run_tests.sh failed)']
    return True, failed


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--mode', default='tests-gnu')
    ap.add_argument('--workdir', type=Path)
    ap.add_argument('--only', nargs='*')
    args = ap.parse_args()
    repo = Path.cwd()
    work = args.workdir or Path(tempfile.mkdtemp(prefix='seth-mutants-'))
    if work.exists():
        shutil.rmtree(work)
    files = run(['git', 'ls-files'], repo).stdout.split() + run(['git', 'ls-files', '--others', '--exclude-standard'],
                                                                repo).stdout.split()
    for f in files:
        dst = work / f
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo / f, dst)
    shutil.copytree(repo / 'src/third_party/PENF', work / 'src/third_party/PENF')
    built, failed = failing_tests(work, args.mode)
    if not built or failed:
        print(f'baseline is not green ({failed}): fix the suite first')
        return 1
    print(f'baseline green in {work}')
    survivors = []
    for name, rel, anchor, old, new in MUTANTS:
        if args.only and name not in args.only:
            continue
        path = work / rel
        original = path.read_text()
        i = original.find(anchor)
        j = original.find(old, i) if i >= 0 else -1
        if i < 0 or j < 0:
            print(f'{name:28s} NOT APPLIED (anchor or text not found)')
            survivors.append(name)
            continue
        path.write_text(original[:j] + new + original[j + len(old):])
        built, failed = failing_tests(work, args.mode)
        path.write_text(original)
        if not built:
            print(f'{name:28s} does not build (invalid mutant)')
            survivors.append(name)
        elif failed:
            print(f'{name:28s} killed by {", ".join(sorted(set(failed)))}')
        else:
            print(f'{name:28s} SURVIVED')
            survivors.append(name)
        sys.stdout.flush()
    print(f'{len(MUTANTS) - len(survivors)} of {len(MUTANTS)} mutants killed' if not args.only else 'done')
    if survivors:
        print('survivors: ' + ', '.join(survivors))
    return 0 if not survivors else 1


if __name__ == '__main__':
    sys.exit(main())
