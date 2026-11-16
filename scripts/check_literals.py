#!/usr/bin/env python3
"""Find the kind-sensitive real literals outside the module that holds them (roadmap issue #4, 6.2).

SETH computes in PENF's build-time kind R_P (R8P by default, R4P with -D_R_P_IS_R4P). A literal whose meaning depends on
the kind belongs to src/lib/seth_kind_parameters.F90, one value per kind behind the same macros. Two classes are
caught here, on any kind suffix (_R_P, _R4P, _R8P, _R16P):

* magnitude below 1e-7 (exponent <= -8): below the epsilon of R4P, a step or a guard of that size vanishes against
  1 (1 - 1e-15 = 1), and below 1.2e-38 it is not a normal R4P number (1e-40 underflows, a GPU flushes it to zero);
* magnitude of 1e38 or more: outside the range of R4P.

A line carrying the marker `kind-ok` (with its reason, e.g. `! kind-ok: a TENO-A cutoff, a value not a step`) is
accepted. Exit status 1 when a literal is found, and every finding is printed as path:line: literal.

Usage: check_literals.py [PATH ...]   (default: src/lib/*.F90 src/lib/*.INC)
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

HOME = "seth_kind_parameters.F90"
MARKER = "kind-ok"
# a real literal with an exponent and a kind suffix: 1.e-40_R_P, 5e-7_R4P, 2.5d-9_R8P
LITERAL = re.compile(
    r"(?<![\w.])(\d+\.?\d*|\.\d+)[eEdD]([+-]?\d+)_(R_P|R4P|R8P|R16P)\b"
)


def findings(path: Path) -> list[str]:
    """The kind-sensitive literals of one file, as 'path:line: literal'."""
    out: list[str] = []
    for n, line in enumerate(path.read_text().splitlines(), start=1):
        if MARKER in line:
            continue
        code = line.split("!", 1)[0]  # a literal in a comment is not code
        for m in LITERAL.finditer(code):
            mantissa, exponent = float(m.group(1)), int(m.group(2))
            if mantissa == 0.0:
                continue
            magnitude = mantissa * 10.0**exponent
            if magnitude < 1.0e-7 or magnitude >= 1.0e38:
                out.append(f"{path}:{n}: {m.group(0)}")
    return out


def main() -> int:
    paths = [Path(p) for p in sys.argv[1:]] or sorted(
        Path("src/lib").glob("*.F90")
    ) + sorted(Path("src/lib").glob("*.INC"))
    found = [f for p in paths if p.name != HOME for f in findings(p)]
    for f in found:
        print(f)
    print(
        f"check_literals: {len(found)} kind-sensitive literal(s) outside {HOME}",
        file=sys.stderr,
    )
    return 1 if found else 0


if __name__ == "__main__":
    sys.exit(main())
