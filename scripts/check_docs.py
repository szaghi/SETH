#!/usr/bin/env python3
"""Check that the documentation covers the code, and that it names nothing the code does not have (GIDEON's scheme).

Usage: check_docs.py [--verbose] [--only CHECK[,CHECK...]]

The pages are the hand-written ones under docs/ (docs/api, written by `formal`, and docs/node_modules excluded), with
the tables they include. Every check prints its count; a failure lists what is missing and the script exits 1.

Coverage (the code -> the pages):
  public     every public name of `seth` (kernels, constants, kinds) is named on a page of docs/reference
  modules    every module of src/lib is named on a page
  tests      every test program of src/tests is named on a page of the verification section
  examples   every program of docs/examples/src is included by a page, and every output too
  modes      every fobos mode is named on a page, or through a family (`static-*-r4`, `tests-{gnu,nvf}`)
  scripts    every script of scripts/ and docs/scripts/ is named on a page
  redirects  every static page of docs/public (a moved page's old URL) leads to a page that exists
  pages      every page is in a sidebar of config.mts; every sidebar entry is a page

Accuracy (the pages -> the code):
  names      every seth_* / SETH_* name on a page is a procedure, constant, module, program or file of the repository,
             or one of their families written with a placeholder (<N>, {3,5,7}, *)
  links      every image, internal link, include and snippet of a page resolves to a file
  anchors    every link to an anchor of a page (/page#anchor, ./page#anchor, #anchor) names one of its headings
"""

from __future__ import annotations

import argparse
import re
import sys
import unicodedata
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / 'docs'
PUBLIC = DOCS / 'public'
LIB = ROOT / 'src' / 'lib'

FAILURES: list[str] = []


def report(check: str, total: int, missing: list[str], verbose: bool) -> None:
    """Print one check's result and record its failures."""
    status = 'ok' if not missing else f'FAIL ({len(set(missing))})'
    print(f'{check:9s} {total:5d} checked  {status}')
    for item in sorted(set(missing)):
        FAILURES.append(f'{check}: {item}')
        print(f'          missing: {item}')
    if verbose and not missing:
        print('          (all covered)')


def pages(section: str | None = None) -> list[Path]:
    """Every hand-written page and included table, optionally of one section (its first path component)."""
    skip = {'api', 'node_modules', '.vitepress', 'public', 'examples'}
    found = [p for p in DOCS.rglob('*.md') if not skip & set(p.relative_to(DOCS).parts) and p.name != 'ford.md']
    if section is not None:
        found = [p for p in found if p.relative_to(DOCS).parts[0] == section]
    return sorted(found)


def text_of(paths: list[Path]) -> str:
    return '\n'.join(p.read_text(encoding='utf-8') for p in paths)


def strip_comments(src: str) -> str:
    return '\n'.join(line.split('!', 1)[0] for line in src.splitlines())


def public_names() -> set[str]:
    """The names `use seth` gives: the public statements of src/lib/seth.F90."""
    src = strip_comments((LIB / 'seth.F90').read_text(encoding='utf-8'))
    names: set[str] = set()
    for stmt in re.findall(r'^\s*public\s*::\s*(.*?)\s*$', src, re.M | re.I):
        names |= {n.strip() for n in stmt.split(',') if n.strip()}
    return names


def library_entities() -> set[str]:
    """Every procedure, module and named constant of src/lib (public or not)."""
    found: set[str] = set()
    for f in sorted(LIB.glob('*.F90')) + sorted(LIB.glob('*.INC')):
        src = strip_comments(f.read_text(encoding='utf-8'))
        found |= set(re.findall(r'^\s*module\s+(\w+)\s*$', src, re.M | re.I))
        found |= set(re.findall(r'^\s*(?:pure\s+)?(?:elemental\s+)?(?:logical\s+)?(?:subroutine|function)\s+(\w+)',
                                src, re.M | re.I))
        for decl in re.findall(r'parameter\s*::\s*(.*)', src, re.I):
            found |= set(re.findall(r'(?:^|,)\s*(\w+)\s*(?:\([^=]*\))?\s*=', decl))
    return found


def modules() -> set[str]:
    names: set[str] = set()
    for f in LIB.glob('*.F90'):
        names |= set(re.findall(r'^\s*module\s+(\w+)\s*$', f.read_text(encoding='utf-8'), re.M | re.I))
    return names


def programs(directory: Path) -> set[str]:
    names: set[str] = set()
    if directory.is_dir():
        for f in sorted(directory.glob('*.F90')):
            names |= set(re.findall(r'^\s*program\s+(\w+)', f.read_text(encoding='utf-8'), re.M | re.I))
    return names


def fobos_modes() -> list[str]:
    text = (ROOT / 'fobos').read_text(encoding='utf-8')
    block = re.search(r'^\[modes\]\n(.*?)(?=^\[)', text, re.M | re.S)
    if block is None:
        return []
    body = re.sub(r'#.*', '', block.group(1)).replace('modes =', ' ')
    return body.split()


def sidebar_links() -> set[str]:
    config = (DOCS / '.vitepress' / 'config.mts').read_text(encoding='utf-8')
    return set(re.findall(r"link:\s*'(/[^']*)'", config))


def page_of_link(link: str) -> Path | None:
    """The page a site link resolves to (without its hash), or None."""
    path = link.split('#', 1)[0].split('?', 1)[0]
    if not path.startswith('/'):
        return None
    candidates = [DOCS / (path.strip('/') + '.md'), DOCS / path.strip('/') / 'index.md']
    if path == '/':
        candidates = [DOCS / 'index.md']
    for c in candidates:
        if c.is_file():
            return c
    return None


def family_regex(name: str) -> re.Pattern[str]:
    """A name written with placeholders, as a pattern: `<N>`, `<p>` any word; `{a,b}` one of them; `*` anything."""
    out, i = [], 0
    for m in re.finditer(r'<[^>\s]+>|\{[\w,]+\}|\*', name):
        out.append(re.escape(name[i:m.start()]))
        tok = m.group(0)
        if tok == '*':
            out.append(r'\w*')
        elif tok.startswith('<'):
            out.append(r'\w+')
        else:
            out.append('(' + '|'.join(re.escape(a) for a in tok[1:-1].split(',')) + ')')
        i = m.end()
    out.append(re.escape(name[i:]))
    return re.compile('^' + ''.join(out) + '$', re.I)


def expand(name: str) -> list[str]:
    """The names a `{a,b}` family stands for, each of which must exist (other placeholders kept)."""
    m = re.search(r'\{([\w,]+)\}', name)
    if m is None:
        return [name]
    return [e for alt in m.group(1).split(',') for e in expand(name[:m.start()] + alt + name[m.end():])]


def check_public(verbose: bool) -> None:
    words = set(re.findall(r'\b\w+\b', text_of(pages('reference'))))
    names = public_names()
    report('public', len(names), [n for n in names if n not in words], verbose)


def check_modules(text: str, verbose: bool) -> None:
    words = set(re.findall(r'\b\w+\b', text))
    names = modules()
    report('modules', len(names), [m for m in names if m not in words], verbose)


def check_tests(verbose: bool) -> None:
    words = set(re.findall(r'\b\w+\b', text_of(pages('verification'))))
    names = programs(ROOT / 'src' / 'tests')
    report('tests', len(names), [n for n in names if n not in words], verbose)


def check_examples(text: str, verbose: bool) -> None:
    ex = DOCS / 'examples'
    srcs = sorted((ex / 'src').glob('*.F90')) if (ex / 'src').is_dir() else []
    missing = []
    for s in srcs:
        if f'examples/snippets/{s.stem}' not in text:
            missing.append(f'docs/examples/src/{s.name} (no page includes it or a region of it)')
    for out in sorted((ex / 'output').glob('*.txt')) if (ex / 'output').is_dir() else []:
        if f'examples/output/{out.name}' not in text:
            missing.append(f'docs/examples/output/{out.name} (no page includes it)')
    report('examples', len(srcs), missing, verbose)


def mode_family_named(mode: str, text: str) -> bool:
    for fam in re.findall(r'`([\w-]*\*[\w-]*)`', text):
        if re.fullmatch(re.escape(fam).replace(r'\*', r'[\w]+'), mode):
            return True
    for fam in re.findall(r'[\w-]*\{[\w,-]+\}[\w-]*', text):
        head, alts, tail = re.match(r'([\w-]*)\{([\w,-]+)\}([\w-]*)', fam).groups()
        if any(mode == f'{head}{a}{tail}' for a in alts.split(',')):
            return True
    return False


def check_modes(text: str, verbose: bool) -> None:
    modes = fobos_modes()
    report('modes', len(modes), [m for m in modes if not re.search(rf'(?<![\w-]){re.escape(m)}(?![\w-])', text)
                                 and not mode_family_named(m, text)], verbose)


def check_scripts(text: str, verbose: bool) -> None:
    scripts = [p for p in (ROOT / 'scripts').iterdir() if p.is_file() and p.suffix in ('.py', '.sh')]
    scripts += [p for p in (ROOT / 'scripts' / 'ci').glob('*') if p.is_file()]
    scripts += [p for p in (DOCS / 'scripts').glob('*.py')]
    missing = [str(p.relative_to(ROOT)) for p in scripts if p.name not in text]
    report('scripts', len(scripts), missing, verbose)


def check_redirects(verbose: bool) -> None:
    stubs = sorted(PUBLIC.rglob('*.html'))
    missing = []
    for stub in stubs:
        rel = stub.relative_to(PUBLIC).as_posix()
        targets = set(re.findall(r'url=/SETH(/[^"\s]*)"', stub.read_text(encoding='utf-8')))
        if len(targets) != 1 or page_of_link(next(iter(targets))) is None:
            missing.append(f'docs/public/{rel}: redirect to {sorted(targets)} (no such page)')
        if (DOCS / rel.replace('.html', '.md')).is_file():
            missing.append(f'docs/public/{rel}: shadows the page {rel.replace(".html", ".md")}')
    report('redirects', len(stubs), missing, verbose)


def check_pages(verbose: bool) -> None:
    all_pages = [p for p in pages() if not p.name.startswith('_') and '_tables' not in p.parts and p.parent != DOCS]
    sidebar = sidebar_links()
    reached = {page_of_link(link) for link in sidebar} - {None}
    missing = [f'{p.relative_to(ROOT)} (not in the sidebar)' for p in all_pages if p not in reached]
    missing += [f'sidebar link {link} (no page)' for link in sidebar if page_of_link(link) is None
                and not link.startswith('/api')]
    report('pages', len(all_pages), missing, verbose)


def check_names(text: str, verbose: bool) -> None:
    known = {n.lower() for n in library_entities() | public_names()}
    for d in ('tests', 'bench'):
        known |= {n.lower() for n in programs(ROOT / 'src' / d)}
        known |= {f.stem.lower() for f in (ROOT / 'src' / d).glob('*.F90')}
    known |= {n.lower() for n in programs(DOCS / 'examples' / 'src')}
    known |= {f.stem.lower() for f in (ROOT / 'src' / 'tests' / 'golden').glob('*')}
    known |= {f.stem.lower() for f in (ROOT / 'scripts').glob('*')} | {f.stem.lower() for f in (DOCS / 'scripts').glob('*')}
    known |= {f.stem.lower() for f in LIB.glob('*.INC')}
    # preprocessor macros of the library, and the environment variables the tests and scripts read
    for f in sorted(LIB.glob('*.F90')) + sorted(LIB.glob('*.INC')):
        known |= {m.lower() for m in re.findall(r'#\s*(?:define|undef|ifdef|ifndef|if\s+defined)\s+(SETH_\w+)',
                                                f.read_text(encoding='utf-8'))}
    for f in sorted((ROOT / 'src' / 'tests').glob('*.F90')) + sorted((ROOT / 'scripts').glob('*')):
        if f.is_file():
            known |= {m.lower() for m in re.findall(r"['\"$](SETH_\w+)", f.read_text(encoding='utf-8'))}
    known |= {'seth'}
    # a placeholder: <N>, <5|6>, <2S-1>, {3,5,7}, or a * not followed by a letter or digit (SETH_TENO*_CT, SETH_RBF_*)
    token = r'\b(?:seth|SETH)_(?:\w|<[^>\s]+>|\{[\w,]+\}|\*(?![A-Za-z0-9]))+'
    found = set(re.findall(token, text))
    missing = []
    for n in sorted(found):
        if n.lower() in known:
            continue
        if any(c in n for c in '<{*') and all(any(family_regex(e).match(k) for k in known) for e in expand(n)):
            continue
        missing.append(n)
    report('names', len(found), missing, verbose)


def check_links(verbose: bool) -> None:
    missing, total = [], 0
    for p in pages():
        src = p.read_text(encoding='utf-8')
        for img in re.findall(r'!\[[^\]]*\]\((/[^)\s]+)\)', src):
            total += 1
            if not (PUBLIC / img.lstrip('/')).is_file():
                missing.append(f'{p.relative_to(ROOT)}: image {img}')
        for media in re.findall(r'(?:src|poster)="(/[^"]+)"', src):
            total += 1
            if not (PUBLIC / media.lstrip('/')).is_file():
                missing.append(f'{p.relative_to(ROOT)}: media {media}')
        for link in re.findall(r'(?<!!)\[[^\]]*\]\((/[^)\s]*)\)', src):
            total += 1
            if link.startswith('/api') or (PUBLIC / link.lstrip('/').split('#')[0]).is_file():
                continue
            if page_of_link(link) is None:
                missing.append(f'{p.relative_to(ROOT)}: link {link}')
        for inc in re.findall(r'<!--@include:\s*(\S+?)\s*-->', src):
            total += 1
            if not (p.parent / inc).resolve().is_file():
                missing.append(f'{p.relative_to(ROOT)}: include {inc}')
        for snip in re.findall(r'<<<\s*@/(\S+?)(?:\{|\s|$)', src):
            total += 1
            if not (DOCS / snip).is_file():
                missing.append(f'{p.relative_to(ROOT)}: snippet {snip}')
    report('links', total, missing, verbose)


def slug(heading: str) -> str:
    """VitePress's anchor of a heading (its slugify, on the heading's text without inline markup)."""
    m = re.search(r'\{#([\w-]+)\}\s*$', heading)
    if m:
        return m.group(1)
    text = re.sub(r'\[([^\]]*)\]\([^)]*\)', r'\1', heading)  # a link: its text
    text = text.replace('`', '').replace('*', '')
    text = re.sub(r'[\u0300-\u036F]', '', unicodedata.normalize('NFKD', text))
    text = re.sub(r'[\s~`!@#$%^&*()\-_+=\[\]{}|\\;:"\'“”‘’<>,.?/]+', '-', text)
    text = re.sub(r'-{2,}', '-', text).strip('-')
    text = re.sub(r'^(\d)', r'_\1', text)
    return text.lower()


def anchors_of(page: Path) -> set[str]:
    found, fence = set(), False
    for line in page.read_text(encoding='utf-8').splitlines():
        if line.startswith('```'):
            fence = not fence
        elif not fence and re.match(r'#{1,6}\s', line):
            found.add(slug(re.sub(r'^#{1,6}\s+', '', line)))
    return found


def check_anchors(verbose: bool) -> None:
    """Every link to an anchor of a hand-written page (/page#anchor, ./page#anchor, #anchor) names a heading of it."""
    missing, total = [], 0
    for p in pages():
        src = p.read_text(encoding='utf-8')
        for path, anchor in re.findall(r'(?<!!)\[[^\]]*\]\(([^)\s#]*)#([^)\s]+)\)', src):
            if path.startswith(('/api', 'http')):
                continue
            total += 1
            if path == '':
                target = p
            elif path.startswith('/'):
                target = page_of_link(path)
            else:
                rel = (p.parent / path).resolve()
                target = next((c for c in (rel.with_suffix('.md'), rel / 'index.md', Path(str(rel) + '.md'))
                               if c.is_file()), None)
            if target is None:
                missing.append(f'{p.relative_to(ROOT)}: {path}#{anchor} (no page)')
            elif anchor not in anchors_of(target):
                missing.append(f'{p.relative_to(ROOT)}: {path}#{anchor} (no such heading)')
    for m in re.findall(r"link:\s*(/[^#\s]+)#([\w-]+)", (DOCS / 'index.md').read_text(encoding='utf-8')):
        total += 1
        target = page_of_link(m[0])
        if target is None or m[1] not in anchors_of(target):
            missing.append(f'docs/index.md: {m[0]}#{m[1]}')
    report('anchors', total, missing, verbose)


CHECKS = ('public', 'modules', 'tests', 'examples', 'modes', 'scripts', 'redirects', 'pages', 'names', 'links',
          'anchors')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split('\n', 1)[0])
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('--only', default=','.join(CHECKS), help='comma-separated checks to run (default: all)')
    args = parser.parse_args()
    only = set(args.only.split(','))
    unknown = only - set(CHECKS)
    if unknown:
        parser.error(f'unknown checks: {sorted(unknown)}')
    text = text_of(pages())
    run = {
        'public': lambda: check_public(args.verbose),
        'modules': lambda: check_modules(text, args.verbose),
        'tests': lambda: check_tests(args.verbose),
        'examples': lambda: check_examples(text, args.verbose),
        'modes': lambda: check_modes(text, args.verbose),
        'scripts': lambda: check_scripts(text, args.verbose),
        'redirects': lambda: check_redirects(args.verbose),
        'pages': lambda: check_pages(args.verbose),
        'names': lambda: check_names(text, args.verbose),
        'links': lambda: check_links(args.verbose),
        'anchors': lambda: check_anchors(args.verbose),
    }
    for name in CHECKS:
        if name in only:
            run[name]()
    if FAILURES:
        print(f'check_docs: {len(FAILURES)} failure(s)')
        return 1
    print('check_docs: ok')
    return 0


if __name__ == '__main__':
    sys.exit(main())
