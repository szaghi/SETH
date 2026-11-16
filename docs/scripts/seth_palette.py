"""The colours of SETH's documentation figures and animations: the Neon Moegi identity applied to data.

Why this module exists: every figure script used to hard-code its colours (ADAM's vermilion and amber accents). The
identity changed (Ōgi in Neon Moegi, docs/.vitepress/theme/custom.css), and the colours now live here, once, for the
figures (make_doc_figures, make_convention_figures, make_performance_figures) and the animations
(make_doc_animations, the cast of make_convention_figures).

The categorical lines were validated with the dataviz skill's validate_palette.js: lightness band, chroma floor,
colour-blind and normal-vision separation of adjacent pairs, contrast against the ground, all passing, light against
washi #f3ede1 and dark against silver ink #0e1116. The identity's neon moegi #a6ff2e (OKLCH L 0.91) is above the dark
band of data marks (0.48-0.67): it stays a UI accent, and dark data lines use #6aa41c.

The static figures are drawn twice by their scripts, with the light colours and, inside figure_theme.theme("dark")
(ADAM's harness, unchanged: it maps the neutral greys and the grounds), with the dark ones (SCHEME_STYLE_DARK,
LINES_DARK, LEFT_D, RIGHT_D). The animations draw each theme with its own colours.
"""

from __future__ import annotations

# identity grounds (custom.css, tokens.json of the identity bundle)
WASHI, WASHI_ALT, INK_LIGHT, MESH_LIGHT = '#f3ede1', '#ebe3d3', '#1d1b18', '#6f6a62'
SILVER, SILVER_ALT, INK_DARK, MESH_DARK = '#0e1116', '#0a0d11', '#e8ecf0', '#8b949e'
MOEGI_LIGHT, MOEGI_DARK_UI = '#4a7a06', '#a6ff2e'  # the identity's shock/seal colour on each ground

# categorical data colours, in fixed order (never cycled): blue, moegi, plum, umber
BLUE, MOEGI, PLUM, UMBER = '#1f63a8', '#4a7a06', '#8a4f9e', '#a8661c'
BLUE_D, MOEGI_D, PLUM_D, UMBER_D = '#3f8fd6', '#6aa41c', '#b06ac4', '#c27a2a'
LINES_LIGHT = (BLUE, MOEGI, PLUM, UMBER)
LINES_DARK = (BLUE_D, MOEGI_D, PLUM_D, UMBER_D)

# one colour and marker per scheme: the family is the hue (WENO-JS blue, WENO-Z moegi, TENO5 plum, TENO6 umber), the
# order or variant a lighter or darker shade of it, the marker the family again (a second encoding)
SCHEME_STYLE = {
    'weno1_js': ('0.62', '.'), 'weno3_js': ('0.38', 'v'),
    'weno5_js': (BLUE, 'o'), 'weno7_js': ('#4f8fcf', 'o'), 'weno9_js': ('#8db6e0', 'o'),
    'wenoz3': ('#2f5004', 'v'), 'wenoz5': (MOEGI, 's'), 'wenoz7': ('#76a52a', 's'), 'wenoz9': ('#a5c76a', 's'),
    'teno5': (PLUM, '^'), 'teno5_opt': ('#b77fb4', '^'), 'teno6': (UMBER, 'D'), 'teno6_opt': ('#d29a5a', 'D'), 'teno6_neno6': ('#6e420e', 'D'),
}

# the same for the dark twins, drawn by the scripts themselves (figure_theme keeps any stroke above 3:1, so it would
# leave the light shades, dim on silver ink): the validated dark lines for the four main schemes, shades around them
SCHEME_STYLE_DARK = {
    'weno1_js': ('0.62', '.'), 'weno3_js': ('0.38', 'v'),
    'weno5_js': (BLUE_D, 'o'), 'weno7_js': ('#6aa8e6', 'o'), 'weno9_js': ('#9cc6f0', 'o'),
    'wenoz3': ('#4c7d12', 'v'), 'wenoz5': (MOEGI_D, 's'), 'wenoz7': ('#8fc23e', 's'), 'wenoz9': ('#b5d77a', 's'),
    'teno5': (PLUM_D, '^'), 'teno5_opt': ('#d0a3cd', '^'), 'teno6': (UMBER_D, 'D'), 'teno6_opt': ('#dba264', 'D'), 'teno6_neno6': ('#9a5c1c', 'D'),
}

# the two biased stencils of the Conventions page: vp / the left state, vm / the right state
LEFT, RIGHT = BLUE, MOEGI
LEFT_D, RIGHT_D = BLUE_D, MOEGI_D
