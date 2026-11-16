# M4 sources: what the kernels added after v0.1.0 rest on

Roadmap: issue #4. Each formula is classed as in `notes/m2-sources.md`: **quoted** (copied from the primary text,
equation or table given), **derived** (computed by `scripts/gen_z_teno_tables.py` in exact rational arithmetic from a
quoted definition, and cross-checked against the paper's printed values), or **unsourced** (not used). From this note
on every source was read in the paper itself (the PDFs of the llm-wiki corpus, `~/papers/mosaic/`), not through a
summary.

## TENO5-opt (step 2.1)

Fu, L., Hu, X.Y., Adams, N.A. (2016). *J. Comput. Phys.* 305:333. doi:10.1016/j.jcp.2015.10.037.

| Item | Status | Where |
|---|---|---|
| The optimised background of order K-1: eta g1 + (1 - eta) g2 | quoted | §4.1, Eq. (31) |
| Five-point sub-schemes g1 = (f_{i-2} - 5 f_{i-1} + 13 f_i + 3 f_{i+1})/12, g2 = (-f_{i-1} + 7 f_i + 7 f_{i+1} - f_{i+2})/12 | quoted, and derived (the generator checks they are the fourth-order reconstructions of their cells) | §4.2, Eq. (32) |
| Six-point sub-schemes (2, -13, 47, 27, -3)/60 and (-3, 27, 47, -13, 2)/60 | quoted, and derived (same check) | §4.3, Eq. (33) |
| d = 11/20, 2/5, 1/20 from eta = 0.2 | derived; equal to Table 3 (eta = 0.20) and Table 5 (TENO5-opt: 0.550, 0.400, 0.050) | Tables 3, 5 |
| TENO6-opt's d = 0.462, 0.3, 0.054, 0.184 from eta = 0.54 | derived; equal to Table 4 and Table 5 (until now quoted second-hand via the llm-wiki) | Tables 4, 5 |
| C_T = 1e-5 for TENO5-opt, q = 6, C = 1 | quoted | Table 5 |
| Candidate face coefficients for k = 4, 5 of the eight-point TENO | quoted, and derived | Table 1 |

Measured: linear-scheme order 4 (generator, exact); 3.99-4.01 at a point (`seth_test_teno`); 4.00 at a first-order
critical point and 2.00 at a second-order one (`seth_test_critical_points`).
