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

## TENO5-AZ (step 2.2)

Tang, S. (2024). *Comput. Phys. Commun.* 299:109158.

| Item | Status | Where |
|---|---|---|
| tau = abs(beta_0 + beta_2 - 2 beta_1), WENO labels (beta_0 on {j-2..j}, beta_1 on {j-1..j+1}, beta_2 on {j..j+2}) | quoted; SETH's TENO labels give TC = (-2, 1, 1) on (S0 central, S1 downwind, S2 upwind) | Eqs. (14), (19) |
| C_T = min(lambda_d, tau), lambda_d = 1e-5 | quoted | Eq. (20) |
| gamma_k = (C + tau/(beta_k + eps))^q with tau of Eq. (19), C = 1, q = 6, eps = 1e-40 | quoted | Eqs. (12), (20) and the text after |
| order of tau on smooth data: 4 | derived (generator, exact) | |
| fifth order at critical points of order 2 and 3 (TENO-Z: 2.26-3.0) | quoted; reproduced: 5.00 at n_cp = 1 and 2 | Tables 1-3 |

Not in the paper: the `rmu` argument. The cutoff compares a pure number with tau, which has the dimension of the data
squared; SETH multiplies tau by rmu^2 (rmu = 1: Tang's form). Measured: the descaled form is bitwise scale invariant
and loses the critical-point property (order 2 at n_cp = 2).

## TENO8 (step 2.4)

Liang, T., Fu, L. (2024b). *Comput. Phys. Commun.* 302:109236. Fu, Hu & Adams (2016), Tables 1-2. Hong, Ye & Ye
(2022), *J. Comput. Phys.*

| Item | Status | Where |
|---|---|---|
| Candidates S0..S5, widths 3, 3, 3, 4, 4, 5 on i-3..i+4, and their face coefficients | quoted and derived (generator: equal) | Fu 2016, Table 1, Eq. (21) |
| Max-order weights (30, 18, 4, 12, 1, 5)/70 | derived; equal to Table 2, K = 8 | Fu 2016, Table 2 |
| tau_K = abs(beta_K - (beta_1 + beta_2 + 4 beta_0)/6), O(dx^6) | quoted | Liang & Fu, Eq. (15) |
| tau_8 written out over 62270208000 | quoted; the eight f_{i+4} terms equal the generator's exact beta_8 | after Eq. (15) |
| Background face coefficients (16 digits) | quoted; the weights reproduce them within 1.6e-15 | Eq. (46) |
| Linear weights d_0..d_5 (16 digits; they sum to 1 + 2.9e-15) | quoted; order 6 checked to 1e-12 | after Eq. (46) |
| C_T = 1e-7 for TENO8 | quoted | Hong et al. 2022, section 4 ("C_T = 10^-7 is adopted for TENO6 and TENO8") |

Measured, not in any paper read: the nonlinear ADR of TENO8 is anti-dissipative (3.2e-3) at theta/pi 0.36-0.45 for
C_T = 1e-7, and no mode of that band grows when advected.
