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

## TENO5-A, TENO8-A (step 2.3)

| Item | Status | Where |
|---|---|---|
| m = 1 - min(1, eta_{i+1/2}/C_r), eta_j = (abs(2 Df_{j+1/2} Df_{j-1/2}) + e)/(Df_{j+1/2}^2 + Df_{j-1/2}^2 + e) | quoted | Fu 2021, Eqs. (17)-(18); Tang 2024, Eqs. (15)-(16); Liang & Fu 2024b, Eqs. (47)-(48) |
| g(m) = (1 - m)^4 (1 + 4 m), beta = a1 - a2 (1 - g), C_T = 10^-floor(beta) | quoted (Tang prints 10^floor(beta): the sign slip noted by the llm-wiki; Fu 2021 and Liang & Fu print the minus) | Fu 2021, Eq. (19); Liang & Fu, Eq. (49) |
| e = 0.9 C_r/(1 - 0.9 C_r) xi^2 | quoted | Tang 2024, Eq. (16); Liang & Fu, Eq. (48) |
| TENO5-A: C_r = 0.24, xi = 1e-3, a1 = 10, a2 = 5, min over j-1..j+1 | quoted | Tang 2024, after Eq. (16) |
| TENO8-A: C_r = 0.23, xi (their e_0) = 1e-3, a1 = 10.5, a2 = 3.5, min over i-1..i+2 | quoted | Liang & Fu 2024b, Eqs. (47)-(49) |
| Fu 2021's own set (C_r = 0.265, e = 1e-40, a1 = 14, a2 = 6.4) | quoted, not used (it is for the 8/10-point TENO-AA) | Fu 2021, after Eq. (19) |

Not in the papers: the rmu argument (xi/rmu), for a scale-invariant sensor.

## Affine-invariant descaler (step 1.1)

Wang, B.-S., Don, W.S. (2022). *Appl. Numer. Math.* 181:630.

| Item | Status | Where |
|---|---|---|
| Ai-weights: the indicators divided by mu^2, mu = mu_0 + (1/(2r-1)) sum abs(f_j - mean f) over the global stencil | quoted | Eqs. (9)-(14) |
| mu_0 = 1e-40 (double precision) | quoted; SETH floors mu at 64 eps mean abs(v) and tiny instead (same purpose, exactly scale invariant) | after Eq. (14) |
| Equivalent implementation: subtract the mean, apply the scale-invariant operator, add the mean back | quoted | Remark 2 |
| The Si descaler (mean abs(f)) fails on a small jump on an offset | quoted (their Fig. 5, WENO-Z); reproduced for WENO5-JS: 18 % overshoot | Section 3.2.1 |

SETH applies the descaler per bias (each state's own global stencil) and averages the two mu, since one rmu serves
both states of a face.

## Monotonicity-preserving bound (step 1.2)

Suresh, A., Huynh, H.T. (1997). *J. Comput. Phys.* 136:83.

| Item | Status | Where |
|---|---|---|
| minmod of k arguments, median(x, y, z) = x + minmod(y - x, z - x) | quoted | Eqs. (2.4)-(2.6) |
| v_UL = v_j + alpha (v_j - v_{j-1}), v_MP = v_j + minmod(v_{j+1} - v_j, alpha (v_j - v_{j-1})) | quoted | Eqs. (2.8), (2.12) |
| d_j = v_{j-1} + v_{j+1} - 2 v_j; d^M4 = minmod(4 d_j - d_{j+1}, 4 d_{j+1} - d_j, d_j, d_{j+1}) | quoted | Eqs. (2.19), (2.27) |
| v_MD = v_AV - d^M4_{j+1/2}/2, v_LC = v_j + (v_j - v_{j-1})/2 + 4/3 d^M4_{j-1/2} | quoted | Eqs. (2.28)-(2.29) |
| v_min, v_max: the intersection of I[v_j, v_{j+1}, v_MD] and I[v_j, v_UL, v_LC]; v <- median(v, v_min, v_max) | quoted | Eqs. (2.24)-(2.26) |
| alpha = 4 for Runge-Kutta, CFL <= 1/(1 + alpha) | quoted | section 2, after Eq. (2.10) |
| skip test with tolerance 1e-10 | quoted; not used (the median is applied everywhere, same result where the test would skip) | Eq. (2.30) |

## Sensors (steps 3.2, 3.3)

Li, Y., Chen, C., Ren, Y.-X. (2022). *J. Comput. Phys.* 448:110770, §3.2. Wan, Y., Xia, Y. (2022), §3.

| Item | Status | Where |
|---|---|---|
| Ren's psi_j = (2ab + e)/(a^2 + b^2 + e), a = abs(f_j - f_{j-1}), b = abs(f_j - f_{j+1}), e = 0.9 psi_c/(1 - 0.9 psi_c) xi^2, xi = 1e-3, psi_c = 0.3, sigma = min(1, psi/psi_c), min over j-1..j+2 | quoted | Eqs. (48)-(53) |
| Fu's detector: four stencils, gamma_k = 1/(beta_k + 1e-4)^6, chi normalised, smooth iff every chi_k > 5e-4; beta_3 on {i+1..i+3}: 1/4 (5 f_{i+1} - 8 f_{i+2} + 3 f_{i+3})^2 + 13/12 (...)^2 | quoted; the four indicators equal the generator's exact Jiang-Shu forms on cell i (checked) | Eqs. (43), (55)-(57) |
| Li, Chen & Ren's: a, b with the second differences added, min over j, j+1, xi = 1e-2, psi_c = 0.4 | quoted | Eq. (60) |
| Normalisation by max(abs(f_max), 1) | quoted; the tests pass it as rmu | Eq. (61) |
| Fig. 9: Fu correct for k < 0.64, the new detector up to k = 2.06; Ren's and Harten's fail at very low k | quoted; reproduced: 0.64, 2.07, 0.04 | §3.2 |
| Test (b): sin(12 pi x) - 2 / sin(24 pi x + pi/2) + 2 | quoted | §2.2 |
| Wan & Xia: eta = tau_5/(e + mean beta), e = 1e-6; smooth eta <= 1/2, non-smooth eta >= 1, M(eta) = 2 - 2 eta between | quoted | §3, Eq. (3.14) |

## Central fluxes and Ducros sensors (steps 4.1, 4.2)

| Item | Status | Where |
|---|---|---|
| Pirozzoli's conservative form F = 2 sum_l a_l sum_{m<l} f(i-m, i-m+l), a_l the central first-derivative coefficients | from ADAM's compute_kep_flux (W) and the llm-wiki (A); verified by the order of every kernel | ADAM adam_flume_hybrid_library.F90; wiki split-convective-forms |
| KEEP-PE internal energy flux <p><u_n>/(gamma - 1), mass flux <rho><u_n>, momentum f_rho <u> | quoted | De Michele et al. 2023, Eq. (15), "the scheme KEEP-PE in [Shima et al.]" |
| Kinetic energy flux f_rho (u_c . u_r)/2, split pressure work (u_n,c p_r + u_n,r p_c)/2 | from ADAM (W); verified by the kinetic-energy identity and the pressure equilibrium (to round-off) | |
| Shima et al. 2021 itself | not obtained (not in the corpus) | |
| Liang & Fu's eight-point a_l (0.8441336457279522, -0.2441336457279522, 0.0570096576929319, -0.006723831837710873) | quoted; order 6 measured | Liang & Fu 2024b, after Eq. (46) |
| Ducros phi = theta^2/(theta^2 + omega^2 + eps) | from ADAM's compute_ducros (W), equal to it bitwise with the override off | |
| De Vanna's psi = max(-div u/sqrt(div^2 u + rot^2 u + eps), 0), eps = (u_inf/delta_0)^2 | quoted | De Vanna et al. 2023a, Eq. (9) |

## Centred WENO (6.3)

ADAM's `weno_reconstruct_centered` (`adam_weno_object.F90`, added in ADAM commit 8a56ff17, 2025-10-20, "add new WENO
schemes (untested)"); no paper is cited in ADAM's code or documentation, and none was found.

| Item | Status | Where |
|---|---|---|
| The upwind part: WENO-JS of order 2S-1 on its candidates | quoted (Jiang & Shu 1996, as for `seth_weno*_js`) | — |
| The central rows on cells 1-S..S: (-1, 7, 7, -1)/12, (1, -8, 37, 37, -8, 1)/60, (-3, 29, -139, 533, 533, -139, 29, -3)/840 | derived: the reconstruction of order 2S from the cell averages (`scripts/gen_z_teno_tables.py`, check 7), equal to ADAM's `initialize_centered_S2..4` | ADAM |
| The blend r = r_c + phi (r_up - r_c), phi = max over k of abs(I_{k+1} - I_k)/(I_{k+1} + I_k), I_k = (v_k - v_{k-1})^2 + eps on cells -1..2 | unsourced: ADAM's own construction, ported as it is | ADAM l.1016-1023 |
| phi in terms of Ren's ratio: with a, b two consecutive differences and eps = 0, abs(b^2 - a^2)/(a^2 + b^2) = sqrt(1 - psi^2), psi = 2 abs(ab)/(a^2 + b^2), Ren et al.'s smoothness ratio (Li, Chen & Ren 2022, Eqs. (48)-(53)) | derived here (an identity: the two squares sum to 1) | — |
| D1: ADAM's left-state weights computed on the unshifted cells | found by the audit (`notes/adam-weno-audit.md`); fixed by construction (the upwind part is `seth_weno*_js`); confirmed by `seth_test_weno_c_adam` | ADAM l.1001-1007 |
| The descaler in the blend's indicators, I_k = (rmu dv)^2 + eps | SETH's choice (rmu = 1 is ADAM's arithmetic) | — |

Measured: orders 4.00/4.07, 6.03/6.00, 7.99/8.19 at a point (`seth_test_weno_c`); at first-order critical points
1, 3, 5 and at second-order ones 2, 2, 4 (`seth_test_critical_points`), the rates of the upwind WENO-JS, because phi
is O(1) where a difference vanishes; WENO8-C overshoots 1/280 at a step in its outer cells (-3, 4), which the blend's
detector (cells -1..2) does not read. ADAM's S = 1 (`weno-c-2`, documented as second order) never reaches the blend
and is WENO1: not ported.
