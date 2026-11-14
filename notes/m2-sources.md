# M2 sources: what the WENO-Z and TENO kernels rest on

Night 2026-11-14. Each formula is classed as **quoted** (copied from a primary text, equation number given),
**derived** (computed by `scripts/gen_z_teno_tables.py` in exact rational arithmetic from a quoted definition, and
cross-checked), or **unsourced** (from memory or secondary citation only: not used, or used and flagged).

Texts read in full (pdftotext of the arXiv PDFs):

- [LFA20] Y. Li, L. Fu, N.A. Adams, *A low-dissipation shock-capturing framework with flexible nonlinear dissipation
  control*, arXiv:2010.13289. Sections 2-3, Table 1.
- [NENO6] *A six-point neuron-based ENO (NENO6) scheme for compressible fluid dynamics*, arXiv:2207.08500. Section 2
  (TENO6-opt).
- [RN16] S. Rathan, G. Naga Raju, *An improved non-linear weights for seventh-order WENO scheme*, arXiv:1611.06755.
  Section 6.2 (WENO-Z).
- [LHYD21] Q. Li, X. Huang, P. Yan, Y. Duan, *Third-order scale-independent WENO-Z scheme achieving optimal order at
  critical points*, arXiv:2107.13225. Section 2.

Not obtained: Fu, Hu & Adams 2016 (JCP 305:333-359), Borges et al. 2008, Castro, Costa & Don 2011 (paywalled).

## TENO (Fu et al.), as restated by [LFA20] and [NENO6]

| Item | Status | Source |
|---|---|---|
| Candidate widths {3,3,3,4,...} for order K (eq. 8) and cells: S0 = {i-1,i,i+1}, S1 = {i,i+1,i+2}, S2 = {i-2,i-1,i}, S3 = {i,...,i+3} | quoted | [LFA20] eq. 8, Fig. 1, Table 1 |
| Candidate face values at x_{i+1/2} (k=0: -1/6, 5/6, 2/6; k=1: 2/6, 5/6, -1/6; k=2: 2/6, -7/6, 11/6; k=3: 3/12, 13/12, -5/12, 1/12) | quoted, and derived (the generator must reproduce them) | [LFA20] Table 1 |
| gamma_k = (C + tau_K/(beta_k + eps))^q, C = 1, q = 6, eps = 1e-40 | quoted | [LFA20] eq. 9, [NENO6] eq. 8 |
| beta_k = sum_{j=1}^{r-1} dx^{2j-1} int_{cell i} (d^j f_k/dx^j)^2 dx | quoted (definition); coefficients derived | [LFA20] eq. 10 |
| tau_K = beta_K - (beta_1 + beta_2 + 4 beta_0)/6, beta_K on the full K-point stencil (K = 6) | quoted (no absolute value; q = 6 is even, so the sign does not matter) | [LFA20] eq. 11, [NENO6] eq. 9 |
| chi_k = gamma_k / sum gamma, delta_k = 0 if chi_k < C_T else 1 | quoted | [LFA20] eqs. 12-13, [NENO6] eqs. 10-11 |
| w_k = d_k delta_k / sum d delta | quoted | [LFA20] eq. 18, [NENO6] eq. 12 |
| TENO6-opt weights d = 0.4294317718960898, 0.1727270875843552, 0.0855682281039113, 0.312272912415645; C_T = 1e-6 | quoted | [NENO6] section 2.2 |
| Max-order linear weights (TENO5: 0.6, 0.3, 0.1 on S0, S1, S2; TENO6: the unique d giving 6th order) | derived | exactness of the linear combination |
| tau_5 for TENO5 = abs(beta_1 - beta_2) (the WENO-Z tau_5 on the upwind and downwind 3-point stencils) | **unsourced** for TENO (Fu 2016 not obtained); its order O(dx^5) is checked exactly by the generator | memory + WENO-Z |
| C_T for TENO5 and for the max-order TENO6 | **unsourced**: no default constant is published by SETH; C_T is an argument | |

Evaluation detail taken in SETH (not a change of the scheme): chi_k = (g_k/g_max)^q / sum (g_j/g_max)^q with
g_k = C + tau/(beta_k + eps). It is the same chi, but it cannot overflow: with eps = 1e-40, gamma_k overflows double
precision once tau/beta exceeds about 1e51 (a jump of 1e10 on piecewise-constant data), which turns chi into NaN and
deletes every candidate.

## WENO-Z

| Item | Status | Source |
|---|---|---|
| alpha_k = d_k (1 + (tau/(beta_k + eps))^p), w = alpha / sum alpha | quoted | [RN16] eq. 38, [LHYD21] eq. 6 |
| p = 1 or 2 (WENO5-Z); WENO3-Z uses p = 1 | quoted | [LHYD21] after eq. 6, Table 1 |
| tau_3 = abs(beta_0 - beta_1) | quoted | [LHYD21] eq. 9 |
| tau_5 = abs(beta_0 - beta_2) | **secondary** (a search-engine summary of W.S. Don's slides; Borges 2008 not obtained); order O(dx^5) checked exactly by the generator | |
| tau_7 = abs(beta_0 - beta_1 - beta_2 + beta_3) (Castro et al., O(dx^7)); abs(beta_0 + 3 beta_1 - 3 beta_2 - beta_3) (Don & Borges, O(dx^8), better convergence) | quoted, with their attributions | [RN16] eqs. 39-40 |
| tau_9 | **unsourced**: WENO-Z9 is not implemented (left to the owner) | |
| beta numbering: beta_0 on the most upwind stencil {i-r+1..i} | consistent with [RN16] eq. 37 (the 4-point beta_3 printed there is ADAM's stencil-0 form mirrored); every tau above is invariant up to sign under the reversal k -> r-1-k, so the mirrored (right-state) reconstruction needs no other numbering | |

SETH uses the Don & Borges tau_7 (the better-converging one) and makes p an argument (1 or 2 recommended).

## Cross-checks the generator runs

1. Every candidate reproduces cell averages of polynomials up to degree w-1 at its face.
2. The face coefficients of the four TENO candidates equal [LFA20] Table 1.
3. The eq.-10 indicators of the 3-, 4- and 5-point upwind stencils equal ADAM's WENO-5/7/9 smoothness tables
   (checked to all digits in M1), divided by 1, 240 and 5040.
4. The max-order linear weights are exact to degree K-1 (K = 5, 6) and reproduce the WENO-5 weights for K = 5.
5. The order in dx of every tau on smooth data, measured in exact arithmetic on a Taylor polynomial of exp: tau_3
   >= 3, tau_5 >= 5, tau_6 >= 6, tau_7 (Don & Borges) >= 8.
