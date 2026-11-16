# Issue #5 sources: the Roe flux with entropy fix of the finite-difference TENO papers

Issue #5 asks whether the flux explains why Fu, Hu & Adams (2016) ran the Woodward–Colella blast waves with TENO5
and TENO6 at their published cutoffs and no positivity safeguard, while SETH's harness loses positivity there. Every
item was read in the paper itself (the PDFs of the llm-wiki corpus, `~/papers/mosaic/`, through `pdftotext`), and is
classed as in `notes/m4-sources.md`: **quoted** (copied from the primary text, place given), or **inferred** (not in
the text; the reasoning is given).

## What the TENO papers say

Fu, L., Hu, X.Y., Adams, N.A. (2016). *J. Comput. Phys.* 305:333. doi:10.1016/j.jcp.2015.10.037.

| Item | Status | Where |
|---|---|---|
| "As a default setting, the Roe average is exploited for characteristic decomposition at the cell face, and the Global Lax–Friedrichs (GLF) is utilized for flux splitting." Third-order SSP Runge–Kutta in 1-D, CFL 0.6 in all computations | quoted | §5, p. 345 |
| Blast waves: N = 400, t = 0.038, reference WENO-JS on N = 2500; "A Roe scheme with entropy-fix is applied for evaluating the numerical flux." No equation, no width, no citation | quoted | §5.3, p. 348 |
| No positivity safeguard named for any case | quoted (absence: §5 and §5.3 read in full) | §5 |

Fu, L. (2023). Review of the high-order TENO schemes. *Arch. Comput. Methods Eng.* Algorithm 1: Roe-averaged
characteristic decomposition, then "the Rusanov flux splitting", f± = (f ± λmax u)/2 with λmax "the maximum
characteristic velocity of the lth component over the entire computational domain" (quoted, Algorithm 1, step 3): the
GLF splitting of Fu 2016's default, without a margin. No blast-wave case.

Liang, T., Fu, L. (2024b). *Comput. Phys. Commun.* 302:109236. The algorithm of TENO-S: characteristic decomposition
with the Roe average, "the Rusanov flux splitting" (quoted, algorithm step 7). No blast-wave case.

So among the three TENO papers only Fu 2016 uses a Roe flux, only on this case, and does not define it.

## What "a Roe scheme with entropy fix" means in this lineage

Jiang, G.-S., Shu, C.-W. (1996). *J. Comput. Phys.* 126:202.

| Item | Status | Where |
|---|---|---|
| f± = (f ± αu)/2, α = max \|f'(u)\| over the whole range: "the global Lax–Friedrichs (LF) flux splitting. For other flux splittings, especially the Roe flux splitting with entropy fix (RF), see [15]"; [15] is Shu & Osher (1989) | quoted | §2, Eq. (2.5) and the text after |
| WENO-X: X = LF, Roe, RF "refers resp. to the global Lax–Friedrichs flux splitting, Roe's flux splitting and Roe's flux splitting with entropy fix" | quoted | §3 |
| The interacting blast waves computed with WENO-RF-3/4/5 (199, 399, 799 points, t = 0.038), reference ENO-RF-5 on 1600 points | quoted | §5, Example 3 |

Shu, C.-W., Osher, S. (1989). *J. Comput. Phys.* 83:32. doi:10.1016/0021-9991(89)90222-2.

| Item | Status | Where |
|---|---|---|
| ENO-Roe: the stencil's first point, hence the upwind side, chosen by the sign of the Roe speed a_{j+1/2} = (f(u_{j+1}) - f(u_j))/(u_{j+1} - u_j); the flux itself is reconstructed (no splitting) | quoted | §2, Eq. (2.4), Algorithm 2.1 |
| ENO-LLF: the split fluxes f ± α_{j+1/2} u, α_{j+1/2} = max \|f'(u)\| between u_j and u_{j+1} | quoted | Algorithm 2.2, Eqs. (2.11)–(2.13) |
| For convex f, α_{j+1/2} = max(\|f'(u_j)\|, \|f'(u_{j+1})\|) | quoted | Remark 2.2(a), Eq. (2.14) |
| ENO-RF (the entropy fix): Algorithm 2.1 "if f'(u) does not change sign between u_j and u_{j+1}; otherwise" Algorithm 2.2 | quoted | Algorithm 2.3 |
| For convex f the fix is needed only where f'(u_j) < 0 < f'(u_{j+1}), "the cells which contain expansion shocks", e.g. the Euler equations | quoted | Remark 2.3 |
| Systems: the algorithms applied in each local characteristic field of the Roe matrix at x_{j+1/2}; for the Euler equations (fields genuinely nonlinear or linearly degenerate) α^(p) = max(\|λ^(p)(u_j)\|, \|λ^(p)(u_{j+1})\|) | quoted | §4, Eqs. (4.1)–(4.4), Algorithm 4.1 |
| No Harten–Hyman width δ: the "entropy fix" is the switch to local Lax–Friedrichs | quoted (Algorithm 2.3 has no parameter) | §2 |

## The form SETH's harness uses (`PATH_ROE`)

Inferred: Fu 2016 names neither the form nor a source; "Roe scheme with entropy-fix" for a finite-difference WENO-type
scheme on this very case is Jiang & Shu's WENO-RF, i.e. Shu & Osher's Algorithm 2.3 in characteristic fields. That is
what `PATH_ROE` implements: at face i+1/2, the Roe-averaged eigenvectors project the stencil's states and fluxes,
w = L q, g = L f; per field k,

- if λ_k(q_i), λ_k(q_{i+1}) and the Roe eigenvalue λ_k(i+1/2) are all positive (all negative), g is reconstructed
  with the left (right) bias, unsplit;
- otherwise (a sign change or a zero, the sonic case), local Lax–Friedrichs: g± = (g ± α w)/2 with
  α = max(|λ_k(q_i)|, |λ_k(q_{i+1})|) (Eq. 4.4), reconstructed as on the splitting path.

Including the Roe eigenvalue in the sign test is a choice SETH makes (inferred): Algorithm 2.1 upwinds by the Roe
speed, Algorithm 2.3 switches on the end states; requiring all three to agree makes the upwind side well defined.
Algorithm 2.3 switches on any sign change, not only on expansions (Remark 2.3 relaxes it for convex fluxes); SETH
follows Algorithm 2.3. The switch `roe_fix = .false.` gives Algorithm 2.1 (ENO-Roe without the fix), for the negative
control on the transonic rarefaction.

Not known, and not recoverable from the paper: whether Fu et al. used Algorithm 2.3 or Remark 2.3's relaxed form, and
their CFL on this case (0.6 is the paper-wide setting).

## Measurements (2026-11-16, gfortran 16.2 -O3)

`seth_vv_roe`, section 5, N = 400 against WENO5-JS (splitting path) on 1600 cells: with the Roe flux and the fix,
without the limiter, at CFL 0.5 and at Fu et al.'s 0.6, WENO5-JS, WENO-Z5 and TENO6 (C_T = 1e-7) stay positive and
TENO5 (C_T = 1e-5) does not; the splitting path at CFL 0.6 keeps only WENO5-JS positive, as at 0.5.

One-off scratch variants of the harness (not kept; the same case and grid), to tell which part of the Roe path does
it:

| Variant | WENO5-JS | WENO-Z5 | TENO5 (1e-5) | TENO6 (1e-7) |
|---|---|---|---|---|
| Roe with the fix (`PATH_ROE`) | positive | positive | no | positive |
| Roe without the fix (`roe_fix = .false.`) | positive | positive | no | positive |
| local Lax–Friedrichs in every field, alpha of Eq. (4.4) | positive | positive | no | no |
| global Lax–Friedrichs, alpha = 1.1 max \|lambda\| (`PATH_SPLIT`) | positive | no | no | no |
| global Lax–Friedrichs, alpha = max \|lambda\| (no margin) | positive | no | no | no |

With C_T = 1e-4 for both TENO schemes: positive on the Roe path, with and without the fix, and with local
Lax–Friedrichs; on the splitting path TENO6 positive and TENO5 not (as `docs/verification/shocks/blast.md` records).

Reading (inferred from the table, not from a source): the entropy fix plays no part in positivity on this case; the
local wave speed of the splitting rescues WENO-Z5; TENO6 at 1e-7 needs the unsplit upwind reconstruction of Roe's
splitting; TENO5 at 1e-5 needs a larger cutoff or the limiter whatever the flux.
