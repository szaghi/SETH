# M5 sources: what the non-polynomial reconstructions rest on

Roadmap: issue #6, step 0.1. Each formula is classed as in `notes/m2-sources.md` and `notes/m4-sources.md`:
**quoted** (copied from the primary text, equation, table and page given), **derived** (computed by this session,
2026-11-17, from a quoted definition, in exact or 400-digit arithmetic; not yet by the repository generator of step
0.4, which must redo it), **figure** (read from a rendered figure, no number asserted from it), or **unsourced**
(the papers are silent). Every source was read in the paper itself: `pdftotext -layout` and `-raw` of the corpus
PDFs (`~/papers/mosaic/weno-teno/`); garbled equations were resolved from the raw glyph order, from glyph
positions (`pdftotext -bbox`) and from page renders, and each such case is said below. The three throwaway checks
(`flags_check.py`, `thinc_check.py`, `rbf_check.py`) and the text dumps are in `.night/current/reads/`, outside the
repository.

Stencil convention, common to the three papers and to SETH's `TENO5_CLO = [-1, 0, -2]`: for the face `i+1/2` with
`df/du > 0`, `S0 = {i-1, i, i+1}` (central), `S1 = {i, i+1, i+2}` (downwind), `S2 = {i-2, i-1, i}` (upwind);
`tau_5 = |beta_1 - beta_2|`, i.e. downwind minus upwind indicator; the negative-speed state follows "by symmetry with
respect to the cell interface". The flags `delta_{k,j}` of a cell `j` are those of the left-biased reconstruction of
face `j+1/2` ("computed by the TENO weighting strategy in j-th cell", Takagi p. 6).

## Takagi, Fu, Wakimura & Xiao 2022 (TENO5-THINC)

*J. Comput. Phys.* 452:110899. doi:10.1016/j.jcp.2021.110899. Read in full (§1–§5, Tables 1–3, the captions of
Figs. 1–25).

| Item | Status | Where |
|---|---|---|
| Candidates `(-f_{i-1} + 5f_i + 2f_{i+1})/6`, `(2f_i + 5f_{i+1} - f_{i+2})/6`, `(2f_{i-2} - 7f_{i-1} + 11f_i)/6`; Jiang–Shu `beta_{k,3}` | quoted; equal to SETH's `TENO5_P` (S) | Eqs. (7), (10), p. 4 |
| `gamma_k = (1 + tau_K/(beta_{k,r} + eps))^q`, `eps = 1e-40`, `q = 6`; `tau_5 = abs(beta_{1,3} - beta_{2,3})` | quoted | Eq. (8) and text after Eq. (10), p. 4 |
| `chi_k = gamma_k / sum gamma`, `delta_k = 0 if chi_k < C_T else 1`, `C_T = 1e-5`; `d = 6/10, 3/10, 1/10` | quoted | Eqs. (11)–(13), p. 4–5 |
| The THINC function, Eq. (15) below, `beta = 1.6` (citing Deng et al. 2018, Comput. Fluids 171:1) | quoted | Eq. (15), p. 5 |
| The face formula, Eq. (16) below | quoted; the fraction `A = sgn(alpha)/T (B-1)/(B+1)` is resolved from the raw glyph order (the layout dump interleaves it) | Eq. (16), p. 5 |
| Eq. (16) is the exact face value of Eq. (15) under the conservation constraint | derived: `A = tanh(alpha beta/2)/tanh(beta/2)`; equivalently `f = f+ + f- (cosh beta - e^{-alpha beta})/sinh beta`; equal to a root-find of `d_i` plus quadrature within 2e-50 on 300 random monotone triples (`thinc_check.py`) | — |
| "The formula ... above is different from that in [5], which may cause the symmetry breaking from floating point errors" | quoted: the reason for the odd-in-`alpha` form | p. 5 |
| Guard: where `(f_{i+1} - f_i)(f_i - f_{i-1}) <= eps`, `eps = 1e-15`, the face is `f_i` | quoted. **The guard originates here**, as part of THINC, not in Liang 2022 | Eq. (17), p. 5 |
| The face is bounded by `f_{i-1}`, `f_{i+1}`; no limiter needed | quoted | p. 5 |
| Stronger: on strictly monotone triples the THINC face lies in `[f_i, f_{i+1}]` | derived: `t(alpha) = (cosh beta - e^{-alpha beta})/sinh beta` is concave with `t(-1) = -1`, `t(1) = 1`, so `t >= alpha`; checked on 300 triples | — |
| `beta = 1.6` "may lead to a positive imaginary part of the corresponding modified wavenumber in specific wavenumber ranges" (SuperBee-like anti-diffusion) | quoted | p. 5 |
| Step profile `f_i = 0 (i<5), 0.5 (i=5), 1 (i>5)` and its flags, Table 1, printed for cells 2..8 only | quoted; equal to Liang 2022 Table 1 on cells 2..8 (copied in the next section); recomputed (`flags_check.py`) | Eq. (18), Table 1, p. 6 |
| `zeta_i`, Eq. (19) below | quoted | Eq. (19), p. 6 |
| Detection, Eq. (20) below: reads `zeta` of cells `i-3..i+3` | quoted | Eq. (20), p. 6 |
| "if ... detected to contain a discontinuity at cell i ..., the cell interface flux `f_{i+1/2}` from TENO is replaced by that given by Eq. (16) **and Eq. (17)** from THINC; the same procedure applies ... `df/du < 0`" | quoted: the override includes the guard | §3.1 step 3, p. 7 |
| All `delta_k = 1`: THINC is never activated, TENO5's linear scheme is restored | quoted | §3.2 item 1, p. 7 |
| Euler setup: characteristic (Roe) projection, Rusanov splitting, SSP-RK3, CFL 0.4 "if not mentioned otherwise", `q = 6`, `C_T = 1e-5` "without tuning case-by-case"; symmetry-preserving technique of Fleischmann, Adami & Adams 2019 [11] | quoted | §4, p. 7 |
| Activation plots per characteristic component and per sign of the split flux, with `eps = 1e-8` (TENO's) "to avoid the potential effects from the machine round-off errors in the region with constant solutions"; Fig. 13: double Mach insensitive to `eps` in {1e-8, 1e-12, 1e-40} | quoted | §4.11, p. 19; Fig. 13 |
| Cost per time step, double Mach: 512×128: WENO5-Z 696, TENO5 810, present 946, TENO6 867 ms; 1024×256: 2802, 3255, 3810, 3330 ms ("16.8% (17.1%) more expensive than TENO5") | quoted | Table 3, p. 12 |

Eqs. (15)–(17), (19)–(20) (quoted):

$$h_i(x) = f_+ + f_-\tanh\big(\beta(X_i - d_i)\big),\quad f_+ = \tfrac{f_{i-1}+f_{i+1}}{2},\ f_- = \tfrac{f_{i+1}-f_{i-1}}{2},\ X_i = \tfrac{x - x_{i-1/2}}{x_{i+1/2}-x_{i-1/2}} - 0.5$$

with `d_i` from `f_i = int h` over the cell (printed without the `1/dx`; read as the cell average in the unit
coordinate `X`, the only reading under which Eq. (16) is exact);

$$\hat f_{i+1/2} = f_+ + f_-\,\frac{A+T}{1+AT},\quad T=\tanh(0.5\beta),\ A=\frac{\operatorname{sgn}(\alpha)}{T}\,\frac{B-1}{B+1},\ B=e^{|\alpha|\beta},\ \alpha=\frac{f_i-f_+}{f_-};$$

$$\hat f_{i+1/2} = f_i \quad\text{where } (f_{i+1}-f_i)(f_i-f_{i-1}) \le \varepsilon = 10^{-15};$$

$$\zeta_i = \begin{cases} 2 & \delta_{2,i}=0 \text{ and } \delta_{1,i}=0\\ 1 & \text{else if } \delta_{1,i}=0 \text{ and } \delta_{2,i}=1\\ -1 & \text{else if } \delta_{2,i}=0 \text{ and } \delta_{1,i}=1\\ 0 & \text{otherwise}\end{cases}$$

$$\text{detect}(i):\ \big(\zeta_{i-3}\in\{1,2\} \lor \zeta_{i-2}\in\{1,2\} \lor \zeta_{i-1}\in\{1,2\}\big) \land \zeta_i=0 \land \big(\zeta_{i+1}\in\{-1,2\} \lor \zeta_{i+2}\in\{-1,2\} \lor \zeta_{i+3}\in\{-1,2\}\big).$$

`alpha` is Liang 2022's `B` of Eq. (31): `alpha = (2 f_i - f_{i+1} - f_{i-1})/(f_{i+1} - f_{i-1})`. The guard makes
`alpha` well defined: `f_- = 0` gives a product `<= 0`, and a product above 1e-15 implies a strictly monotone triple,
`|alpha| < 1`.

**One-dimensional cases** (§4.1–4.5, quoted):

| Case | Setup | Parameters |
|---|---|---|
| Gaussian pulse, Eq. (22) | `exp(-300 (x - 0.5)^2)`, [0, 1], periodic, t = 1 | Table 2: h = 0.01, 0.005, 0.0025, 0.00125, 0.000625; `L_inf` 4.63E-3, 1.65E-4, 5.25E-06, 1.65E-07, 5.24E-09; orders —, 4.81, 4.98, 4.99, 4.97; identical columns for upwind-5 and the present scheme |
| Multiple waves (Jiang–Shu), Eqs. (24)–(26) | [0, 2] | N = 200, t = 200, SSP-RK3 (Fig. 2); N = 400, t = 2000, five-stage fourth-order SSP-RK (Fig. 3) |
| Sod, Eq. (27); Lax, Eq. (28) | [0, 1] | N = 100; t = 0.2 and 0.14 |
| Shu–Osher, Eq. (29) | Eq. (29) is written on [-5, 5] (`sin(5x)`, shock at x = -4) but the text says "the computational domain is [0,10]": an inconsistency of the paper | N = 200, t = 1.8, reference WENO5-JS N = 2000; CFL 0.4 and 0.6 compared (Fig. 7) |
| Blast waves, Eq. (30) | [0, 1] | N = 400, "Roe flux splitting with entropy fix", CFL 0.35, t = 0.038, reference WENO5-JS N = 32000. Fig. 9: `C_T = 1e-4` converges monotonically on N = 400..3200, `1e-5` leaves "a small overshoot ... near the strong contact" at very high resolution; the default 1e-5 is kept |

## Liang, Xiao, Shyy & Fu 2022 (TENO5-DR)

*J. Comput. Phys.* 467:111465. Read in full, §4.3 onwards now included.

| Item | Status | Where |
|---|---|---|
| Background `(2, -13, 47, 27, -3)/60`; candidates and indicators as Takagi; `d = 6/10, 3/10, 1/10` | quoted | Eqs. (7), (8), (12), p. 3–4 |
| `gamma_k = (C + tau_K/(beta_{k,r} + eps))^q`, `eps = 1e-40`, `q = 6`, `C = 1`; `tau_5 = abs(beta_{1,3} - beta_{2,3})`; `C_T = 1e-5` | quoted | Eq. (10), text after Eq. (12), Eq. (14), p. 4–5 |
| Table 1 (Eq. (16) profile), cells 0..10 | quoted; recomputed with constant continuation: equal on all 11 cells | Table 1, p. 5 |
| `zeta`, Eq. (17) | quoted; identical to Takagi Eq. (19) | p. 5 |
| Eq. (18) | quoted; Takagi's Eq. (20) verbatim | p. 6 |
| Eqs. (19)–(21) | quoted; Takagi's Eqs. (15)–(16) verbatim (Eq. (20) printed as `f_i = int h(x) dx`) | p. 6 |
| Profiles Eqs. (22)–(23) with `lambda_1 = 0.1`, `lambda_2 = 0.9`; Tables 2–3, "herein, C_T = 1e-5" | quoted; recomputed, equal on all 11 cells, for the profiles and their mirror images in value | p. 6–7 |
| In all four profiles cell 5 has `delta_1 = 0, delta_2 = 1`, `zeta = 1`, "which do not satisfy the condition defined by Eq. (18)" | quoted | p. 7 |
| Eq. (24) below; reads `zeta` of `i-3..i+2` | quoted; "or" between the two conditions and inside the rows read from the raw text, unambiguous | p. 7 |
| Static tests (a)–(e), below | quoted, two parts garbled (said there) | p. 7–8 |
| On (a), (b): "neither the discontinuity-detecting indicator of the TENO5-THINC scheme nor the present improved one has been triggered" | quoted | p. 8 |
| Non-monotone test Eq. (25) and central flux Eq. (26) | quoted | p. 10 |
| Local coordinates `x_{i-1} = -2.5`, `x_i = 0`, `x_{i+1} = 2.5` (`tanh(±2.5) ≈ ±1`), Eqs. (27)–(30) | quoted | p. 11 |
| Eq. (31) below, `A = 1.25 beta` | quoted | p. 11 |
| Step 1: Eq. (24) satisfied → Eq. (31) with `A = 0.91`; Step 2: Eq. (32) below "also satisfied" → `A = 1.23` | quoted; "and" inside the rows of Eq. (32) read from the raw text | p. 11 |
| "**all formulas** for approximating the cell interface flux in this section cannot be deployed in the regions where `(f_{i+1} - f_i)(f_i - f_{i-1}) <= eps`, `eps = 1e-15`. Instead ... `f_{i+1/2} = f_i`" | quoted (Eq. (33)); Takagi's Eq. (17) | p. 11 |
| Eq. (31) = `m + h tanh(A + atanh B)`, `m`, `h` mean and half difference of `f_{i±1}`, `B = (2f_i - f_{i+1} - f_{i-1})/(f_{i+1} - f_{i-1})` (Takagi's `alpha`) | derived (tanh addition); numerically to 1e-51 | — |
| `A = 0.91, 1.23` correspond to `beta = 0.728, 0.984` in Eq. (29) (not THINC's 1.6); `tanh(0.91) = 0.72113225407669986`, `tanh(1.23) = 0.84257932565892954` | derived | — |
| DR faces lie in `[f_i, f_{i+1}]` on strictly monotone triples (`tanh(A + atanh B) >= B` for `A > 0`) | derived; checked on 300 triples | — |
| Eq. (25) true ⇒ `f_i` strictly outside `[min, max]` of `f_{i±1}` ⇒ the product of Eq. (33) is negative ⇒ Eq. (26) never executes under the guard as written | derived; 0 counterexamples in 2·10^5 random triples (`thinc_check.py`) | — |
| Each condition of Eq. (32) implies the matching condition of Eq. (24); Eq. (32) reads only `i-2..i+2` | derived (by inspection) | — |
| Table 4, below; time step "reduced so that the time integration error can be ignored" (no value) | quoted | p. 12 |
| Euler setup: Roe-average characteristic projection, Rusanov splitting, SSP-RK3, CFL 0.4 "unless otherwise mentioned", "all parameters ... set as the recommended values"; symmetry-preserving techniques [51] Fleischmann et al. 2019 and [52] Wakimura, Takagi & Xiao 2022 | quoted | §4, p. 12 |
| No table after Table 4; no cost data; no ghost cells, non-uniform grids or basis of a neighbour's flags | quoted by absence (full text keyword-searched) | — |

Tables 1–3 (quoted; cells 0..10; each row recomputed bit for bit by `flags_check.py`):

| Profile `f_0..f_10` | `delta_0` | `delta_1` | `delta_2` | `zeta` (derived) |
|---|---|---|---|---|
| Table 1: 0 0 0 0 0 0.5 1 1 1 1 1 | 1 1 1 1 0 1 0 1 1 1 1 | 1 1 1 0 0 1 1 1 1 1 1 | 1 1 1 1 1 1 0 0 1 1 1 | 0 0 0 1 1 0 -1 -1 0 0 0 |
| Table 2: 0 0 0 0 0 0 1 1 1 1 1, and 1 1 1 1 1 1 0 0 0 0 0 | 1 1 1 1 1 0 0 1 1 1 1 | 1 1 1 1 0 0 1 1 1 1 1 | 1 1 1 1 1 1 0 0 1 1 1 | 0 0 0 0 1 1 -1 -1 0 0 0 |
| Table 3: 0 0 0 0 0 0.1 1 1 1 1 1, and 1 1 1 1 1 0.9 0 0 0 0 0 | 1 1 1 1 0 0 0 1 1 1 1 | 1 1 1 0 0 0 1 1 1 1 1 | 1 1 1 1 1 1 0 0 1 1 1 | 0 0 0 1 1 1 -1 -1 0 0 0 |

At the centred jump cell of Table 1, `beta_1 = beta_2 = 5/6` exactly, so `tau_5 = 0`, `chi_k = 1/3` and `delta =
(1,1,1)` (derived). Detector levels on cells 0..10 (derived, from the printed tables and from the recomputed flags):
Table 1 → cell 5 at level 2, Table 2 → cell 5 at level 1, Table 3 → cell 5 at level 2, every other cell at level 0;
Takagi's Eq. (20) flags cell 5 of Table 1 only. Reversing the data in space swaps `delta_1` and `delta_2` and keeps
`delta_0` (derived, on Table 3's profile).

Eqs. (24), (31), (32) (quoted):

$$\underbrace{\big(\zeta_{i-3}\in\{1,2\}\lor\zeta_{i-2}\in\{1,2\}\lor\zeta_{i-1}\in\{1,2\}\big)\land\zeta_i=0\land\big(\zeta_{i+1}\in\{-1,2\}\lor\zeta_{i+2}\in\{-1,2\}\big)}_{\text{condition I}}\ \lor\ \underbrace{\big(\zeta_{i-2}\in\{1,2\}\lor\zeta_{i-1}\in\{1,2\}\big)\land\big(\zeta_i=1\land\zeta_{i+1}=-1\big)\land\zeta_{i+2}\in\{-1,2\}}_{\text{condition II}}$$

$$\hat f_{i+1/2} = \frac{f_{i+1}+f_{i-1}}{2} + \frac{f_{i+1}-f_{i-1}}{2}\times\frac{\tanh(A) + \frac{2f_i-(f_{i+1}+f_{i-1})}{f_{i+1}-f_{i-1}}}{1+\tanh(A)\,\frac{2f_i-(f_{i+1}+f_{i-1})}{f_{i+1}-f_{i-1}}}$$

$$\big(\zeta_{i-2}\in\{1,2\}\land\zeta_{i-1}\in\{1,2\}\big)\land\zeta_i=0\land\big(\zeta_{i+1}\in\{-1,2\}\land\zeta_{i+2}\in\{-1,2\}\big)\ \lor\ \big(\zeta_{i-2}\in\{1,2\}\land\zeta_{i-1}\in\{1,2\}\big)\land\big(\zeta_i=1\land\zeta_{i+1}=-1\big)\land\zeta_{i+2}\in\{-1,2\}$$

Eq. (25): `|f_i - f_{i-1}| > |f_{i+1} - f_{i-1}|` or `|f_{i+1} - f_i| > |f_{i+1} - f_{i-1}|` → Eq. (26)
`(f_i + f_{i+1})/2`.

**Static tests** (§3.1, quoted), all on x in [0, 2], "assigned data distributions" (point values; where the nodes sit
and how the ends are treated is not stated, **unsourced**):

- (a) `sin(10 pi x)`, (b) `sin(24 pi x)`, N = 120 (scaled wavenumbers pi/6 and 2 pi/5, about 12 and 5 points per
  wavelength).
- (c) Jiang–Shu's multiple waves, `G(x, beta, z) = exp(-beta (x - z)^2)`, `F(x, alpha, a)` printed with a radical
  whose glyph is lost in the extraction (the figure shows the half-ellipse, i.e. Jiang–Shu's
  `sqrt(max(1 - alpha^2 (x - a)^2, 0))`), `a = 0.5, z = -0.7, theta = 0.005, alpha = 10, beta = log 2/(36 theta^2)`.
- (d) `sin(8.0x)` on [0, 0.4]; −1.0 on (0.4, 0.5]; 1.0 on (0.5, 0.6]; 0.0 on (0.9, 1.0]; 1.0 on (1.0, 1.1]; 0.5 on
  (1.1, 1.2]; `e^x - 5.0` on (1.5, 1.8]; `sin(9.42x)` otherwise. The layout dump misaligns values and intervals; the
  pairing is from the glyph rows (`pdftotext -bbox`, page 8).
- (e) 1.25 on (0.27, 0.29], (0.63, 0.66], (1.14, 1.17], (1.53, 1.56]; −1.25 on [0.25, 0.27], (0.6, 0.63], (1.5, 1.53];
  `sin(2 pi x e^{x+1})` otherwise. The superscript `x+1` and the trailing `x` were read from glyph positions; the
  figure confirms it (`f(0) = 0`, first maximum near x = 0.08).
- N: the text says "Case (c) and case (d) are computed on uniform grids with N = 100 and N = 400, respectively"; the
  captions of Figs. 3 and 4 show each case at both N = 100 and N = 400. The captions decide. (e): N = 200 and 400.

What the figures show (figure; pages rendered at 300–600 dpi, positions read against the axis):

| Case | DR detector ("present") | THINC detector (TENO5-THINC) |
|---|---|---|
| (c), N = 100 and 400 | x ≈ 0.6, 0.8 (the square wave's jumps) and ≈ 1.1 (the triangle's peak: a jump of the derivative) | N = 100: several cells at 1.05–1.15; N = 400: 1.1 only; the square wave missed |
| (d), N = 100 and 400 | all nine jumps (0.4, 0.5, 0.6, 0.9, 1.0, 1.1, 1.2, 1.5, 1.8) | one activation near 0.5 at N = 100; none at N = 400 |
| (e), N = 400 | x ≈ 0.245, 0.27, 0.287, 0.60, 0.63, 0.66, 1.135, 1.50, 1.53, **1.93**; no activation at the jumps at **1.17** (1.25 → 1.00) and **1.56** (1.25 → 0.90), the two smallest of the eleven; 1.93 is in the smooth sine at 3.6 points per wavelength | 0.27, 1.93 |
| (e), N = 200 | activations at the jumps up to 1.53, plus about seven at 1.58–1.87 in the sine beyond the last jump, at 2.0–2.9 points per wavelength | 0.27, 0.58, 1.13, 1.47, and about seven at 1.58–1.87 |

Table 4 (quoted; the three columns are identical):

| N | fifth-order linear `L_inf` | order | TENO5 `L_inf` | order | present `L_inf` | order |
|---|---|---|---|---|---|---|
| 80 | 1.40E-2 | * | 1.40E-2 | * | 1.40E-2 | * |
| 160 | 4.94E-4 | 4.826 | 4.94E-4 | 4.826 | 4.94E-4 | 4.826 |
| 320 | 1.59E-5 | 4.955 | 1.59E-5 | 4.955 | 1.59E-5 | 4.955 |
| 640 | 5.02E-07 | 4.988 | 5.02E-07 | 4.988 | 5.02E-07 | 4.988 |
| 1280 | 1.57E-08 | 4.997 | 1.57E-08 | 4.997 | 1.57E-08 | 4.997 |

**Cases of §4.1 onwards** (quoted; there is no table after Table 4, and no cost data):

| § | Case | Grid, time | Changes of `C_T`, CFL, flux |
|---|---|---|---|
| 4.1 | Gaussian, Eq. (35), [0, 1], periodic | N = 80..1280, t = 1 | time step reduced, value not stated |
| 4.2 | Sod, Eq. (36); Lax, Eq. (37) | N = 100; t = 0.2, 0.14 | — |
| 4.3 | Blast waves, Eq. (38), reflective at x = 0 and 1 | N = 400, t = 0.038 (Fig. 9); reference WENO5-JS N = 2500 | Roe with entropy fix; "**C_T is particularly adjusted to 10^-4 for the TENO5 scheme** to obtain a stable result" (the present scheme's `C_T` there is not stated) |
| 4.4 | Shu–Osher, Eq. (39): `(3.857, 2.629, 10.333)` on [0, 1), `(1 + 0.2 sin(5x), 0, 1)` on [1, 10] | N = 200, t = 1.8; reference N = 2000 | — (note `sin(5x)`, not `sin(5(x - 5))`: a phase-shifted variant of the standard case) |
| 4.5 | Double Mach, Eq. (40), [0, 4]×[0, 1] | 360×90, 512×128, 1024×256; t = 0.2 | CFL 0.2, 0.4, 0.6 compared (Fig. 15) |
| 4.6 | 2-D Riemann config. 3, Eq. (41); config. 7, Eq. (42) | 512², 1024² (t = 0.3); 1024² (t = 0.25) | — |
| 4.7 | Shock–vortex, Eq. (43), [0, 2]×[0, 1] | 300×150, t = 0.6; reference 2000×1000 | — |
| 4.8 | Single-material triple point, Eq. (44), [0, 7]×[0, 3] | 350×150, t = 5.0 | — |
| 4.9 | Implosion, Eq. (45), periodic [−0.3, 0.3]² | 800×800, t = 0.8 | — |
| 4.10 | Rayleigh–Taylor, Eq. (46) | 64×256, 128×512, t = 1.95 | Rusanov |
| 4.11 | Astrophysical jets, Mach 80 and 2000 | 800×400 | CFL 0.35, positivity-preserving flux limiter |

## Liang & Fu 2024 (Gaussian-RBF TENO)

*J. Comput. Phys.* 497:112618. Re-read in full for the items below.

| Item | Status | Where |
|---|---|---|
| TENO5 as above: `gamma_k = (1 + tau_5/(beta_{k,r} + eps))^q`, `eps = 1e-40`, `q = 6`; `tau_5 = abs(beta_{1,3} - beta_{2,3})`; `C_T = 1e-5` "if not mentioned otherwise"; `d = 0.6, 0.3, 0.1` | quoted (in Eq. (11) the square of `beta_{2,3}`'s first term is displaced in the extraction) | Eqs. (8)–(14), p. 3–4 |
| The six stencil combinations, below | quoted | §3.1, p. 5 |
| Fig. 2, six panels (non-smooth location → target stencil): none → `i-2..i+2`; between `x_{i+1}` and `x_{i+2}` → `i-2..i+1`; between `x_{i-2}` and `x_{i-1}` → `i-1..i+2`; both → `i-1..i+1`; at `x_{i-1/2}` → `i..i+2`; at `x_{i+1/2}` → `i-2..i`. The stencil "utilizes as many smooth data as possible while not crossing the non-smooth region" | figure (rendered); the principle quoted | Fig. 2, p. 5 |
| The conditions overlap: `(1,0,1)` satisfies those of `f^2, f^4, f^6`; `(1,1,0)` those of `f^3, f^4, f^5`; `(0,1,1)` those of `f^5, f^6` and matches no panel of Fig. 2 | derived (by inspection) | — |
| Gaussian kernel `phi = exp(-lambda^2 x^2)`; RBF of the primitive `H` through the faces of the stencil, `R'(x)` the flux | quoted: `lambda` multiplies a length, so it has dimension 1/length | Eqs. (16)–(24), p. 5–6 |
| "Applying the Taylor expansion ... and truncating the expansion up to degree `lambda^2 dx^2` **for brevity**, the explicit candidate stencil fluxes ... are given as" Eq. (25) | quoted. Whether the authors' code used the truncated rows or the full Gaussian is not stated (**unsourced**) | p. 6 |
| Eq. (25), six rows, below | quoted; every coefficient re-derived (next row) | p. 7 |
| Each row is the `O(s^2)` truncation, `s = lambda dx`, of the Gaussian-RBF face flux on its stencil | derived: 400-digit solve at `s = 1e-30, 2e-30`; all printed coefficients match to 1e-20, including the absent `f_{i-2}` term in `R'_6`'s `s^2` part (`rbf_check.py`) | — |
| Every `s^2` row sums to zero; at `lambda = 0`: `R'_1` = upwind-5, `R'_2 = (-1, 7, 7, -1)/12 = q0/2 + q1/2`, `R'_3 = (1, -5, 13, 3)/12 = 3q0/4 + q2/4`, `R'_4..R'_6 = q0, q1, q2` | derived. `R'_2`, `R'_3` at `lambda = 0` are TENO5-opt's `g2`, `g1` (Fu, Hu & Adams 2016, Eq. (32); `notes/m4-sources.md`) | — |
| **Label swap**: `R'_2` sits on `i-1..i+2` (the stencil of `f^3`, `delta_2 = 0`), `R'_3` on `i-2..i+1` (the stencil of `f^2`, `delta_1 = 0`). Coding `R'_k` for `f^k` by index would put the rejected downwind point `i+2` into the `(1,0,1)` face | derived from Eq. (25) against §3.1 and Fig. 2: an inconsistency of the paper | — |
| Eq. (28) cell averages in Taylor form; Eq. (29) truncation errors, below | quoted | p. 7 |
| Eq. (29) from the printed Eq. (25) | derived symbolically (`rbf_check.py`): the printed leading terms of all six rows match. `R'_1`'s next term is `dx^6 (28 lambda^2 h^(4) + h^(6))/140`; the paper prints only its `lambda^2` part (`12/60 = 1/5`, matches) | — |
| Sensor Eqs. (30)–(31), `c_r = 0.4`, `eps = 1e-40` in numerator and denominator; Eq. (32), `alpha_1 = 3.0`, `alpha_2 = 3.0` "determined by numerical experiments"; `m ≈ 0` → `lambda = alpha_1`, `m ≈ 1` → `alpha_1 + alpha_2` | quoted | p. 8 |
| "as the computational mesh is progressively refined, i.e., `dx → 0`, the present Gaussian RBF reconstruction will degenerate to the classical polynomial reconstruction" | quoted: the paper's own statement that `lambda` is fixed while `dx` shrinks | p. 8 |
| ADR of Fig. 3, "conducted on 128 uniformly distributed points": `lambda`, domain and `dx` not stated, so the `s` of the curves is unknown | quoted; **unsourced** `s` | Fig. 3, p. 8 |
| Eq. (33) = Liang 2022 Eq. (17); Eq. (34) = Eq. (24); Eq. (35) = Eq. (31) with the steepness 0.91 (its symbol is lost in the extraction; Liang 2022 calls it `A`); Eq. (36) = Eq. (32), 1.23 "as the same value in [48]"; guard 1e-15 → `f_i` for Eq. (35); no central flux | quoted; compared term by term: identical | §3.3, p. 8–9 |
| Setup: Rusanov splitting, Roe-average characteristic projection, SSP-RK3 Eq. (37), CFL 0.4; symmetry-preserving technique [59] (Fleischmann et al. 2019) | quoted | §4, p. 9 |
| Table 1, double Mach, total cost: 400×100: WENO5-Z 36.51, TENO5 40.04, TENO6 45.00, present 55.23 s; 1000×250: 593.16, 630.89, 717.20, 895.52 s | quoted; ratios to TENO5 1.379 and 1.419 (derived) | Table 1, p. 13 |
| Optimal-shape-parameter RBF schemes: "the numerical division between higher-order derivatives ... is unstable ... especially ... at the critical points" | quoted | §1, item (iii), p. 2 |

The stencil combinations of §3.1 (quoted):
`f^1 | delta_0 delta_1 delta_2 = 1`, `f^2 | delta_0 delta_2 = 1, delta_1 = 0`, `f^3 | delta_0 delta_1 = 1, delta_2 = 0`,
`f^4 | delta_1 delta_2 = 0, delta_0 = 1`, `f^5 | delta_0 delta_2 = 0, delta_1 = 1`, `f^6 | delta_0 delta_1 = 0, delta_2 = 1`.

Eq. (25) (quoted; `s^2 = lambda^2 dx^2`):

$$\begin{aligned}
R'_1 &= -\tfrac{1}{60}\big(3f_{i+2}-27f_{i+1}-47f_i+13f_{i-1}-2f_{i-2} + s^2(18f_{i+2}-42f_{i+1}+18f_i+18f_{i-1}-12f_{i-2})\big)\\
R'_2 &= -\tfrac{1}{12}\big(f_{i+2}-7f_{i+1}-7f_i+f_{i-1} + s^2(4f_{i+2}-4f_{i+1}-4f_i+4f_{i-1})\big)\\
R'_3 &= \tfrac{1}{12}\big(3f_{i+1}+13f_i-5f_{i-1}+f_{i-2} + s^2(9f_{i+1}-15f_i+3f_{i-1}+3f_{i-2})\big)\\
R'_4 &= \tfrac{1}{6}\big(2f_{i+1}+5f_i-f_{i-1} + s^2(4f_{i+1}-2f_i-2f_{i-1})\big)\\
R'_5 &= -\tfrac{1}{6}\big(f_{i+2}-5f_{i+1}-2f_i + s^2(2f_{i+2}+2f_{i+1}-4f_i)\big)\\
R'_6 &= -\tfrac{1}{6}\big(-11f_i+7f_{i-1}-2f_{i-2} + s^2(18f_i-18f_{i-1})\big)
\end{aligned}$$

Eq. (29) (quoted; `h` and its derivatives at `x_{i+1/2}`):

$$\begin{aligned}
R'_1 &= -\tfrac{1}{60}\big(-60h + dx^5(30h'''\lambda^2 + h^{(5)}) - 12h^{(4)}dx^6\lambda^2\big) + O(dx^7) &
R'_2 &= -\tfrac{1}{30}\big(-30h + dx^4(20h''\lambda^2 + h^{(4)})\big) + O(dx^6)\\
R'_3 &= -\tfrac{1}{120}\big(-120h + dx^4(-120h''\lambda^2 - 6h^{(4)})\big) + O(dx^5) &
R'_4 &= \tfrac{1}{360}\big(360h + dx^3(360h'\lambda^2 + 30h''')\big) + O(dx^4)\\
R'_5 &= -\tfrac{1}{360}\big(-360h + dx^3(360h'\lambda^2 + 30h''')\big) + O(dx^4) &
R'_6 &= -\tfrac{1}{120}\big(-120h + dx^3(360h'\lambda^2 + 30h''')\big) + O(dx^4)
\end{aligned}$$

From `R'_1` (derived here, the plan's §1.3 reproduced): `u_t + u_x = dx^5 (lambda^2/2 u_xxxx + 1/60 u_xxxxxx) + ...`;
the imaginary part of the modified wavenumber is `xi^4 (s^2/2 - xi^2/60)`, positive for `xi^2 < 30 s^2`, with
maximum `(200/3) s^6` at `xi^2 = 20 s^2`.

Eqs. (30)–(32) (quoted):

$$m = 1-\min(1,\theta_i/c_r),\quad \theta_i=\min(\eta_{i-1},\eta_i,\eta_{i+1}),\quad \eta_i=\frac{|2(f_i-f_{i-1})(f_{i+1}-f_i)|+\varepsilon}{(f_i-f_{i-1})^2+(f_{i+1}-f_i)^2+\varepsilon},\quad c_r = 0.4,\ \varepsilon=10^{-40}$$

$$\lambda=\alpha_1+\alpha_2\big(1-g(m)\big),\quad g(m)=(1-m)^4(1+4m),\quad \alpha_1=\alpha_2=3.$$

**Cases** (§4, quoted): Gaussian pulse (Fig. 4, 100 cells and a refinement, "time step ... small enough"); Sod and
Lax, 100 points; blast waves, 400 points, reflective, Roe with entropy fix, t = 0.038, "**C_T is particularly
adjusted to 10^-4 for the TENO5 scheme**", reference WENO5-JS N = 10000; Shu–Osher, Eq. (43),
`1 + 0.2 sin(5(x - 5))` on [1, 10], domain [0, 10], N = 200, t = 1.8; double Mach 400×100 and 1000×250; Rayleigh–
Taylor 64×256; Kelvin–Helmholtz, Eq. (46), [0, 1]², periodic, 400×400 (Fig. 13 caption), t = 0.5; 2-D Riemann
configurations 3 (t = 0.3), 8 (t = 0.25), 11 (t = 0.3), 800×800. No table of errors; Figs. 5–8 are Sod, Lax, blast,
Shu–Osher.

## Other findings

- **Symmetry.** All three papers state that their low-dissipation schemes break the symmetry of symmetric flows
  through round-off and apply the symmetry-preserving technique of Fleischmann, Adami & Adams (2019) (Liang 2022 also
  Wakimura, Takagi & Xiao 2022). Step 4.3 asserts symmetry without that technique: expect a failure that is not a bug.
- **Shu–Osher has three variants**: Takagi's Eq. (29) on [−5, 5] with a text saying [0, 10]; Liang 2022's `sin(5x)`
  on [1, 10]; Liang & Fu 2024's `sin(5(x - 5))` on [1, 10] (the standard case shifted). Figure comparisons (step
  3.4, E5) must use the variant of the paper compared with.
- **`C_T` on the blast waves.** Both Liang papers run the comparison TENO5 at `C_T = 1e-4`; the composite schemes'
  `C_T` there is not stated (unsourced). Takagi runs TENO5-THINC at 1e-5 and shows 1e-4 converges monotonically.
- **Flags near round-off.** With `eps = 1e-40` the flags on constant data polluted by round-off are arbitrary; Takagi
  switches to `eps = 1e-8` for the activation plots. The detector inherits whatever `eps` the caller passes.
- **THINC's cost** (Takagi, Table 3): 1.17 times TENO5 per time step, whole solver, double Mach. A second published
  anchor for E1 beside Liang & Fu's 1.38–1.42 for RBF+DR.
- **Long-time contact data** for E2 exist in Takagi: the square wave of the multiple-wave case at N = 200, t = 200
  (100 periods) and N = 400, t = 2000 (Figs. 2–3, figure only).

## Plan claims checked

| Claim (issue #6) | Plan tag | Verdict | Evidence |
|---|---|---|---|
| §1.2: THINC is Takagi's Eq. 15, `h = f+ + f- tanh(beta (X - d))`, `beta = 1.6` | P1 | confirmed | Takagi Eq. (15), p. 5 |
| §1.2, step 1.2: the THINC indicator is "Eq. 18" / "Takagi's Eq. 18" | P2 | **corrected**: Takagi's own numbers are Eq. (19) (`zeta`) and Eq. (20) (detection); Takagi's Eq. (18) is the step profile. Liang 2022 Eqs. (17)–(18) restate them verbatim. Now P1 | Takagi p. 6; Liang 2022 p. 5–6 |
| §1.2: DR is Liang 2022 Eqs. 17, 24, 31–33 | P1 | confirmed; the and/or structure of Eqs. (24), (32) re-read in the raw text | Liang 2022 p. 5, 7, 11 |
| §1.2, step 3.1: ST from Liang & Fu §3.1, Fig. 2; `(1,0,1)` → `i-2..i+1`, `(1,1,0)` → `i-1..i+2` | P1 | confirmed by Fig. 2 (rendered). **Caveat**: Eq. (25)'s labels `R'_2`, `R'_3` are swapped against §3.1's `f^2`, `f^3`; rows must be chosen by stencil, never by index (the plan does so) | Fig. 2; Eq. (25) |
| §1.2, step 3.2: RBF from Eqs. 25, 30–32 | P1 | confirmed; every coefficient of Eq. (25) re-derived from the Gaussian RBF | `rbf_check.py` |
| Eqs. 33–36 of Liang & Fu are Liang 2022's 17, 24, 31, 32 unchanged | P1 | confirmed term by term | Liang & Fu p. 8–9 |
| The central flux, Eqs. 25–26, never executes | D | confirmed as derivation: Eq. (25) implies a negative product, and the guard covers "all formulas ... in this section"; 0 counterexamples in 2·10^5 triples. Liang & Fu omit it (P1) | Liang 2022 p. 11 |
| Eq. 25 contains no exponential: the truncation at `lambda^2 dx^2` | P1 | confirmed, with a caveat: "for brevity" — whether the authors' code used the truncated rows is not stated (unsourced). Step 0.4(f) should report both | Liang & Fu p. 6 |
| Every `lambda^2` row sums to zero; `lambda = 0` gives the polynomial flux | D | confirmed (derived, exact) | `rbf_check.py` |
| `R2 = q0/2 + q1/2`, `R3 = 3q0/4 + q2/4` at `lambda = 0` | D | confirmed; they are TENO5-opt's `g2`, `g1` | Fu et al. 2016 Eq. (32) |
| The detector reads `zeta` of `i-3..i+2` (DR), `i-3..i+3` (THINC); data footprint `-5:6`; 6 ghost cells | D | confirmed as reading (Eq. (32) alone reads `i-2..i+2`); the footprint follows: left state `i-5..i+4` (DR) or `i-5..i+5` (THINC), right state mirrored, union `i-5..i+6`. The halo stays D until N6 | Takagi Eq. (20); Liang 2022 Eq. (24) |
| Criterion calibrated on plain TENO5, `tau_5 = |beta_1 - beta_2|`, `C_T = 1e-5`; centred jump cell `(1,1,1)` because `tau_5` vanishes | P1 + D | confirmed: `beta_1 = beta_2 = 5/6`, `tau_5 = 0`. Add to the calibration: `eps = 1e-40`, `q = 6`, `C = 1`, `d = (0.6, 0.3, 0.1)` (all three papers) | Takagi p. 4–6; Liang 2022 p. 4–6 |
| Papers silent on per-kernel cost, ghost cells, non-uniform grids, the basis of a neighbour's flags | U | confirmed for all three; Takagi gives whole-solver cost per step only; activations are per characteristic component and flux sign, the basis of the neighbours' `zeta` unstated | — |
| Liang & Fu Table 1: 1.38, 1.42 times TENO5 | P1 | confirmed (1.379, 1.419) | Table 1 |
| §1.3.1 length-unit dependence | D | supported by quoted text (`lambda` multiplies a length; fixed `alpha`; "degenerates to polynomial as `dx → 0`"); still D until E5 | Eqs. (21), (32), p. 8 |
| §1.3.2 long-wave anti-dissipation, growth `dx^5 k^4 (lambda^2/2 - k^2/60)` | D | the input, Eq. (29), confirmed exactly; the derivation reproduced; still D until E6 | Eq. (29); `rbf_check.py` |
| Step 1.1: bit `k` is candidate `k`, `S0` central, `S1` downwind, `S2` upwind | P1/S | confirmed | all three papers, Fig. 1; `TENO5_CLO` |
| Step 1.1 spec 2: the digits of Tables 1–3 | P1 | confirmed and recomputed; constant continuation reproduces cells 0..10, not only 2..8 (Takagi prints 2..8) | `flags_check.py` |
| Step 1.1 spec 3: mirroring swaps bits 1 and 2 | — | confirmed numerically (Table 3's profile) | `flags_check.py` |
| Step 1.2: `zeta` by Eq. 17 | P1 | confirmed | Liang 2022 Eq. (17) |
| Step 1.2 frozen spec: Table 1 → level 2, Table 2 → 1, Table 3 → 2, all else 0; THINC flags Table 1, misses 2 and 3 | D | confirmed, from the printed tables and from recomputed flags | `flags_check.py`; Liang 2022 p. 7 |
| Step 1.3: (a), (b) at N = 120; (c), (d) at 100 and 400; (e) at 200 and 400, on [0, 2] | P1 | confirmed from the captions (the text's "respectively" contradicts them); (c)'s `F` has a square root; (e) is `sin(2 pi x e^{x+1})`; node placement and ends unsourced | Liang 2022 p. 7–8, Figs. 3–5 |
| Step 1.3: zero activations on (a), (b) | P1 | confirmed | p. 8 |
| Step 1.3 frozen spec: "every jump flagged by the DR detector; every flagged cell within one cell of a jump of the function or of its derivative" | plan's own | **contradicted by the paper's Fig. 5 on (e)**: the jumps at 1.17 and 1.56 are not flagged at either N, and the detector fires at x ≈ 1.93 (N = 400) and at about seven points in 1.58–1.87 (N = 200) in the under-resolved sine. Holds by figure on (c) (including the kink at 1.1) and (d). Correction: assert locality and completeness on (c), (d), and on (e) only for x < 1.5; report (e)'s counts beside Fig. 5 | Fig. 5 at 600 dpi |
| Step 2.1: replace the state by Takagi's face value where the level is 1 | P1 | **incomplete**: Takagi replaces the face by Eq. (16) *and* Eq. (17); the guard `(f_{i+1} - f_i)(f_i - f_{i-1}) <= 1e-15 → f_i` belongs to THINC too (it is Takagi's, inherited by Liang 2022), and it is needed: `alpha` is 0/0 when `f_{i+1} = f_{i-1}`. Spec: add the guard (with `rmu`, as step 2.2) and its G10 mutant | Takagi §2.3, §3.1 step 3 |
| Step 2.1: `B = exp(|alpha| beta)`, `(B-1)/(B+1) = tanh(|alpha| beta/2)`, one `tanh`, no `exp` | D | confirmed: `A = tanh(alpha beta/2)/tanh(beta/2)`; Eq. (16) is exact for the conservation constraint (2e-50). Keep the form odd in `alpha`: Takagi chose it for bitwise mirror symmetry | `thinc_check.py`; Takagi p. 5 |
| N3: "an overridden state lies between the two neighbours' values" | — | needs the interval named: the papers guarantee `[f_{i-1}, f_{i+1}]` (P1); `[f_i, f_{i+1}]` also holds, for THINC and DR, on strictly monotone triples (derived) | Takagi p. 5; `thinc_check.py` |
| Step 2.2: Eq. 31, `A = 0.91` (level 1), `1.23` (level 2), first order by Eq. 33 | P1 | confirmed; level 2 = Eq. (32) true (it implies Eq. (24)) | Liang 2022 p. 11 |
| Step 2.2 specs 2–4: `m + h tanh(A + atanh B)`; linear data → `f_i + d tanh(A)`; denominator positive on monotone triples | D | confirmed (derived; 1e-51) | `thinc_check.py` |
| Step 2.2 spec 6: Table 4's digits; time step not stated | P1 | confirmed; the paper says only "reduced so that the time integration error can be ignored" | Liang 2022 Table 4, p. 12 |
| Step 2.3 early exit: `zeta_i` = 2 or −1 rules detection out; two neighbour masks suffice on smooth data | D | first part confirmed (both detectors need `zeta_i` = 0 or 1); the second holds for DR (`zeta_{i+1}`, `zeta_{i+2}`) but THINC needs three (`i+1..i+3`) | Eqs. (20), (24) |
| Step 3.1: the six conditions overlap, `(0,1,1)` matches two | P1, D | confirmed; the overlap is wider (`(1,0,1)` and `(1,1,0)` each match three). The plan's precedence follows the quoted principle "as many smooth data as possible"; `(0,1,1)` remains a gap (D4) | §3.1, Fig. 2 |
| Step 3.2: `c_r = 0.4`, minimum over `i-1..i+1`, regulariser 1e-40, `g(m)`, `alpha_1 = alpha_2 = 3` | P1 | confirmed | Eqs. (30)–(32) |
| Step 3.2 spec 3: Eq. 29's coefficients | P1 | confirmed exactly from the printed Eq. (25) | `rbf_check.py` |
| §6 table of `s`: Shu–Osher 200 on [0, 10], Sod/Lax 100, double Mach 400×100, KH 400×400, Riemann 800² | P1 | confirmed; missing published grids: blast 400 (`dx` 0.0025), double Mach 1000×250 (0.004), Rayleigh–Taylor 64×256 (0.0039) | Liang & Fu §4 |
| D7: peak `≈ 67 s^6` at `xi^2 = 20 s^2`, validated range near `s = 0.12` | D | arithmetic confirmed (`200/3 s^6`); still D. The paper's own ADR (Fig. 3) cannot be compared: its `s` is not stated | — |
| Step 3.4: Liang & Fu print no table of results; Figs. 5–8 are Sod, Lax, blast, Shu–Osher | P1 | confirmed (Table 1 is cost only) | Liang & Fu §4 |
| §8: optimal-shape-parameter RBF-WENO unstable at critical points per Liang & Fu's introduction | P1 | confirmed | Liang & Fu p. 2 |
| G7 on the blast waves "TENO5's bounds" | — | new information: the papers' TENO5 runs at `C_T = 1e-4` there; the composite's `C_T` is unsourced | Liang 2022 §4.3; Liang & Fu §4.3 |
