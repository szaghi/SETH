# Schemes

SETH computes one thing: the two values of a field at a face $x_{i+1/2}$, one from the left and one from the right,
from the values around it, without oscillating where the field jumps. This page sets out where those face values are
needed, how a weighted ENO scheme builds them, and then each family SETH implements, with its formulas. The first
three sections are the setting, which the caller's solver implements; from [Smoothness
indicators](#smoothness-indicators) on, every formula is implemented as written. Where each one comes from, and
what SETH derives from it, is in the repository's `notes/m2-sources.md`.

## Notation

| Symbol | Meaning |
|---|---|
| $I_j = [x_{j-1/2}, x_{j+1/2}]$, $h$ | cell $j$, of width $h$; $x_j$ its centre, $x_{j+1/2}$ the face between cells $j$ and $j+1$ |
| $\bar u_j = \frac1h \int_{I_j} u\,dx$ | cell average |
| $u_j = u(x_j)$ | point value |
| $u^-_{i+1/2}$, $u^+_{i+1/2}$ | the left and the right state at the face (`rp`, `rm`) |
| $S$ | the number of candidate stencils of WENO-$(2S-1)$; each candidate has $S$ cells |
| $q_k$ | the value at the face of the polynomial of candidate $k$ |
| $d_k$ | the linear (optimal) weight of candidate $k$ |
| $\beta_k$ | the smoothness indicator of candidate $k$ |
| $\omega_k$ | the nonlinear weight of candidate $k$ |
| $\tau$ | the global smoothness measure of WENO-Z and TENO |
| $\varepsilon$ | the regulariser of the weights (`eps`) |

## Where face values are needed

A conservation law $\partial_t u + \partial_x f(u) = 0$ is discretised in conservative form,

$$ \frac{d}{dt}\,U_i = -\frac{1}{h}\left(\hat f_{i+1/2} - \hat f_{i-1/2}\right), $$

and the whole accuracy and robustness of the scheme is in the numerical flux $\hat f_{i+1/2}$, which needs values at
the faces. Three formulations get them in three ways, and SETH serves all three.

**Finite volumes.** $U_i = \bar u_i$, the cell averages. The face states $u^\mp_{i+1/2}$ are **reconstructed** from
the averages around the face, and a Riemann solver gives $\hat f_{i+1/2} = F(u^-_{i+1/2}, u^+_{i+1/2})$. In one
dimension this is high order as it stands; in several, the face flux also needs a quadrature over the face.

**Conservative finite differences** (Shu & Osher 1989). $U_i = u_i$, the point values. If $\mathcal H$ is defined by
$f(u(x)) = \frac1h \int_{x-h/2}^{x+h/2} \mathcal H(\xi)\,d\xi$, then $\partial_x f(u(x_i)) = \big(\mathcal
H(x_{i+1/2}) - \mathcal H(x_{i-1/2})\big)/h$ exactly. The point values $f(u_j)$ are therefore the *cell averages*
of $\mathcal H$, and the same reconstruction, applied to $f(u_j)$ in place of $\bar u_j$, gives $\hat f_{i+1/2}
\approx \mathcal H(x_{i+1/2})$ to the order of the reconstruction, dimension by dimension. Upwinding comes from a
flux splitting $f = f^+ + f^-$ with $df^+/du \ge 0 \ge df^-/du$, for instance Lax–Friedrichs, $f^\pm = \tfrac12(f
\pm \alpha u)$ with $\alpha \ge \max|f'(u)|$:

$$ \hat f_{i+1/2} = \underbrace{\big(f^+\big)^-_{i+1/2}}_{\text{left state of } f^+} +
\underbrace{\big(f^-\big)^+_{i+1/2}}_{\text{right state of } f^-}, $$

which in SETH is `vp` $= f^+$, `vm` $= f^-$, $\hat f$ = `rp + rm`. For systems the splitting is done in local
characteristic fields.

**The alternative finite-difference form** (Shu & Osher 1988; Jiang, Shu & Zhang 2013). $U_i = u_i$ again, but the
face states are **interpolated** from the point values $u_j$, any Riemann solver is applied at the face, and even
derivatives of the flux correct the difference to high order:

$$ \hat f_{i+1/2} = F\!\left(u^-_{i+1/2}, u^+_{i+1/2}\right) - \frac{h^2}{24}\,\partial_x^2 f\Big|_{i+1/2}
+ \frac{7h^4}{5760}\,\partial_x^4 f\Big|_{i+1/2} - \dots $$

It admits any Riemann solver and keeps free streams on curvilinear grids, which the splitting form does not; the
weighted compact nonlinear schemes (WCNS) interpolate the same way.

| Formulation | Data | Kernel | Kernel output |
|---|---|---|---|
| finite volume | cell averages $\bar u_j$ | reconstruction: `seth_<scheme>` | $u^\mp_{i+1/2}$, then a Riemann solver |
| conservative FD | $f^\pm(u_j)$ at the points | reconstruction: `seth_<scheme>` | $\hat f^{\pm}$, summed to $\hat f_{i+1/2}$ |
| alternative FD, WCNS | point values $u_j$ | interpolation: `seth_<scheme>_interp` | $u^\mp_{i+1/2}$, then a Riemann solver and the corrections |

The [verification harness](/verification/shocks/sod#models-and-numerics) runs both finite-difference paths with
every scheme.

## Candidate stencils and linear weights

For the left state at $x_{i+1/2}$, WENO-$(2S-1)$ uses the $2S-1$ cells $i-S+1, \dots, i+S-1$, split into
$S$ **candidate stencils** of $S$ cells,

$$ S_k = \{\, i-S+1+k, \dots, i+k \,\},\qquad k = 0, \dots, S-1. $$

On each candidate, the polynomial $p_k$ of degree $S-1$ is fitted to the data, and its value at the face is a fixed
linear combination of the data:

- **reconstruction** (cell averages): $\frac1h \int_{I_j} p_k\,dx = \bar u_j$ for $j \in S_k$;
- **interpolation** (point values): $p_k(x_j) = u_j$ for $j \in S_k$;

$$ q_k = p_k(x_{i+1/2}) = \sum_{l=0}^{S-1} c_{kl}\,v_{i-S+1+k+l} = u(x_{i+1/2}) + O(h^S) $$

for smooth data, $v$ the averages or the point values. The polynomial of degree $2S-2$ on the whole stencil gives the
value of order $2S-1$, and it is itself a combination of the candidates, with positive **linear weights**:

$$ P(x_{i+1/2}) = \sum_{k=0}^{S-1} d_k\, q_k = u(x_{i+1/2}) + O(h^{2S-1}),\qquad \sum_k d_k = 1. $$

For WENO5 ($S = 3$, cells $i-2..i+2$):

| | reconstruction (cell averages) | interpolation (point values) |
|---|---|---|
| $q_0$, cells $i-2..i$ | $\frac16(2\bar u_{i-2} - 7\bar u_{i-1} + 11\bar u_i)$ | $\frac18(3u_{i-2} - 10u_{i-1} + 15u_i)$ |
| $q_1$, cells $i-1..i+1$ | $\frac16(-\bar u_{i-1} + 5\bar u_i + 2\bar u_{i+1})$ | $\frac18(-u_{i-1} + 6u_i + 3u_{i+1})$ |
| $q_2$, cells $i..i+2$ | $\frac16(2\bar u_i + 5\bar u_{i+1} - \bar u_{i+2})$ | $\frac18(3u_i + 6u_{i+1} - u_{i+2})$ |
| $d_0, d_1, d_2$ | $\frac1{10}, \frac6{10}, \frac3{10}$ | $\frac1{16}, \frac{10}{16}, \frac5{16}$ |

The right state $u^+_{i+1/2}$ is the same construction on the stencil mirrored about the face. SETH takes it from
`vm` and mirrors internally ([Conventions](/guide/conventions), with figures of the stencils of every scheme).

The linear combination $\sum d_k q_k$ is a fixed-stencil scheme: of high order, and oscillatory at a discontinuity
(an overshoot of $O(1)$ that does not shrink with $h$).

## Nonlinear weights

A weighted ENO scheme replaces the linear weights by **nonlinear** ones, computed from the data at every face:

$$ u^-_{i+1/2} = \sum_{k=0}^{S-1} \omega_k\, q_k,\qquad \omega_k \ge 0,\qquad \sum_k \omega_k = 1. $$

Two requirements fix them:

- **Accuracy.** Writing $u^-_{i+1/2} = \sum_k d_k q_k + \sum_k (\omega_k - d_k)(q_k - u)$, the scheme keeps order
  $2S-1$ on smooth data if $\omega_k = d_k + O(h^{S-1})$ (Jiang & Shu 1996).
- **ENO property.** A candidate whose stencil crosses a discontinuity must get a weight that vanishes as $h \to 0$,
  so that the value is built from the smooth candidates only.

Both rest on a **smoothness indicator** $\beta_k$ per candidate: $O(h^2)$ on smooth data, $O(1)$ across a jump. The
families differ in how $\beta_k$ becomes $\omega_k$:

| Family | Weights | On smooth data |
|---|---|---|
| WENO-JS (Jiang & Shu 1996) | $\omega_k \propto d_k\,/\,(\varepsilon + \beta_k)^{w}$ | $\omega_k \to d_k$ as $h \to 0$; slower at critical points |
| WENO-Z (Borges et al. 2008) | $\omega_k \propto d_k\,[1 + (\tau/(\beta_k + \varepsilon))^{p}]$, $\tau$ a higher-order combination of the $\beta_k$ | closer to $d_k$: less dissipation, better at critical points |
| TENO (Fu, Hu & Adams 2016) | a cutoff: candidate $k$ keeps $d_k$ or gets 0 | exactly $d_k$: the linear scheme, unchanged |

The next sections give each in full.

## Smoothness indicators

For a candidate of $w$ cells reconstructed by the polynomial $p_k$ of degree $w-1$ (Jiang & Shu 1996):

$$ \beta_k = \sum_{j=1}^{w-1} h^{2j-1} \int_{x_{i-1/2}}^{x_{i+1/2}} \left(\frac{d^j p_k}{dx^j}\right)^2 dx. $$

For WENO5 this gives the familiar $\beta_0 = \frac{13}{12}(\bar u_{i-2} - 2\bar u_{i-1} + \bar u_i)^2 +
\frac14(\bar u_{i-2} - 4\bar u_{i-1} + 3\bar u_i)^2$, and so on.

SETH derives every indicator from this definition in exact rational arithmetic. WENO-Z and TENO evaluate it as a sum
of squares, $\beta_k = \sum_i D_i \big(\sum_a L_{ai} u_a\big)^2$ with $D_i \ge 0$ (an exact $LDL^T$), non-negative in
floating point by construction. WENO-JS keeps ADAM's expanded form, whose WENO-7 and WENO-9 indicators are
Balsara & Shu's (2000) integer forms: 240 and 5040 times the normalised ones (`WENO7_JS_D_SCALE`,
`WENO9_JS_D_SCALE`).

## WENO-JS: WENO-1/3/5/7/9

$S$ candidates of $S$ cells (order $2S-1$), weights

$$ \alpha_k = \frac{d_k}{(\varepsilon + \beta_k\,\mu^{-2})^{w}},\qquad \omega_k = \frac{\alpha_k}{\sum_j \alpha_j},
\qquad w = S\ (S \le 4),\ w = 4\ (S = 5), $$

with ADAM's exponent $w$ and $\varepsilon = 10^{-6}$ (`SETH_WENO_JS_EPS`). $\mu^{-1}$ is the argument `rmu`: `1` for
the Jiang–Shu weights, or the inverse of a magnitude of the data for scale-invariant weights (Don et al. 2022;
`seth_weno_descaling` gives ADAM's generic one). WENO-1 is the first-order upwind value, $u^-_{i+1/2} = v_i$.

- **Scale invariance.** With `rmu = 1` the weights depend on the scale of the data (an absolute $\varepsilon$): data of
  magnitude $10^{-4}$ is reconstructed differently from the same data at $1$. With the descaler the result scales
  bitwise under $v \to 2^k v$.
- **Affine invariance.** The mean $|v|$ fails for small fluctuations on a large offset, the usual case of a pressure or
  a density: $\mu$ is then the offset, $\varepsilon$ swamps the indicators, the weights become the linear ones and the
  scheme oscillates. `seth_weno_descaling_ai` takes $\mu$ = the mean deviation from the stencil's mean (Wang & Don 2022,
  Eq. (13)), so that the weights ignore the offset too. Measured on a step of height $10^{-7}$ on the offset 1: WENO5-JS
  overshoots by 18% of the jump with the mean-$|v|$ descaler and by $1.6\cdot10^{-6}$ (gfortran) to $3\cdot10^{-4}$
  (nvfortran, fused multiply-adds) with the affine-invariant one.
  The descaler alone depends on that round-off, and for WENO7/9-JS it is not enough (0.2% to 5% overshoot by compiler):
  their expanded integer-form indicators cancel the offset only to its round-off, which $\mu^{-2}$ amplifies.
  Subtracting the stencil's mean from the data first and adding it back to the states (Wang & Don's Remark 2) gives no
  overshoot at all, with every compiler: that is the recommended use. The cost is one more pass over the
  stencil: WENO5-JS with it 1.336 ns per face on the RTX 4070, with the mean-$|v|$ descaler 1.105.
- **Range.** $(\varepsilon + \beta)^{w}$ overflows double precision near $|v| \sim 10^{77}$ (WENO-3), $10^{51}$
  (WENO-5), $10^{38}$ (WENO-7/9); with the descaler the limit is $|v| \sim 10^{150}$ ($\beta$ itself).
- **Critical points.** Where $u' = 0$ the $\beta_k$ differ at relative $O(h)$, $\omega_k - d_k$ is larger than
  $O(h^{S-1})$ and WENO5-JS drops to third order there. WENO-Z and TENO were designed against this
  ([Critical points](/verification/accuracy/critical-points)).

The `_wr` kernels also return the **weights ratio** $\min_k \omega_k/d_k \in [0, 1]$ of each state: 1 on smooth data,
toward 0 near a discontinuity.

## WENO-Z: WENO-Z3/5/7/9

The candidates and linear weights of WENO-JS, the weights of Borges et al. (2008):

$$ \alpha_k = d_k\left[1 + \left(\frac{\tau}{\beta_k + \varepsilon}\right)^{p}\right],\qquad
\omega_k = \frac{\alpha_k}{\sum_j \alpha_j}. $$

$\tau$ is a combination of the $\beta_k$ in which their common leading term cancels, so that $\tau = O(h^{S+1})$ or
higher on smooth data while each $\beta_k = O(h^2)$: $\tau/\beta_k \to 0$ and the weights approach $d_k$ faster than
WENO-JS's. Across a jump, $\tau = O(1)$ and the smooth candidates, with small $\beta_k$, take the weight.

| Scheme | $\tau$ | Source |
|---|---|---|
| WENO-Z3 | $\lvert\beta_0 - \beta_1\rvert$ | Don & Borges; Li et al. 2021 |
| WENO-Z5 | $\lvert\beta_0 - \beta_2\rvert$ | Borges et al. 2008 |
| WENO-Z7 | $\lvert\beta_0 + 3\beta_1 - 3\beta_2 - \beta_3\rvert$, $O(h^7)$ | Don & Borges, as quoted by Rathan & Naga Raju 2016 |
| WENO-Z9 | $\lvert\beta_0 - \beta_4\rvert$, $O(h^7)$ | Castro et al.'s rule for odd $r$, as restated by Li et al. 2022 |

$p$ = 1 (Borges) or 2 (fifth order at first-order critical points); $\varepsilon = 10^{-40}$ is the usual choice.
What SETH's verification runs, and recommends, is measured on its shock cases ($L_1(\rho)$ against the reference or
the exact solution; the overshoot behind the Burgers shock; limited faces on the blast wave):

| Scheme | $p$ | Why | $p = 1$ against $p = 2$ |
|---|---|---|---|
| WENO-Z3, WENO-Z5 | 1 | Borges et al.'s choice; the low-dissipation claims of the [Sod](/verification/shocks/sod) and [Lax](/verification/shocks/lax) cases are asserted with it | — |
| WENO-Z7 | 1 | $p = 2$ is more robust but more dissipative everywhere | composite profile $L_1$ 0.0479 / 0.0510, Shu–Osher 0.142 / 0.180, Sod 0.00242 / 0.00280, blast 0.074 / 0.090 with 1826 / 23 limited faces, Burgers overshoot $6\cdot10^{-3}$ / $9\cdot10^{-7}$ |
| WENO-Z9 | 2 | with $p = 1$ it is not robust at shocks | Burgers overshoot 0.28 / $9\cdot10^{-4}$, composite profile within $[-0.055, 1.056]$ / $[-8\cdot10^{-5}, 1.0005]$, blast 6294 / 585 limited faces |

Castro et al. (2011) recommend $p = r - 1$ for robustness at high order and report $p \le 2$ unstable on the double
Mach reflection for $r > 4$ (llm-wiki high-order-capturing-schemes `methods/weno-z.md`). SETH's cases are
one-dimensional, so that claim is neither confirmed nor refuted here. The orders of $\tau$ are measured exactly by the table generator (Castro et al.'s
$\tau_7 = \lvert\beta_0-\beta_1-\beta_2+\beta_3\rvert$ is $O(h^6)$, as Li et al. 2022 state).

## TENO: TENO5, TENO5-opt, TENO5-AZ, TENO5-A, TENO6, TENO6-opt, TENO6-NENO6, TENO8, TENO8-A

TENO (Fu, Hu & Adams 2016) changes two things. Its candidates have **incremental widths**: three-cell candidates
around the face and, for TENO6, a fourth one of four cells that reaches downwind, so that the full stencil can carry
a linear scheme of a chosen order and spectral behaviour. TENO5's linear weights give exactly WENO5's linear scheme;
TENO6's give a sixth-order one, and TENO5-opt, TENO6-opt and NENO6 trade order for lower dispersion and dissipation. And its
weights are a **cutoff**: each candidate either keeps its linear weight or is dropped, so on smooth data the scheme
*is* its linear scheme, with that scheme's spectral properties, not an approximation of it.

Candidates: $S_0 = \{i-1, i, i+1\}$, $S_1 = \{i, i+1, i+2\}$, $S_2 = \{i-2, i-1, i\}$, and for TENO6 the downwind
$S_3 = \{i, \dots, i+3\}$.

$$ \gamma_k = \left(1 + \frac{\tau_K}{\beta_k + \varepsilon}\right)^{6},\quad
   \chi_k = \frac{\gamma_k}{\sum_j \gamma_j},\quad
   \delta_k = \begin{cases} 0 & \chi_k < C_T \\ 1 & \text{otherwise} \end{cases},\quad
   \omega_k = \frac{d_k \delta_k}{\sum_j d_j \delta_j} $$

with $\tau_5 = \lvert\beta_1 - \beta_2\rvert$ and $\tau_6 = \beta_6 - \tfrac16(\beta_1 + \beta_2 + 4\beta_0)$
($\beta_6$ on the six-point stencil), $\varepsilon = 10^{-40}$ (`SETH_TENO_EPS`). The large power separates the
scales: a candidate across a discontinuity has $\chi_k$ many orders below the others and falls under the cutoff
$C_T$, while smooth candidates stay far above it. SETH evaluates $\chi$ as $(g_k/g_{\max})^6 / \sum (g_j/g_{\max})^6$,
$g = 1 + \tau/(\beta+\varepsilon)$: the same values, without the overflow of $\gamma$ on large jumps.

| Kernel | $d_k$ | Order | $C_T$ constant | Source |
|---|---|---|---|---|
| `seth_teno5` | 0.6, 0.3, 0.1 | 5 | `SETH_TENO5_CT` = $10^{-5}$ | Fu et al. 2016, Table 5 |
| `seth_teno5_opt` | 0.55, 0.4, 0.05 | 4 | `SETH_TENO5_OPT_CT` = $10^{-5}$ | Fu et al. 2016, Tables 3 and 5 |
| `seth_teno5_az` | 0.6, 0.3, 0.1 | 5 | `SETH_TENO5_AZ_CT` = $10^{-5}$, the ceiling $\lambda_d$ | Tang 2024, Eqs. (19)–(20) |
| `seth_teno6` | 0.45, 0.3, 0.05, 0.2 | 6 | `SETH_TENO6_CT` = $10^{-7}$ | Fu et al. 2016, Table 5 |
| `seth_teno6_opt` | 0.462, 0.3, 0.054, 0.184 | 5 | `SETH_TENO6_OPT_CT` = $10^{-7}$ | Fu et al. 2016, Table 5 |
| `seth_teno6_neno6` | 0.42943, 0.17273, 0.08557, 0.31227 | 4 | `SETH_TENO6_NENO6_CT` = $10^{-6}$ | the NENO6 paper |
| `seth_teno8` | 0.41309, 0.21931, 0.06459, 0.14392, 0.02747, 0.13162 | 6 | `SETH_TENO8_CT` = $10^{-7}$ | Liang & Fu 2024b; $C_T$: Hong et al. 2022 |

**The optimised weights.** Fu et al. (2016, §4, Eq. (31)) build the background of order $K-1$ as
$\eta\, g^{(1)} + (1-\eta)\, g^{(2)}$, two schemes of order $K-1$ on the $K$-point stencil: for $K = 5$ the upwind
$\tfrac1{12}(f_{i-2} - 5f_{i-1} + 13f_i + 3f_{i+1})$ and the central $\tfrac1{12}(-f_{i-1} + 7f_i + 7f_{i+1} - f_{i+2})$
(Eq. (32)); for $K = 6$ two fifth-order ones that share the dispersion of the sixth-order central scheme (Eq. (33)).
$\eta = 0.4$ gives back the upwind fifth-order scheme; smaller $\eta$ moves toward the central one and removes
dissipation. The table generator solves for the candidates' weights that reproduce that combination: $\eta = 0.2$
gives exactly $\tfrac{11}{20}, \tfrac25, \tfrac1{20}$ and $\eta = 0.54$ exactly $0.462, 0.3, 0.054, 0.184$, the
values of Fu et al.'s Tables 3 to 5, which it asserts. TENO5-opt's linear scheme is fourth order, measured 3.99–4.01,
and stays fourth order at a first-order critical point.

**TENO5-AZ: a cutoff that vanishes on smooth data.** Tang (2024) keeps TENO5's candidates and weights and changes
two things: the global measure, $\tau = \lvert\beta_0 + \beta_2 - 2\beta_1\rvert$ in WENO labels (upwind, central,
downwind), $\lvert\beta_1 + \beta_2 - 2\beta_0\rvert$ in TENO's, measured $O(h^4)$ by the generator, and the cutoff,
$C_T = \min(\lambda_d, \tau)$ with $\lambda_d = 10^{-5}$ (Eqs. (19)–(20)). Where the data is smooth $\tau \to 0$, so
$C_T \to 0$ and no candidate is dropped: the scheme keeps fifth order at critical points of any order, where TENO5
falls to second. SETH measures 5.00 at first- and second-order critical points (TENO5: 5.00 and 2.00), as Tang's
Tables 1–3 report.

The price is dimensional: $\tau$ scales with the square of the data, a pure number does not, so the cutoff depends
on the units. `seth_teno5_az` takes an inverse descaler `rmu` and compares $\tau\,\mu^{-2}$: `rmu = 1` is Tang's
scheme, which is what the verification runs (it is asserted not scale invariant, like plain WENO-JS);
`rmu = seth_weno_descaling(3, 1, vp, vm)` makes it bitwise invariant under scaling by powers of 2, and also takes the
property away, since $\tau/|v|^2$ does not vanish at a critical point: order 2 at a second-order one, measured. No
choice of `rmu` keeps both.

**TENO8: eight points, a sixth-order dissipative background.** Six candidates on the cells $i-3 \dots i+4$: the
three of TENO5, $S_3 = \{i, \dots, i+3\}$, $S_4 = \{i-3, \dots, i\}$ and $S_5 = \{i, \dots, i+4\}$ (Fu et al.
2016, Table 1), every one of them with an end at $x_i$ and at least one upwind point; $\tau_8 = \lvert\beta_8 -
\tfrac16(\beta_1 + \beta_2 + 4\beta_0)\rvert$, $O(h^6)$, with $\beta_8$ on the whole stencil. The linear weights
are Liang & Fu's (2024b, §3.2): a sixth-order scheme whose dispersion is optimised and whose dissipation is an
explicit eighth-order term ($b_0 = 0.01$), dissipative at every wavenumber (measured: Im $\tilde\xi \le 5\cdot10^{-17}$,
resolution 0.547 against 0.42 for WENO7-JS). The weights are printed to 16 digits from a numerical optimisation;
the generator checks that they reproduce the paper's face coefficients (Eq. (46)) within $1.6\cdot10^{-15}$, that
the linear scheme is sixth order, that its own $\beta_8$ equals the paper's printed $\tau_8$ term by term, and that
the max-order weights of these candidates are Fu et al.'s $(30, 18, 4, 12, 1, 5)/70$ (Table 2).

Its nonlinear scheme is **not** free of the anti-dissipation of TENO6: the [ADR](/verification/spectral/adr) is
positive, up to $3.2\cdot10^{-3}$, for $\vartheta/\pi = 0.36$–$0.45$, where the cutoff starts to drop candidates, and
more so for smaller $C_T$ ($2.2\cdot10^{-3}$ at $10^{-5}$, $9.6\cdot10^{-3}$ at $10^{-10}$). The background's
dissipation does not reach that band. Advected for four periods, no mode of the band grows (its amplitude never
exceeds the initial one): the positive ADR is a one-step measure, not a growth, for TENO6 and TENO8 alike.

**TENO5-A and TENO8-A: the cutoff from a sensor.** Fu, Hu & Adams (2019) let the cutoff follow the local
smoothness (Fu 2021, Eqs. (17)–(19); Liang & Fu 2024b, Eqs. (47)–(49)):

$$ \eta_j = \frac{\lvert 2\Delta_{j+1/2}\Delta_{j-1/2}\rvert + e}{\Delta_{j+1/2}^2 + \Delta_{j-1/2}^2 + e},\quad
   m = 1 - \min\!\Big(1, \frac{\min_j \eta_j}{C_r}\Big),\quad g = (1-m)^4 (1 + 4m),\quad
   C_T = 10^{-\lfloor a_1 - a_2 (1 - g) \rfloor}, $$

$\Delta_{j+1/2} = f_{j+1} - f_j$, $e = \frac{0.9 C_r}{1 - 0.9 C_r}\,\xi^2$, the minimum over the nodes $i-1 \dots i+1$
(TENO5-A) or $i-1 \dots i+2$ (TENO8-A). Smooth monotone data has $\eta \to 1$, $m = 0$ and the smallest cutoff
$10^{-a_1}$; a jump has $\eta \to 0$, $m = 1$ and $10^{-(a_1 - a_2)}$.

| Kernel | Background | $C_r$, $\xi$, $a_1$, $a_2$ | $C_T$ range | Source |
|---|---|---|---|---|
| `seth_teno5_a` | TENO5 | 0.24, $10^{-3}$, 10, 5 | $10^{-10}$ … $10^{-5}$ | Tang 2024, Eqs. (15)–(16); Huang, Liang & Fu 2023, Eqs. (18)–(20) |
| `seth_teno8_a` | TENO8 | 0.23, $10^{-3}$, 10.5, 3.5 | $10^{-10}$ … $10^{-7}$ | Liang & Fu 2024b, Eqs. (47)–(49) |

With $a_2 = 0$ either kernel is the fixed-cutoff one at $10^{-\lfloor a_1 \rfloor}$, bitwise (asserted). What SETH
measures:

- **The sensor reads extrema as jumps.** At a crest one difference vanishes and $\eta \to 0$, so on a sine wave
  TENO5-A keeps $C_T = 10^{-5}$ and its [ADR](/verification/spectral/adr) is TENO5's to the last digit; at a
  second-order critical point it falls to second order like TENO5. Li, Chen & Ren (2022) make the same observation
  of Ren's sensor.
- **TENO8-A keeps sixth order at a second-order critical point** (TENO8: 2), its eight-point stencil and four-node
  minimum leaving the cutoff at its smallest there; and its ADR is anti-dissipative like TENO8's, more so (Im up to
  $9.6\cdot10^{-3}$ at $\vartheta/\pi = 0.38$–$0.51$, the band of a $10^{-10}$ cutoff); no mode of that band grows.
- **Not scale invariant:** $e$ is absolute. Both kernels take an `rmu` (as TENO5-AZ) that scales $\xi$ by the data's
  magnitude; `rmu = 1` is the papers' form and what the verification runs.
- **Cost** (RTX 4070): TENO5-A 1.41 times TENO5, over the budget of 1.25 set before coding (three divisions and a
  power per state on a kernel that does little else); TENO8-A 1.16 times TENO8.

"TENO6-opt" denotes different weights in different papers; SETH's is the founding paper's. Lowering $C_T$ keeps more
candidates (less dissipation, less robustness); the defaults are the papers'.

## Centred WENO: WENO4/6/8-C

ADAM's centred WENO (`WENO_C_4/6/8`), ported with its defect D1 fixed. With $S$ candidates, on the big stencil of
cells $i+1-S..i+S$ around the face, it blends two reconstructions of the same face value:

$$
u^-_{i+1/2} = r_c + \phi\,(r_{\rm up} - r_c),\qquad
\phi = \max_{k=1,2} \frac{\lvert I_k - I_{k-1}\rvert}{I_k + I_{k-1}},\qquad
I_k = \big(r_\mu (\bar u_{i+k} - \bar u_{i+k-1})\big)^2 + \varepsilon,\ k = 0, 1, 2,
$$

with $r_{\rm up}$ the WENO-JS of order $2S-1$ (`seth_weno{3,5,7}_js`, on its own left-biased cells) and $r_c$ the
central reconstruction of order $2S$ on all $2S$ cells:

| Kernel | $S$ | Upwind part | Central row (cells $i+1-S..i+S$) | Order |
|---|---|---|---|---|
| `seth_weno4_c` | 2 | WENO3-JS | $\frac1{12}(-1, 7, 7, -1)$ | 4 |
| `seth_weno6_c` | 3 | WENO5-JS | $\frac1{60}(1, -8, 37, 37, -8, 1)$ | 6 |
| `seth_weno8_c` | 4 | WENO7-JS | $\frac1{840}(-3, 29, -139, 533, 533, -139, 29, -3)$ | 8 |

On smooth data $\phi = O(h)$ and the error is $O(h^{2S})$: measured 4.00, 6.03, 7.99 at a point
([test programs](/verification/unit-tests)). Where a difference of the cells $i-1..i+2$ vanishes (a critical
point), $\phi = O(1)$ and the state is the upwind one: the orders there are WENO-JS's, 1, 3, 5 at first-order
critical points. The right state is the same on the mirrored stencil. $r_\mu$ (`rmu`) is the descaler of the
upwind weights and of the $I_k$ alike: $r_\mu = 1$ gives ADAM's indicators. $\varepsilon > 0$ (`eps`,
`SETH_WENO_JS_EPS`) regularises both. In R8P $\phi$ takes one division: the larger ratio is chosen by cross
products, $|I_1 - I_0|(I_2 + I_1) \ge |I_2 - I_1|(I_1 + I_0)$, and only it is divided; the value is ADAM's except
where the two ratios tie to rounding (then the other, a few ulps away), and it halves the blend's divisions
([Performance](/guide/performance)). In R4P the products would overflow at large jumps, and $\phi$ keeps ADAM's two
divisions.

What it is, and what it is not:

- **D1 fixed.** ADAM computes the left state's weights on the unshifted cells, one cell off (the audit's D1). SETH's
  upwind part is its WENO-JS kernel, which takes the face-offset stencil, so the fix is by construction;
  `seth_test_weno_c_adam` shows SETH equal to ADAM's routine with that one line corrected (to $6\cdot10^{-16}$ of the
  data), and ADAM's as it is different on the left state for most random stencils.
- **The detector reads cells $i-1..i+2$ for every $S$.** WENO8-C's outer cells $i-3$ and $i+4$ are outside it: a step
  there leaves $\phi = 0$, and the central row's edge coefficient $-\tfrac3{840}$ passes, an overshoot of $\tfrac1{280}$
  of the jump. WENO6-C's edge coefficient $\tfrac1{60}$ is positive: a precursor of the jump inside the data's
  range, not an overshoot.
- **Behind a shock, smooth regions converge slowly.** After the Burgers shock forms, WENO6-C and WENO8-C leave a
  grid-scale packet upstream of it, and their smooth-region error converges at first order or below, about
  $3\cdot10^{-6}$ at 640 cells where their upwind parts reach $10^{-10}$ and below
  ([S6](/verification/shocks/burgers)). The blend damps a small disturbance on a smooth flow at $\phi$ times the
  upwind rate, and $\phi$ there is the flow's, $O(h)$: what a shock leaves behind is damped at a rate that does not
  grow with the resolution. WENO8-C also loses positivity on the blast wave
  ([S5](/verification/shocks/blast)).
- **No paper.** The blend is ADAM's own construction; its pair term is $\sqrt{1-\psi^2}$ of Ren et al.'s smoothness
  ratio $\psi = 2\lvert ab\rvert/(a^2+b^2)$ (for $\varepsilon = 0$, an identity): a continuous form of the switch of the
  hybrid schemes ([Hybrids](#hybrids)). ADAM's $S = 1$ (`weno-c-2`) never reaches the blend and is WENO1: not ported.

## Interpolation with WENO-Z and TENO

`seth_wenoz<N>_interp` and `seth_teno<5|6>_interp` interpolate the face value from point values at the cell centres,
for the alternative finite-difference formulation and WCNS-type schemes ([above](#where-face-values-are-needed)).
Only two things change from the reconstruction kernels: the candidates' face coefficients, those of the Lagrange
interpolant of their points, and the max-order linear weights. Everything else is shared, the body included: the
indicators, $\tau$, the mirroring and the cutoff.

| Scheme | Interpolation linear weights | Order | Source |
|---|---|---|---|
| WENO-Z3/5/7/9 | those of WENO-JS interpolation: 1/4, 3/4; 1/16, 5/8, 5/16; ... | 3, 5, 7, 9 | equal to ADAM's `initialize_interpolation_S*` |
| TENO5 | 5/8, 5/16, 1/16 | 5 | Hiejima (2022), Chu, Kurganov & Na (2022) |
| TENO6 | 15/32, 5/16, 1/32, 3/16 | 6 | derived (positive; no published set) |

The indicators are the reconstruction ones applied to point values. That choice keeps every $\tau$ order (measured
exactly on a generic polynomial): integrating the derivatives of the point-value interpolant instead would drop
WENO-Z7's $\tau$ from 7 to 5 and TENO6's from 6 to 4 (`notes/interpolation-plan.md`). Verification: every page of
the [Verification](/verification/) section covers the interpolations, and the [123](/verification/shocks/test123)
and [blast](/verification/shocks/blast) cases show that the WENO-Z and TENO interpolations need a positivity
safeguard near vacuum and on strong shocks (`seth_positivity_euler`, [API](/guide/api#positivity-limiter)).

## Monotonicity-preserving bound

`seth_mp_bound` bounds the two states of a face, after any reconstruction, by the interval of Suresh & Huynh (1997,
Eqs. (2.8)–(2.29)), from the five values $v_{j-2} \dots v_{j+2}$ around the upwind cell:

$$ d_k = v_{k-1} - 2v_k + v_{k+1},\quad
   d^{M4}_{j\pm1/2} = \operatorname{minmod}(4d_j - d_{j\pm1},\ 4d_{j\pm1} - d_j,\ d_j,\ d_{j\pm1}), $$
$$ v^{UL} = v_j + \alpha(v_j - v_{j-1}),\quad v^{MD} = \tfrac12(v_j + v_{j+1}) - \tfrac12 d^{M4}_{j+1/2},\quad
   v^{LC} = v_j + \tfrac12(v_j - v_{j-1}) + \tfrac43 d^{M4}_{j-1/2}, $$
$$ v^{\min} = \max\!\big(\min(v_j, v_{j+1}, v^{MD}),\ \min(v_j, v^{UL}, v^{LC})\big),\quad
   v^{\max} = \min\!\big(\max(v_j, v_{j+1}, v^{MD}),\ \max(v_j, v^{UL}, v^{LC})\big), $$

and the state becomes the median of itself, $v^{\min}$ and $v^{\max}$. On monotone data the interval is
$[v_j, v_{j+1}]$; near a smooth extremum the curvature terms widen it so that an accurate state is left alone.
$\alpha = 4$ (`SETH_MP_ALPHA`), Suresh & Huynh's value for Runge–Kutta stepping, whose monotonicity proof needs CFL
$\le 1/(1+\alpha) = 0.2$. Suresh & Huynh skip the bound where the state already lies in a first-order interval, with
an absolute tolerance $10^{-10}$ (Eq. (2.30)); SETH always applies the median, which there returns the state
unchanged: no branch, and scale invariant.

Measured (`seth_test_mp`): on monotone data the bounded WENO5-JS states stay in $[v_0, v_1]$; the steps that WENO9-JS,
WENO-Z9 with $p = 1$ and TENO6 reconstruct end in $[0, 1]$ exactly; on $\sin^4(\pi x)$, the test where Balsara & Shu
(2000) report order loss for MPWENO, the bound never acts on WENO5-JS's states for $N = 20 \dots 1280$ (order 5.3).
Cost on the RTX 4070: 0.42 ns per face, 0.43 times WENO5-JS (over the budget of 0.35 set before coding).

## Hybrids

A **hybrid** scheme uses a cheap linear reconstruction where a sensor judges the data smooth, and a nonlinear one
(WENO or TENO) elsewhere:

$$ u^-_{i+1/2} = \begin{cases} \sum_k d_k q_k \ \text{(or a central, compact or energy-preserving value)} &
\text{sensor: smooth} \\ \sum_k \omega_k q_k & \text{sensor: not smooth} \end{cases} $$

Two motives recur: less numerical dissipation where the flow is smooth, which matters for turbulence, and less cost,
since the linear branch needs no weights and, for systems, no characteristic decomposition. Their weak point is the
sensor and its threshold. Examples are Pirozzoli's (2002) compact–WENO hybrid, switched on density differences, and
the central–WENO and central–TENO hybrids of the turbulence literature.

**In SETH, today.** The pieces of a hybrid, not yet a hybrid of its own (roadmap: issue #4, Phase 3):

- **the linear kernels**, the smooth branch: `seth_linear3/5/7/9` and their `_interp` twins (the upwind backgrounds of
  WENO-(2S-1)), `seth_linear6[_interp]` (TENO6's) and `seth_linear8` (TENO8's), each the background $\sum_k d_k q_k$
  collapsed by the table generator to one row of face coefficients. They equal TENO5, TENO6, TENO8 and the TENO5/6
  interpolations on smooth data to $4.4\cdot10^{-16}$ of the data, reach their design orders, and cost one dot product
  per state: 0.095 ns per face for `seth_linear5` on the RTX 4070, a tenth of WENO5-JS. A linear background is not a
  central scheme: the upwind ones keep their dissipation; the gain is the cost of the weights and, for systems, of the
  characteristic projection the caller can skip where the sensor says smooth;
- **the sensors**, the switch: `seth_sensor_teno` (Fu's TENO-based detector), `seth_sensor_ren` (Ren et al.),
  `seth_sensor_lcr` (Li, Chen & Ren 2022) and `seth_sensor_tau` (Wan & Xia 2022), each on the six values around a
  face, returning $s \in [0, 1]$ with $s = 1$ where the published test says smooth. On a sine they judge smooth up to
  $k = 0.64$, $0.04$, $2.07$ and $1.36$ ([H1](/verification/hybrids/sensors): the paper's 0.64 and 2.06 reproduced).
  Measured further on H1: Li, Chen & Ren's never flags a critical point, and is the most selective in flows: 6 % of
  Shu–Osher's faces nonlinear, where TENO6's cost hybrid stays within 1.044 of the pure scheme's $L_1$ (Fu's detector
  1.020 at 14 %, Wan & Xia's 1.006 at 27 %; TENO5's upwind background within 1.01 with these three); none is scale
  invariant by itself, and $r_\mu = 1/(\max f - \min f)$ over the field makes all four so. **Recommended** (from
  these one-dimensional measurements, published thresholds, nothing tuned; turbulent flows may revise it): Li, Chen &
  Ren's sensor, for an upwind or a nearly non-dissipative background (TENO6; by extension TENO8 and a central flux,
  not measured); Wan & Xia's where the last few per cent of accuracy on a nearly non-dissipative background is worth
  four times the nonlinear faces; never Ren's (it flags extrema at every resolution); and always $r_\mu = 1/(\max f - \min f)$, one reduction over the field;
- **the cost**, measured at kernel level on the RTX 4070: a sensor 0.30 ns per face, a linear kernel 0.095, WENO-Z5
  0.985. A hybrid that calls the sensor on both biases and the linear kernel where both say smooth costs 0.67 on a
  mostly smooth field, and breaks even with WENO-Z5 at about 30% of faces flagged (60% with one sensor call per face,
  as a finite-volume caller with `vp = vm` needs). On a GPU a warp with one flagged face runs both branches. The
  larger saving of a hybrid, the characteristic projection skipped where smooth, is the caller's and is not measured
  here;
- **the central fluxes**, the non-dissipative branch: `seth_keep_pe2/4/6/8` and `seth_keep_pe8_opt`, the KEEP-PE
  flux (Shima et al. 2021, as written by De Michele et al. 2023, Eq. (15): the kinetic energy flux
  $\overline{\rho u_n}\,(\mathbf u_l\cdot\mathbf u_r)/2$, the internal energy flux $\bar p\,\bar u_n/(\gamma-1)$, the
  split pressure work) in Pirozzoli's conservative form of order $2L$. Because a kinetic-energy-preserving flux above
  second order is a sum of two-point products of the state, these take the primitive state of $2L$ cells and return
  the flux vector: a second interface. The sixth-order one is ADAM's `compute_kep_flux`, bitwise
  ([H2](/verification/hybrids/central));
- **the shock sensors** of the central branch, `seth_ducros` (Ducros et al. 1999) and `seth_ducros_compression` (De
  Vanna et al. 2023a), from the velocity gradient the caller computes;
- the linear weights and candidate coefficients of every scheme, as named constants (`seth_weno_js_tables`,
  `seth_z_teno_tables`);
- the weights ratios of the `_wr` kernels, $\min_k \omega_k/d_k$, a sensor already computed by the nonlinear branch;
- TENO itself, whose cutoff returns the linear scheme exactly on smooth data: a hybrid inside the weights, with the
  scale-separated $\chi_k$ as its sensor.

## The kernels

| Family | Kernels | Candidates | Order | Parameter |
|---|---|---|---|---|
| WENO-JS | `seth_weno{1,3,5,7,9}_js` (+ `_wr`, `_interp`, `_interp_wr`) | $S = 1..5$, width $S$ | $2S-1$ | `rmu` |
| WENO-Z | `seth_wenoz{3,5,7,9}` (+ `_interp`) | $S = 2..5$, width $S$ | $2S-1$ | `p` |
| TENO | `seth_teno5`, `seth_teno6` (+ `_interp`), `seth_teno5_opt`, `seth_teno5_az`, `seth_teno6_opt`, `seth_teno6_neno6`, `seth_teno8`, `seth_teno5_a`, `seth_teno8_a` | 3 to 6, incremental | 5, 6, 8-point 6 (5-opt 4, 6-opt 5, NENO6 4) | `ct` (TENO5-AZ: `ct`, `rmu`; -A: `cr`, `xi`, `a1`, `a2`, `rmu`) |
| centred WENO | `seth_weno{4,6,8}_c` | WENO-JS ($S = 2..4$) blended with the central row on $2S$ cells | $2S$ | `eps`, `rmu` |
| positivity | `seth_positivity_euler`, `seth_euler_admissible` | | | |

Arguments and constants: [API](/guide/api). Stencil layout: [Conventions](/guide/conventions).
