# Schemes

Every formula below is implemented as written; where it comes from, and what is derived from it, is in the
repository's `notes/m2-sources.md`. $\beta_k$ are the Jiang–Shu smoothness indicators of the candidate stencils,
$d_k$ the linear (optimal) weights, $q_k$ the candidate values at the face.

## Jiang–Shu smoothness indicators

For a candidate of $w$ cells reconstructed by the polynomial $p_k$ of degree $w-1$ (Jiang & Shu 1996):

$$ \beta_k = \sum_{j=1}^{w-1} \Delta x^{2j-1} \int_{x_{i-1/2}}^{x_{i+1/2}} \left(\frac{d^j p_k}{dx^j}\right)^2 dx. $$

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
`seth_weno_descaling` gives ADAM's generic one).

- **Scale invariance.** With `rmu = 1` the weights depend on the scale of the data (an absolute $\varepsilon$): data of
  magnitude $10^{-4}$ is reconstructed differently from the same data at $1$. With the descaler the result scales
  bitwise under $v \to 2^k v$.
- **Range.** $(\varepsilon + \beta)^{w}$ overflows double precision near $|v| \sim 10^{77}$ (WENO-3), $10^{51}$
  (WENO-5), $10^{38}$ (WENO-7/9); with the descaler the limit is $|v| \sim 10^{150}$ ($\beta$ itself).

**Interpolation** (`_interp`): face values from point values (Jiang, Shu & Zhang 2013), other candidate coefficients
and linear weights, the same indicators.

## WENO-Z: WENO-Z3/5/7/9

The candidates and linear weights of WENO-JS, the weights of Borges et al. (2008):

$$ \alpha_k = d_k\left[1 + \left(\frac{\tau}{\beta_k + \varepsilon}\right)^{p}\right] $$

| Scheme | $\tau$ | Source |
|---|---|---|
| WENO-Z3 | $\lvert\beta_0 - \beta_1\rvert$ | Don & Borges; Li et al. 2021 |
| WENO-Z5 | $\lvert\beta_0 - \beta_2\rvert$ | Borges et al. 2008 |
| WENO-Z7 | $\lvert\beta_0 + 3\beta_1 - 3\beta_2 - \beta_3\rvert$, $O(\Delta x^7)$ | Don & Borges, as quoted by Rathan & Naga Raju 2016 |
| WENO-Z9 | $\lvert\beta_0 - \beta_4\rvert$, $O(\Delta x^7)$ | Castro et al.'s rule for odd $r$, as restated by Li et al. 2022 |

$p$ = 1 (Borges) or 2 (fifth order at first-order critical points; recommended for higher orders); $\varepsilon =
10^{-40}$ is the usual choice. The orders of $\tau$ are measured exactly by the table generator (Castro et al.'s
$\tau_7 = \lvert\beta_0-\beta_1-\beta_2+\beta_3\rvert$ is $O(\Delta x^6)$, as Li et al. 2022 state).

## TENO: TENO5, TENO6, TENO6-opt, TENO6-NENO6

Candidates of incremental width (Fu, Hu & Adams 2016): $S_0 = \{i-1, i, i+1\}$, $S_1 = \{i, i+1, i+2\}$,
$S_2 = \{i-2, i-1, i\}$, and for TENO6 the downwind $S_3 = \{i, \dots, i+3\}$.

$$ \gamma_k = \left(1 + \frac{\tau_K}{\beta_k + \varepsilon}\right)^{6},\quad
   \chi_k = \frac{\gamma_k}{\sum_j \gamma_j},\quad
   \delta_k = \begin{cases} 0 & \chi_k < C_T \\ 1 & \text{otherwise} \end{cases},\quad
   \omega_k = \frac{d_k \delta_k}{\sum_j d_j \delta_j} $$

with $\tau_5 = \lvert\beta_1 - \beta_2\rvert$ and $\tau_6 = \beta_6 - \tfrac16(\beta_1 + \beta_2 + 4\beta_0)$
($\beta_6$ on the six-point stencil), $\varepsilon = 10^{-40}$ (`SETH_TENO_EPS`). SETH evaluates $\chi$ as
$(g_k/g_{\max})^6 / \sum (g_j/g_{\max})^6$, $g = 1 + \tau/(\beta+\varepsilon)$: the same values, without the overflow of
$\gamma$ on large jumps.

| Kernel | $d_k$ | Order | $C_T$ constant | Source |
|---|---|---|---|---|
| `seth_teno5` | 0.6, 0.3, 0.1 | 5 | `SETH_TENO5_CT` = $10^{-5}$ | Fu et al. 2016, Table 5 |
| `seth_teno6` | 0.45, 0.3, 0.05, 0.2 | 6 | `SETH_TENO6_CT` = $10^{-7}$ | Fu et al. 2016, Table 5 |
| `seth_teno6_opt` | 0.462, 0.3, 0.054, 0.184 | 5 | `SETH_TENO6_OPT_CT` = $10^{-7}$ | Fu et al. 2016, Table 5 |
| `seth_teno6_neno6` | 0.42943, 0.17273, 0.08557, 0.31227 | 4 | `SETH_TENO6_NENO6_CT` = $10^{-6}$ | the NENO6 paper |

"TENO6-opt" denotes different weights in different papers; SETH's is the founding paper's.
