# Gaussian-RBF fluxes

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| NP6 | Non-polynomial | P, C | CPU (gfortran, nvfortran, ifx), device | the band of growth at four $s$, E6; animations | `seth_test_nonpoly`, `seth_vv_adr`, `seth_vv_nonpoly`, `seth_test_tables` |

## Case and scope

Liang & Fu (2024) give each row of stencil targeting ([NP5](./targeting)) a Gaussian radial-basis correction:
the face value is $(P_m + s^2 C_m)\cdot f$, $P_m$ the polynomial row of the mask $m$, $C_m$ its $s^2$ row, the
expansion of the Gaussian RBF truncated at $s^2 = \lambda^2\Delta x^2$ as the paper prints it (Eq. 25). The shape
parameter comes from a sensor (Eqs. 30–32):

$$\eta_j = \frac{\lvert 2\Delta_- \Delta_+\rvert + \varepsilon}{\Delta_-^2 + \Delta_+^2 + \varepsilon},\quad
m = 1 - \min\Big(1, \frac{\min(\eta_{i-1}, \eta_i, \eta_{i+1})}{c_r}\Big),\quad g = (1-m)^4(1+4m),\quad
s = s_1 + s_2\,(1 - g)$$

with $c_r = 0.4$, $\varepsilon = 10^{-40}$ and, in the paper, $s_k = \alpha_k \Delta x$, $\alpha_1 = \alpha_2 = 3$. The
caller passes $s_1$, $s_2$ (owner decision D5): no grid enters the library, and the results depend on the length
unit. This case verifies the kernel `seth_teno5_rbf` (step 3.2 of [issue #6](https://github.com/szaghi/SETH/issues/6))
and measures the long-wave growth that the plan derived from the paper's own truncation error (E6, gate G6).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Kernel | the body of `seth_teno5` with `SETH_TENO_ROWS` and `SETH_TENO_RBF`: TENO5's cutoff, the sensor on each state's stencil, the row $P + s^2 C$ of its mask | [`seth_nonpoly.F90`](https://github.com/szaghi/SETH/blob/main/src/lib/seth_nonpoly.F90) |
| Rows | re-derived from the Gaussian RBF in 400-digit arithmetic, equal to Eq. 25 once its labels $R'_2$, $R'_3$ are swapped to match §3.1's stencils; every $C$ row sums to zero; the leading errors equal Eq. 29 | [`gen_nonpoly_tables.py`](https://github.com/szaghi/SETH/blob/main/scripts/gen_nonpoly_tables.py) |
| Specs 15–16 | $s_1 = s_2 = 0$ is `seth_teno5_st` bitwise; order 5.0 ± 0.1 with $\alpha$ fixed and $\Delta x \to 0$; the leading error of each of the paper's six rows (the masks 7, 3, 5, 1, 2, 4) with $\lambda = 3$ equals Eq. 29 within 1 % by Richardson on $\Delta x$, e.g. the five-point row's $-\Delta x^5(\lambda^2 g'''/2 + g^{(5)}/60)$ (measured: 0.007–0.04 %); every $C$ row zero outside its stencil and summing to zero (`seth_test_tables`); the kernel equals a reference written in the test (the sensor, the row) on $10^5$ stencils with the sensor active | [`seth_test_nonpoly.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_test_nonpoly.F90) |
| E6 | $\sin 2\pi x$ advected 100 periods with $\lambda$ frozen at 3 and 6, $N = 20$ and 40, SSP-RK3 at CFL 0.4 and 0.05: the growth rate of the amplitude against the $\lambda$ term of Eq. 29, beside TENO5's; asserted inside the validated range: at most twice the prediction, amplitude below twice the initial | [`seth_vv_nonpoly.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_nonpoly.F90) |
| Spec 18 | the linear five-point row with $s$ fixed: $\operatorname{Im}\tilde\theta = -\operatorname{Re}\sum_m c_m (e^{i\theta m} - e^{i\theta(m-1)})$ on $\theta \in (0, \pi]$, at the four $s^2$ the plan froze; asserted at the two smaller: the band of growth within 10 % of $\sqrt{30}\,s$, its peak within a factor 2 of $(200/3)\,s^6$ | |

| G6 | the nonlinear ADR of `seth_teno5_rbf` (`seth_vv_adr`, $N = 200$, unit amplitude, $s = 3\,$`rbf_dx`) at the four $s$ against its linear row: inside the validated range ($s = 0.03, 0.06$) the band's edge within 10 % of $\sqrt{30}\,s$ and its peak within twice $(200/3)\,s^6$, asserted; outside it reported | [`seth_vv_adr.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_adr.F90) |

## Program

```bash
./exe/seth_test_nonpoly
scripts/run_vv.sh --only seth_test_nonpoly
```

## Results

The band of growth of the linear RBF scheme against the prediction of the plan (§1.3: from Eq. 29,
$\operatorname{Im}\tilde\theta = \theta^4 (s^2/2 - \theta^2/60)$ at leading order):

<!--@include: ../_tables/nonpoly-adr.md-->

The nonlinear scheme against its linear row (G6, `seth_vv_adr`):

![The imaginary part of the modified wavenumber of TENO5-RBF at s = 0.03, 0.06, 0.15, 0.3, against its linear row](/verification/nonpoly-adr-rbf.png)

<div class="seth-casts">
<figure class="seth-cast seth-cast-wide" style="--ar: 2.25">
<video class="light-only" autoplay loop muted playsinline poster="/verification/casts/nonpoly-adr-morph-light.png" aria-label="The ADR of the linear RBF row as s grows"><source src="/verification/casts/nonpoly-adr-morph-light.webm" type="video/webm"><source src="/verification/casts/nonpoly-adr-morph-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/verification/casts/nonpoly-adr-morph-dark.png" aria-label="The ADR of the linear RBF row as s grows"><source src="/verification/casts/nonpoly-adr-morph-dark.webm" type="video/webm"><source src="/verification/casts/nonpoly-adr-morph-dark.mp4" type="video/mp4"></video>
<figcaption>The modified wavenumber of the linear five-point RBF row as $s = \lambda\Delta x$ grows from 0 to 0.3: the band of growth widens with $s$ (dotted: $\sqrt{30}\,s$).</figcaption>
</figure>
</div>

The growth in time (E6):

![The measured growth rate of a sine under TENO5-RBF against the prediction of Eq. 29](/verification/nonpoly-e6.png)

<!--@include: ../_tables/nonpoly-e6.md-->

<div class="seth-casts">
<figure class="seth-cast seth-cast-wide" style="--ar: 2">
<video class="light-only" autoplay loop muted playsinline poster="/verification/casts/nonpoly-lambda-light.png" aria-label="lambda along Shu-Osher"><source src="/verification/casts/nonpoly-lambda-light.webm" type="video/webm"><source src="/verification/casts/nonpoly-lambda-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/verification/casts/nonpoly-lambda-dark.png" aria-label="lambda along Shu-Osher"><source src="/verification/casts/nonpoly-lambda-dark.webm" type="video/webm"><source src="/verification/casts/nonpoly-lambda-dark.mp4" type="video/mp4"></video>
<figcaption>Shu–Osher under TENO5-RBF, N = 200; under it, $\lambda = s/\Delta x$ of the sensor (Eqs. 30–32) on the density: 3 where the flow is smooth, up to 6 at the shock and at the crests of the entropy waves.</figcaption>
</figure>
</div>

## Findings

- **The Gaussian-RBF correction is anti-dissipative at long waves**, exactly as its own truncation error says: the
  term $+\tfrac{\lambda^2}{2}\Delta x^5\,\partial_x^4 u$ of Eq. 29 grows every wave with $\theta^2 < 30 s^2$. The band and
  its peak match the leading-order prediction at small $s$ (the plan's derivation, tier D until now, confirmed);
  at larger $s$ the higher terms reduce both.
- **The growth is real in time** (E6): where TENO5 damps a resolved sine, the RBF scheme grows it at about the
  predicted rate when the time step is small (CFL 0.05); at $s = 0.3$ the amplitude nearly triples in 100 periods.
  At CFL 0.4, the paper's, the dissipation of SSP-RK3 hides part of it (a likely reason the paper's runs show no
  growth; not tested against them).
  Inside the validated range ($s = 0.075$) the measured rate is the predicted one.
- **The validated range is $s \le 0.12$**: there the peak stays below TENO6's $2\cdot10^{-4}$, the bound decision D7
  set. The paper's Shu–Osher run ($\Delta x = 0.05$, $s = 0.15$ to $0.30$) is outside it; its 2-D runs ($s \le 0.06$)
  are inside, where $s^2 C$ changes the five-point row's coefficients by at most $2.5\cdot10^{-3}$ (double Mach at
  $400\times100$) and by $2\cdot10^{-4}$ or less on the finer grids.
- **Outside the validated range the nonlinear scheme is more anti-dissipative than its row** (G6): its sensor reads
  the crests of a unit wave as non-smooth and raises $s$ toward $s_1 + s_2$. At $s = 0.15$, the $s$ of the paper's
  Shu–Osher run, the peak is $1.6\cdot10^{-2}$, 22 times the linear prediction, and the band reaches
  $\theta = 1.26$; at $s = 0.3$ the peak is 0.198. Inside the range the nonlinear band is the row's: edges 0.157
  and 0.314 against 0.164 and 0.329, peaks within 11 % of the prediction.
- The results depend on the length unit: $s = \alpha\,\Delta x$ with $\alpha$ fixed. The same flow on $[0, 10]$ and
  on $[0, 1]$ gets corrections 100 times apart ([NP8](./units): 0.47 in density on Shu–Osher).
