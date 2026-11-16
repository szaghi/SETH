# Approximate dispersion relation

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| F2 | Spectral | P | CPU (gfortran, nvfortran, ifx) | yes | `seth_vv_adr` |

## Case and scope

A nonlinear scheme has no modified wavenumber, but its response to a single Fourier mode defines an approximate one,
the ADR (Pirozzoli 2006). The TENO and WENO papers use it to show how close a nonlinear scheme stays to its linear
background, and where it departs from it (llm-wiki weno-teno `definitions/approximate-dispersion-relation.md`).

What the case catches: a TENO cutoff that discards candidates on resolved waves; nonlinear weights that do not
approach the optimal ones at low wavenumbers; anti-dissipation (Im $\tilde\vartheta > 0$), which amplifies a mode.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Signal | $u_j = \cos(\vartheta_n j)$, $\vartheta_n = 2\pi n/N$, $n = 1 … N/2 - 1$, $N = 200$, periodic, $h = 1$, unit amplitude | [`seth_vv_adr.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_adr.F90) |
| Operator | $D_j = F_{j+1/2} - F_{j-1/2}$ of $u_t + u_x = 0$, $F$ the left-biased state | `adr` |
| ADR | $G = \frac{2}{N}\sum_j D_j e^{-i\vartheta_n j}$, $\tilde\vartheta = -iG$: Re the dispersion (ideal $\vartheta$), Im the dissipation (ideal 0, negative dissipative); the sign convention of the wiki page and of [F1](./linear-spectra) | |
| Background | the linear scheme's $\tilde\vartheta$ from the tables (as F1) | `background` |
| Parameters | each scheme's defaults: WENO-JS $\varepsilon = 10^{-6}$ (not scale invariant: its ADR changes with the amplitude), WENO-Z and TENO $\varepsilon = 10^{-40}$ | |

Asserted:

1. TENO5 and TENO6 equal their linear background to $10^{-12}$ for every $\vartheta \le \pi/4$;
2. every scheme of order $\ge 5$ is within $10^{-3}\vartheta$ of its background for $\vartheta \le 0.3$;
3. Im $\tilde\vartheta \le 10^{-12}$ everywhere for every scheme with a dissipative background, except TENO6, whose
   measured anti-dissipative band must lie beyond $\vartheta = 0.4\pi$ and stay below $10^{-3}$.

## Program

```bash
./exe/seth_vv_adr
scripts/run_vv.sh --only seth_vv_adr
```

## Results

![Approximate dispersion relation](/verification/adr.png)

*Dispersion and dissipation of every nonlinear reconstruction; dotted: the linear backgrounds of TENO5 and TENO6.*

## Findings

- TENO5 and TENO6 are their linear schemes, to round-off, up to the wavenumber where the cutoff starts to discard
  candidates: $\vartheta/\pi \approx 0.57$ for TENO5 and $0.46$ for TENO6. WENO-Z5 departs from its background
  immediately but stays within $2\cdot10^{-6}\vartheta$ at low wavenumbers; WENO5-JS within $2\cdot10^{-4}\vartheta$.
- The third-order weights (WENO3-JS, WENO-Z3) stay far from their background even at low wavenumbers at unit
  amplitude, a known property of 3rd-order WENO: printed, not asserted.
- The centred WENO's background is its central row of order $2S$, which does not dissipate. WENO6-C and WENO8-C
  approach it at low wavenumbers, within $5\cdot10^{-5}\vartheta$ and $3\cdot10^{-6}\vartheta$ (asserted, as for
  every scheme of order $\ge 5$); WENO4-C, whose upwind part is third order, stays within $4.5\cdot10^{-3}\vartheta$
  (printed). None of the three is anti-dissipative at any wavenumber (Im $\tilde\vartheta \le -2\cdot10^{-15}$).
  This analysis advects one wave of unit amplitude with no background, which the blend's detector sees in full; a
  small disturbance on a smooth flow it barely sees, and damps far less ([S6](../shocks/burgers)).
- **TENO6 is anti-dissipative**, Im $\tilde\vartheta$ up to $2\cdot10^{-4}$, in a narrow band
  $\vartheta/\pi = 0.46$–$0.49$ just where its cutoff starts to act. Modes of that wavelength (about $4h$) are the
  packets TENO6 leaves undamped behind the [Burgers](../shocks/burgers) shock. The optimised TENO6 weights, with the
  same candidates, show no anti-dissipation.
