# Linear backgrounds

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| F1 | Spectral | P | CPU (gfortran, nvfortran, ifx) | yes | `seth_test_spectral` |

## Case and scope

On smooth, well-resolved data every nonlinear scheme returns, or approaches, its linear background, the combination
of its candidates with the optimal weights. The modified wavenumber of that background (Pirozzoli 2002; Martín et al.
2006; llm-wiki high-order-capturing-schemes `definitions/modified-wavenumber.md`) says how it disperses and
dissipates each Fourier mode. It is computed analytically from the compiled tables, so a wrong coefficient shows as a
wrong curve.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Operator | $D u_i = (F_{i+1/2} - F_{i-1/2})/h$, $F_{i+1/2} = \sum_m c_m u_{i+m}$, $c_m$ the background's face coefficients | [`seth_test_spectral.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_test_spectral.F90) |
| Modified wavenumber | $\tilde\xi = -i\sum_m c_m e^{im\xi}(1 - e^{-i\xi})$; Re $\tilde\xi$ the dispersion, Im $\tilde\xi$ the dissipation (negative: dissipative) | `xit` |
| Asserted | Im $\tilde\xi \le 10^{-14}$ on $(0, \pi]$ for the max-order backgrounds; low-wavenumber error slope $\ge o + 1 - 0.3$; resolution (largest $\xi/\pi$ with $|{\rm Re}\,\tilde\xi - \xi|/\xi < 1\%$) 0.35 for the 5th-order and 0.42 for the 7th-order WENO background, within 0.01 (Martín et al. 2006, via llm-wiki `methods/bandwidth-optimized-weno.md`) | |

The interpolation tables are not a derivative operator and are left out (their accuracy is in [A1](../accuracy/order-of-accuracy)).

## Program

```bash
./exe/seth_test_spectral
scripts/run_vv.sh --only seth_test_spectral
```

## Results

![Modified wavenumber of the linear backgrounds](/verification/spectral.png)

*Dispersion and dissipation of every linear background (WENO-Z shares WENO-JS's).*

## Findings

- No max-order background is anti-dissipative, and every low-wavenumber slope is $o + 1$.
- Martín et al.'s resolution efficiencies are reproduced: 0.35 (5th order) and 0.42 (7th order).
- The optimised TENO6 weights trade one order (Fu et al.) or two (NENO6) for a wider resolved band: their curves
  follow $\tilde\xi = \xi$ farther before they bend.
