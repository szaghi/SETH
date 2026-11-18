# Length units

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| NP8 | Non-polynomial | P | CPU (gfortran, nvfortran, ifx) | the two profiles | `seth_vv_nonpoly` |

## Case and scope

The Gaussian-RBF kernel's shape parameter is $s = \lambda\,\Delta x$ with $\lambda = \alpha_1 + \alpha_2(1 - g)$ a pure
number (Liang & Fu 2024, Eq. 32): $\lambda$ multiplies a length, so the same flow written in other units gets another
$s$, and another result. TENO5 has no length scale. Experiment E5 of [issue #6](https://github.com/szaghi/SETH/issues/6)
(step 3.3) measures the difference.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Case | Shu–Osher on $[-5, 5]$, $N = 200$, $t = 1.8$, and the same flow with $x$ and $t$ divided by 8 (a power of 2, so that every operation of the harness scales exactly; with 10 the rounding of $\Delta x/10$ moves a TENO cutoff and the profile by $10^{-4}$) | [`seth_vv_nonpoly.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_nonpoly.F90) |
| Asserted | TENO5's density profiles equal bitwise | |
| Reported | TENO5-RBF ($s = 3\,\Delta x$): the largest density difference | |

## Program

```bash
./exe/seth_vv_nonpoly
```

## Results

<!--@include: ../_tables/nonpoly-units.md-->

![Shu-Osher at N = 200 on [-5, 5] and on the domain scaled by 1/8: TENO5 and TENO5-RBF](/verification/nonpoly-units.png)

## Findings

- **TENO5 does not depend on the length unit**: the two runs are equal bit for bit.
- **The Gaussian-RBF scheme does, strongly**: the same Shu–Osher flow differs by 0.47 in density between the two
  units, because $s$ goes from 0.15 to 0.019. A caller must choose $s_1$, $s_2$ knowing the unit of its grid; the
  validated range ($s \le 0.12$, [NP6](./rbf)) is stated in $s$, not in $\alpha$.
