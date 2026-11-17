# Flags

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| NP1 | Non-polynomial | P, C | CPU (gfortran, nvfortran, ifx), device | yes | `seth_test_nonpoly` |

## Case and scope

Every non-polynomial reconstruction of the TENO line starts from the same input: the flags of TENO5's cutoff,
$\delta_k = 1$ where candidate $k$ survives, $\delta_k = 0$ where it is cut. TENO5 computes them at every face and
discards them; the detectors of TENO5-DR and TENO5-THINC read them at the neighbouring cells (Liang, Xiao, Shyy & Fu
2022, §3). `seth_teno5_mask` returns them as one integer per cell, the **mask**
$m = \delta_0 + 2\delta_1 + 4\delta_2$, from the five values centred on the cell. This case verifies that the mask
is TENO5's own cutoff and that it reproduces the flags the paper prints.

| Bit | Candidate | Cells | Role |
|---|---|---|---|
| $\delta_0$ | $S_0$ | $\{-1, 0, 1\}$ | central |
| $\delta_1$ | $S_1$ | $\{0, 1, 2\}$ | downwind |
| $\delta_2$ | $S_2$ | $\{-2, -1, 0\}$ | upwind |

SETH's candidate order is the papers'. The flag of a cell is the cutoff of the left-biased TENO5 reconstruction at
its right face; for the other bias the caller mirrors the data. The criterion is plain TENO5's:
$\tau_5 = \lvert\beta_1 - \beta_2\rvert$, $\chi_k = \gamma_k / \sum\gamma$, $\gamma_k = (1 + \tau_5/(\beta_k +
\varepsilon))^6$, $\delta_k = [\chi_k \ge C_T]$, $C_T = 10^{-5}$, $\varepsilon = 10^{-40}$ (R4P: $10^{-30}$).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Kernel | the body of `seth_teno5` (`seth_z_teno_body.INC`) compiled with `SETH_TENO_MASK`: the indicators, $\tau_5$ and the cutoff, no candidate values, no weights; the preprocessed code of every other kernel unchanged, token for token | [`seth_nonpoly.F90`](https://github.com/szaghi/SETH/blob/main/src/lib/seth_nonpoly.F90) |
| Reference | a twin of `seth_teno5` in test code, the same body with `SETH_TENO_FLAGS`, returning the kernel's states and the masks of its two states | [`seth_test_teno_flags.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_test_teno_flags.F90) |
| Spec 1 | $10^6$ random stencil pairs (smooth cubics, independent values, jumps of height $10^{-6}$ to $10^2$, kinks): the twin's states equal `seth_teno5`'s bitwise, and the mask of $v_p(-2{:}2)$ and of the mirrored $v_m(1-m)$ equal the twin's two masks | [`seth_test_nonpoly.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_test_nonpoly.F90) |
| Spec 2 | Liang et al. 2022, Tables 1–3: the profiles on cells 0..10, constant continuation beyond the ends, the flags of cells 2..8 asserted digit for digit; profiles 2 and 3 also mirrored in value | |
| Spec 3 | $v(m) \to v(-m)$ swaps $\delta_1$ and $\delta_2$ and leaves $\delta_0$, on $10^5$ random stencils, exactly | |
| Spec 4 | constant data and a jump of $10^{30}$ (R4P: $10^{16}$): the mask lies in 1..7 | |
| Device | `seth_test_device`, case 51: both masks of 4096 stencils on the GPU equal the host's | [`seth_test_device.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_test_device.F90) |

The bounds were frozen in issue #6 (step 1.1) before the kernel was written.

## Program

```bash
./exe/seth_test_nonpoly
scripts/run_vv.sh --only seth_test_nonpoly
```

## Results

![The flags of the three profiles of Liang et al. 2022, Tables 1–3, computed by seth_teno5_mask](/verification/nonpoly-flags.png)

The flags of cells 2..8, as `seth_teno5_mask` returns them, beside the profile (the paper prints cells 0..10; the
outer two of each side see the constant continuation):

<!--@include: ../_tables/nonpoly-flags.md-->

The masks seen on the random stencils of spec 1 and on constant data:

<!--@include: ../_tables/nonpoly-masks.md-->

## Findings

- The mask is TENO5's cutoff, bit for bit: on $10^6$ random stencils the flag kernel and the kernel's own cutoff agree
  everywhere, because they are the same code.
- Liang et al.'s Tables 1–3 are reproduced digit for digit, and the mirrored profiles give the same flags.
- **Constant data does not give the mask 7.** On a constant that is not exact in binary, the indicators are
  round-off, of the order of $(\epsilon_{\rm mach} c)^2$, far above $\varepsilon = 10^{-40}$; $\tau_5 / \beta_k$ is
  then of order one and the cutoff drops candidates by the accident of rounding. Which candidates depends on the
  compiler: with gfortran 65 % of 1000 constants in $[-10, 10]$ give a mask other than 7 (the table above), among
  them the mask 1, which the detectors read as $\zeta = 2$; with nvfortran `-fast` every one of them gives the mask 5
  ($\zeta = 1$). In R4P ($\varepsilon = 10^{-30}$) even data equal to 1 gives the mask 5: the indicators of flat data
  are about $(6\cdot10^{-8} f)^2$, and Tables 1–3 are reproduced in R4P only with $\varepsilon \ge 10^{-14}$ (the test
  uses $10^{-12}$ there). TENO5's state is unaffected (every candidate reproduces a constant), and so is the state of
  the overrides, whose guards return the cell value on flat data (Takagi et al., Eq. 17; Liang et al., Eq. 33); a
  count of the detector's activations in uniform regions is not. Takagi et al. plot activations with
  $\varepsilon = 10^{-8}$ for this reason (§4.11). The detectors' pages report activations with this in view.
