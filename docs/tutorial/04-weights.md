---
title: 4. The weights
---

# 4. The weights

WENO's nonlinear weights compare smoothness indicators $\beta_k$, quadratic in the data, with a regulariser
$\varepsilon$ that is a pure number. So the weights depend on the **magnitude** of the data, unless the indicators are
rescaled. That is what the argument `rmu` does: the kernel uses $\beta_k\,r_\mu^2$, and $r_\mu = 1/\mu$ with $\mu$ a
magnitude of the data makes the weights independent of it.

## A small jump

The step of chapter 1, multiplied by a scale, with Jiang and Shu's weights (`rmu = 1`) and with the descaler
`seth_weno_descaling`, $\mu$ the mean $|v|$ of the stencils:

<<< @/examples/snippets/t04_weights-scale.F90{fortran}

## A small jump on a large offset

A pressure or a density is often a small fluctuation on a large mean. There the mean $|v|$ is the offset, and the
descaler sees no jump. The remedy (Wang & Don 2022, Remark 2): subtract the stencil's mean, call the kernel with the
affine-invariant descaler, and add the mean back:

<<< @/examples/snippets/t04_weights-offset.F90{fortran}

## Running it

<<< @/examples/output/t04_weights.txt{text}

- With `rmu = 1` the step of height 1 is reconstructed without overshoot, but at $10^{-4}$ the regulariser
  $\varepsilon = 10^{-6}$ swamps the indicators: the weights become the linear ones and the state overshoots by 18%
  of the jump, as the linear scheme of chapter 1 does.
- With the descaler the result does not depend on the scale.
- On the offset, the mean-$|v|$ descaler fails the same way; subtracting the mean and using the affine-invariant
  descaler gives no overshoot.

WENO-Z and TENO take no `rmu` (TENO5-AZ and the adaptive TENO do): their weights are ratios of indicators. WENO-Z's
exponent `p` is the other choice to make: 1 for WENO-Z3/5/7, **2 for WENO-Z9**, which is not robust at shocks with 1
([WENO-Z](/reference/weno-z)).

::: tip What you learned
The weights see the data through `eps` and `rmu`. For data of arbitrary magnitude, pass a descaler; for data with a
large offset, subtract the stencil's mean first. Reference: [Descalers](/reference/weno-js#descalers).
:::

<details><summary>The whole program</summary>

<<< @/examples/snippets/t04_weights.F90{fortran}

</details>
