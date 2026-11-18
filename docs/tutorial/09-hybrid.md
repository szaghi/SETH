---
title: 9. A hybrid
---

# 9. A hybrid

Most of a flow is smooth, where a linear scheme is as good as a nonlinear one and much cheaper: no smoothness
indicators, no weights, and for a system no characteristic projection. A **hybrid** asks a sensor at each face and
calls the nonlinear kernel only where the sensor says the data is not smooth.

## The switch

The sensor of Li, Chen & Ren (2022), `seth_sensor_lcr`, on the six cells around the face, once for each bias; the
linear scheme where both say smooth ($s = 1$), TENO5 elsewhere:

<<< @/examples/snippets/t09_hybrid-hybrid.F90{fortran}

`rmu`, the inverse of the field's range ($\max u - \min u$, one reduction per stage), makes the sensor's decisions
independent of the scale and offset of the data: its threshold `SETH_SENSOR_LCR_PSIC` is a pure number.

## Running it

The advection of chapter 3, one period on 200 cells:

<<< @/examples/output/t09_hybrid.txt{text}

The hybrid calls TENO5 at 6.8% of the faces, around the edges of the square wave and the steep flanks of the bump,
and the linear scheme everywhere else. Its error is TENO5's to 1%; its oscillations are those of the faces where the
linear scheme sits next to a nonlinear face, $3\cdot10^{-4}$, against the linear scheme's $7.5\cdot10^{-2}$.

Which sensor, and which threshold: the measurements of [H1](/verification/hybrids/sensors) recommend Li, Chen & Ren's
for an upwind background. For a non-dissipative branch, the central KEEP-PE fluxes with the Ducros shock sensors
([Central fluxes](/reference/central)).

::: tip What you learned
A hybrid is a branch per face: a sensor, the linear kernel, the nonlinear kernel. SETH supplies the three pieces; the
switch is a line of yours. Reference: [Sensors](/reference/sensors), [Linear schemes](/reference/linear).
:::

<details><summary>The whole program</summary>

<<< @/examples/snippets/t09_hybrid.F90{fortran}

</details>

Next: [10. Non-polynomial TENO](./10-nonpoly).
