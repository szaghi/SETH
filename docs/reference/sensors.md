---
title: Sensors
---

# Sensors

The switch of a hybrid: on the six values around a face, a smoothness $s \in [0, 1]$, with $s = 1$ where the
published binary test says smooth (the linear branch may be used); below 1, the published continuous form, usable as
a blending weight. Module `seth_sensors`.

```fortran
call seth_sensor_teno    (v, rmu, thr, s)        ! v(-2:3); thr = SETH_SENSOR_TENO_CT
call seth_sensor_teno_eps(v, rmu, thr, eps, s)   ! eps = SETH_SENSOR_TENO_EPS_1D or SETH_SENSOR_TENO_EPS_MD
call seth_sensor_ren     (v, rmu, thr, s)        ! thr = SETH_SENSOR_REN_PSIC
call seth_sensor_lcr     (v, rmu, thr, s)        ! thr = SETH_SENSOR_LCR_PSIC
call seth_sensor_tau     (v, rmu, thr, s)        ! thr = SETH_SENSOR_TAU_ETA
```

| Argument | Declaration | Meaning |
|---|---|---|
| `v` | `real(R_P), intent(in) :: v(-2:3)` | the values around the face |
| `rmu` | `real(R_P), intent(in)` | multiplies the data: 1, or `1/(max f - min f)` over the field (recommended) |
| `thr` | `real(R_P), intent(in)` | the published threshold, below |
| `eps` | `real(R_P), intent(in)` | the regulariser of Fu's detector |
| `s` | `real(R_P), intent(out)` | smoothness, 1: smooth |

| Sensor | Threshold constant | Value | Source |
|---|---|---|---|
| `seth_sensor_teno` | `SETH_SENSOR_TENO_CT` | $5\cdot10^{-4}$ | Fu 2019, Eq. 4.4 ($\varepsilon$ = `SETH_SENSOR_TENO_EPS_1D`) |
| `seth_sensor_teno_eps` | `SETH_SENSOR_TENO_CT`; `SETH_SENSOR_TENO_EPS_1D`, `SETH_SENSOR_TENO_EPS_MD` | $5\cdot10^{-4}$; $10^{-4}$ (1-D), $10^{-3}$ (several dimensions) | Fu 2019 |
| `seth_sensor_ren` | `SETH_SENSOR_REN_PSIC` | 0.3 | Ren et al., as Li, Chen & Ren 2022, Eq. 52 |
| `seth_sensor_lcr` | `SETH_SENSOR_LCR_PSIC` | 0.4 | Li, Chen & Ren 2022, Eq. 60 |
| `seth_sensor_tau` | `SETH_SENSOR_TAU_ETA` | 0.5 | Wan & Xia 2022 |

- **Scale.** The regularisers are absolute. With `rmu = 1/(max f - min f)`, the range of the field the caller reduces
  once, every sensor is invariant to the scale and offset of the data; a per-stencil magnitude is not, and under an
  offset it changes decisions.
- **Which one.** From the one-dimensional measurements of [H1](/verification/hybrids/sensors): Li, Chen & Ren's
  (`seth_sensor_lcr`) for an upwind or nearly non-dissipative background; Wan & Xia's where the last per cent of
  accuracy is worth four times the nonlinear faces; never Ren's, which flags extrema at every resolution.

Two more sensors are already computed by the nonlinear kernels: the weights ratios `wrp`, `wrm` of the WENO-JS
`_wr` kernels ([WENO-JS](./weno-js)), and TENO's cutoff itself. API: [`seth_sensors`](/api/src/lib/seth_sensors).
