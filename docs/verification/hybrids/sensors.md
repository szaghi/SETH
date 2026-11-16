# Sensors

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| H1 | Hybrids | P | CPU (gfortran, nvfortran, ifx), device (one case each) | yes | `seth_test_sensors` |

## Case and scope

A hybrid scheme is only as good as its sensor: the function of the stencil that says where the cheap linear branch
may replace the nonlinear one. Each source in the literature asserts that its own is "more universal", and few
compare them on the same data. Li, Chen & Ren (2022, §3.2) did, on two static tests; this case runs SETH's four
sensor kernels on the same two tests and asserts the paper's numbers.

| Kernel | Sensor | Published threshold |
|---|---|---|
| `seth_sensor_teno` | Fu's TENO-based detector: four three-point stencils, $\gamma_k = (\beta_k + 10^{-4})^{-6}$, smooth iff every normalised $\chi_k > C_T^*$ (Li, Chen & Ren, Eqs. (55)–(57)) | $C_T^* = 5\cdot10^{-4}$ |
| `seth_sensor_ren` | Ren et al.'s ratio of the first differences, the minimum over four nodes (Eqs. (48)–(53)) | $\psi_c = 0.3$ |
| `seth_sensor_lcr` | Li, Chen & Ren's: Ren's ratio with second differences added (Eq. (60)) | $\psi_c = 0.4$ |
| `seth_sensor_tau` | Wan & Xia's (2022) $\eta = \tau_5/(\varepsilon + \bar\beta)$ on WENO5's candidates (§3) | $\eta \le 1/2$ smooth, $\ge 1$ not |

Every kernel returns $s \in [0, 1]$, $s = 1$ exactly where the published binary test says smooth.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Test 1 | $f_j = \sin(k j)$, $j = 0 \dots 191$, $k = 0.01 \dots 3$; the fraction of faces with $s = 1$ (their Fig. 9) | [`seth_test_sensors.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_test_sensors.F90) |
| Test 2 | $\sin(12\pi x) - 2$ for $x < 0$, $\sin(24\pi x + \pi/2) + 2$ for $x \ge 0$, 192 cells on $[-1, 1]$ (their test (b), Fig. 10) | |
| Normalisation | the stencil divided by $\max(\lvert f\rvert_{\max}, 1)$ (their Eq. (61)) | `rmu` |
| Asserted | test 1: Fu's detector smooth up to $k = 0.64 \pm 0.05$, Li, Chen & Ren's up to $2.06 \pm 0.05$, Ren's failing below 0.3; test 2: Li, Chen & Ren's flags the jump and nothing farther than $3h$ from it; constant data smooth, a unit step anywhere in the stencil not | |

## Program

```bash
./exe/seth_test_sensors
scripts/run_vv.sh --only seth_test_sensors
```

## Results

![Sensors on a sine sweep and on a jump between two sines](/verification/sensors.png)
*Left: the fraction of faces judged smooth on a sine of wavenumber $k$. Right: $1 - s$ along test 2 (each sensor
offset), over the data.*

| Sensor | smooth up to $k$ (paper) | test 2: faces flagged, of which farther than $3h$ from the jump |
|---|---|---|
| TENO-based (Fu) | 0.64 (0.64) | 73, 68 |
| Ren et al. | 0.04 ("very low wavenumbers") | 154, 147 |
| Li, Chen & Ren | 2.07 (2.06) | 5, 0 |
| Wan & Xia | 1.36 (not in the paper) | 4, 0 |

## Findings

- **The paper's numbers are reproduced**: Fu's detector judges a sine smooth up to $k = 0.64$, Li, Chen & Ren's up
  to $2.07$ against the printed $2.06$; Ren's flags critical points from the lowest wavenumbers.
- On test 2, Fu's detector flags every extremum of the faster wave ($k = \pi/4 > 0.64$) and Ren's every extremum of
  both; Li, Chen & Ren's and Wan & Xia's flag the jump alone.
- Wan & Xia's detector, not in Li, Chen & Ren's comparison, sits between the two: smooth up to $k = 1.36$, clean on
  test 2. Its $\tau_5$ vanishes on smooth data, so it is a WENO-Z indicator used as a switch.
- The thresholds are the papers', fixed: no case-by-case tuning. Every regulariser is absolute, so the sensors are
  not scale invariant; the normalisation $\max(\lvert f\rvert, 1)$ is the paper's choice, passed as `rmu`.
- What this case does not show: how the sensors behave in a flow with shocklets and turbulence, which no static test
  can. That comparison is the roadmap's Phase 5, on ADAM's flows.
