# Sensors

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| H1 | Hybrids | P | CPU (gfortran, nvfortran, ifx), device (one case each) | yes | `seth_test_sensors`, `seth_vv_sensors`, `seth_test_sensor_fields` |

## Case and scope

A hybrid scheme is only as good as its sensor: the function of the stencil that says where the cheap linear branch
may replace the nonlinear one. Each source in the literature asserts that its own is "more universal", and few
compare them on the same data. Li, Chen & Ren (2022, §3.2) did, on two static tests; this case runs SETH's four
sensor kernels on the same two tests and asserts the paper's numbers, then carries the comparison further (the sensor
study of the roadmap, [issue #4](https://github.com/szaghi/SETH/issues/4), Phase 5): critical points under refinement,
scaled and offset data, two flows with shocks where each sensor drives a cost hybrid, and the Ducros sensors on
analytic two-dimensional fields.

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

| Critical points | $\cos(\pi x)$ (order 1) and $\sin^3(\pi x)$ (order 2) around $x = 0$, $m = 2 \dots 1024$ cells per unit length; is any of the five faces nearest the point flagged? | `seth_test_sensors` |
| Scaled data | test 2 as $a f + b$, $a = 10^{-6} \dots 10^{6}$, $b = 0, 100$; three normalisations: none ($r_\mu = 1$), the paper's per-stencil $1/\max(\lvert v\rvert, 1)$, the field's range $1/(\max f - \min f)$; the fraction of faces whose decision changes | |
| Flows | Shu–Osher and the Woodward–Colella blast, 400 cells, the [1-D harness](../shocks/sod#models-and-numerics) (blast with the positivity limiter), against WENO5-JS on 16 times the cells; the sensors on the density as the harness's hybrid reads it; the cost hybrid (`HYB_COST`: the linear background where the sensor, dilated by two faces, says smooth) of TENO5 and TENO6 driven by each | [`seth_vv_sensors.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_sensors.F90) |
| Ducros | sampled 2-D fields, $N \times N$ on $[-1, 1]^2$: the stencil sensors on the density along $x$ and $y$ (a cell flagged when a face of it is), against `seth_ducros` and `seth_ducros_compression` on the central-difference velocity gradient with ADAM's policy (threshold 0.2, $\varepsilon = 10^{-8}(c/h)^2$, no density-jump override) | [`seth_test_sensor_fields.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_test_sensor_fields.F90) |

## Program

```bash
./exe/seth_test_sensors
./exe/seth_vv_sensors
./exe/seth_test_sensor_fields
scripts/run_vv.sh --only seth_test_sensors --only seth_vv_sensors --only seth_test_sensor_fields
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

### Critical points

Faces flagged near a critical point, by resolution ($m$ cells per unit length; a period of $\cos(\pi x)$ is $2m$
cells):

<!--@include: ../_tables/sensors-cp.md-->

### Scaled and offset data

The fraction of the faces of test 2 whose decision changes when the data become $a f + b$:

<!--@include: ../_tables/sensors-scaling.md-->

### In flows with shocks

![Faces flagged by each sensor in x and t, Shu–Osher and the blast waves](/verification/sensor-flows.png)
*The faces each sensor flags (black) on the density of the pure TENO5 runs, 60 frames. Ren's flags extrema of the
pre-shock sine ($x \approx \pm 2.2, \pm 2.8$) and the post-shock entropy waves cell by cell; Fu's TENO detector and
Wan & Xia's flag the entropy waves as bands; Li, Chen & Ren's flags the shocks and little else.*

Averaged over the 60 frames, the fraction of faces flagged on Shu–Osher is 4.1 % (TENO-based), 8.9 % (Ren), 1.7 %
(Li, Chen & Ren), 6.2 % (Wan & Xia); on the blast 6.7 %, 4.9 %, 4.6 %, 8.6 %. The cost hybrid driven by each:

<!--@include: ../_tables/sensor-flows.md-->

### Against Ducros on analytic fields

<!--@include: ../_tables/sensor-fields.md-->

## Findings

- **The paper's numbers are reproduced**: Fu's detector judges a sine smooth up to $k = 0.64$, Li, Chen & Ren's up
  to $2.07$ against the printed $2.06$; Ren's flags critical points from the lowest wavenumbers.
- On test 2, Fu's detector flags every extremum of the faster wave ($k = \pi/4 > 0.64$) and Ren's every extremum of
  both; Li, Chen & Ren's and Wan & Xia's flag the jump alone.
- Wan & Xia's detector, not in Li, Chen & Ren's comparison, sits between the two: smooth up to $k = 1.36$, clean on
  test 2. Its $\tau_5$ vanishes on smooth data, so it is a WENO-Z indicator used as a switch.
- The thresholds are the papers', fixed: no case-by-case tuning. Every regulariser is absolute, so the sensors are
  not scale invariant; the normalisation $\max(\lvert f\rvert, 1)$ is the paper's choice, passed as `rmu`.
- **At critical points Li, Chen & Ren's sensor is the only one that never flags**, at any resolution and for both
  orders; the others leave a first-order point alone from 4 (Wan & Xia) to 128 (Ren) cells per unit length, a
  second-order one from 64 to 128. A sensor that flags a critical point sends a smooth extremum to the nonlinear
  scheme: harmless for accuracy, lost for cost.
- **No sensor is scale invariant by itself, and the paper's normalisation does not make it so.** On data of amplitude
  $10^{-6}$ the absolute regularisers change 39 % (TENO-based) and 82 % (Ren) of the decisions, 2–3 % for Li, Chen &
  Ren's and Wan & Xia's. The per-stencil $1/\max(\lvert v\rvert, 1)$ only rescales data above 1, which were already
  unaffected, and makes an offset matter: unit data on a mean of 100 change 36 % of the TENO detector's decisions.
  Normalising by the field's range ($r_\mu = 1/(\max f - \min f)$, a reduction the caller computes once) makes every
  sensor invariant to scale and offset (asserted).
- **In the flows, the sensor decides whether TENO6's cost hybrid works.** With Li, Chen & Ren's sensor, the most
  selective, TENO6's hybrid is 1.81 times the pure scheme's $L_1$ on Shu–Osher: the post-shock entropy waves pass as
  smooth and go to TENO6's nearly central background ([H3](./harness)). With Fu's TENO detector it is 1.04, with Wan &
  Xia's 1.00, at 14 % and 28 % of the faces nonlinear instead of 6 %. TENO5's upwind background tolerates every sensor
  (0.995–1.07), and on the blast every hybrid is within 0.94–1.11 of its pure scheme. Selectivity and safety pull
  apart: the sensor that flags least is the one that can starve a nearly non-dissipative background. The guide's
  recommendation follows from it ([Schemes](/guide/schemes#hybrids)): Li, Chen & Ren's for an upwind background,
  Fu's TENO detector or Wan & Xia's for a nearly non-dissipative one, never Ren's, always the range normalisation.
- **Ducros sees what the velocity says, and only that.** It never flags the Taylor–Green vortex at any resolution
  (the discrete divergence vanishes); there Fu's and Wan & Xia's detectors flag the under-resolved field ($N \le 16$
  and $N = 8$), Li, Chen & Ren's flags nothing at any $N$, and Ren's flags the extrema at every $N$. A shock inside a vortex escapes it once the vorticity exceeds twice the shock's
  dilatation (half the shock cells missed at $\omega_0 = 4\lvert\theta\rvert$; $\phi > 0.2$ iff
  $\omega < 2\lvert\theta\rvert$), and a contact at uniform velocity is invisible to it, the case ADAM's
  density-jump override exists for; the stencil sensors on the density catch both. The classic $\phi$ flags a steep
  expansion, the compression form $\psi$ does not.
- What this case does not show: how the sensors behave in turbulence, where the vorticity is everywhere and the
  smooth scales are under-resolved by design. That comparison needs ADAM's flows (Taylor–Green at $Re = 1600$, the
  channel, the boundary layer; ADAM issue 119), with the flagged fraction and the statistics side by side.
