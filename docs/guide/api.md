# API

`use seth` gives every kernel, the constants and the kinds `I4P`, `I8P`, `R4P`, `R8P` and `R_P` (from PENF). All
arrays are explicit-shape with the scheme's bounds; every kernel is `pure subroutine`, `!$acc routine seq`,
`!$omp declare target`.

**Kind.** Every real argument is `real(R_P)`, PENF's working kind, fixed when the library is built: R8P by default,
R4P with the `*-r4` build modes (`-D_R_P_IS_R4P`, which PENF and SETH must both see; a build where they disagree stops
with a division by zero in a constant expression). One precision per build: no generic interfaces, the kernels are
compiled once for that kind. The literals whose meaning depends on the kind (regularisers, back-off steps) live in
`seth_kind_parameters`, one set per kind. In R4P, `SETH_TENO_EPS` is $10^{-30}$ (the papers' $10^{-40}$ is below the
smallest normal R4P number, and a GPU flushes it to zero): pass it as `eps` to WENO-Z and TENO in either kind.

The R4P build computes the WENO-JS, WENO-Z and TENO weights as ratios to the smallest indicator: the same weights,
with every factor at most 1. The R8P forms overflow R4P's range (3.4e38): $(\tau/\varepsilon)^2$ on a unit step
(WENO-Z with $p = 2$), or $1/\beta^q$ of WENO-JS ($q$ = 2 to 4 with the order) once the data passes about $10^4$
(WENO9) to $10^9$ (WENO3). Two limits remain, both of
R4P itself:

- the indicators are quadratic in the data, so data above about $10^{19}$ overflow them, less by the square root of
  the indicator's coefficients: about $3\cdot10^{15}$ for WENO9-JS, whose integer coefficients reach $2.5\cdot10^6$
  (the indicators are computed on the raw data, before any descaler);
- the expanded WENO-JS indicator cancels to about $10^{-7} v^2$, so data with a large mean needs the descaler `rmu`.

Two more things hold in R4P only, both measured by the test suite:

- the affine-invariant descaler `seth_weno_descaling_ai` alone does not keep WENO-JS ENO on a small step over a large
  offset: the expanded indicators cancel the offset to $\varepsilon\,\text{offset}^2$, which $\mu^{-2}$ then makes
  $O(1)$ (an overshoot of 0.25 of a step of $2^{-13}$ on 1, as without a descaler). Subtract the stencil's mean first
  (Wang & Don), which stays ENO in R4P;
- a scheme of order 6 or more reaches R4P's round-off within the first refinements of a smooth problem (errors near
  $10^{-5}$ for a derivative on grids of 20 to 640 points): its design order cannot be observed in R4P, only that it
  converged as far as R4P allows.

The test suite passes in R4P (`fobis build --mode tests-gnu-r4 && bash scripts/run_tests.sh --exe-dir exe-r4`, 33
programs: the comparisons with ADAM are R8P only) on gfortran, nvfortran and ifx, against expectations measured in R4P
([test programs](/verification/unit-tests#single-precision)), every case held to its bound. What R4P costs and buys,
measured: [Precision](/guide/precision).

## WENO-JS

```fortran
call seth_weno<N>_js         (vp, vm, eps, rmu, rp, rm)              ! N = 1, 3, 5, 7, 9
call seth_weno<N>_js_wr      (vp, vm, eps, rmu, rp, rm, wrp, wrm)
call seth_weno<N>_js_interp  (vp, vm, eps, rmu, rp, rm)
call seth_weno<N>_js_interp_wr(vp, vm, eps, rmu, rp, rm, wrp, wrm)
```

| Dummy | Type | Meaning |
|---|---|---|
| `vp`, `vm` | `real(R_P), intent(in) :: v(1-S:S)` | left-/right-biased stencils, $S = (N+1)/2$ |
| `eps` | `real(R_P), intent(in)` | regulariser; `SETH_WENO_JS_EPS` = $10^{-6}$ (ADAM's) |
| `rmu` | `real(R_P), intent(in)` | inverse descaler: `1` (Jiang–Shu), `seth_weno_descaling(...)` or `seth_weno_descaling_ai(...)` |
| `rp`, `rm` | `real(R_P), intent(out)` | left and right states at the face |
| `wrp`, `wrm` | `real(R_P), intent(out)` | weights ratios $\min_k \omega_k/d_k \in [0,1]$: 1 on smooth data, toward 0 near a jump |

```fortran
rmu = seth_weno_descaling(S, sigma, vp, vm)   ! sigma = 0: rmu = 1 exactly; sigma = 1: 1/mean|v|
rmu = seth_weno_descaling_ai(S, vp, vm)       ! affine invariant: 1/mean|v - mean v| (Wang & Don 2022)
```

## WENO-Z

```fortran
call seth_wenoz<N>       (vp, vm, eps, p, rp, rm)   ! N = 3, 5, 7, 9; vp, vm(1-S:S); p integer, 1 or 2
call seth_wenoz<N>_interp(vp, vm, eps, p, rp, rm)   ! the same from point values at the cell centres
```

## Centred WENO

```fortran
call seth_weno4_c(vp, vm, eps, rmu, rp, rm)   ! vp, vm(-1:2); WENO3-JS blended with the 4th-order central row
call seth_weno6_c(vp, vm, eps, rmu, rp, rm)   ! vp, vm(-2:3); WENO5-JS, 6th order
call seth_weno8_c(vp, vm, eps, rmu, rp, rm)   ! vp, vm(-3:4); WENO7-JS, 8th order
```

`eps` > 0 (`SETH_WENO_JS_EPS`) and `rmu` (1, or `seth_weno_descaling`) as for WENO-JS: both go to the upwind weights
and to the blend's indicators. WENO8-C is not ENO for a jump in its two outer cells ([Schemes](/guide/schemes#centred-weno-weno4-6-8-c)).

## TENO

```fortran
call seth_teno5      (vp, vm, eps, ct, rp, rm)   ! vp, vm(-2:3)
call seth_teno5_opt  (vp, vm, eps, ct, rp, rm)
call seth_teno5_az   (vp, vm, eps, ct, rmu, rp, rm) ! ct = lambda_d; rmu = 1 (Tang 2024) or the descaler
call seth_teno6      (vp, vm, eps, ct, rp, rm)
call seth_teno6_opt  (vp, vm, eps, ct, rp, rm)
call seth_teno6_neno6(vp, vm, eps, ct, rp, rm)
call seth_teno8     (vp, vm, eps, ct, rp, rm)   ! vp, vm(-3:4)
call seth_teno5_a   (vp, vm, eps, cr, xi, a1, a2, rmu, rp, rm) ! adaptive cutoff; SETH_TENO5_A_CR, _XI, _A1, _A2
call seth_teno8_a   (vp, vm, eps, cr, xi, a1, a2, rmu, rp, rm) ! vp, vm(-3:4); SETH_TENO8_A_*
call seth_teno5_interp(vp, vm, eps, ct, rp, rm)   ! interpolation from point values at the cell centres
call seth_teno6_interp(vp, vm, eps, ct, rp, rm)
```

The interpolation kernels take the same arguments as the reconstruction ones; TENO5-opt, TENO6-opt and TENO6-NENO6 have none
(their weights are tuned for the reconstruction).

`eps = SETH_TENO_EPS` ($10^{-40}$); `ct` the cutoff: `SETH_TENO5_CT`, `SETH_TENO5_OPT_CT`, `SETH_TENO5_AZ_CT` (the ceiling of the adaptive one), `SETH_TENO6_CT`, `SETH_TENO6_OPT_CT`,
`SETH_TENO6_NENO6_CT`, `SETH_TENO8_CT` (see [Schemes](./schemes)).

## Non-polynomial TENO: flags, detectors, overrides

```fortran
call seth_teno5_dr   (vp, vm, eps, ct, t1, t2, rmu, rp, rm)  ! vp, vm(-5:6): TENO5-DR (Liang et al. 2022), one call per face
call seth_teno5_thinc(vp, vm, eps, ct, beta, rmu, rp, rm)    ! vp, vm(-5:6): TENO5-THINC (Takagi et al. 2022)
call seth_teno5_st   (vp, vm, eps, ct, rp, rm)                ! vp, vm(-2:3): stencil targeting (Liang & Fu 2024)
call seth_teno5_rbf  (vp, vm, eps, ct, cr, s1, s2, rp, rm)    ! vp, vm(-2:3): Gaussian-RBF; s1 = alpha_1 dx, s2 = alpha_2 dx
call seth_teno5_rbf_dr(vp, vm, eps, ct, cr, s1, s2, t1, t2, rmu, rp, rm) ! vp, vm(-5:6): Liang & Fu 2024's complete scheme
! the pieces, for a caller building its own composite (a stored-mask, two-pass variant)
call seth_teno5_mask(v, eps, ct, m)                 ! v(-2:2) centred on a cell: m = delta_0 + 2 delta_1 + 4 delta_2
z = seth_teno5_zeta(m)                              ! zeta of a mask (Liang et al. 2022, Eq. 17): 2, 1, -1 or 0
call seth_dr_detect(mp, mm, lp, lm)                 ! mp, mm(-3:3): the masks around the face; lp, lm in 0..2
call seth_thinc_detect(mp, mm, lp, lm)              ! lp, lm in 0..1
call seth_dr_override(vp, vm, lp, lm, t1, t2, rmu, rp, rm)    ! vp, vm(-1:2); rp, rm inout
call seth_thinc_override(vp, vm, lp, lm, beta, rmu, rp, rm)
```

`eps = SETH_TENO_EPS`, `ct = SETH_TENO5_CT` (the papers' TENO5: the detectors are calibrated on its flags);
`t1, t2 = SETH_DR_T1, SETH_DR_T2` ($\tanh 0.91$, $\tanh 1.23$); `beta = SETH_THINC_BETA` (1.6); `rmu` scales the
first-order guard of the overrides, $(f_{i+1}-f_i)(f_i-f_{i-1})\,r_\mu^2 \le$ `SETH_NONPOLY_GUARD` ($10^{-15}$) gives
the cell value ($r_\mu = 1$: the papers' form). The composite kernels read **twelve cells**, `-5..6`: a solver needs 6
ghost cells and projects all twelve into the face's characteristic basis ([Conventions](./conventions)). The masks of
the neighbouring cells are recomputed per face; `mp(c)` is the mask of cell `c` of the left-biased data, `mm(c)` of
cell `c` of the mirrored right-biased data `vm(1-m)`. `cr = SETH_RBF_CR` (0.4), `s1 = SETH_RBF_ALPHA1*dx`,
`s2 = SETH_RBF_ALPHA2*dx` (both 3): the RBF kernel's results depend on the length unit, and its validated range is
$s \le 0.12$. Verification: [Non-polynomial](/verification/nonpolynomial/).

## Linear kernels

```fortran
call seth_linear<N>       (vp, vm, rp, rm)   ! N = 3, 5, 7, 9 (WENO backgrounds), 6 (TENO6), 8 (TENO8, bounds -3:4)
call seth_linear<N>_interp(vp, vm, rp, rm)   ! N = 3, 5, 7, 9, 6: from point values
```

The linear background of a scheme, one row of face coefficients (`LIN<N>_C` in `seth_z_teno_tables`): the smooth
branch of a hybrid. Bounds as the scheme's (`1-S:S`; TENO6 and linear6 `-2:3`).

## Sensors

```fortran
call seth_sensor_teno(v, rmu, thr, s)   ! v(-2:3); thr = SETH_SENSOR_TENO_CT  (5e-4)
call seth_sensor_teno_eps(v, rmu, thr, eps, s) ! eps = SETH_SENSOR_TENO_EPS_1D (1e-4) or _MD (1e-3)
call seth_sensor_ren (v, rmu, thr, s)   !          thr = SETH_SENSOR_REN_PSIC (0.3)
call seth_sensor_lcr (v, rmu, thr, s)   !          thr = SETH_SENSOR_LCR_PSIC (0.4)
call seth_sensor_tau (v, rmu, thr, s)   !          thr = SETH_SENSOR_TAU_ETA  (0.5)
```

`s` in $[0, 1]$, `s = 1` where the published binary test says smooth (the linear branch may be used); below 1 the
published continuous form, usable as a blending weight. `rmu` multiplies the data (1, or Li, Chen & Ren's
`1/max(maxval(abs(v)), 1)`): the regularisers are absolute. With `rmu = 1/(max f - min f)`, the range of the field the
caller reduces once, every sensor is invariant to the scale and offset of the data; the per-stencil form is not, and
under an offset it changes decisions. See [H1](/verification/hybrids/sensors).

`seth_sensor_teno`'s regulariser is $\varepsilon$ = `SETH_SENSOR_TENO_EPS_1D` = $10^{-4}$, Fu's (2019) one-dimensional
value; `seth_sensor_teno_eps` takes it as an argument, and Fu uses `SETH_SENSOR_TENO_EPS_MD` = $10^{-3}$ in several
dimensions. $\varepsilon$ is absolute: with `rmu = 1` it sets the smallest jump flagged, about 0.005 of the data's
units in one dimension and 0.017 in several.

## Central fluxes and shock sensors

```fortran
call seth_keep_pe6(d, gamma, w, f)   ! w(5,-2:3): rho, u, v, w, p of the cells around the face; f(5): the flux
                                     ! seth_keep_pe2 (w(5,0:1)), _pe4 (-1:2), _pe8 and _pe8_opt (-3:4)
call seth_ducros(eps, g, phi)              ! g(3,3) = d u_a/d x_b; phi = theta^2/(theta^2 + |omega|^2 + eps)
call seth_ducros_compression(eps, g, psi)  ! psi = max(-theta/sqrt(theta^2 + |omega|^2 + eps), 0)
```

`d` the normal direction (1 to 3), `gamma` the ratio of specific heats. See [H2](/verification/hybrids/central).

## Monotonicity-preserving bound

```fortran
call seth_mp_bound(vp, vm, alpha, rp, rm)   ! vp, vm(-2:3): the data; rp, rm: the states, bounded in place
```

`alpha = SETH_MP_ALPHA` (4). Call it after any reconstruction kernel, on the five values around each state's upwind
cell (pass `vp(-2:3)`, `vm(-2:3)` of a wider stencil). See [Schemes](./schemes#monotonicity-preserving-bound).

## Positivity limiter

```fortran
call seth_positivity_euler(nv, gamma, lam, eps_rho, eps_p, ql, qr, f_high, f_low, f, theta)
ok = seth_euler_admissible(nv, gamma, q, eps_rho, eps_p)
```

The positivity-preserving flux limiter of Hu, Adams & Shu (2013) at one face of an Euler (ideal-gas) solver:
`f = theta f_high + (1 - theta) f_low`, `theta` the largest value in $[0, 1]$ for which the one-sided states
`ql - lam f` and `qr + lam f` keep $\rho \ge$ `eps_rho` and $p \ge$ `eps_p`. Every argument is plain data:

| Argument | Meaning |
|---|---|
| `nv` | 3, 4 or 5: the state is $(\rho, \rho u_1, …, \rho u_d, E)$ in $d = $ `nv` $- 2$ dimensions |
| `gamma` | ratio of specific heats |
| `lam` | $2\Delta t/h$ in 1-D; in $d$ dimensions $2\Delta t/(w_k h_k)$ for the faces of direction $k$, the partition weights $w_k$ summing to 1 |
| `eps_rho`, `eps_p` | the floors; Hu, Adams & Shu take $\min(10^{-13}$, the initial minimum$)$; `SETH_POSITIVITY_EPS` $= 10^{-13}$ |
| `ql`, `qr` | the states of the cells left and right of the face |
| `f_high`, `f_low` | the high-order face flux and a first-order positivity-preserving one (Lax–Friedrichs with $\alpha \ge \max(\vert u_n\vert  + c)$ is, for `lam` $\alpha \le 1$) |
| `f`, `theta` | the limited flux and its factor (1: `f_high` unchanged) |

If even `f_low` leaves a one-sided state below the floors (`lam` too large, or an inadmissible cell), the kernel
returns `theta = 0` and `f = f_low`; `seth_euler_admissible` lets the caller check its cells. The kernel is `pure`,
`acc routine seq` and `omp declare target`, like the reconstruction kernels. See the [positivity
verification](/verification/contracts/positivity).

- **`nv` at most `SETH_POSITIVITY_NV_MAX` = 5.** The kernel's local arrays have that constant dimension, since an array
  sized by a dummy argument is allocated from the device heap at every call ([measured](/guide/performance#device-kernel-rules-arrays-and-loop-extents):
  333× slower). For a larger `nv` it returns `theta = 0` and `f = f_low`.
- **Round-off.** The floors hold for the one-sided states as the kernel evaluates them. A caller that evaluates
  `ql - lam f` differently (a fused multiply-add) gets states that differ by $\epsilon(|q| + |\lambda f|)$ and pressures
  that differ by about that times $1 + |u|^2$: choose the floors above it.

## Tables

`seth_weno_js_tables` and `seth_z_teno_tables` expose every coefficient as a named constant (`WENO5_JS_A`,
`TENO6_D`, ...), generated by `scripts/gen_weno_js_tables.py` and `scripts/gen_z_teno_tables.py`; the descriptors of
WENO-Z and TENO give the candidate cells (`*_CLO`, `*_CW`), face coefficients (`*_P`), indicator factors (`*_L`,
`*_DQ`, `*_LG`, `*_DG`), linear weights (`*_D`) and $\tau$ coefficients (`*_TC`, `*_TG`).
