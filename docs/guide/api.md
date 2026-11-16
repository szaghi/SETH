# API

`use seth` gives every kernel, the constants and the kinds `I4P`, `I8P`, `R4P`, `R8P` (from PENF). All arrays are
explicit-shape with the scheme's bounds; every kernel is `pure subroutine`, `!$acc routine seq`,
`!$omp declare target`.

## WENO-JS

```fortran
call seth_weno<N>_js         (vp, vm, eps, rmu, rp, rm)              ! N = 1, 3, 5, 7, 9
call seth_weno<N>_js_wr      (vp, vm, eps, rmu, rp, rm, wrp, wrm)
call seth_weno<N>_js_interp  (vp, vm, eps, rmu, rp, rm)
call seth_weno<N>_js_interp_wr(vp, vm, eps, rmu, rp, rm, wrp, wrm)
```

| Dummy | Type | Meaning |
|---|---|---|
| `vp`, `vm` | `real(R8P), intent(in) :: v(1-S:S)` | left-/right-biased stencils, $S = (N+1)/2$ |
| `eps` | `real(R8P), intent(in)` | regulariser; `SETH_WENO_JS_EPS` = $10^{-6}$ (ADAM's) |
| `rmu` | `real(R8P), intent(in)` | inverse descaler: `1` (Jiang–Shu), `seth_weno_descaling(...)` or `seth_weno_descaling_ai(...)` |
| `rp`, `rm` | `real(R8P), intent(out)` | left and right states at the face |
| `wrp`, `wrm` | `real(R8P), intent(out)` | weights ratios $\min_k \omega_k/d_k \in [0,1]$: 1 on smooth data, toward 0 near a jump |

```fortran
rmu = seth_weno_descaling(S, sigma, vp, vm)   ! sigma = 0: rmu = 1 exactly; sigma = 1: 1/mean|v|
rmu = seth_weno_descaling_ai(S, vp, vm)       ! affine invariant: 1/mean|v - mean v| (Wang & Don 2022)
```

## WENO-Z

```fortran
call seth_wenoz<N>       (vp, vm, eps, p, rp, rm)   ! N = 3, 5, 7, 9; vp, vm(1-S:S); p integer, 1 or 2
call seth_wenoz<N>_interp(vp, vm, eps, p, rp, rm)   ! the same from point values at the cell centres
```

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
call seth_sensor_ren (v, rmu, thr, s)   !          thr = SETH_SENSOR_REN_PSIC (0.3)
call seth_sensor_lcr (v, rmu, thr, s)   !          thr = SETH_SENSOR_LCR_PSIC (0.4)
call seth_sensor_tau (v, rmu, thr, s)   !          thr = SETH_SENSOR_TAU_ETA  (0.5)
```

`s` in $[0, 1]$, `s = 1` where the published binary test says smooth (the linear branch may be used); below 1 the
published continuous form, usable as a blending weight. `rmu` multiplies the data (1, or Li, Chen & Ren's
`1/max(maxval(abs(v)), 1)`): the regularisers are absolute. With `rmu = 1/(max f - min f)`, the range of the field the
caller reduces once, every sensor is invariant to the scale and offset of the data; the per-stencil form is not, and
under an offset it changes decisions. See [H1](/verification/hybrids/sensors).

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
