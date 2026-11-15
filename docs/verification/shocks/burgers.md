# Burgers equation

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| S6 | Shock capturing | P | CPU (gfortran, nvfortran, ifx; debug and optimised) | yes, animation | `seth_vv_burgers` |

## Case and scope

The inviscid Burgers equation from $u_0 = 0.5 + \sin\pi x$ on $[-1, 1]$, periodic (Shu 1998, §4): the classic WENO
accuracy test, here with the exact solution before and after the shock forms. Before ($t = 0.5/\pi$) the solution is
smooth, with first-order critical points, and every scheme must converge, TENO at its design order. After
($t = 1.5/\pi$) the shock moves at speed 0.5. The error over the domain converges at first order at best, the error
farther than 0.1 from the shock should stay high order, and the overshoot must stay bounded.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Scheme | conservative finite differences, global Lax–Friedrichs splitting $f^\pm = (f \pm \alpha u)/2$, $\alpha = 1.1\max|u|$, SSP-RK3 | [`seth_vv_scalar.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_scalar.F90) |
| Time step | before the shock $\Delta t = 0.4\,h\,(h/h_0)^{(o-3)/3}$, $h_0 = 2/40$: the $O(\Delta t^3)$ error of RK3 falls with the $O(h^o)$ spatial one (Henrick et al. use $\Delta t \sim h^{5/3}$ for 5th order) | [`seth_vv_burgers.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_burgers.F90) |
| Exact solution | in the frame moving at 0.5, $v = u - 0.5$ solves Burgers from $\sin\pi y$; the characteristics cross first at $y = \pm 1$ at $t = 1/\pi$ and by odd symmetry the shock stays there; $v(y) = \sin\pi\xi$, $\xi$ the root of $\xi + t\sin\pi\xi = y$ on the rising branch $[0, \xi_m]$, $\xi_m = \arccos(-1/\pi t)/\pi$ (the characteristics beyond it have entered the shock); bisection to $10^{-15}$ | `exact`, `v` |

Why $\alpha = 1.1\max|u|$: with $\alpha = \max|u|$ exactly, $\partial_x f^- = (u - \alpha)u_x/2$ vanishes to second
order at the maximum of $u$, a second-order critical point where TENO falls to second order. Measured before the
shock without the margin, TENO6 converged at 5.7, 5.0, 4.7, 4.2, 4.1 from 80 to 1280 cells; with the margin, at 5.7,
5.8, 6.0, 6.0.

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[-1, 1]$, periodic |
| Grids | $N = 40 … 640$ (report mode), $40 … 320$ (test run); profiles and animation at $N = 80$ |
| Times | $t = 0.5/\pi$ (smooth), $1.5/\pi$ (shocked) |

## Program

```bash
./exe/seth_vv_burgers
scripts/run_vv.sh --only seth_vv_burgers
```

## Results

![Burgers: convergence before and after the shock](/verification/burgers.png)

*Left: $L_1$ before the shock. Centre: $L_1$ farther than 0.1 from the shock, after it. Right: the pointwise error
after the shock at $N = 80$.*

<!--@include: ../_tables/burgers.md-->

<div class="seth-casts">
<figure class="seth-cast seth-cast-wide" style="--ar: 2.25">
<video class="light-only" autoplay loop muted playsinline poster="/verification/casts/burgers-light.png" aria-label="Burgers, shock formation"><source src="/verification/casts/burgers-light.webm" type="video/webm"><source src="/verification/casts/burgers-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/verification/casts/burgers-dark.png" aria-label="Burgers, shock formation"><source src="/verification/casts/burgers-dark.webm" type="video/webm"><source src="/verification/casts/burgers-dark.mp4" type="video/mp4"></video>
<figcaption>From $t = 0$ to $1.5/\pi$, 80 cells: the exact solution, WENO5-JS and TENO6.</figcaption>
</figure>
</div>

## Findings

- Before the shock every scheme converges, TENO5 and TENO6 at 5 and 6 (asserted), once the splitting has its margin.
- After the shock, the error away from it falls at high order for WENO-JS, WENO-Z and TENO5. **TENO6 leaves a
  sign-alternating grid-scale packet** in the smooth region, far from the shock: about $1.6\cdot10^{-3}$ at 80
  cells, around $x \approx 0.69$. Its smooth-region error then converges at about second order only
  ($8\cdot10^{-6}$ at 640 cells, against $3\cdot10^{-11}$ for TENO5 and WENO-Z5). TENO6-opt and TENO6-NENO6, with the
  same candidates, show the same pollution, weaker ($\sim10^{-6}$). The exact
  characteristics all run into the shock, so the packet is numerical: high-wavenumber dispersive errors leave the
  shock with an upstream group velocity, and TENO6 does not damp them. The [ADR](../spectral/adr) shows why: TENO6 is
  anti-dissipative in a narrow band near $\vartheta/\pi \approx 0.47$, a wavelength of about $4h$, where its cutoff
  starts to act.
- The overshoot after the shock stays below $10^{-4}$ for most schemes and below $10^{-2}$ for all but one
  (asserted). WENO-Z7 reaches $6\cdot10^{-3}$.
- **WENO-Z9 overshoots by 0.28 behind the shock** (printed, exempted from the bound; to be fixed:
  [issue #1](https://github.com/szaghi/SETH/issues/1)). Its
  $\tau_9 = |\beta_0 - \beta_4|$ (Castro et al.'s rule for odd $r$, as restated by Li et al. 2022b; the choice was
  already recorded as unsourced for WENO-Z9) is small when the two outer candidates straddle a smeared shock
  symmetrically. The weights then stay close to the linear ones across the shock. The same weakness makes WENO-Z9
  lose positivity in the first step of the [blast wave](./blast) and on the [123 problem](./test123).
