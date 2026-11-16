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
| Scheme | conservative finite differences, global Lax–Friedrichs splitting $f^\pm = (f \pm \alpha u)/2$, $\alpha = 1.1\max\vert u\vert$, SSP-RK3 | [`seth_vv_scalar.F90`](https://github.com/szaghi/SETH/blob/main/src/tests/seth_vv_scalar.F90) |
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
- **WENO6-C and WENO8-C leave a packet of the same kind**, sign-alternating, upstream of the shock: at 80 cells
  $1.9\cdot10^{-5}$ and $1.8\cdot10^{-4}$ at most, against $1.0\cdot10^{-5}$ for WENO7-JS. Their smooth-region error
  then converges at first order or below on the finest pairs ($2.6\cdot10^{-6}$ and $3.6\cdot10^{-6}$ at 640 cells,
  against $10^{-10}$ and $3\cdot10^{-12}$ for WENO5-JS and WENO7-JS, their upwind parts), where before the shock both
  converge at 5 to 8. Neither is anti-dissipative ([ADR](../spectral/adr)); the cause is the blend's damping, measured
  with $\phi$ traced and varied in a copy of this harness (note
  [`weno-c-burgers-packet.md`](https://github.com/szaghi/SETH/blob/main/notes/weno-c-burgers-packet.md)):
  - every scheme's error grows in the last steepening before the shock forms; the centred schemes leave about 5 times
    their upwind parts' there, as a grid-scale packet (about 80% of neighbouring cells alternate in sign), at
    wavenumbers where the central rows' group velocity is reversed ($-1.2$ to $-2.6$ times the exact one): it moves
    into the smooth region against the flow;
  - the upwind parts then damp their error 23 and 12 times; the centred ones do not. The detector sees a small packet
    only as a small change of the background's differences: where the packet is, $\phi \approx 10^{-2}$, 1 to 3.5
    times its value on the exact solution. A small packet is damped at $\phi$ times the upwind rate, and the
    background's $\phi$ is $O(h)$ while that rate is $O(1/h)$: a damping that does not grow with the resolution. On
    linear advection a packet of amplitude $10^{-5}$ on $\sin(\pi x)$ keeps 30% of itself after 20 cells of travel
    under WENO6-C, 0.16% under WENO5-JS;
  - refuted: the upwind part ($\phi = 1$: $8\cdot10^{-11}$ and $4\cdot10^{-12}$), the regulariser as the cause
    ($\varepsilon = 10^{-40}$ in the blend changes $\phi$ there by less than 15%; the error falls 7.5 and 3.9 times), and
    the detector's reach (the whole stencil: at most 1.6 times lower).

  WENO4-C stays at WENO3-JS's level ($6\cdot10^{-7}$).
- The overshoot after the shock stays below $10^{-4}$ for most schemes and below $10^{-2}$ for all but one (asserted).
  WENO-Z7 reaches $6\cdot10^{-3}$, WENO-Z9 $9\cdot10^{-4}$.
- **TENO5-AZ overshoots behind the shock**: $3.0\cdot10^{-2}$ at 160 cells, $3.5\cdot10^{-2}$ at 320,
  $4.6\cdot10^{-3}$ at 640 (asserted below $5\cdot10^{-2}$, printed), where TENO5 stays below $3\cdot10^{-5}$. Its
  cutoff, $\min(10^{-5}, \tau)$ with $\tau = \lvert\beta_0 + \beta_2 - 2\beta_1\rvert$ (Tang 2024), drops to zero where
  the smeared shock lies across the stencil so that the upwind and downwind indicators balance the central one; then
  no candidate is cut. The descaled cutoff gives the same overshoot to three digits (the data is of order 1). The
  same scheme keeps fifth order at critical points where TENO5 falls to second
  ([critical points](../accuracy/critical-points)): the two properties come from the same vanishing cutoff.
- **WENO-Z9 needs the exponent $p = 2$.** With $p = 1$ it overshoots by 0.19 to 0.28 behind the shock on every
  grid, leaves $[-0.05, 1.05]$ on the [composite profile](./advection) and needs 6294 limited faces on the
  [blast wave](./blast); with $p = 2$, which every case of this section now runs, the overshoot is at most
  $9\cdot10^{-4}$, the composite profile stays within $[-8\cdot10^{-5}, 1.0005]$ and the blast wave needs 585
  limited faces. The convergence before the shock is the same for both. The cause is the exponent, not
  $\tau_9 = |\beta_0 - \beta_4|$: Castro et al. (2011) use that $\tau_9$ and report low exponents unstable at
  high order (llm-wiki high-order-capturing-schemes `methods/weno-z.md`). Measured on these one-dimensional cases
  only ([issue #1](https://github.com/szaghi/SETH/issues/1)); $p = 3$ and 4 bring no further gain here and cost
  accuracy (blast $L_1(\rho)$ 0.077 at $p = 2$, 0.083 at $p = 4$).
