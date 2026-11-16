# The centred WENO's packet behind the Burgers shock (issue #4, 6.3)

2026-11-16. Experiment `notes/experiments/weno-c-packet/` (`weno_c_packet.F90`, `run.sh`; recorded output
`output-2026-11-16.txt`). Every number below is from that output, gfortran 16.2, R8P.

## Question

`seth_vv_burgers` (S6, report mode): after the shock forms, WENO6-C and WENO8-C leave a sign-alternating packet in the
smooth region; their L1 farther than 0.1 from the shock is 2.6e-6 and 3.6e-6 at 640 cells, converging at first order
or below on the finest pairs, where their upwind parts WENO5-JS and WENO7-JS reach 1e-10 and 3e-12. Their ADR (F2) has
no anti-dissipation, so TENO6's explanation does not apply. Why?

The centred WENO is r = r_c + phi (r_up - r_c): r_c the central row of order 2S, r_up WENO-JS of order 2S-1, phi the
largest |I_k - I_{k-1}|/(I_k + I_{k-1}) over cells -1..2, I_k = (Delta v)^2 + eps. The harness runs it with the
registry's eps = 1e-6 and rmu = 1 on the Lax-Friedrichs-split fluxes f+-.

## Method

The experiment reproduces the harness (same residual, splitting, RK3, time step, exact solution) with phi in variants:
`adam` (the library kernel), `copy` (the same composed in the program), `up` (phi = 1), `cen` (phi = 0), `eps40`
(eps = 1e-40 in the blend only), `wide` (the detector over every pair of cells 1-S..S), `wide40` (both).

- `copy` equals the library kernel bitwise on 3 x 100000 random stencils.
- `adam` reproduces `seth_vv_burgers` (run in report mode by `run.sh`): L1 before, after, away from the shock and the
  overshoot within 1e-8 relative + 1e-14 (largest: 0.40 of that bound; the smooth-region L1 to 6.5e-10 relative). Not
  bitwise: the two programs are compiled apart and differ in the last bits, which the shock amplifies; the largest
  relative difference, 8.8e-4, is on a pre-shock L1 at the round-off floor (9e-13).

## Results

L1 farther than 0.1 from the shock, t = 1.5/pi, 640 cells:

| phi | WENO6-C | WENO8-C |
|---|---|---|
| adam (the library) | 2.6e-6 | 3.6e-6 |
| up (phi = 1) | 8.5e-11 | 3.9e-12 |
| cen (phi = 0) | 9.9e-2 | 1.1e-1 |
| eps40 | 3.4e-7 | 9.2e-7 |
| wide | 2.2e-6 | 2.2e-6 |
| wide40 | 1.9e-7 | 3.2e-7 |

WENO4-C is not polluted (6.4e-7, against 1.2e-6 for its WENO3-JS part).

**Where it comes from** (N = 160, the error farther than 0.1 from the shock, or from where it forms):

| t pi | WENO6-C | its WENO5-JS (phi = 1) | WENO8-C | its WENO7-JS |
|---|---|---|---|---|
| 0.80 | 6.6e-7 | 1.6e-6 | 1.0e-8 | 2.9e-8 |
| 0.90 | 1.3e-6 | 2.1e-6 | 1.2e-6 | 6.6e-8 |
| 1.00 (the shock forms) | 1.2e-5 | 2.7e-6 | 1.7e-5 | 3.1e-6 |
| 1.50 | 5.5e-6 | 1.2e-7 | 3.0e-5 | 2.7e-7 |

The centred schemes are more accurate than their upwind parts until t = 0.85/pi. In the last steepening before the
shock forms (0.85/pi to 1/pi) every error grows: WENO6-C's 15 times, WENO5-JS's 1.5 times, WENO8-C's 120 times,
WENO7-JS's 97 times; when the shock forms the centred schemes are about 5 times above their upwind parts. The
difference that decides is afterwards: the upwind parts damp theirs by 23 and 12 times by t = 1.5/pi, the centred ones
keep theirs (WENO6-C) or let it grow (WENO8-C), and end 48 and 110 times above.

**What the packet is.** About 80 % of neighbouring cells in it alternate in sign (101 of 125 at 640 cells for WENO6-C,
128 of 152 for WENO8-C): theta/pi about 0.8 to 1. There the central rows' group velocity is reversed (FD flux form,
relative to the exact one): CEN6 -1.2 at 0.8 pi, -2.2 at pi; CEN8 -1.2 and -2.6. Such modes move away from where they
are made, against the flow, into the smooth region.

**phi where the packet is** (the final state of `adam`; mean over the cells above 0.1 of the largest error):

| | WENO6-C, 160 | WENO6-C, 640 | WENO8-C, 160 | WENO8-C, 640 |
|---|---|---|---|---|
| phi+ (mean; max) | 0.018; 0.096 | 0.013; 0.047 | 0.050; 0.24 | 0.013; 0.051 |
| phi+ of the exact solution, same faces (mean) | 0.016 | 0.0047 | 0.021 | 0.0038 |
| phi+ with eps = 1e-40 (mean) | 0.018 | 0.014 | 0.051 | 0.014 |
| (Delta f+)^2 / eps (mean) | 322 | 19 | 302 | 20 |

The detector sees the packet only as a small change of the background's differences: phi is 1.1 to 3.5 times the
exact solution's, and stays at 1e-2. eps is not what holds it there: the indicators are 6 to 400 times eps, and eps = 1e-40
changes phi by less than 15 %.

**The damping of a small packet** (linear advection u_t + u_x = 0, N = 160, u0 = sin(pi x) + A (-1)^i
exp(-((x + 0.5)/0.1)^2), the grid-scale part of u(A) - u(0) left after t = 0.25, 20 cells of travel, over A):

| A | WENO6-C, background | alone | phi = 1, background | phi = 0 |
|---|---|---|---|---|
| 1e-1 | 4.5e-3 | 1.1e-2 | 7.4e-4 | 0.98 |
| 1e-3 | 0.13 | 0.34 | 4.5e-2 | 0.98 |
| 1e-5 | 0.30 | 0.98 | 1.6e-3 | 0.98 |
| 1e-7 | 0.30 | 0.98 | 1.7e-5 | 0.98 |

(WENO8-C alike: 0.31 and 0.96 at A = 1e-7; its WENO7-JS part 4.4e-5.) A large packet is damped. A small one on a
smooth background keeps 30 % whatever its amplitude: the damping is the background's phi times the upwind scheme's.
Alone, a packet below sqrt(eps) = 1e-3 is not damped at all (0.98, as phi = 0): its differences are below eps.

## Conclusion

The pollution is made where the shock forms and kept because the blend barely damps it:

1. **Generation.** In the last steepening (t = 0.85/pi to 1/pi) the front is under-resolved and every scheme's error
   grows; the centred schemes leave about 5 times their upwind parts' when the shock forms, as a grid-scale packet.
2. **Transport.** At those wavenumbers the central rows' group velocity is reversed: the packet moves into the smooth
   region against the flow.
3. **Weak damping.** A small packet on a smooth background is damped at phi times the upwind rate, and phi is that of
   the background, O(h) (its mean over the smooth region falls from 0.017 at 160 cells to 0.0052 at 640). The upwind
   rate at a fixed wavenumber is O(1/h), so the blend damps the packet at a rate that does not grow with the
   resolution: the error falls at first order or below. The upwind WENO-JS damps the same packet at O(1/h).

Refuted: the upwind part (phi = 1 is clean); eps as the cause (phi in the packet is the same with eps = 1e-40; eps
matters only for a packet with no background); the detector's reach as the cause (`wide` changes the level by
1.6 times at most). TENO6's anti-dissipation does not occur.

The ADR (F2) cannot show this: it advects one wave of unit amplitude with no background, which the detector sees in
full (phi = O(1) at high wavenumbers).

## Decision (owner, 2026-11-16): accepted as documented

SETH keeps ADAM's scheme unchanged; the packet is a documented property of the centred WENO (S6, A5, the guide).
Whether it matters in turbulence (grid-scale energy that is never damped) is ADAM's to verify on its flows, under the
SETH/ADAM verification split.

Measured and not adopted (each a change of ADAM's scheme):

- eps = 1e-40 in the blend: 3.4e-7 and 9.2e-7 (7.5 and 3.9 times lower), about second order.
- eps = 1e-40 and the widened detector: 1.9e-7 and 3.2e-7 (13 and 11 times lower), about second order; it costs the
  smooth flow (before the shock, 640 cells: WENO6-C 1.8e-11 -> 6.2e-11) and was not run on the rest of the suite.

No variant reaches the upwind parts' level: the mechanism (damping proportional to the background's phi) is the
blend's construction. The principled alternative, a floor on phi (equivalently a dissipative smooth branch), costs
one order (2S -> 2S-1) and would be a new scheme rather than ADAM's; not tried.
