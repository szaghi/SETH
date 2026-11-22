---
title: WENO-6 point interpolation
---

# WENO-6 point interpolation

The value at an off-node point of a line of six equispaced point values, non-oscillatory: a WENO-Z combination of
the cubic interpolants of the 4-point sub-stencils that contain the target's interval. Unlike every other SETH kernel
it does not compute the two states at a face: it takes one array of point values and a target anywhere between the
nodes, on a quarter of the node spacing. ADAM's 2:1 AMR seam transfer is this kernel (ADAM issue 68, C3 and D2): the
fine-to-coarse restriction at the midpoints, the coarse-to-fine fill at the quarters. Module `seth_point_weno6`, its
weights in `seth_point_weno6_tables` (written by `scripts/gen_weno6_tables.py`).

## Calls

```fortran
v = seth_weno6_point(u, t4)     ! real(R_P) :: u(1:6), v;  integer(I4P) :: t4 (3..25)
call seth_weno6_jumps(u, jump)  ! real(R_P) :: u(1:6);     logical :: jump(1:5)
```

| Argument | Meaning |
|---|---|
| `u(1:6)` | the values at six equispaced nodes, node $k$ at $x = k$ |
| `t4` | the target, at $x = t_4/4$ node spacings from node 0: $t_4 = 3..25$ |
| `v` | the interpolated value |
| `jump(i)` | interval $i$ (nodes $i$, $i+1$) is a jump |

| $t_4$ | Target | Candidates | Order on smooth data |
|---|---|---|---|
| odd: 3, 5, …, 25 | the quarter positions $x = 3/4 … 25/4$ | the sub-stencils containing the target's interval | 6 at $x = 13/4$, $15/4$; 5 at $9/4$, $11/4$, $17/4$, $19/4$; 4 elsewhere |
| $2 \bmod 4$: 6, 10, 14, 18, 22 | the midpoints between nodes 1–2, …, 5–6 | idem | 6 at $x = 14/4$; 5 at $10/4$, $18/4$; 4 at $6/4$, $22/4$ |
| $0 \bmod 4$: 4, 8, …, 24 | node $t_4/4$ | none | its value, unchanged |

Other values of `t4` are outside the contract: no check is made in the kernel (it is device code), and a
bounds-checked build stops on the table index. At $t_4 = 3$ and $25$, a quarter beyond the end nodes, the one
candidate is the one-sided cubic of the end interval. The order 4 toward the ends is the price of never
extrapolating across a jump.

## The scheme

Sub-stencil $s$ is the nodes $s..s+3$ ($s = 1, 2, 3$); its cubic interpolant at the target is $q_s$, with the weights
`WENO6_POINT_W(:,s,t4)`. The target's interval is $m = \mathrm{clamp}(\lfloor t_4/4 \rfloor, 1, 5)$.

1. **Jumps** (`seth_weno6_jumps`): interval $i$ is a jump when $|u_{i+1} - u_i|$ exceeds `SETH_WENO6_JUMP_RATIO`
   $= 2$ times the larger of its neighbours' undivided differences (at the ends, the inner neighbour and its linear
   extrapolation). Resolved smooth data never trip it: the neighbours' differences match the interval's up to
   $O(h^3)$.
2. **Candidates**: the sub-stencils with a positive linear weight $d_s$ (`WENO6_POINT_D(s,t4)`, positive exactly for
   the sub-stencils containing the interval $m$) and no jump. The $d_s$ make $\sum_s d_s q_s$ the Lagrange interpolant
   of the candidates' union, and all are positive.
3. **No candidate** (a jump inside the interval $m$): the troubled-cell fallback, $\tfrac12(u_m + u_{m+1})$ at a
   midpoint, the linear interpolation on the interval clamped to $[\min(u_m, u_{m+1}), \max(u_m, u_{m+1})]$ at a
   quarter.
4. **One candidate**: its cubic.
5. **Several**: WENO-Z weights (Borges, Carmona, Costa & Don 2008, $p = 2$),

$$ \alpha_s = d_s\left(1 + \left(\frac{\tau}{\beta_s + \varepsilon}\right)^2\right), \qquad
   v = \frac{\sum_s \alpha_s q_s}{\sum_s \alpha_s}, $$

on the Jiang–Shu indicator of each cubic over its central interval (Jiang & Shu 1996),

$$ \beta_s = c^2 + \tfrac{13}{12} a^2 + \tfrac{781}{720} t^2, \quad
   c = u_{s+2} - u_{s+1},\; a = \tfrac12(u_{s+3} - u_{s+2} - u_{s+1} + u_s),\;
   t = u_{s+3} - 3u_{s+2} + 3u_{s+1} - u_s, $$

with $\tau$ the square of the highest undivided difference of the candidates' union ($O(h^{8})$ to $O(h^{12})$ on
smooth data, ADAM's choice so that the weights reach the linear ones fast; the orders at critical points are not
measured here) and $\varepsilon$ = `SETH_WENO6_EPS`.

## Constants

| Constant | Value | Meaning |
|---|---|---|
| `SETH_WENO6_JUMP_RATIO` | 2 | the jump detector's ratio (ADAM's) |
| `SETH_WENO6_EPS` | $10^{-40}$ (R8P), $10^{-30}$ (R4P) | the guard of the weights against $0/0$ on constant data; per kind, in `seth_kind_parameters` ([Kinds and constants](./kinds)) |

In R4P no $\varepsilon$ small enough to be negligible keeps $(\tau/(\beta_s+\varepsilon))^2$ finite next to a flat
sub-stencil ($\beta_s = 0$: with $u = 1, 1, 1, 1, 1.001, 1.002$ it is $\infty$ and the weights NaN; without the scaling
below, 6 of the 23 targets return NaN). The R4P kernel therefore scales every
$\alpha_s$ by $(b_{\min}/\tau)^2$, $b_{\min} = \min_s(\beta_s + \varepsilon)$, whenever $\tau > b_{\min}$:
$\alpha_s = d_s\left((b_{\min}/\tau)^2 + (b_{\min}/(\beta_s+\varepsilon))^2\right)$, each factor in $[0, 1]$ and the
normalised weights unchanged. The R8P kernel keeps ADAM's formula, operation for operation.

## Equivalence with ADAM

In R8P, `seth_weno6_point(u, 4*code + 6)` is ADAM's `seam_restriction_weno6_1d(u, code)` ($code = 0..4$) and
`seth_weno6_point(u, 4*p + 2*sub - 3)` its `seam_fill_weno6_1d(u, sub, p)` ($p = 1..6$, $sub = 1, 2$), bit for bit with
gfortran, debug and optimised; `seth_weno6_jumps` is `seam_weno6_jumps`. The tables are derived from the definition
by the generator, which also asserts that they equal ADAM's four, and each entry is written as ADAM writes it (an
integer over 16, 128 or 320), so that it folds to the same double ([C1](/verification/contracts/adam-equivalence)).
The properties: `seth_test_weno6_point` ([Test programs](/verification/unit-tests)). API:
[`seth_point_weno6`](/api/src/lib/seth_point_weno6).
