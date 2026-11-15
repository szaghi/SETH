# Positivity limiter

| ID | Group | Type | Backends | Figure | Programs |
|---|---|---|---|---|---|
| C7 | Contracts | C | CPU (gfortran, nvfortran, ifx; debug and optimised), GPU | no | `seth_test_positivity`, `seth_test_device`; in use: `seth_vv_test123`, `seth_vv_blast` |

## Case and scope

`seth_positivity_euler` is the positivity-preserving flux limiter of Hu, Adams & Shu (2013) at one face ([API](/guide/api#positivity-limiter);
llm-wiki high-order-capturing-schemes `methods/positivity-preserving-flux-limiter.md`). It must return a flux whose
one-sided states are admissible, with the largest admissible blending factor, and leave admissible fluxes untouched.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Limited flux | $F = \theta F_{\rm high} + (1 - \theta)F_{\rm low}$, the one-sided states $q_l - \lambda F$, $q_r + \lambda F$ above the floors | [`seth_positivity.F90`](https://github.com/szaghi/SETH/blob/main/src/lib/seth_positivity.F90) |
| $\theta$ | closed form: a ratio for the density, the smallest root of $\rho E - \tfrac12|m|^2 - \varepsilon_p\rho/(\gamma - 1)$ for the pressure (the admissible set is convex, so the admissible factors form an interval $[0, \theta^*]$) | `segment_bound` |
| Rounding guard | the one-sided states of the returned flux are re-checked; a state a few ulp below its floor lowers $\theta$ by relative steps $10^{-15}, 2\cdot10^{-15}, …$ | `seth_positivity_euler` |
| Layout | $q = (\rho, \rho u_1, …, \rho u_d, E)$, `nv` $= 3, 4, 5$; ideal gas | |
| Device rule | no local array sized by `nv` (dimension `SETH_POSITIVITY_NV_MAX` = 5, loops to `nv`): automatic arrays made the first version 333× slower on the GPU ([performance](/guide/performance#device-kernel-rules-arrays-and-loop-extents)) | |

## Setup

`seth_test_positivity` draws random cell pairs (fixed seed): densities and pressures log-uniform in $[10^{-8}, 1]$,
velocities up to 10 (near vacuum: kinetic energy far above the internal one), `nv` $= 3, 4, 5$. The low-order flux
is Lax–Friedrichs with the maximal speed of the pair and $\lambda = 1/\alpha$ (CFL $\tfrac12$, the limit of Hu, Adams
& Shu); the high-order flux is that flux plus a random perturbation of up to $10\times$ the fluxes' magnitude.
Asserted, every check with the test's own admissibility function so that a defect of the library's cannot hide itself:

| Check | Tolerance |
|---|---|
| $\theta \in [0, 1]$ and $F = \theta F_{\rm high} + (1 - \theta)F_{\rm low}$ | round-off (measured 0) |
| the limited one-sided states admissible ($\rho, p \ge 10^{-13}$) | the round-off of the state as the test evaluates it (see below) |
| $\theta$ maximal: within $10^{-9}$ of the largest admissible factor found by 60 bisection steps | measured $\le 8\cdot10^{-12}$ |
| $\theta = 1$ and $F = F_{\rm high}$ bitwise when $F_{\rm high}$ is admissible | bitwise |
| $\theta = 0$ and $F = F_{\rm low}$ when even $F_{\rm low}$ fails ($\lambda$ ten times too large) | bitwise |
| the Lax–Friedrichs flux admissible at CFL $\tfrac12$ (the premise, Hu, Adams & Shu eq. 13) | exact |
| `seth_euler_admissible` agrees with the test's own check | exact |

The kernel guarantees its floors on the one-sided states as it evaluates them. The test evaluates them again, and
with nvfortran `-fast` the two evaluations differ (a fused multiply-add) by $\epsilon(|q| + |\lambda f|)$; near vacuum
the pressure of the test's state then falls up to $8\cdot10^{-16}$ below the floor, and in one case of 60 000 (density
$10^{-7}$ at $|u| = 400$, next to a cell $10^9$ times denser) to $-10^{-12}$. The check allows that round-off
(`adm_ro`); the [API](/guide/api#positivity-limiter) says what it means for the choice of the floors.

About 26 000 limited cases per run (8 300–9 100 per layout); seeds 1006 (default), 7, 99, 12345 checked.

`seth_test_device` runs the kernel inside an OpenACC loop on 4096 pairs per layout (`nv` = 3, 5) against the host:
$\theta$ and $F$ within $10^{-10}$ (the closed form's square root and the multiply-adds may round differently),
every device flux admissible (to the same round-off). On the WSL RTX 4070, 2026-11-15: worst difference $4.1\cdot10^{-14}$, 2228 and 1866 of
4096 bitwise.

## In use

The verification harness limits with the kernel ([Sod](../shocks/sod#models-and-numerics)). On the
[123 problem](../shocks/test123) and the [blast wave](../shocks/blast) every scheme, which without a limiter often loses
positivity, completes the run; on smooth flows ([Euler order](../accuracy/euler-order)) and on Sod and Lax the limiter
never acts (asserted).

## Mutation testing

Four mutants of the kernel (the pressure quadratic, the density bound, the pressure floor, a one-sided state) are
each caught ([C6](./mutation)).
