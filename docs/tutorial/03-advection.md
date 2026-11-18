---
title: 3. A conservation law
---

# 3. A conservation law

The linear advection $u_t + (a u)_x = 0$, $a = 1$, on the periodic unit interval, by finite volumes: the cell
averages change by the difference of the fluxes at their two faces,

$$ \frac{d\bar u_i}{dt} = -\frac{F_{i+1/2} - F_{i-1/2}}{h}. $$

## Face fluxes

The face flux comes from a **flux splitting**: Lax–Friedrichs writes $f = f^+ + f^-$, with
$f^\pm = \tfrac12(f \pm \alpha u)$ and $\alpha \ge \max|f'|$, so that $f^+$ travels right and $f^-$ left. $f^+$ is
reconstructed from the left (`vp`) and $f^-$ from the right (`vm`), and the face flux is `rp + rm`:

<<< @/examples/snippets/t03_advection-residual.F90{fortran}

Here $f' = a$ is constant and $\alpha = |a|$: $f^- = 0$ and $f^+ = au$. With a nonlinear flux both parts
are nonzero, and $\alpha$ a little above $\max|f'|$ avoids a loss of order where $|f'|$ peaks
([S6](/verification/shocks/burgers)).

## Time stepping

The three-stage strong-stability-preserving Runge–Kutta scheme of Shu and Osher, with CFL 0.4:

<<< @/examples/snippets/t03_advection-rk3.F90{fortran}

## Running it

The initial profile is a square wave on $[0.2, 0.4]$ and a Gaussian bump at 0.7; after one period the exact solution
is the initial one.

<<< @/examples/output/t03_advection.txt{text}

- The **linear** scheme oscillates at the edges of the square wave: 7.5% of the jump below 0 and above 1.
- **WENO5-JS** keeps the over- and undershoots to $3\cdot10^{-4}$, at the price of more dissipation.
- **TENO5** is both sharper (the smallest $L_1$ error) and nearly free of oscillations: on smooth parts it is its
  linear scheme exactly, and at the edges it drops the candidates that cross them.

::: tip What you learned
A finite-volume solver is a flux splitting, a face loop calling the kernel on $f^+$ and $f^-$, and a Runge–Kutta
scheme. SETH supplies the face values; the splitting, the fluxes and the time stepping are yours.
:::

<details><summary>The whole program</summary>

<<< @/examples/snippets/t03_advection.F90{fortran}

</details>

Next: [4. The weights](./04-weights).
