---
title: 2. A line of faces
---

# 2. A line of faces

A grid of $N$ cells on the periodic interval $[0, 1]$, the cell averages of $\sin 2\pi x$, and every face
reconstructed. The error of the left states against the exact point values measures the order of accuracy.

## Ghost cells

The face $x_{i+1/2}$ needs the cells $i-2 \dots i+3$. For the faces at the ends of the grid those cells lie outside
it: the program extends the array by `NG = 3` ghost cells on each side and fills them periodically.

<<< @/examples/snippets/t02_line-fill.F90{fortran}

The number of ghost cells is the scheme's: 3 for WENO5 and TENO5, $S$ for WENO-$(2S-1)$, 4 for TENO8, 6 for the
composite non-polynomial kernels ([Reference](/reference/#every-kernel-at-a-glance)).

## The face loop

The kernel goes inside your own loop. The slice `u(i-2:i+3)` is passed to a dummy declared `v(-2:3)`: inside the
kernel, cell 0 is cell $i$.

<<< @/examples/snippets/t02_line-faces.F90{fortran}

- WENO-Z and TENO take `SETH_TENO_EPS` ($10^{-40}$), not WENO-JS's $10^{-6}$.
- `1` is WENO-Z's exponent $p$; `SETH_TENO5_CT` is TENO5's cutoff.

## Running it

<<< @/examples/output/t02_line.txt{text}

All three schemes converge at fifth order. WENO5-JS's error is ten times the others': at the two extrema of the sine
its weights stray from the linear ones by more than the design allows (a critical point). WENO-Z and TENO were
built against this ([A2](/verification/accuracy/critical-points)).

::: tip What you learned
A field is a loop of kernel calls over its faces, with as many ghost cells as the scheme reads. Reference:
[WENO-Z](/reference/weno-z), [TENO](/reference/teno).
:::

<details><summary>The whole program</summary>

<<< @/examples/snippets/t02_line.F90{fortran}

</details>

Next: [3. A conservation law](./03-advection).
