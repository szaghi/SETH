# Conventions

## One face, two stencils

Every SETH kernel computes the two states at **one face** from **two stencils**. Cells are indexed by their offset
from the face $x_{i+1/2}$: cell $0$ is the cell left of it, cell $1$ the cell right of it.

![One face of WENO5: the offsets, vp and vm with the cells each reads, the candidate stencils and their linear weights](/guide/conv-face.png)

*WENO5 ($S = 3$): both arrays are `(-2:3)`. The left state `rp` reads `vp` on cells $-2..2$ through three candidate
stencils; the right state `rm` reads `vm` on cells $-1..3$ through their mirrors. The linear weights $d_k$ are the
ones the nonlinear weights tend to on smooth data.*

| Argument | Meaning |
|---|---|
| `vp(lo:hi)` | values reconstructed with a **left** bias: the cell averages for the left state, or $f^+$ of a flux splitting |
| `vm(lo:hi)` | values reconstructed with a **right** bias: the cell averages for the right state, or $f^-$ |
| `rp` | the left state at the face, from `vp` |
| `rm` | the right state at the face, from `vm` |

Both arrays use the **same indexing, unshifted**; the kernel mirrors `vm` internally (cell $m$ of the right-biased
reconstruction is cell $1-m$ of the left-biased one). The bounds are the scheme's full stencil, symmetric about the
face: `1-S:S` for WENO-$(2S-1)$, `-2:3` for TENO5 and TENO6.

- To reconstruct the two states of one field, pass the same array twice: `call seth_weno5_js(u, u, ...)`.
- For a flux-split scheme pass $f^+$ and $f^-$; the numerical flux is `rp + rm`.
- With finite differences (point values $f_i$ treated as the averages of the flux function), the same kernels give
  the flux at the face; `seth_test_advection` and `seth_test_critical_points` use them so.

### One face at a time

<div class="seth-casts">
<figure class="seth-cast seth-cast-wide" style="--ar: 1.565">
<video class="light-only" autoplay loop muted playsinline poster="/guide/casts/conv-sweep-light.png" aria-label="WENO5 sweeping a field with a jump"><source src="/guide/casts/conv-sweep-light.webm" type="video/webm"><source src="/guide/casts/conv-sweep-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/guide/casts/conv-sweep-dark.png" aria-label="WENO5 sweeping a field with a jump"><source src="/guide/casts/conv-sweep-dark.webm" type="video/webm"><source src="/guide/casts/conv-sweep-dark.mp4" type="video/mp4"></video>
<figcaption>WENO5-JS sweeping the faces of 32 cell averages, a smooth wave with a jump. At each face, the stencils of
<code>vp</code> (blue, upper band) and <code>vm</code> (green, lower band), the candidate polynomials (drawn thicker as their weight
grows), and the two states ◀ <code>rp</code> and ▶ <code>rm</code>; below, the nonlinear weights against the linear
ones. On the smooth wave the weights stay near <i>d<sub>k</sub></i> and the two states agree; next to the jump each state puts its
weight on the candidate that does not cross it.</figcaption>
</figure>
</div>

## Every scheme's stencil

![The bounds of vp and vm and the cells each reads, for every scheme; TENO's candidate stencils](/guide/conv-schemes.png)

*Read from the coefficient tables (`seth_weno_js_tables`, `seth_z_teno_tables`). WENO-JS and WENO-Z read the same
cells at equal order. TENO5 reads the cells of WENO5 through the same three candidates, numbered as Fu, Hu & Adams
number them; TENO6 adds a fourth, four cells wide, that reaches cell $3$, so both its arrays are read in full.*

## Mirror symmetry

![A profile vp and its mirror vm(m) = vp(1-m): the left state of one equals the right state of the other](/guide/conv-mirror.png)

With `vm(m) = vp(1-m)` the right state equals the left one. WENO-Z and TENO compute both states with the same
arithmetic on mirrored data, so the equality is **bitwise**; WENO-JS uses ADAM's mirrored tables, summed in another
order, and agrees to round-off.

## No setup, no state

The coefficients are compile-time constants. A kernel needs no table, no object and no device copy: call it inside
your own loops, on the host or on the device. Every kernel is `pure`, `!$acc routine seq` and `!$omp declare target`.

## From ADAM's rows

ADAM's kernels take `v(1:2, 1-S:S-1)`: row 2 the left-biased stencil, row 1 the right-biased one shifted one cell to
the left, and return the left state in `vr(2)`, the right one in `vr(1)`.

![SETH's vp and vm against ADAM's rows v(1:2, -2:2) for S = 3](/guide/conv-adam-rows.png)

The mapping is

$$ v(2,m) = v_p(m),\qquad v(1,m) = v_m(m+1),\qquad m = 1-S,\dots,S-1, $$

and SETH's WENO-JS kernels return ADAM's values bitwise (`seth_test_weno_js_adam`).

---

*The figures and the animation are drawn by `docs/scripts/make_convention_figures.py` (light and dark, with the figure
harness of the ADAM and FLUME documentation). The stencil map reads the tables of the Fortran sources; the animation
runs WENO5-JS with SETH's coefficients, and its states equal `seth_weno5_js` to $3.3\cdot10^{-16}$.*
