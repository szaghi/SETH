# Conventions

## One face, two stencils

Every SETH kernel computes the two states at **one face** from **two stencils**. Cells are indexed by their offset
from the face $x_{i+1/2}$: cell $0$ is the cell left of it, cell $1$ the cell right of it.

```
      -2     -1      0  |   1      2      3
...|------|------|------|------|------|------|...
                      face i+1/2
```

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

## Mirror symmetry

With `vm(m) = vp(1-m)` the right state equals the left one. WENO-Z and TENO compute both states with the same
arithmetic on mirrored data, so the equality is **bitwise**; WENO-JS uses ADAM's mirrored tables, summed in another
order, and agrees to round-off.

## No setup, no state

The coefficients are compile-time constants. A kernel needs no table, no object and no device copy: call it inside
your own loops, on the host or on the device. Every kernel is `pure`, `!$acc routine seq` and `!$omp declare target`.

## From ADAM's rows

ADAM's kernels take `v(1:2, 1-S:S-1)`: row 2 the left-biased stencil, row 1 the right-biased one shifted one cell to
the left, and return the left state in `vr(2)`, the right one in `vr(1)`. The mapping is

$$ v(2,m) = v_p(m),\qquad v(1,m) = v_m(m+1),\qquad m = 1-S,\dots,S-1, $$

and SETH's WENO-JS kernels return ADAM's values bitwise (`seth_test_weno_js_adam`).
