---
title: Point values or cell averages
---

# Point values or cell averages

A reconstruction kernel takes **cell averages** and returns the point value at the face: the finite-volume setting,
and the finite-difference one where the point values of a split flux play the role of averages (Shu & Osher). An
interpolation kernel takes **point values** at the cell centres and returns the point value at the face: WCNS, and the
alternative finite-difference formulation (Jiang, Shu & Zhang 2013). Given the wrong data, a kernel loses its order.

<<< @/examples/snippets/cookbook_schemes-interp.F90{fortran}

<<< @/examples/output/cookbook-interp.txt{text}

The reconstruction kernel on point values converges at second order: the difference between a point value and a
cell average is $h^2 f''/24$. The interpolation kernel gives fifth order. Every reconstruction kernel of WENO-JS,
WENO-Z, TENO5, TENO6 and the linear schemes has an `_interp` twin ([Reference](/reference/)).

<details><summary>The whole program</summary>

<<< @/examples/snippets/cookbook_schemes.F90{fortran}

</details>
