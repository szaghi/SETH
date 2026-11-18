---
title: Single precision
---

# Single precision

SETH computes in one real kind per build, `R_P`. Write your code in `R_P` too, and the same source builds in either
precision:

<<< @/examples/snippets/cookbook_kinds-kinds.F90{fortran}

The documentation's examples are built in double precision:

<<< @/examples/output/cookbook-kinds.txt{text}

For single precision, build the library and your program with the `*-r4` modes, which define `_R_P_IS_R4P` for PENF
and SETH alike:

```bash
fobis build --mode static-gnu-r4      # also static-nvf-r4, static-intel-r4, shared-*-r4
```

What changes: `SETH_TENO_EPS` becomes $10^{-30}$ (the papers' $10^{-40}$ is not a normal R4P number); the weights are
computed as ratios to the smallest indicator, so they do not overflow; data with a large mean needs the descaler, and
a mean subtracted first for the affine-invariant one; schemes of order 6 or more reach R4P's round-off within the
first refinements of a smooth problem. What it costs and buys, measured: [Precision](/guide/precision).
