---
title: Cookbook
---

# Cookbook

Short recipes for specific tasks. Most are a recipe of one program in `docs/examples/src/`, chosen on the command
line; the code and output on each page are that program and what it printed when the documentation was built. The
first and the last recipes are references, with the numbers taken from the verification and the benchmarks.

| Recipe | Task |
|---|---|
| [Choose a scheme](./choose) | which kernel for which data, order, robustness and cost |
| [Point values or cell averages](./point-values) | reconstruction or interpolation: the kernel must match the data |
| [Flag the shocks](./shock-flags) | the weights ratios of the `_wr` kernels, a sensor already computed |
| [Bound a reconstruction](./mp-bound) | the monotonicity-preserving bound after any kernel, inactive on smooth extrema |
| [A central branch for the Euler equations](./central) | the KEEP-PE flux and the Ducros sensors |
| [Single precision](./precision) | the R4P build, and what changes with it |
| [Time a kernel on a GPU](./performance) | benchmark discipline, and the costs measured |

The recipes that the tutorial already covers: a descaler for data of any magnitude ([4](/tutorial/04-weights)), a
system by characteristics ([7](/tutorial/07-euler)), positivity ([8](/tutorial/08-positivity)), a hybrid
([9](/tutorial/09-hybrid)), a GPU loop ([11](/tutorial/11-gpu)).
