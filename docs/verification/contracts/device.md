# Device

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| C5 | Contracts | C | GPU (nvfortran OpenACC, `tests-nvf-oac`); host fallback elsewhere | no | `seth_test_device` |

## Case and scope

Every kernel is `pure`, `!$acc routine seq` and `!$omp declare target`, meant to run inside a client's device loop.
The test runs each one inside an `!$acc parallel loop` over 4096 deterministic stencils (smooth cubics, jumps, an
integer LCG), then the same loop on the host, and compares. Built without OpenACC the directives are ignored and the
comparison is trivial; built with `tests-nvf-oac` it is the device test.

## Setup

| Item | Value |
|---|---|
| Cases | 63 cases: the 20 WENO-JS entry points, WENO-JS with the descaler called on the device, WENO-Z3/5/7/9, TENO5/6/6-opt/6-NENO6, the six WENO-Z and TENO interpolations (35, the v0.1.0 set); TENO5-opt, TENO5-AZ, TENO8, TENO5-A, TENO8-A, the affine-invariant descaler, the MP bound, the linear kernels, two sensors, the centred WENO4/6/8-C, Fu's sensor eps (36–50); the non-polynomial kernels of issue #6: the mask, TENO5-DR, TENO5-THINC, TENO5-ST, TENO5-RBF, TENO5-RBF+DR (51–56), the two detectors with $\zeta$ and the two overrides (57–59, each asserted to act on some stencils); the stored-mask pair, `seth_teno5_face_masks` and the three `*_stored` kernels with their masks computed on the device (60–63); and the positivity limiter on 4096 cell pairs for `nv` = 3 and 5 ([C7](./positivity)) |
| Tolerance | states within $10^{-13}$ of the stencil's magnitude, weights ratios within $10^{-9}$ (the device and the host compile the kernel separately and contract multiply-adds differently: the tolerances of the golden test) |
| Run | `fobis build --mode tests-nvf-oac --varset cc89`, `CUDA_VISIBLE_DEVICES=1 ./exe/seth_test_device`, on an idle GPU |

## Results

On the WSL RTX 4070 (GPU 1), nvfortran 26.9, 2026-11-15: 35 of 35 within tolerance, 26 bitwise on all 4096 stencils,
worst state difference $7.9\cdot10^{-16}$ of the data (`notes/device-gpu1-2026-11-15.txt`). The whole suite also
passes in that build (25 of 25).
