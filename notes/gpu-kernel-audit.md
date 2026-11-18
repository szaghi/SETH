# Device-kernel audit: automatic arrays and loop extents

Night 2026-11-15 (owner's request). Two rules of ADAM's device kernels, checked on every routine SETH offloads:

1. **No automatic local arrays** in a routine called on the device (a local whose bounds depend on a dummy argument,
   e.g. `q(1-S:S)` with `S` a dummy). Over-dimension it with a constant (`S_MAX`).
2. **Loop extents known at compile time** (`do i=1, S` with `S` a named constant).

Program: `src/bench/seth_bench_device_arrays.F90`, modes `bench-arrays-nvf-oac`, `bench-arrays-nvf`,
`bench-arrays-gnu`. Raw outputs: `notes/gpu-audit-data/`.

## Result

- **The WENO-JS, WENO-Z and TENO kernels already follow both rules**, by construction: one kernel per scheme, every
  bound and extent a named constant of that kernel. Nothing to change there (the M3 benchmark already had them
  1.26–1.42× faster than ADAM's kernels, `m3-gpu-benchmark.md`).
- **`seth_positivity_euler` broke rule 1** (written the night before, commit ea4959b): four automatic arrays
  `a0..b1(nv)` and array expressions over `nv` as arguments. On the GPU it was **333× slower** than it should be
  (899 against 2.70 ns per face, nv = 3; 278× for nv = 5), and with the default device heap a kernel with automatic
  arrays can **fail to launch** (`DEVICE FORTRAN AUTO ALLOCATION FAILED`: seen on the WENO-5 copy with automatic
  arrays at 2^24 faces; the limiter's own threshold was not searched). Fixed in this night: locals dimensioned
  `SETH_POSITIVITY_NV_MAX = 5`, explicit loops to `nv`, API unchanged.
- **Rule 2 is worth 4–7 %**, not orders of magnitude: WENO-5 with over-dimensioned arrays and run-time loop extents
  is 7 % slower than the library kernel. For the limiter, loops to the dummy `nv` cost nothing measurable against
  `nv` constant (≤ 1 %), so no specialised 1-D/2-D/3-D kernels are needed.
- **`seth_weno_descaling(S, ..)`** has a run-time loop extent and no local array: 0.8 % against a constant-`S` copy.
  Left as it is.

## Static audit

Every `!$acc routine seq` of `src/lib`:

| Routine(s) | Local arrays | Loop extents | Verdict |
|---|---|---|---|
| `seth_weno{1,3,5,7,9}_js` (+ `_wr`, `_interp`, `_interp_wr`): `seth_weno_js_body.INC` | `v(1:2,1-S:S-1)`, `q`, `is`, `al`, `w(1:2,0:S-1)`: `S` a named constant of each kernel | `S`, constant | follows both rules |
| `seth_wenoz{3,5,7,9}`, `seth_teno5/6/6_opt/6_neno6` (+ `_interp`): `seth_z_teno_body.INC` | `st(LO:UMAX)`, `q`, `beta`, `wgt`, `g(0:NC-1)`: constants | `NC`, `LO`, `UMAX`, `GW` constants; the inner loops run to `CW(k)`, an element of a constant array | follows both rules (see notes) |
| `seth_weno_descaling(S, ..)` | none | `S`, a dummy | rule 2 not followed; 0.8 % |
| `seth_positivity_euler(nv, ..)`, `segment_bound`, `seth_euler_admissible` | **before:** `a0, a1, b0, b1(nv)`, `d(nv)` automatic; `ql - lam*f` temporaries; `dot_product` over `2:nv-1`. **now:** dimension 5 | `nv`, a dummy | rule 1 was broken, fixed; rule 2 not followed, ≤ 1 % |
| `seth_weno{4,6,8}_c`, `blend` (`seth_weno_c.F90`, 2026-11-16) | `C(CEN<2S>_LO:CEN<2S>_HI)`, a parameter; scalars only otherwise | `CEN<2S>_LO`, `_HI`, constants; the upwind part is `seth_weno{3,5,7}_js` (above) | follows both rules; no array expression passed (the benchmark's WENO4-C copies its stencil into `(-1:2)` privates); the R8P blend (one division, 2026-11-16) adds four scalars and a logical, no array |

| `seth_nonpoly.F90` (issue #6, 2026-11-17; audited 2026-11-18): `seth_teno5_mask`, `seth_teno5_zeta`, `seth_dr_detect`, `seth_thinc_detect`, `seth_dr_override`, `seth_thinc_override`, `seth_teno5_dr`, `seth_teno5_thinc`, `seth_teno5_st`, `seth_teno5_rbf`, `seth_teno5_rbf_dr`, and the private `side_level`, `cell_mask`, `dr_face`, `thinc_face`, `dr_level`, `thinc_level` | `a(-2:3)`, `b(-2:3)`, `w(-5:5)`, `c(-1:2)`, `d(-1:2)`, `mk(-3:3)`, `u(-2:2)`, `z(-3:3)`: constant bounds; the tables `ROWP`, `ROWC` parameters | constant, but **`side_level` ran its downwind loop to `merge(3, 2, thinc)`**, a run-time extent: now `do c=1, 3` with an `exit` at 3 for DR (bitwise the same results); no `optional` dummy; the composites pass contiguous sections (`vp(-1:2)`) to explicit-shape dummies, no array expression | follows both rules since 2026-11-18. `nvfortran -acc -gpu=cc89 -Minfo=accel`: every routine "Generating acc routine seq / NVIDIA GPU code", no device allocation; `side_level` reports "reference argument passing prevents parallelization: c" (a seq routine: informational). Output in `gpu-audit-data/minfo-nonpoly-2026-11-18.txt` |

The dummy arrays (`vp(1-S:S)`, `ql(nv)`) are explicit-shape dummies, not automatic arrays: they are the caller's
storage. In the callers (the benchmark, the device test) the stencil copies are `private` arrays of constant bounds.

Notes on the Z/TENO body, not measured tonight:

- the inner loops `do a=0, CW(k)-1` have an extent that is a constant array element indexed by the outer loop
  variable: known at compile time, but only after the outer loop is unrolled;
- WENO-Z selects on the dummy `p` (1 or 2) once per candidate.

Both are candidates for a later measurement; `seth_wenoz5` already runs at the speed of `seth_weno5_js`
(1.02 against 0.97 ns per face, M3).

## Measurements

RTX 4070 (WSL GPU 1, idle), nvfortran 26.9 `-O4 -fast -acc -gpu=cc89`, 2^24 faces (2^22 for the limiter), median of
20 passes after 2 warm-up passes, `!$acc wait` + `system_clock`. `S` and `nv` are read from `volatile` variables, so
the compiler cannot fold them. The WENO variants are copies of the library body with the bounds and extents as
macros (`seth_bench_weno_js_dyn_body.INC`); each equals the library kernel on every face.

### WENO-JS 5 (ns per face)

| Variant | Local arrays | Loop extents | ns/face | vs library |
|---|---|---|---|---|
| `seth_weno5_js` (library) | constant | constant | 0.97 | 1.00 |
| `weno5_js_smax` | constant, `S_MAX = 5` | run-time `S` | 1.04 | 1.07 |
| `weno5_js_auto_cloop` | automatic | constant | 553 | 569 |
| `weno5_js_auto` | automatic | run-time `S` | 528 | 543 |
| `seth_weno5_js` + `seth_weno_descaling(S, ..)` | | run-time `S` | 1.086 | 1.00 |
| `seth_weno5_js` + descaler with `S` constant | | constant | 1.077 | 0.99 |

The automatic-array rows need `NVCOMPILER_ACC_CUDA_HEAPSIZE`: with the default heap the launch fails. 256 MiB, 1 GiB
and 4 GiB give the same times (553, 555, 555 ns): the cost is the allocation at every call, not the heap's size.

### Positivity limiter (ns per face; 34 % of the faces limited for nv = 3, 45 % for nv = 5)

| Variant | nv = 3 | nv = 5 |
|---|---|---|
| first version: automatic arrays, array expressions (ea4959b, kept as `seth_bench_positivity_auto`) | 899 | 979 |
| **library now: locals of dimension 5, loops to `nv`** | **2.70** | **3.52** |
| bench copy with `nv` constant (both rules) | 2.67 | 3.53 |
| speed-up of the fix | 333× | 278× |

Three runs of the fixed library agree within 0.4 % (`bench-gpu1-fixed.txt`).

### CPU, one core

One i9-13900KF core, 2^22 faces (2^20 for the limiter), median of 10 passes; ns per face and ratio to the library:

| Variant | nvfortran `-O4 -fast` | gfortran `-O3` |
|---|---|---|
| `seth_weno5_js` (library) | 28.8 | 23.5 |
| WENO-5, `S_MAX` arrays, run-time loops | 44.9 (1.56×) | 41.9 (1.78×) |
| WENO-5, automatic arrays, constant loops | 50.1 (1.74×) | 31.5 (1.34×) |
| WENO-5, automatic arrays, run-time loops | 72.9 (2.53×) | 96.5 (4.10×) |
| descaler, `S` dummy / constant | 39.4 / 39.1 | 23.9 / 24.0 |
| limiter nv = 3: first version (automatic) | 94.1 (2.34×) | 66.4 (2.28×) |
| limiter nv = 3: **library now** | 40.1 | 29.1 |
| limiter nv = 3: `nv` constant | 36.9 (0.92×) | 26.0 (0.89×) |
| limiter nv = 5: first version (automatic) | 115.4 (1.90×) | 94.4 (2.14×) |
| limiter nv = 5: **library now** | 60.8 | 44.1 |
| limiter nv = 5: `nv` constant | 43.9 (0.72×) | 42.9 (0.97×) |

On a CPU core the two rules weigh differently from the GPU: automatic arrays cost 1.3–1.7× (a heap or stack
allocation per call, not a device-heap one), and run-time loop extents 1.6–1.8×, since the compiler can no longer
unroll and vectorise the fixed-size loops; together 2.5–4.1×. The fix of the limiter is worth 1.9–2.3× on the CPU. A
limiter specialised for `nv` would gain a further 3–11 % with gfortran and 8–28 % with nvfortran on the CPU, and
nothing on the GPU: left to the owner (it adds kernels to the API).

In the nvfortran CPU run the first version and the library differ by 9·10⁻¹⁰ in θ on one face of 2^20 (not investigated
further; the closed form is ill-conditioned near a double root of the pressure quadratic, which is the likely cause;
the benchmark's tolerance is 10⁻⁸); on the GPU and with gfortran they agree to 10⁻¹⁰.

## What changed in SETH

- `src/lib/seth_positivity.F90`: no automatic arrays, no array temporaries; `SETH_POSITIVITY_NV_MAX = 5` (public);
  `nv > 5` returns the low-order flux. Same algorithm and API.
- The device test had passed with the first version because it launches 4096 faces: it checks results, not speed, and
  the default heap is enough for that many.
- Round-off found on the way: nvfortran `-fast` evaluates the kernel's loop and the tests' array expression of the
  one-sided states differently (fused multiply-add). Near vacuum that difference exceeds the floor 1e-13, so the
  tests now allow the round-off of the state (`adm_ro`), and the kernel's header states what the floors guarantee.

## Rule for new kernels

A routine with `!$acc routine` declares no local array whose bound depends on a dummy argument and passes no array
expression of such a size as an argument. Where the size is a property of the scheme, it is a named constant of the
kernel; where it is the caller's (`nv`), the locals are over-dimensioned by a named constant and the loops run to the
dummy. `seth_bench_device_arrays` measures what the alternative costs.
