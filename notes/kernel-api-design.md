# SETH kernel API: design proposal

Status: decisions 1-4 of section 8 taken by the owner (2026-11-14: all as recommended); M1 implemented. Input:
`notes/adam-weno-audit.md`.

## 1. Requirements

- **R1, no dynamic dispatch in the hot path**: the scheme is chosen once, outside the loop over faces.
- **R2, one source for CPU and GPU**: every kernel `pure`, `!$acc routine seq`, `!$omp declare target`; no copy per
  backend.
- **R3, agnostic signatures**: kinds, plain arrays and scalars only; no grid, no solver, no configuration format, no
  singletons, no device-memory library in the core.
- **R4, open to TENO and the newer schemes**: candidate stencils of different widths, downwind candidates, a family
  of weight functions (JS, Z, M, TENO cutoff), without touching existing kernels.
- **R5, extraction verifiable against ADAM**: SETH's WENO-JS reproduces ADAM's to round-off on the same stencils.

## 2. The central decision: coefficients are compile-time constants

ADAM passes the tables (`weno_a`, `weno_p`, `weno_d`) to every kernel call because `S` is a run-time value. That one
choice causes most of the coupling and cost found in the audit:

- the tables must be mirrored to the device (`weno_fnl_object`, FUNDAL, `dev_assign_to_device`);
- assumed-shape table dummies need descriptors in device code (why the per-width FNL variants were written);
- `S_MAX`-sized locals, a run-time `wexp`, and loops whose trip counts the compiler cannot see.

**Proposal: each scheme is a kernel whose coefficients are `parameter` arrays.** The compiler folds them, unrolls the
constant-trip loops and keeps everything in registers. There is no device data at all, so no device mirror object and
no FUNDAL dependency in SETH. The descriptor of a scheme (its candidate stencils, widths, offsets, optimal weights,
polynomial and smoothness coefficients) is a set of named constants in the scheme's module, and the kernel body is
written once, as an include, over those constants. This is ADAM's `adam_fnl_weno_width_body.INC` idea taken to its
end: there the width is constant but the tables are still run-time arrays.

Run-time order selection (ROR, the IB reduction) stays possible. The caller hoists `select case (S)` out of the loop
where it can; where it cannot (a per-cell order), a `select case` inside the kernel call is a uniform branch per
cell, with no divergence when neighbouring threads share the order. SETH provides the per-scheme kernels; the policy
stays in the solver.

## 3. Stencil convention (one face, two biases)

ADAM's kernel natively returns the two faces of the centre cell (row 1: u+ at i-1/2, row 2: u- at i+1/2); a
face-based caller passes row 1 shifted one cell to the left so both rows give the two states at one face (audit D3:
a convention every caller must reproduce, explained nowhere in one place). Proposal: **index cells
by their offset from the face.** For the face between cells 0 and 1 (`i+1/2`):

```
      -2     -1      0  |   1      2      3
...|------|------|------|------|------|------|...
                       face
```

- `vp(lo:hi)`: the values reconstructed with a **left** bias (e.g. f+ of a flux splitting, or the cell averages for
  the left state);
- `vm(lo:hi)`: the values reconstructed with a **right** bias (f−, or the averages for the right state), **same
  indexing, not shifted**; the kernel mirrors internally (`vm(1-m)` plays the role of `vp(m)`);
- `lo:hi` is the scheme's full stencil: `1-S:S` for WENO-(2S-1) with both biases in one array, `-2:3` for TENO6 and
  WENO-5 alike.

Result: `rp` (left state at the face, from `vp`) and `rm` (right state at the face, from `vm`). One stencil array
per bias, with no shifted copy. A caller that reconstructs a single field `u` for both states passes `vp = vm = u`;
a flux-split caller passes `f+` and `f−`.

Mapping to ADAM, for the extraction tests and an ADAM-side adapter: `v(2, m) = vp(m)` and `v(1, m) = vm(m+1)` for
`m = 1-S..S-1`; ADAM's `vr(2)`/`vr(1)` are `rp`/`rm`. Checked on the WENO-5 tables: row 2 is the left-biased
reconstruction at i+1/2 (`a(2,:,3) = [0.1, 0.6, 0.3]`, stencil 0 = {i-2,i-1,i} with 1/3, -7/6, 11/6); row 1, unshifted,
is u+ at i-1/2, so shifted by one cell it is u+ (the right state) at i+1/2.

## 4. Kernel signatures

One module per family; one kernel per scheme; every kernel has the same shape:

```fortran
pure subroutine seth_weno5_js(vp, vm, eps, rmu, rp, rm)
!$acc routine seq
real(R8P), intent(in)  :: vp(-2:3), vm(-2:3) !< left- / right-biased stencil values, cells by offset from the face
real(R8P), intent(in)  :: eps                !< smoothness regulariser
real(R8P), intent(in)  :: rmu                !< inverse descaler, shared by both biases (1: plain weights; ADAM's weno_rmu)
real(R8P), intent(out) :: rp, rm             !< left and right states at the face
```

- **Explicit-shape, constant bounds**: no descriptors, no `S` argument; the declared bounds are the scheme's stencil.
- **Optional outputs are separate kernels, not `optional` dummies** (an `optional` in device code is a run-time test
  plus an ABI cost): `seth_weno5_js_wr(..., rp, rm, wrp, wrm)` returns ADAM's weights ratio too, as its `_wratio_`
  routines do.
- **The regulariser is an argument** (audit D5), with ADAM's value as a public constant (`SETH_WENO_JS_EPS =
  1e-6_R8P`). It acts on the indicators as the tables scale them: WENO-7 and WENO-9 keep Balsara & Shu's integer forms,
  240 and 5040 times the normalised ones (`WENO7_JS_D_SCALE`, `WENO9_JS_D_SCALE`), so the same eps is effectively
  240 and 5040 times smaller there. Found by the table generator; ADAM behaves the same.
- **The descaler stays the caller's** (`rmu`, one per face as in ADAM and FLUME, which computes it from the magnitude
  of its own state); SETH exports ADAM's generic descaler as `seth_weno_descaling`.
- **Two biases per call** keep ADAM's arithmetic shape (the two reconstructions share loops), which is what R5's
  round-off comparison needs.

Interpolation (point values to face values, Jiang-Shu-Zhang 2013) uses other polynomial coefficients and the same
smoothness ones: separate kernels `seth_weno5_js_interp`, generated from the same body with the other constants.

## 5. The scheme descriptor (named constants per scheme)

```fortran
integer,   parameter :: NC = 4                         ! candidate stencils
integer,   parameter :: LO = -2, HI = 3                ! full stencil, cell offsets from the face (left bias)
integer,   parameter :: CLO(NC) = [-1, -2,  0,  0]     ! first cell of each candidate (TENO6 example, to verify)
integer,   parameter :: CW (NC) = [ 3,  3,  3,  4]     ! width of each candidate
real(R8P), parameter :: D  (NC) = [...]                ! optimal (linear) weights
real(R8P), parameter :: P  (MAXW, NC) = reshape(...)   ! face-value coefficients of each candidate
real(R8P), parameter :: B  (MAXW, MAXW, NC) = ...      ! smoothness-indicator quadratic forms
```

WENO-(2S-1) is the special case `NC = S`, `CW = S`, `CLO = 1-S..0`. TENO6 (Fu, Hu & Adams 2016, JCP 305:333-359):
three 3-point candidates and one 4-point downwind candidate. **The offsets above are from memory and must be checked
against the paper's tables before coding.** The weight function is the family's: JS `alpha = D/(eps + beta)^q`; Z
adds the global `tau`; TENO computes `gamma = (C + tau/(beta + eps))^q`, normalises, cuts at `C_T`, and renormalises
`D` over the kept candidates.

## 6. Layout of the code

```
src/lib/
  seth.F90                         re-exports kinds and every family module
  seth_kinds.F90                   R8P (R4P later: see 8)
  seth_weno_js.F90                 WENO-JS 1/3/5/7/9: constants + kernels (include of the bodies)
  seth_weno_js_body.INC            reconstruction body over the named constants (JS weights)
  seth_weno_z.F90, seth_teno.F90   next families, same pattern, own weight body
  seth_descaling.F90               weno_descaling (generic descaler helper)
src/tests/
  seth_test_weno_js_adam.F90       round-off agreement with ADAM's kernels (copied in tests/ as reference)
  seth_test_weno_js_order.F90      convergence order on smooth data, per scheme
  seth_test_weno_js_eno.F90        no spurious extrema across a step (bounded overshoot)
```

No setup object is needed for the core: the schemes have no state. A convenience `seth_weno_tables` (the old
coefficient arrays, built by pure functions, without `mpih`, files, ROR or IB fields) can be offered to clients that
keep ADAM's run-time-table kernels during the migration.

## 7. Milestones

| M | Content | Exit criterion |
|---|---|---|
| M1 | WENO-JS 1/3/5/7/9 reconstruction and interpolation, `_wr` variants, descaling helper | round-off agreement with ADAM on random stencils (all S, both table sets); order tests; gfortran, nvfortran (CPU and `-acc` build), ifx |
| M2 | WENO-Z, TENO5, TENO6 | order tests; TENO cutoff behaviour on a step; coefficients checked against the papers |
| M3 | ADAM migration: ADAM calls SETH's kernels, `weno_fnl_object` table mirroring removed for the schemes that no longer need it | ADAM's FLUME/PRISM regression goldens unchanged (to round-off) |
| M4 | precision genericity (R4P), centered WENO (audit D1 fixed), TENO8 | per-feature tests |

## 8. Decisions (owner, 2026-11-14: 1 compile-time, 2 face-offset arrays, 3 R8P now, 4 fix and port in M4)

1. **Compile-time coefficients (section 2)**: accept, or keep run-time tables passed to every call as ADAM does?
2. **Stencil convention (section 3)**: unshifted face-offset indexing with internal mirroring, or keep ADAM's shifted
   `v(1:2,·)` for drop-in compatibility?
3. **Precision**: R8P only until M4, or kind-generic from M1 (each body instantiated per kind, names suffixed)?
4. **Centered WENO**: fix D1 and port it in M4, or drop it (no caller in ADAM)?

## 9. M1 as implemented (2026-11-14)

- `src/lib/seth_weno_js_tables.F90`, generated by `scripts/gen_weno_js_tables.py` from ADAM's initialisers
  (expressions verbatim, every table checked in exact rational arithmetic); `seth_weno_js.F90`: 20 kernels
  (`seth_weno{1,3,5,7,9}_js[_interp][_wr]`), one body `seth_weno_js_body.INC`; `seth_descaling.F90`.
- `seth_test_weno_js_adam`: against ADAM's kernels copied verbatim (`src/tests/seth_adam_reference.F90`), 1.4 M
  results: bitwise equal with gfortran (debug and -O3) and nvfortran debug; nvfortran -fast 91 % bitwise, ifx 98.5 %,
  worst difference 2.0e-15 of the data magnitude (fused multiply-add contractions).
- `seth_test_weno_js_order`: design order 2S-1 on exp(x), every scheme, both states, reconstruction and interpolation.
- `seth_test_weno_js_eno`: worst overshoot across a unit jump 2.5e-13.
- The OpenACC build (`tests-nvf-oac`) compiles every kernel as `acc routine seq`; not run (GPU reserved for ADAM).
