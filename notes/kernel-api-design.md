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

## 10. M2 as implemented (night 2026-11-14)

- `scripts/gen_z_teno_tables.py` derives the WENO-Z3/5/7, TENO5, TENO6 and TENO6-opt descriptors (candidate cells and
  face coefficients, Jiang-Shu indicators incl. the 6-point global one, linear weights, tau coefficients) in exact
  arithmetic and cross-checks them (TENO Table 1, ADAM's indicators, exactness, tau orders): `seth_z_teno_tables`.
- The section-5 descriptor became the actual layout: `NC`, `LO`, `UMAX`, `CLO`, `CW`, `CP`, `CU` (indicators as
  upper-triangular forms), `D`, `TC` (tau as a linear combination of the beta_k), `TG`/`CUG` (the global indicator
  of TENO6). One body, `seth_z_teno_body.INC`, serves both families; the right state is the left-biased
  reconstruction of the mirrored data.
- Signatures: WENO-Z `(vp, vm, eps, p, rp, rm)`; TENO `(vp, vm, eps, ct, rp, rm)`. No `rmu`: both weight families
  depend on tau/beta, already scale-invariant up to eps.
- Measured: orders Z3 2.9-3.0, Z5 5.0, Z7 7.0, TENO5 5.0, TENO6 6.0, TENO6-opt 4.0; TENO returns its linear scheme on
  smooth data to round-off; no overshoot across a jump; a 1e30 jump stays finite (normalised chi).
- Not done, left to the owner: WENO-Z9 (tau_9 unsourced), C_T defaults for TENO5/TENO6, TENO8 (M4).

## 11. Composite non-polynomial kernels (issue #6, 2026-11-17)

Input: issue #6 (the plan, decisions D1–D8 taken by the owner on 2026-11-16) and `notes/m5-sources.md` (the primary
reads of step 0.1). This section **supersedes the row "TENO-THINC: not self-contained per face" of issue #4, section
8** (owner decision D1): the composite schemes are implemented, at the cost measured in step 2.4.

### 11.1 Three pieces

TENO5-THINC (Takagi et al. 2022), TENO5-DR (Liang, Xiao, Shyy & Fu 2022) and the complete scheme of Liang & Fu 2024
(stencil targeting, Gaussian-RBF fluxes, DR on top) share one structure: a **flag** per cell, a **detector** that
reads the flags of the neighbouring cells, an **override** that replaces the face state where the detector fires.

| Piece | Kernel | Input | Output |
|---|---|---|---|
| Flag | `seth_teno5_mask(v, eps, ct, m)` | `v(-2:2)`, the five values centred on one cell, upwind to the left | `m = delta_0 + 2 delta_1 + 4 delta_2`, TENO5's cutoff outcome (D2: the mask, not `zeta`) |
| Detector | `seth_dr_detect(mp, mm, lp, lm)`, `seth_thinc_detect(mp, mm, lp, lm)` | the masks of the cells around the face, by offset from it, both biases, `mp, mm(-3:3)` | integer levels `lp`, `lm` (DR: 0, 1 by Liang 2022 Eq. 24, 2 by Eq. 32; THINC: 0, 1 by Takagi 2022 Eq. 20) |
| Override | `seth_dr_override(vp, vm, lp, lm, t1, t2, rmu, rp, rm)`, `seth_thinc_override(vp, vm, lp, lm, beta, rmu, rp, rm)` | `vp, vm(-1:2)`, the levels, the states | the states, `inout`, replaced where the level is positive |

The bit `k` of a mask is candidate `k` of `seth_teno5` in SETH's order, which for TENO5 is the papers': `S0 = {-1,0,1}`
central, `S1 = {0,1,2}` downwind, `S2 = {-2,-1,0}` upwind (`TENO5_CLO = [-1, 0, -2]`). The flag of a cell is the
cutoff of the left-biased TENO5 reconstruction at its right face, computed from the five values centred on it; the
detector reads `zeta` (Liang 2022, Eq. 17), a function of the mask, so the mask carries one bit more (`delta_0`) and
Tables 1–3 of Liang 2022 test it digit for digit.

`seth_teno5_mask` is **the cutoff of `seth_teno5` and nothing else**: the same body (`seth_z_teno_body.INC`) compiled
with a mask switch that skips the candidate values and the weights, so the flags cannot drift from the kernel's. The
criterion is calibrated on plain TENO5 (`tau_5 = |beta_1 - beta_2|`, `C_T = 1e-5` in the papers): no detector is
offered on the flags of TENO5-AZ, TENO5-A, TENO6 or TENO8 (issue #6, section 8).

### 11.2 Cells by offset from the face, mirrored for the right state

As every SETH kernel (section 3), the composite kernels index cells by their offset from the face `i+1/2` (cell 0 on
its left) and take two arrays, `vp` for the left-biased and `vm` for the right-biased values, unshifted. The left
state reads the masks of cells `-3..2` (DR; THINC `-3..3`), each from `vp(c-2:c+2)`, so `vp(-5:4)` (THINC `vp(-5:5)`).
The right state is the left-biased construction on the mirrored data `w(m) = vm(1-m)`, so it reads `vm(-3:6)`
(THINC `vm(-4:6)`: the mirror of `-5:4` is `-3:6`, of `-5:5` is `-4:6`). With both biases the footprint is **`-5:6`, twelve cells**. The
detectors take the masks in the same offsets (`mp(-3:3)` of the left state, `mm(-3:3)` of the mirrored right state),
so one integer routine serves both states.

### 11.3 The default: one self-contained kernel per face (D3, D6)

```fortran
pure subroutine seth_teno5_dr(vp, vm, eps, ct, t1, t2, rmu, rp, rm)
!$acc routine seq
real(R_P), intent(in)  :: vp(-5:6), vm(-5:6) !< left- / right-biased values, cells by offset from the face
real(R_P), intent(in)  :: eps, ct            !< TENO5's regulariser and cutoff (the papers: C_T = 1e-5)
real(R_P), intent(in)  :: t1, t2             !< tanh of the two steepness constants (SETH_DR_T1, SETH_DR_T2)
real(R_P), intent(in)  :: rmu                !< inverse descaler of the first-order guard (1: the paper's form)
real(R_P), intent(out) :: rp, rm             !< left and right states at the face
```

`seth_teno5_thinc(vp, vm, eps, ct, beta, rmu, rp, rm)` has the same shape. Each call computes TENO5's two states, the
masks of the neighbouring cells by `seth_teno5_mask`, the detector and the override: no state is kept between faces,
the shape of every other SETH kernel (R1–R3 hold unchanged). The price is the recomputation of the neighbours' masks:
up to six (THINC seven) mask evaluations per state, each about 0.7 of a TENO5 state (budget of step 1.1), against
the owner's budget of 3 times `seth_teno5` per face (step 2.4, E1).

**Early exit, results unchanged.** The face's own cell decides first. If its `zeta` is 2 or −1 neither DR condition
can hold (Eqs. 24 and 32 need it 0 or 1, THINC's Eq. 20 needs 0: `notes/m5-sources.md`) and no neighbour is
evaluated; on smooth data the masks of the downwind cells `1, 2` (DR) or `1, 2, 3` (THINC) rule detection out when
none has `zeta` in `{-1, 2}`. The order of evaluation changes the cost only: the result is
the same as with every mask computed, which the two-pass reference of the tests (N1) checks bitwise.

**Basis.** For a system the caller projects all twelve cells into the characteristic basis of the face being
reconstructed (as it projects the five or six cells of TENO5 today), so the masks of the neighbours are computed in
that face's basis. Tables 1–3 and every scalar statement of the papers hold per characteristic field. The papers are
silent on the basis of a neighbour's flags (issue #6, section 1.2, tier U); this choice is the one that keeps the
kernel self-contained.

**Ghost cells.** A solver needs **6 ghost cells** for these kernels (3 for TENO5, 5 for WENO-9). The test harness
(`seth_vv_euler`, `NG = 6`) and the registry stencils (`seth_test_schemes`, `VLO:VHI = -5:6`) carry the footprint
since step 0.3; gate N6 (experiment E4) proves the width: two blocks with 6 ghost cells equal one block bitwise, and a
harness mutant with 5 differs.

### 11.4 The second variant: stored masks, two passes

If the default exceeds 3 times `seth_teno5` per face (E1), a second variant is implemented: a first pass computes
each cell's mask once (in the basis of that cell's own face) and stores it, a second pass reads the stored masks and
reconstructs. It costs one mask per cell instead of six to seven per state, but the caller keeps an integer array
and its halo, and on a system the masks are computed in another basis than the face's (the papers are silent on
which is meant, tier U). The two variants coincide on scalar problems (N1). The pieces of section 11.1 are exported
for this variant and for callers building their own composite.

**Shipped as a kernel pair (issue #6, step 2.6, owner decision D12, 2026-11-18).** E1 measured THINC over the budget
everywhere and DR and RBF+DR over it with nvfortran; the owner kept the self-contained kernels as the default and
asked for this variant as library kernels too:

```fortran
pure subroutine seth_teno5_face_masks(vp, vm, eps, ct, ml, mr)   ! pass 1, vp, vm(-2:3)
pure subroutine seth_teno5_dr_stored(vp, vm, mp, mm, eps, ct, t1, t2, rmu, rp, rm)                  ! pass 2
pure subroutine seth_teno5_thinc_stored(vp, vm, mp, mm, eps, ct, beta, rmu, rp, rm)
pure subroutine seth_teno5_rbf_dr_stored(vp, vm, mp, mm, eps, ct, cr, s1, s2, t1, t2, rmu, rp, rm)
```

Pass 1 runs over the faces, not the cells: at face `i+1/2`, on TENO5's six cells already projected into that face's
basis, it returns the left-state mask of cell `i` and the mirrored right-state mask of cell `i+1`. So "each mask in
the basis of its own face" has one meaning: the left mask of a cell in the basis of the face to its right, the right
mask in the basis of the face to its left. The caller stores `msl(i) = ml`, `msr(i+1) = mr` and passes
`mp(c) = msl(i+c)`, `mm(c) = msr(i+1-c)`, `c = -3..3`, to pass 2, which reads TENO5's stencil `-2:3` only. Pass 1
over the faces `-3..N+3` of an `N`-cell block covers the masks the faces `0..N` read. No state is kept in the library
(R3): the arrays and their halo are the caller's.

### 11.5 Parameters and kinds

- The steepness constants enter as `t1 = tanh(0.91)`, `t2 = tanh(1.23)` (`SETH_DR_T1`, `SETH_DR_T2`), literals
  printed by the generator for each kind and compared with the intrinsic to one ulp (`seth_test_tables`).
- THINC's `beta = 1.6` (`SETH_THINC_BETA`, Takagi 2022, Eq. 15); the face formula (Eq. 16) needs one `tanh` and no
  `exp`: `A = tanh(alpha beta/2)/tanh(beta/2)` (derived in `notes/m5-sources.md`; the generator of step 0.4 proves
  it). The override includes Takagi's first-order guard, Eq. 17 (the same as Liang 2022, Eq. 33): `alpha` is 0/0
  where `f_{i+1} = f_{i-1}`.
- The first-order guard `(f_{i+1}-f_i)(f_i-f_{i-1}) <= 1e-15` of both overrides (Takagi Eq. 17, Liang Eq. 33) is
  dimensional; the kernel multiplies each difference by
  `rmu` (as `seth_mp_bound`), `rmu = 1` reproducing the paper. In R4P the threshold is the same number (it is above
  `tiny` in both kinds).
- The RBF kernel takes `s1 = alpha_1 dx`, `s2 = alpha_2 dx` from the caller (D5): no grid enters the library. Its
  results depend on the length unit; the guide states the validated range of `s` (D7).

### 11.6 What stays with the caller (R3)

Which variable carries the flags, where stored masks live (variant 11.4), the halo exchange of six cells, the
projection of twelve cells per face, `dx` for the RBF kernel. SETH ships kernels; a solver's choices are not made in
the library.
