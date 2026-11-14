# Audit: ADAM WENO machinery (common CPU + FNL), extraction readiness for SETH

ADAM commit `302f093e`, 2026-11-14. Files: `src/lib/common/adam_weno_object.F90` (1126 l.),
`src/lib/fnl/adam_fnl_weno_object.F90`, `adam_fnl_weno_kernels.F90`, `adam_fnl_weno_width_kernels.F90`,
`adam_fnl_weno_width_body.INC`. Callers checked: FLUME FNL face kernels, `src/tests/flume/unit/test_flume_weno_*`.

## 0. Scope correction

**ADAM contains no TENO.** `grep -i teno src/lib` is empty. What exists: upwind WENO 1/3/5/7/9 (S = 1..5),
centered WENO 2/4/6/8 tables, Jiang-Shu weights plus a scale-invariant descaler (Don et al. 2022), WENO
interpolation tables (Jiang, Shu & Zhang 2013), a weights-ratio smoothness sensor, ROR (recursive order reduction)
and IB order reduction bookkeeping.

## 1. Verdict on API agnosticism

The code has three layers, and they differ sharply:

| Layer | Signatures | Module dependencies | Agnostic? |
|---|---|---|---|
| **Kernels**: the non-TBP `weno_*` of the common module, `adam_fnl_weno_kernels`, `adam_fnl_weno_width_kernels` | `S`, plain arrays `weno_a(1:,0:,1:)`, `weno_p(1:,0:,0:,1:)`, `weno_d(0:,0:,0:,1:)`, scalars `weno_zeps`, `weno_rmu`, stencil `v(1:2,1-S:S-1)`, result `vr(1:2)` (+ `wr(1:2)`) | `penf` kinds, `S_MAX`/`S_MAX_M1`, `weno_weights_exponent`; **but obtained through `use adam_weno_object`**, which pulls in `adam_mpih_global` (`mpih`) and FiNeR | **Yes in the signatures, no in the module graph** |
| **Setup object** `weno_object` | `initialize(file_parameters, scheme, nb, ngc, ni, nj, nk)`, `load_from_file(file_ini)` | `mpih` singleton (messages, `error_stop`, `description`), FiNeR `file_ini` | **No** |
| **FNL helper** `weno_fnl_object` | `initialize(weno: type(weno_object))` | `mpih_fnl` singleton, FUNDAL device pointers | **No** |

So the hot path is already in the right shape. It was written for device kernels ("Non TBP procedures are also
provided for use without class object in device backend kernels", l.41), and the FNL interpolation test drives it
with bare arrays copied out of the object. The ADAM coupling sits in the setup and device-mirroring layers, plus one
import edge in the kernel modules.

## 2. ADAM-specific coupling to remove (setup layer)

1. **Singletons**: `use adam_mpih_global, only : mpih` (l.44): `print_message`, `error_stop`, rank prefixes in
   `description` (l.150-193). FNL: `mpih_fnl`.
2. **Configuration format**: `load_from_file(type(file_ini))` with the `[weno]` section and its option names
   (l.337-396). This belongs in the application, not in the library.
3. **Grid-shaped state inside the scheme object**: `initialize` takes `nb, ngc, ni, nj, nk` only to allocate
   `cell_scheme(nb, 1-ngc:ni+ngc, ..., 3)` and `ror_stats(...)` (l.298-300). These are per-cell solver fields that
   use ADAM's 5-D block layout; they have nothing to do with the scheme.
4. **Solver policy inside the scheme object**: ROR (`ror_number`, `ror_schemes`, `ror_threshold`, `ror_vars_number`,
   `ror_ivar` = indices of the *solver's* variables, `enable_ror_stats`) and IB reduction (`ib_reduction_extent`,
   `ib_reduced_order`). The order-reduction *mechanism* is generic: the tables hold every S up to the scheme's, and a
   kernel call takes S. The *policy* (which variable, which threshold, near which solid) is the solver's.
5. **Device mirroring tied to FUNDAL and to the ADAM object**: `weno_fnl_object` copies the tables plus the four
   ADAM fields above.

## 3. Defects found (verify before relying on them; none touches the FLUME production path)

- **D1, `weno_reconstruct_centered`: the left weights are computed on the wrong stencil** (l.1001-1007). The left
  (f=1) values are shifted into `vu` for the polynomials, but the weights are computed from the *unshifted*
  `v(1:2,1-S:-1+S)`, so the smoothness indicators of the left reconstruction use the right-biased cells, one cell
  off. Latent: **the routine has no caller in `src/`**, the `WENO_C_*` tables (`c`) are built but used by nothing,
  and there is no device twin. Do not port it as-is.
- **D2, `ror_schemes` filled with an uninitialised `S`** (l.368-372). `initialize` calls `load_from_file` *before*
  `select case` sets `S`, so `ror_schemes = self%S` (`ror_number = 0`) reads an undefined integer, and
  `self%S = ror_schemes(1)` (`ror_number > 0`) is overwritten right after. No check that `ror_schemes(r) <= S`
  (the tables are allocated `1:S`). The `scheme=` path never allocates `ror_schemes`/`ror_ivar`, and
  `weno_fnl_object%initialize` passes them unconditionally to `dev_assign_to_device`. Only NASTO (unmaintained) reads
  `ror_schemes`.
- **D3, withdrawn after checking the tables (2026-11-14): not a contradiction, an unexplained convention.** The
  kernel natively returns *both faces of the centre cell* from that cell's stencil: row 1 gives u+ at i-1/2
  (`a(1,:,3) = [0.3, 0.6, 0.1]`, stencil 0 = {i-2,i-1,i} with -1/6, 5/6, 1/3), row 2 gives u- at i+1/2. That is the
  FNL comment ("1 => i-1/2, 2 => i+1/2"). The header (l.6-39) describes the *caller's* use: row 1 is passed shifted
  one cell to the left, so both rows land on the same face (the two states at i+1/2), which is what FLUME does. Both
  are right; neither says how it relates to the other.
- **D4, a bitwise-identity claim that the toolchain does not guarantee** (l.951-953). The host spells out
  `d**wexp` as products, "bitwise the same" as the binary exponentiation. The device kernels still use
  `(…)**wexp`, and nvfortran's lowering of an integer power plus GPU FMA contraction make host = device bitwise a
  hope, not a property. Harmless numerically; it matters only if a test asserts bitwise CPU/GPU agreement.
- **D5, `zeps = 1e-6` hard-wired** (l.288). It cannot be set from the API or the file.
- **D6, fragile initialisation order.** `initialize_upwind_S1` assigns the *whole* `a`, `p`, `d` arrays to 1
  (l.451-457), and the S2..S5 routines overwrite their slices. Correct only because S1 is always called first.

## 4. Structural findings that matter for SETH

1. **Three copies of one algorithm**: the common non-TBP routines, the generic FNL kernels and the per-width
   `.INC` body, plus NVF and GMP twins. The common routines are `pure` but carry no `!$acc routine seq` (only
   `weno_descaling` and `weno_weights_exponent` do, l.1099/1120), so the FNL backend had to copy them. One set of
   `pure` kernels with `!$acc routine seq` / `!$omp declare target`, compiled once, can serve both the CPU and the
   device. The per-width variants (S a named constant, constant-bound locals, no calls) are a real performance win
   (no descriptors, no scratch memory on AMD offload); keep that idea as the main form.
2. **The table layout cannot express TENO.** `a(2,0:S-1,1:S)`, `p(2,0:S-1,0:S-1,1:S)`, `d(0:S-1,0:S-1,0:S-1,1:S)`
   assume S candidate stencils, all of width S. TENO uses candidate stencils of *different* widths (TENO6: three
   3-point stencils plus one 4-point; TENO8: 3,3,3,4,5) and a downwind stencil, and replaces the weight formula with a
   cutoff. SETH needs a stencil-set description (count, offset, width per candidate) instead of a single S.
3. **The weight formula is hard-coded in the kernel.** Only JS is implemented, plus scale invariance through
   `weno_rmu`. WENO-Z/M or the TENO cutoff would mean editing the kernel. Following the no-dispatch-in-the-loop rule,
   give each weight family its own kernel (per width), and let the caller pick the routine once, outside the loop.
4. **Good choices to keep**: the descaler `weno_rmu` is passed in by the caller (FLUME computes it from its own
   magnitude: library mechanism, caller policy); the weights ratio `wr` as a smoothness sensor; constant-bound locals
   sized by `S_MAX` (the comment measures 10% of a CPU run lost to heap automatic arrays without it); the
   interpolation tables sharing layout and `d` with the reconstruction ones.
5. **Precision is fixed to `R8P`**: no kind genericity. Decide in SETH whether R4P matters (on consumer GPUs the
   FP64:FP32 ratio is 1:64).
6. **One stencil per call (scalar API)**: `v(2,2S-1) -> vr(2)` per variable per face; the caller owns the loops and
   the data layout. That is what keeps it agnostic of ADAM's 5-D arrays; keep it, and add array wrappers for CPU
   callers only if needed.

## 5. Extraction plan (proposal)

| SETH piece | From ADAM | Change |
|---|---|---|
| `seth_weno_kernels` (pure, `acc routine seq`, `omp declare target`) | the common non-TBP routines + the FNL width body | one copy; per-width variants as the main form; generic-S variant kept; document the row/face convention (D3); `S_MAX` and `weights_exponent` local to SETH |
| `seth_weno_tables` | `initialize_upwind_S*`, `initialize_interpolation_S*` | pure builders returning plain arrays (no `mpih`, no file); no ROR/IB/grid fields; `zeps` a parameter (D5); no order dependence (D6) |
| (optional) `seth_weno_object` | `weno_object` minus §2 | a plain holder of the tables, no type-bound dispatch in the hot path; errors through a status argument, not `error_stop` via `mpih` |
| device mirroring | `weno_fnl_object` | left to the client (ADAM keeps it with FUNDAL), or a thin optional FUNDAL module: SETH's core needs only `!$acc routine` kernels taking device-resident arrays |
| centered WENO | `weno_reconstruct_centered` | fix D1 and add a test, or drop it (unused) |
| stays in ADAM | `load_from_file`, ROR/IB policy, `cell_scheme`, `ror_stats`, `description` with rank prefixes | — |
| new | — | stencil-set descriptor and weight-family kernels for TENO, WENO-Z/M |
