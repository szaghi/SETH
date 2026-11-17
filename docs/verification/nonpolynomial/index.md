# Non-polynomial reconstructions

The schemes of the TENO line that leave the polynomial candidates: Takagi et al.'s THINC step function (TENO5-THINC,
2022), Liang, Xiao, Shyy & Fu's discontinuity-resolving interpolation (TENO5-DR, 2022) and Liang & Fu's stencil
targeting with Gaussian-RBF fluxes (2024). They share one structure, each piece verified on its own before the
composite: a flag per cell (TENO5's cutoff), a detector that reads the flags of the neighbouring cells, an override
of the face state where the detector fires. The plan, its gates and its decisions are
[issue #6](https://github.com/szaghi/SETH/issues/6); the kernel design is `notes/kernel-api-design.md`, section 11.

The composite kernels read twelve cells per face (`-5:6`): a solver using them needs **6 ghost cells**. The test
harness carries that footprint since 2026-11-17, with every earlier result unchanged.

| Case | Type | What it verifies |
|---|---|---|
| [NP1 Flags](./flags) | P | The 3-bit mask of TENO5's cutoff per cell, the input of every detector: equal to the cutoff inside `seth_teno5` on random stencils, and Liang, Xiao, Shyy & Fu's (2022) Tables 1–3 reproduced digit for digit. |
| [NP2 Detectors](./detectors) | P | The two detectors on the masks of the neighbouring cells: the levels of the three profiles of Tables 1–3, and the published static tests (a)–(e) of Liang et al. 2022, §3.1, with zero activations on the two sines. |
| [NP3 THINC override](./thinc) | P | Takagi et al.'s (2022) step function at the flagged faces: bitwise inert where not flagged, bounded, the cell value recovered, and the self-contained kernel equal to a two-pass reference. |
| [NP4 DR override](./dr) | P | The two-steepness interpolation of Liang et al. 2022 (Eq. 31): the Gaussian pulse of their Table 4, contacts on Sod and Lax, the cost against three times TENO5 (E1), contact spreading (E2), steepening of smooth gradients (E3), the halo width (E4). |
| [NP5 Stencil targeting](./targeting) | P | TENO5 whose masks select one of six stencils (Liang & Fu 2024, §3.1): order 5, order 4 beside a jump, the rows equal to the pair weights, and how often the mask (0,1,1) occurs (E7). |
| [NP6 Gaussian-RBF fluxes](./rbf) | P | The six RBF flux rows of Liang & Fu 2024 (Eq. 25): rows summing to zero, truncation errors of Eq. 29, order 5 as dx goes to 0, the band of anti-dissipation predicted at long waves (E6) and the validated range of s. |
| [NP7 Ablation](./ablation) | P | Seven schemes on every 1-D case (TENO5, ST, RBF, THINC, DR, ST+DR, RBF+DR): L1 error, overshoot, contact width, activations, limited faces; which piece buys what. |
| [NP8 Length units](./units) | P | Shu–Osher on [0, 10] and rescaled to [0, 1] (E5): TENO5 equal to round-off, the RBF scheme's dependence on the length unit measured, and the s of every published case. |

NP1 is implemented; the other seven are planned.
