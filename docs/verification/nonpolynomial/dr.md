# DR override

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| NP4 | Non-polynomial | P | — | — | `seth_test_nonpoly`, `seth_vv_nonpoly` |

::: warning Planned
This case is planned, not implemented: [issue #6](https://github.com/szaghi/SETH/issues/6), step 2.2, 2.3, 2.4, 2.5. The kernels
(`seth_dr_override`, `seth_teno5_dr`) do not exist yet; the page will carry the case, the program, the results and the findings when they do.
:::

## Case and scope

The two-steepness interpolation of Liang et al. 2022 (Eq. 31): the Gaussian pulse of their Table 4, contacts on Sod and Lax, the cost against three times TENO5 (E1), contact spreading (E2), steepening of smooth gradients (E3), the halo width (E4).
