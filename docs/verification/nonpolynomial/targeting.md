# Stencil targeting

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| NP5 | Non-polynomial | P | — | — | `seth_test_nonpoly` |

::: warning Planned
This case is planned, not implemented: [issue #6](https://github.com/szaghi/SETH/issues/6), step 3.1, 3.3. The kernels
(`seth_teno5_st`) do not exist yet; the page will carry the case, the program, the results and the findings when they do.
:::

## Case and scope

TENO5 whose masks select one of six stencils (Liang & Fu 2024, §3.1): order 5, order 4 beside a jump, the rows equal to the pair weights, and how often the mask (0,1,1) occurs (E7).
