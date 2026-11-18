| case | scheme | $L_1(\rho)$ | overshoot | cells across the contact | detector activations | faces limited |
|---|---|---|---|---|---|---|
| Sod | TENO5 | 0.00273 | 2.498e-16 | 5 | 0 | 0 |
| Sod | TENO5-ST | 0.002742 | 8.327e-17 | 5 | 0 | 0 |
| Sod | TENO5-RBF | 0.002742 | 8.327e-17 | 5 | 0 | 0 |
| Sod | TENO5-THINC | 0.001857 | 1.388e-16 | 3 | 38472 | 0 |
| Sod | TENO5-DR | 0.001386 | 2.220e-16 | 3 | 30778 | 0 |
| Sod | TENO5-ST+DR | 0.001351 | 8.327e-17 | 3 | 27192 | 0 |
| Sod | TENO5-RBF+DR | 0.001341 | 5.551e-17 | 3 | 27269 | 0 |
| Lax | TENO5 | 0.0092 | 0.00402 | 4 | 0 | 0 |
| Lax | TENO5-ST | 0.00919 | 0.003525 | 4 | 0 | 0 |
| Lax | TENO5-RBF | 0.009193 | 0.00359 | 4 | 0 | 0 |
| Lax | TENO5-THINC | 0.006667 | 0.008094 | 3 | 46501 | 0 |
| Lax | TENO5-DR | 0.005389 | 0.005532 | 2 | 37599 | 0 |
| Lax | TENO5-ST+DR | 0.006109 | 0.002556 | 3 | 37456 | 0 |
| Lax | TENO5-RBF+DR | 0.005982 | 0.001921 | 3 | 39697 | 0 |
| 123 | TENO5 | 0.004992 | 0.003146 | — | 0 | 76 |
| 123 | TENO5-ST | 0.004916 | 0.003946 | — | 0 | 72 |
| 123 | TENO5-RBF | 0.004914 | 0.003944 | — | 0 | 72 |
| 123 | TENO5-THINC | 0.003899 | 0.004668 | — | 12857 | 154 |
| 123 | TENO5-DR | 0.004838 | 0.005387 | — | 9620 | 120 |
| 123 | TENO5-ST+DR | 0.004791 | 0.005945 | — | 11684 | 120 |
| 123 | TENO5-RBF+DR | 0.004791 | 0.00595 | — | 14563 | 120 |
| Shu–Osher | TENO5 | 0.572 | — | — | 0 | 0 |
| Shu–Osher | TENO5-ST | 0.5655 | — | — | 0 | 0 |
| Shu–Osher | TENO5-RBF | 0.2934 | — | — | 0 | 0 |
| Shu–Osher | TENO5-THINC | 0.5644 | — | — | 74626 | 0 |
| Shu–Osher | TENO5-DR | 0.5508 | — | — | 60325 | 0 |
| Shu–Osher | TENO5-ST+DR | 0.5519 | — | — | 57540 | 0 |
| Shu–Osher | TENO5-RBF+DR | 0.2488 | — | — | 57251 | 0 |
| blast waves | TENO5 | 0.08622 | — | — | 0 | 514 |
| blast waves | TENO5-ST | 0.08531 | — | — | 0 | 538 |
| blast waves | TENO5-RBF | 0.08532 | — | — | 0 | 502 |
| blast waves | TENO5-THINC | 0.0902 | — | — | 246752 | 67 |
| blast waves | TENO5-DR | 0.08703 | — | — | 181616 | 108 |
| blast waves | TENO5-ST+DR | 0.08605 | — | — | 126937 | 93 |
| blast waves | TENO5-RBF+DR | 0.08632 | — | — | 152282 | 74 |

Sod, Lax, 123: $N = 200$ against the exact solution (the overshoot beyond its extrema; 123 with the positivity limiter); Shu–Osher: $N = 200$ on $[-5, 5]$ against WENO5-JS on 3200 cells; blast waves: $N = 400$, the limiter on, against WENO5-JS on 6400 cells. Activations: states whose detector fires, per characteristic field, summed over the run (DR's detector for DR, ST+DR, RBF+DR; THINC's for THINC). The RBF kernels with $s = 3\,\Delta x$ (Liang & Fu's $\alpha = 3$): 0.15 on Shu–Osher, outside the validated range $s \le 0.12$.
