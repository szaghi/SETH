| scheme | $L_1(\rho)$ Sod | $L_1(\rho)$ Lax | $L_1(\rho)$ Shu–Osher | $L_1(\rho)$ blast waves | blast: faces limited |
|---|---|---|---|---|---|
| TENO5 | 0.00273 | 0.0092 | 0.572 | 0.08622 | 514 |
| TENO5-ST | 0.002742 | 0.00919 | 0.5655 | 0.08531 | 538 |
| TENO5-RBF | 0.002742 | 0.009193 | 0.2934 | 0.08532 | 502 |
| TENO5-THINC | 0.001857 | 0.006667 | 0.5644 | 0.0902 | 67 |
| TENO5-DR | 0.001386 | 0.005389 | 0.5508 | 0.08703 | 108 |
| TENO5-ST+DR | 0.001351 | 0.006109 | 0.5519 | 0.08605 | 93 |
| TENO5-RBF+DR | 0.001341 | 0.005982 | 0.2488 | 0.08632 | 74 |

Sod, Lax: $N = 200$ against the exact solution; Shu–Osher: $N = 200$ on $[-5, 5]$ against WENO5-JS on 3200 cells; blast waves: $N = 400$, the positivity limiter on, against WENO5-JS on 6400 cells. The RBF kernels with $s = 3\,\Delta x$ (Liang & Fu's $\alpha = 3$): 0.15 on Shu–Osher, outside the validated range $s \le 0.12$.
