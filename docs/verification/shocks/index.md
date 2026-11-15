# Shock capturing

The kernels inside the two ways a finite-difference code uses them, on the canonical problems by which shock-capturing
schemes are judged. The solvers are test harnesses (`src/tests/seth_vv_*.F90`), not library API: SETH ships kernels,
and these programs show what the kernels do in a complete scheme.

- **The splitting path** uses the reconstruction kernels: global Lax–Friedrichs flux splitting on characteristic
  variables (Roe-averaged eigenvectors at each face), $g^\pm$ reconstructed with the left and right biases, as Jiang &
  Shu (1996).
- **The interpolation path** uses the interpolation kernels: the alternative finite-difference formulation (Shu &
  Osher 1988; Jiang, Shu & Zhang 2013). Characteristic variables are interpolated to the faces, an HLLC Riemann solver
  gives the flux, and the flux is corrected to high order with its even derivatives.

Both run SSP-RK3 at CFL 0.5. A positivity-preserving flux limiter (Hu, Adams & Shu 2013) can be switched on. Each
stress test runs every scheme without it, recording whether density and pressure stay positive, and with it,
counting the faces it limits. The harness is described on the [Sod](./sod#models-and-numerics) page.

| Case | Type | What it verifies |
|---|---|---|
| [S1 Sod shock tube](./sod) | P | $L_1(\rho)$ against the exact solution and its convergence; WENO-Z5 and TENO5 below WENO5-JS; the limiter inactive |
| [S2 Lax shock tube](./lax) | P | the same on a strong contact behind a shock |
| [S3 123 problem](./test123) | P | near vacuum: positivity of every scheme with and without the limiter |
| [S4 Shu–Osher](./shu-osher) | P | entropy waves behind a Mach 3 shock against a 16× finer reference |
| [S5 Woodward–Colella blast](./blast) | P | colliding blast waves: which schemes need the limiter, and how often |
| [S6 Burgers equation](./burgers) | P | exact solution with a moving shock: design order before, the smooth region after |
| [S7 Composite profile](./advection) | P | long-time advection of smooth, discontinuous and kinked waves: dissipation, boundedness |

![Positivity: limited faces per scheme](/verification/positivity.png)

*How often the positivity limiter acts on the two stress tests (log scale); hatched bars: the scheme loses positivity
without it.*
