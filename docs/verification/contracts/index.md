# Contracts

What must hold exactly or bit for bit, and what protects it. These pages document their setup in tables; most have no
figure, since a plot of a difference that is zero shows nothing.

| Case | Type | What it verifies |
|---|---|---|
| [C1 Equivalence with ADAM](./adam-equivalence) | C | WENO-JS bitwise equal to ADAM's kernels, both table sets, Jiang–Shu and scale-invariant weights |
| [C2 Properties and invariances](./properties) | C | polynomial exactness, constants, translation, mirror, power-of-2 scaling, robustness to $10^{\pm 200}$, ENO across jumps |
| [C3 Coefficient tables](./tables) | C | derived in exact arithmetic, cross-checked with papers and ADAM, re-verified as compiled |
| [C4 Golden regression](./golden) | C | 8192 values of every kernel on exact binary inputs |
| [C5 Device](./device) | C | every kernel on the GPU against the host |
| [C6 Mutation testing](./mutation) | C | 25 deliberate defects, each caught |
| [C7 Positivity limiter](./positivity) | C | the limited flux admissible, its factor maximal, admissible fluxes untouched; on the GPU against the host |
