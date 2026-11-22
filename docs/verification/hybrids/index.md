# Hybrids

The pieces of a hybrid scheme, verified one by one before any hybrid is: the sensors that decide where the linear
branch may be used, and (later) the junction of a linear and a nonlinear flux.

| Case | Type | What it verifies |
|---|---|---|
| [H3 Hybrids in a solver](./harness) | P | the cost hybrid and the KEEP-PE hybrid in the 1-D Euler harness: smooth, junction, conservation, shocks |
| [H2 Central fluxes](./central) | P, C | the KEEP-PE fluxes of order 2 to 8: order, kinetic-energy identity, pressure equilibrium, bitwise equal to ADAM's former one; the Ducros sensors |
| [H1 Sensors](./sensors) | P | four published smoothness sensors against Li, Chen & Ren's (2022) static tests: the wavenumbers up to which each judges a sine smooth, and which faces of a jump between two sines each flags |
