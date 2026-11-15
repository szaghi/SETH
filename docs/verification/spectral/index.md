# Spectral properties

How well each scheme resolves waves: the modified wavenumber of every linear background, computed from the compiled
tables, and the approximate dispersion relation of every nonlinear scheme, measured by applying the kernel to single
Fourier modes.

| Case | Type | What it verifies |
|---|---|---|
| [F1 Linear backgrounds](./linear-spectra) | P | no anti-dissipation, low-wavenumber error slope $o + 1$, Martín et al.'s resolution efficiencies |
| [F2 Approximate dispersion relation](./adr) | P | TENO equals its linear background for resolved waves; every scheme of order $\ge 5$ approaches its background at low wavenumbers; dissipation never positive except TENO6's measured band |
