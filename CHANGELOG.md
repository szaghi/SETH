# Changelog

All notable changes to this project are documented here.
Versions follow [Semantic Versioning](https://semver.org/).
Format follows [Keep a Changelog](https://keepachangelog.com/).

## [0.2.0] — 2026-11-17
### Added
- **positivity**: Positivity-preserving flux limiter for the Euler equations

- **license**: Offer SETH under GPL-3.0, BSD-2, BSD-3 or MIT

- **teno**: TENO5-opt, its weights derived from Fu et al.'s optimisation

- **teno**: TENO5-AZ, the adaptive cutoff that vanishes on smooth data

- **teno**: TENO8, eight points on Liang & Fu's sixth-order background

- **teno**: TENO5-A and TENO8-A, the cutoff from a discontinuity sensor

- **descaling**: The affine-invariant descaler of Wang & Don

- **mp**: The monotonicity-preserving bound of Suresh & Huynh

- **linear**: The linear backgrounds as kernels, the smooth branch of a hybrid

- **sensors**: Four published smoothness sensors, checked on Li, Chen & Ren's tests

- **kep**: KEEP-PE central fluxes and the Ducros sensors, the hybrid's other half

- **kinds**: The single-precision literal set, and a lint for kind-sensitive literals

- **weno-c**: The centred WENO of orders 4, 6 and 8, ADAM's with its defect D1 fixed

- **sensors**: Fu's TENO detector with its regulariser explicit, and Fu's two values exported


### Changed
- **kinds**: The kernels compute in PENF's R_P, the build chooses the kind


### Documentation
- **positivity**: The positivity limiter in the API and the verification

- **mutation**: 25 of 25 mutants killed with the positivity limiter mutants

- **home**: Landing page with the README sections; fobis dependency

- **style**: New identity, Ōgi mark in the Neon Moegi palette

- **performance**: Device-kernel rules, measured; the limiter's limits

- **mutation**: Record of the run after the limiter rewrite

- **conventions**: Figures and an animation of the stencil convention

- **schemes**: The background of reconstruction, interpolation, hybrids

- **performance**: Figures of the benchmarks, drawn from their outputs

- **verification**: One verification page, reached from the guide too

- **figures**: The Neon Moegi palette in every figure and animation script

- **figures**: Re-render the figures and animations in the Neon Moegi palette

- **verification**: WENO-Z9 runs at p = 2, pages and tables redrawn

- **schemes**: The WENO-Z exponent per order, measured

- **blast**: What the published runs used, and what the TENO cutoff does

- **verification**: Tables, figures and animations regenerated with the new schemes

- **performance**: The cost of every kernel added after v0.1.0, GPU and CPU

- **notes**: What the TENO papers' Roe flux with entropy fix is, from the sources

- **blast**: The flux explains most of the published positivity, not all

- **sensors**: The sensor study, critical points to Ducros, and what each sensor does to a hybrid

- **sensors**: Recommend a sensor per background (decision D2)

- **api**: The single-precision build, its weights and its limits

- **tests**: The test suite per kind, and its R4P build

- **kinds**: What the single-precision build verifies, and what it changes

- **kinds**: The precision page, what single precision changes, where it stops, what it buys

- **weno-c**: The centred WENO in the guide and the verification catalog, A5, and what the flows show

- **verification**: What SETH verifies, and what ADAM verifies

- **notes**: Why the centred WENO pollutes Burgers' smooth region after the shock, phi traced

- **weno-c**: The cause of the centred WENO's packet behind the Burgers shock

- **notes**: The centred WENO's Burgers packet accepted as documented

- **weno-c**: WENO8-C kept as ported, for smooth data or with the positivity limiter

- **perf**: WENO4-C's budget miss accepted, the budget kept as frozen

- **sources**: Balsara & Shu, Huang, Liang & Fu and Fu's review read and checked; TENO-A's original dated 2019

- **sources**: The originals of TENO-A, the TENO detector, KEEP-PE and the Ducros sensor read; Fu's multi-D eps

- **teno**: SETH's TENO5-A is TENO-A in Tang's sense, not Fu, Hu & Adams's original scheme


### Fixed
- **docs**: Run the dev server with the es2022 target

- **docs**: Maths with absolute values broke the tables that held it

- **tests**: The device case of the MP bound, and the Ai tolerance under FMA

- **kernels**: Single-precision weights that do not overflow

- **tests**: A fixed inflow for Shu-Osher's hybrids, whose central branch drifted the transmissive inflow


### Performance
- **positivity**: No automatic arrays in the limiter (333x on the GPU)

- **bench**: SETH's kernels timed in R4P against R8P, on the RTX 4070 and one CPU core

- **weno-c**: Mutants and benchmark of the centred WENO, two cost budgets missed

- **weno-c**: One division per state in the centred WENO's blend (R8P)

- **weno-c**: The one-division blend measured, WENO6-C within its budget, WENO4-C still over


## [0.1.0] — 2026-11-15
### Added
- **weno-js**: WENO-1/3/5/7/9 reconstruction and interpolation kernels

- **tables**: WENO-Z and TENO scheme descriptors derived in exact arithmetic

- **wenoz**: WENO-Z3/5/7 reconstruction kernels

- **teno**: TENO5, TENO6 and TENO6-opt reconstruction kernels

- **bench**: OpenACC benchmark of SETH's WENO kernels against ADAM's

- **tables**: Interpolation descriptors for WENO-Z and TENO

- **kernels**: WENO-Z and TENO interpolation kernels


### Documentation
- **notes**: Audit of ADAM's WENO machinery and the SETH kernel API design

- **notes**: Sources of the WENO-Z and TENO formulas for M2

- **m2**: Document the WENO-Z and TENO kernels, add their fpm tests

- **notes**: M3 phase A benchmark results, SETH faster than ADAM on GPU

- **guide**: Conventions, schemes, API, verification, performance, sources

- **verification**: Regression, device, convergence, mutation testing

- **harness**: Figures and animations with the ADAM/FLUME figure harness

- **verification**: A verification section structured like FLUME's

- **verification**: Final gate evidence of the night

- **style**: Adopt the ADAM sumi identity; link the open issues


### Fixed
- **teno**: TENO6-opt is Fu et al.'s scheme; add WENO-Z9 and sourced C_T

- **z-teno**: Indicators as sums of squares, no NaN on constant data


## [0.0.1] — 2026-11-14
### Added
- Scaffold SETH, the WENO/TENO reconstruction library for HPC



