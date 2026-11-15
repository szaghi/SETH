# Changelog

All notable changes to this project are documented here.
Versions follow [Semantic Versioning](https://semver.org/).
Format follows [Keep a Changelog](https://keepachangelog.com/).

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



