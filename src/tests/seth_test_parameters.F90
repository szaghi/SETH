!< SETH tests: the literals of the tests whose meaning depends on the working kind.
module seth_test_parameters
!< SETH tests: the literals of the tests whose meaning depends on the working kind.
!<
!< The tests compute in R_P, the kind the library is built with (seth_kind_parameters). A tolerance, a regulariser or a
!< step that is right for one kind is wrong for another: 1e-14 is round-off in R8P and below the epsilon of R4P
!< (1.2e-7), 1e-40 is normal in R8P and underflows R4P. Each such literal of the tests lives here, one set per kind
!< behind the macros of seth_kind_parameters, never in a test program (roadmap issue #4, 6.2, phase 3); the library's
!< own are in seth_kind_parameters. scripts/check_literals.py fails on any one left outside the two modules.
!<
!< The R8P set is the values the suite was verified with: each name replaces the literal it was, unchanged.
!<
!< The R4P set is provisional until phase 4 records the R4P expectations (measured, with a margin). A round-off
!< tolerance keeps its size in ulps: the R8P value times epsilon(R4P)/epsilon(R8P) = 2^29 = 5.4e8, rounded up to one
!< digit. A quantity that is not a round-off (a regulariser, a step, a coordinate, a bound of the limiter) has its own
!< derivation next to it.
use seth, only : R_P
implicit none
private
public :: EPS_ZT, EPS_FIELD_SENSOR, COMPLEX_STEP, T_SLACK
public :: TOL_RIEMANN, TOL_KINK, TOL_BISECT, ORDER_FLOOR, ORDER_FLOOR_EULER
public :: TOL_TENO_LINEAR, TOL_FLUX, TOL_KE_PEP, TOL_DEVICE_HOST, TOL_DUCROS, TOL_DUCROS_SQRT
public :: TOL_DEVICE_STATE, TOL_DEVICE_WR, TOL_DEVICE_FLUX, POSITIVITY_FLOOR, TOL_GOLDEN, TOL_GOLDEN_WR
public :: GOLDEN_FILE, GOLDEN_SCALE_EXP0, GOLDEN_SCALE_DEXP
public :: BIG_JUMP, ROBUST_EXP_MIN, ROBUST_EXP_MAX, ROBUST_EXP_MAX_JS
public :: TOL_BOUNDS, TOL_THETA, TOL_PROPERTIES, TOL_X, TOL_ANTIDISSIPATION, TOL_TABLES, TOL_ADR, TOL_CONSERVATION

#if defined _R_P_IS_R4P
! R4P: epsilon 1.2e-7, smallest normal 1.2e-38. Provisional (phase 4 measures): see the module's notes.
character(*), parameter :: GOLDEN_FILE = 'src/tests/golden/seth_golden_r4_v1.txt' !< The R4P golden values.
integer,      parameter :: GOLDEN_SCALE_EXP0 = -60 !< The golden file's scaled stencils: data times 2^-60 ..
integer,      parameter :: GOLDEN_SCALE_DEXP = 10  !< 2^50 (12 of them): normal, and their indicators (squares, times
                                                   !< coefficients up to ~1e2) stay inside R4P's range (1e-38 .. 3e38).
! data magnitudes: the indicators square the data, so R4P's range ends near sqrt(3.4e38) = 1.8e19 over the indicator
! coefficients (up to ~1e3): data up to 1e16. Plain WENO-JS's power of the indicator no longer overflows first in R4P
! (its weights are ratios to the smallest indicator there), so it has the same bound.
real(R_P), parameter :: BIG_JUMP          = 1.e16_R_P !< The TENO test's large jump: at the edge of R4P's range.
real(R_P), parameter :: ROBUST_EXP_MIN    = -30._R_P  !< Robustness magnitudes from 1e-30 (normal; squares underflow) ..
real(R_P), parameter :: ROBUST_EXP_MAX    = 16._R_P   !< .. 1e16, every scheme ..
real(R_P), parameter :: ROBUST_EXP_MAX_JS = 16._R_P   !< .. and plain WENO-JS.
! regularisers and steps
real(R_P), parameter :: EPS_ZT           = 1.e-30_R_P !< As SETH_TENO_EPS in R4P: 1e-40 is not a normal R4P number.
real(R_P), parameter :: EPS_FIELD_SENSOR = 1.e-8_R_P  !< As in R8P: it only keeps 0/0 out of the sensor (it is added to
                                                      !< terms of order 1/h^2), and 1e-8/h^2 is normal.
real(R_P), parameter :: COMPLEX_STEP     = 1.e-30_R_P !< As in R8P: normal in R4P, and its square underflows to 0 as it
                                                      !< should (the step's second-order terms vanish).
real(R_P), parameter :: T_SLACK          = 6.e-6_R_P  !< Round-off: 1e-14 in ulps.
! reference solutions
real(R_P), parameter :: TOL_RIEMANN = 6.e-6_R_P  !< Round-off: 1e-14 in ulps.
real(R_P), parameter :: TOL_KINK    = 6.e-6_R_P  !< Round-off: 1e-14 in ulps.
real(R_P), parameter :: TOL_BISECT  = 2.5e-7_R_P !< A bracket on [0, 1] cannot be narrower than the spacing of R4P
                                                 !< there (6e-8 below 1): two ulps of 1.
! orders: errors below the floor are round-off, no order is measured from them
real(R_P), parameter :: ORDER_FLOOR       = 6.e-3_R_P !< Round-off: 1e-11 in ulps.
real(R_P), parameter :: ORDER_FLOOR_EULER = 6.e-2_R_P !< Round-off: 1e-10 in ulps.
! identities and comparisons, each a test's own
real(R_P), parameter :: TOL_TENO_LINEAR     = 6.e-5_R_P !< Round-off: 1e-13 in ulps.
real(R_P), parameter :: TOL_FLUX            = 6.e-6_R_P !< Round-off: 1e-14 in ulps.
real(R_P), parameter :: TOL_KE_PEP          = 6.e-4_R_P !< Round-off: 1e-12 in ulps.
real(R_P), parameter :: TOL_DEVICE_HOST     = 6.e-5_R_P !< Round-off: 1e-13 in ulps.
real(R_P), parameter :: TOL_DUCROS          = 6.e-7_R_P !< Round-off: 1e-15 in ulps.
real(R_P), parameter :: TOL_DUCROS_SQRT     = 6.e-6_R_P !< Round-off: 1e-14 in ulps.
real(R_P), parameter :: TOL_DEVICE_STATE    = 6.e-5_R_P !< Round-off: 1e-13 in ulps.
real(R_P), parameter :: TOL_DEVICE_WR       = 1.e-3_R_P !< Not carried in ulps (1e-9 would become 0.5 of a ratio in
                                                        !< [0, 1]): R8P's 1e-9 is 1e4 times its state tolerance; the
                                                        !< same factor over R4P's would exceed 0.1, so 1e-3 is taken.
real(R_P), parameter :: TOL_DEVICE_FLUX     = 1.e-3_R_P !< As TOL_DEVICE_WR: R8P's 1e-10 carried in ulps is 0.05.
real(R_P), parameter :: POSITIVITY_FLOOR    = 1.e-13_R_P !< As SETH_POSITIVITY_EPS, which keeps 1e-13 in R4P.
real(R_P), parameter :: TOL_GOLDEN          = 6.e-5_R_P !< Round-off: 1e-13 in ulps.
real(R_P), parameter :: TOL_GOLDEN_WR       = 1.e-3_R_P !< As TOL_DEVICE_WR.
real(R_P), parameter :: TOL_BOUNDS          = 6.e-4_R_P !< Round-off: 1e-12 in ulps.
real(R_P), parameter :: TOL_THETA           = 1.e-5_R_P !< The limiter backs theta off by SETH_POSITIVITY_BACKOFF (5e-7
                                                        !< in R4P), doubled per retry: a few retries below the maximum.
real(R_P), parameter :: TOL_PROPERTIES      = 6.e-3_R_P !< Round-off: 1e-11 in ulps.
real(R_P), parameter :: TOL_X               = 1.e-6_R_P !< Coordinates of order 1, computed in one or two operations: a
                                                        !< few ulps of 1 (1e-12 in ulps, 5e-4, would exceed a cell of 2/N
                                                        !< for N > 4000).
real(R_P), parameter :: TOL_ANTIDISSIPATION = 6.e-6_R_P !< Round-off: 1e-14 in ulps.
real(R_P), parameter :: TOL_TABLES          = 6.e-4_R_P !< Round-off: 1e-12 in ulps.
real(R_P), parameter :: TOL_ADR             = 6.e-4_R_P !< Round-off: 1e-12 in ulps.
real(R_P), parameter :: TOL_CONSERVATION    = 6.e-5_R_P !< Round-off: 1e-13 in ulps.
#else
! R8P: the values the suite was verified with.
character(*), parameter :: GOLDEN_FILE = 'src/tests/golden/seth_golden_v6.txt' !< The R8P golden values.
integer,      parameter :: GOLDEN_SCALE_EXP0 = -498 !< The golden file's scaled stencils: data times 2^-498 .. 2^96
integer,      parameter :: GOLDEN_SCALE_DEXP = 54    !< (12 of them), about 1e-150 .. 1e29.
! data magnitudes
real(R_P), parameter :: BIG_JUMP          = 1.e30_R_P !< The TENO test's large jump (tau/eps would overflow without the
                                                      !< normalisation by g_max).
real(R_P), parameter :: ROBUST_EXP_MIN    = -200._R_P !< Robustness magnitudes from 1e-200 ..
real(R_P), parameter :: ROBUST_EXP_MAX    = 140._R_P  !< .. 1e140, every scheme (the squares bound them at 1e154) ..
real(R_P), parameter :: ROBUST_EXP_MAX_JS = 30._R_P   !< .. 1e30, plain WENO-JS ((eps + IS)^wexp overflows near 1e38).
! regularisers and steps
real(R_P), parameter :: EPS_ZT           = 1.e-40_R_P !< eps of WENO-Z and TENO in the tests, the papers' value.
real(R_P), parameter :: EPS_FIELD_SENSOR = 1.e-8_R_P  !< Ducros' regulariser of the field sensors, before /h^2.
real(R_P), parameter :: COMPLEX_STEP     = 1.e-30_R_P !< Imaginary step of the complex-step derivative (exact to
                                                      !< round-off for any step, as long as it is normal).
real(R_P), parameter :: T_SLACK          = 1.e-14_R_P !< Relative slack of the end time of a run, so that rounding
                                                      !< of t does not add a step: a few ulps of 1.
! reference solutions
real(R_P), parameter :: TOL_RIEMANN = 1.e-14_R_P !< Relative convergence of the exact Riemann solver's Newton loop.
real(R_P), parameter :: TOL_KINK    = 1.e-14_R_P !< Distance from the kink of the Burgers solution at which the
                                                 !< exact solution takes its limit value.
real(R_P), parameter :: TOL_BISECT  = 1.e-15_R_P !< Bracket width ending the bisection of the Burgers characteristics.
! orders: errors below the floor are round-off, no order is measured from them
real(R_P), parameter :: ORDER_FLOOR       = 1.e-11_R_P !< Scalar runs (central fluxes, Burgers).
real(R_P), parameter :: ORDER_FLOOR_EULER = 1.e-10_R_P !< Euler runs (WENO and Roe paths).
! identities and comparisons, each a test's own
real(R_P), parameter :: TOL_TENO_LINEAR     = 1.e-13_R_P !< TENO equals its linear scheme on a smooth Gaussian.
real(R_P), parameter :: TOL_FLUX            = 1.e-14_R_P !< A central flux of constant states is the Euler flux.
real(R_P), parameter :: TOL_KE_PEP          = 1.e-12_R_P !< Kinetic-energy identity and pressure equilibrium.
real(R_P), parameter :: TOL_DEVICE_HOST     = 1.e-13_R_P !< seth_keep_pe6 on the device against the host.
real(R_P), parameter :: TOL_DUCROS          = 1.e-15_R_P !< Ducros' sensor at its limit values (1 or 0).
real(R_P), parameter :: TOL_DUCROS_SQRT     = 1.e-14_R_P !< psi = sqrt(phi) of the compression sensor.
real(R_P), parameter :: TOL_DEVICE_STATE    = 1.e-13_R_P !< Device against host, a state over its data's magnitude.
real(R_P), parameter :: TOL_DEVICE_WR       = 1.e-9_R_P  !< Device against host, a weights ratio (in [0, 1]).
real(R_P), parameter :: TOL_DEVICE_FLUX     = 1.e-10_R_P !< Device against host, the positivity limiter's theta and flux.
real(R_P), parameter :: POSITIVITY_FLOOR    = 1.e-13_R_P !< The floor of density and pressure the limiter keeps
                                                         !< (seth_positivity's SETH_POSITIVITY_EPS).
real(R_P), parameter :: TOL_GOLDEN          = 1.e-13_R_P !< A state against the golden file, over the stencil magnitude.
real(R_P), parameter :: TOL_GOLDEN_WR       = 1.e-9_R_P  !< A weights ratio against the golden file.
real(R_P), parameter :: TOL_BOUNDS          = 1.e-12_R_P !< Excursion of a state beyond the data's bounds at a step.
real(R_P), parameter :: TOL_THETA           = 1.e-9_R_P  !< The limiter's theta below the bisected maximal one.
real(R_P), parameter :: TOL_PROPERTIES      = 1.e-11_R_P !< Shift invariance and mirror symmetry of a reconstruction.
real(R_P), parameter :: TOL_X               = 1.e-12_R_P !< Coordinates of faces compared.
real(R_P), parameter :: TOL_ANTIDISSIPATION = 1.e-14_R_P !< Largest imaginary part of a dissipative modified wavenumber.
real(R_P), parameter :: TOL_TABLES          = 1.e-12_R_P !< Tables against the exact rationals of the generators.
real(R_P), parameter :: TOL_ADR             = 1.e-12_R_P !< Imaginary parts and growth in the approximate dispersion
                                                         !< relation.
real(R_P), parameter :: TOL_CONSERVATION    = 1.e-13_R_P !< Change of mass, momentum and energy over 50 steps.
#endif
endmodule seth_test_parameters
