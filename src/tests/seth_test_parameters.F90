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
use seth, only : R_P
implicit none
private
public :: EPS_ZT, EPS_FIELD_SENSOR, COMPLEX_STEP, T_SLACK
public :: TOL_RIEMANN, TOL_KINK, TOL_BISECT, ORDER_FLOOR, ORDER_FLOOR_EULER
public :: TOL_TENO_LINEAR, TOL_FLUX, TOL_KE_PEP, TOL_DEVICE_HOST, TOL_DUCROS, TOL_DUCROS_SQRT
public :: TOL_DEVICE_STATE, TOL_DEVICE_WR, TOL_DEVICE_FLUX, POSITIVITY_FLOOR, TOL_GOLDEN, TOL_GOLDEN_WR
public :: TOL_BOUNDS, TOL_THETA, TOL_PROPERTIES, TOL_X, TOL_ANTIDISSIPATION, TOL_TABLES, TOL_ADR, TOL_CONSERVATION

#if defined _R_P_IS_R4P
#error "SETH tests: the R4P set of seth_test_parameters is not defined yet (issue #4, 6.2, phase 3)"
#else
! R8P: the values the suite was verified with.
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
