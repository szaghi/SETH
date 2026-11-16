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
public :: ORDERS_AT_FLOOR_OK, CONV_FLOOR_OP, CONV_ORDER_LOSS, INTERP_FLOOR, INTERP_CP_ORDER_LOSS, CP2_AZ_TANG
public :: HENRICK_TOL, POINT_FLOOR, POINT_ORDER_LOSS, MP_NG, SPECTRAL_FLOOR, SPECTRAL_SLOPE_LOSS
public :: TAU_FLOOR, TAU_ORDER_LOSS, CENTRAL_FLOOR_N, BURGERS_FLOOR_N, EULER_FLOOR_N2
public :: POS_RANGE_EXP, SMALL_STEP, AI_STEP_TOL, SENSOR_SCALES, SENSOR_FLOW_RATIO, SENSOR_FLOW_EXEMPT
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
! orders (phase 4: measured on gfortran 16.2, tests-gnu-r4). R4P's round-off plateau sits ~1e-5 for a derivative on
! the grids of the tests, so a scheme of order 6 or more reaches it between the first two grids: no pair of grids is
! above round-off, and the order cannot be observed. Such a scheme converged at least as fast as R4P allows (a
! degraded one would stay above the floor), which the tests accept in R4P. The coarse pairs that remain are
! pre-asymptotic: the orders lost there are recorded per test, the measured loss plus 5 %.
logical,   parameter :: ORDERS_AT_FLOOR_OK = .true. !< No pair above round-off: accepted (R4P), not a failure.
real(R_P), parameter :: CONV_FLOOR_OP    = 20._R_P  !< seth_test_convergence: the operator's floor, in eps/h: 5 times the
                                                    !< measured R4P plateau (4 eps/h at N = 160 .. 640).
real(R_P), parameter :: CONV_ORDER_LOSS  = 0.55_R_P !< seth_test_convergence: order lost on the pairs above the floor
                                                    !< (WENO-Z5 advection 4.0 at N = 20 -> 40, against 4.5): 0.5 + 5 %.
real(R_P), parameter :: INTERP_FLOOR     = 20._R_P  !< seth_test_interpolation: the floor of face errors, in eps: 5 times
                                                    !< the measured R4P plateau (5e-7, 4 eps).
real(R_P), parameter :: INTERP_CP_ORDER_LOSS = 2.65_R_P !< seth_test_interpolation, first-order critical points: the
                                                    !< pairs R4P can observe (h = 0.4 .. 0.05) have R8P's errors, and
                                                    !< they are pre-asymptotic (R8P reaches the design order on pairs
                                                    !< below R4P's round-off). Worst: WENO-Z7 (p = 2), minus state, 4.0
                                                    !< at h = 0.4 -> 0.2 against 6.5: 2.5 + 5 %. In R4P the check guards
                                                    !< against a collapse at critical points, not the design order.
logical,   parameter :: CP2_AZ_TANG = .false. !< seth_test_critical_points: TENO5-AZ at the second-order point is not
                                              !< held to Tang's 5. Its error falls from 4e-4 (h = 0.05) to 2e-10 (h =
                                              !< 0.025) as the cutoff switches, below R4P's round-off: the one pair R4P
                                              !< observes is the coarse one, 2.0 in R8P too. The generic floor 1.7 of
                                              !< the second-order point still applies.
real(R_P), parameter :: HENRICK_TOL = 1.7_R_P !< seth_test_critical_points: WENO-JS against Henrick's rates. The
                                              !< floor (10 eps max|f|/h) leaves R4P the coarse pairs, which are
                                              !< pre-asymptotic (R8P has the same errors there): WENO7-JS at the
                                              !< second-order point 5.61 at h = 0.1 -> 0.05, against 4: 1.61 + 5 %.
real(R_P), parameter :: POINT_FLOOR = 20._R_P  !< The pointwise order tests (seth_test_weno_js_order, _wenoz, _teno,
                                              !< _linear): the floor of a face error, in eps (|g| < 2): 5 times R4P's
                                              !< plateau of a face value (4 eps, seth_test_interpolation).
real(R_P), parameter :: POINT_ORDER_LOSS = 0.8_R_P !< The pointwise order tests: order lost on the pairs above the
                                                !< floor, which are pre-asymptotic (R8P has the same errors there).
                                                !< Worst: WENO5-JS interpolation, minus state, 3.74 at the coarsest
                                                !< pair, against 4.5: 0.76 + 5 %.
integer,   parameter :: MP_NG = 3 !< seth_test_mp: the sin^4 grids N = 20, 40, 80. The cell averages are differences of
                                 !< the primitive (~0.4): they lose eps 0.4/h, and from N = 160 the L1 error grows as 1/h
                                 !< (1.38e-6 at N = 80, 1.66e-6 at 160, 1.3e-5 at 1280), noise the bound then reacts to.
real(R_P), parameter :: SPECTRAL_FLOOR = 20._R_P !< seth_test_spectral: the round-off floor of |xi~ - xi| at xi, in
                                               !< eps xi (R4P's plateau ~1e-8 = 0.1 eps at xi = 0.1): below it the
                                               !< slope moves to a pair twice as large.
real(R_P), parameter :: SPECTRAL_SLOPE_LOSS = 0.22_R_P !< seth_test_spectral: TENO8's slope on (0.4, 0.8) is 6.50
                                               !< against 6.7, the same in R8P on that pair (pre-asymptotic there):
                                               !< 0.20 + 5 %.
real(R_P), parameter :: TAU_FLOOR = 1.e-6_R_P !< seth_test_tables: tau cancels to h^order out of indicators of order
                                            !< h^2, so in R4P its values at h <= 0.05 are round-off (1e-8 .. 4e-8
                                            !< for WENO-Z9, TENO8 on gfortran, larger on nvfortran and ifx, and
                                            !< WENO-Z9's is inexact already at h = 0.1, 6.2e-8 against 9.1e-8):
                                            !< the order is taken on the finest pair above 1e-6.
real(R_P), parameter :: TAU_ORDER_LOSS = 0.39_R_P !< seth_test_tables: above the floor, the fifth-order taus keep only
                                            !< the pair (0.2, 0.1): 4.33 against 4.7, 4.32 in R8P on that pair
                                            !< (pre-asymptotic): 0.37 + 5 %.
real(R_P), parameter :: CENTRAL_FLOOR_N = 4.8e-6_R_P !< seth_test_central: the derivative's round-off grows as 1/h =
                                                  !< n: plateau 7.5 eps n (KEEP-PE6, n = 256: 2.3e-4), floor 40 eps n.
real(R_P), parameter :: BURGERS_FLOOR_N = 1.e-7_R_P !< seth_vv_burgers: the smooth L1 error's round-off grows with N
                                                  !< (accumulated over N steps): plateau 1.9e-8 N (3e-6 at N = 160),
                                                  !< floor 5 times that.
real(R_P), parameter :: EULER_FLOOR_N2 = 3.3e-8_R_P !< seth_vv_euler_order, seth_vv_roe: the density wave's round-off,
                                                  !< accumulated over the steps of the time-accurate dt (their number
                                                  !< grows faster than N): 4.2e-5 at N = 80, 1.35e-4 at 160, 4.7e-4 at
                                                  !< 320 (TENO5), about 6.6e-9 N^2; floor 5 times that.
real(R_P), parameter :: POS_RANGE_EXP = 5._R_P !< seth_test_positivity: densities and pressures from 1e-5 to 1. With
                                             !< the R8P range (1e-8) a one-sided state q - lam f next to an O(1)
                                             !< cell rounds to a negative pressure (its round-off, eps 50, exceeds
                                             !< p): the Lax-Friedrichs premise fails at round-off, and it still does
                                             !< from 1e-6 (measured, 4 pairs); 1e-5 is clean.
real(R_P), parameter :: SMALL_STEP = 1.220703125e-4_R_P !< seth_test_properties: the small step on the offset 1, 2^-13
                                                     !< (1024 ulps of 1; R8P's 1e-7 is one ulp of R4P), its square
                                                     !< still far below the WENO-JS eps (the mean-|v| check needs it).
real(R_P), parameter :: AI_STEP_TOL = 0.26_R_P !< seth_test_properties: the affine-invariant WENO-5 alone overshoots the
                                             !< small step by 0.248 of it (measured; 0.27 without the descaler): the
                                             !< expanded indicators cancel the offset to eps offset^2, which rmu^2 ~
                                             !< 1/step^2 makes O(1) in R4P, for any step whose square is below eps.
                                             !< Recorded + 5 %, so that a change shows. Subtracting the stencil's
                                             !< mean first stays ENO in R4P (the check at 1e-5 holds in both kinds).
real(R_P), parameter :: SENSOR_SCALES(5) = [1.e-3_R_P, 1.e-2_R_P, 1._R_P, 1.e3_R_P, 1.e6_R_P] !< seth_test_sensors:
                                             !< the scales a of a f + b. On the offset b = 100 (ulp 7.6e-6), R8P's
                                             !< a = 1e-6 leaves f below R4P's resolution, and 1e-5, 1e-4 (13 ulps)
                                             !< still change the range-normalised decisions (measured); 1e-3 does not.
real(R_P), parameter :: SENSOR_FLOW_RATIO = 3._R_P !< seth_vv_sensors: a cost hybrid's L1 at most 3 times the pure
                                             !< scheme's, as in R8P, but for the case below.
character(*), parameter :: SENSOR_FLOW_EXEMPT = 'shu-osher teno6 lcr' !< seth_vv_sensors: not held to the bound in R4P
                                             !< (printed): its L1 is 3.1 (gfortran), 4.3 (nvfortran), 9.3 (ifx)
                                             !< times the pure scheme's, 2.1 in R8P, the flagged fraction alike: an
                                             !< R4P instability of that hybrid, compiler-dependent, not analysed
                                             !< (left to the owner). Every other case keeps the bound 3.
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
real(R_P), parameter :: ORDER_FLOOR       = 0._R_P    !< Not carried in ulps (it is an absolute error, and R4P's grows
                                                      !< with N): the per-N floors below (CENTRAL_FLOOR_N, ...) hold it.
real(R_P), parameter :: ORDER_FLOOR_EULER = 0._R_P    !< Not carried in ulps (absolute, and growing with N in R4P):
                                                      !< EULER_FLOOR_N2 holds it.
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
real(R_P), parameter :: TOL_GOLDEN_WR       = 6.5e-3_R_P !< Measured: the golden file is gfortran's; nvfortran's weights
                                                         !< ratios differ by up to 6.18e-3 (fused multiply-adds), + 5 %.
real(R_P), parameter :: TOL_BOUNDS          = 6.e-4_R_P !< Round-off: 1e-12 in ulps.
real(R_P), parameter :: TOL_THETA           = 1.35e-3_R_P !< Measured: theta below the bisected maximum by up to 1.28e-3
                                                         !< (nvfortran, nv = 3; gfortran 3.99e-4; data from 1e-5:
                                                         !< POS_RANGE_EXP): theta solves a
                                                         !< quadratic pressure condition whose
                                                         !< root R4P's rounding moves where it is nearly tangent (an
                                                         !< optimality gap, not a positivity one): + 5 %.
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
! orders: the test's own floors, no order loss, no acceptance at round-off
logical,   parameter :: ORDERS_AT_FLOOR_OK = .false. !< No pair above round-off is a failure.
real(R_P), parameter :: CONV_FLOOR_OP    = 1.e3_R_P  !< seth_test_convergence: the operator's floor, in eps/h.
real(R_P), parameter :: CONV_ORDER_LOSS  = 0._R_P    !< seth_test_convergence: no order lost.
real(R_P), parameter :: INTERP_FLOOR     = 1.e2_R_P  !< seth_test_interpolation: the floor of face errors, in eps.
real(R_P), parameter :: INTERP_CP_ORDER_LOSS = 0._R_P !< seth_test_interpolation: no order lost at critical points.
logical,   parameter :: CP2_AZ_TANG = .true. !< seth_test_critical_points: TENO5-AZ held to Tang's 5 at the second-order point.
real(R_P), parameter :: HENRICK_TOL = 0.5_R_P !< seth_test_critical_points: WENO-JS within 0.5 of Henrick's rates.
real(R_P), parameter :: POINT_FLOOR = 2.e2_R_P !< The pointwise order tests: the floor of a face error, in eps.
real(R_P), parameter :: POINT_ORDER_LOSS = 0._R_P !< The pointwise order tests: no order lost.
integer,   parameter :: MP_NG = 7 !< seth_test_mp: the sin^4 grids N = 20 .. 1280.
real(R_P), parameter :: SPECTRAL_FLOOR = 0._R_P !< seth_test_spectral: no floor, the slope between 0.1 and 0.2.
real(R_P), parameter :: SPECTRAL_SLOPE_LOSS = 0._R_P !< seth_test_spectral: no slope lost.
real(R_P), parameter :: TAU_FLOOR = 0._R_P !< seth_test_tables: no floor, tau's order on h = 0.05, 0.025.
real(R_P), parameter :: TAU_ORDER_LOSS = 0._R_P !< seth_test_tables: no order lost.
real(R_P), parameter :: CENTRAL_FLOOR_N = 0._R_P !< seth_test_central: the absolute floor ORDER_FLOOR only.
real(R_P), parameter :: BURGERS_FLOOR_N = 0._R_P !< seth_vv_burgers: the absolute floor ORDER_FLOOR only.
real(R_P), parameter :: EULER_FLOOR_N2 = 0._R_P !< seth_vv_euler_order, seth_vv_roe: the absolute floor only.
real(R_P), parameter :: POS_RANGE_EXP = 8._R_P !< seth_test_positivity: densities and pressures from 1e-8 to 1.
real(R_P), parameter :: SMALL_STEP = 1.e-7_R_P !< seth_test_properties: the small step on the offset 1.
real(R_P), parameter :: AI_STEP_TOL = 1.e-3_R_P !< seth_test_properties: the affine-invariant WENO-5's overshoot on it (1.6e-6
                                              !< gfortran, 3.1e-4 nvfortran with fused multiply-adds).
real(R_P), parameter :: SENSOR_SCALES(5) = [1.e-6_R_P, 1.e-3_R_P, 1._R_P, 1.e3_R_P, 1.e6_R_P] !< seth_test_sensors: the
                                              !< scales a of a f + b.
real(R_P), parameter :: SENSOR_FLOW_RATIO = 3._R_P !< seth_vv_sensors: a cost hybrid's L1 at most 3 times the pure scheme's.
character(*), parameter :: SENSOR_FLOW_EXEMPT = '' !< seth_vv_sensors: every case held to the bound.
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
