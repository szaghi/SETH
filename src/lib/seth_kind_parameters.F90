!< SETH, the working real kind and the literals that depend on it.
module seth_kind_parameters
!< SETH, the working real kind and the literals that depend on it.
!<
!< Every kernel computes in PENF's R_P, chosen when SETH (and PENF, compiled from its sources by the same build mode) is
!< built: -D_R_P_IS_R4P, -D_R_P_IS_R8P, or nothing (R8P, PENF's default). One precision per build, no generic
!< interfaces: the kernels are compiled once, for that kind.
!<
!< A literal whose meaning depends on the kind lives here, one set per kind behind the same macros, never in a kernel:
!< a regulariser below the smallest normal number of a kind underflows (1e-40 in R4P), a step below its epsilon does
!< nothing (1 - 1e-15 = 1 in R4P). Each value states why it is what it is. Roadmap issue #4, 6.2 (the specification is
!< in its comment of 2026-11-16).
!<
!< A compile-time check ties the macro to PENF's R_P: a build that defines it for SETH but not for PENF (or the
!< reverse) stops with a division by zero in a constant expression.
use penf, only : R_P, R4P, R8P
implicit none
private
public :: R_P
public :: SETH_TENO_EPS, SETH_POSITIVITY_BACKOFF, SETH_NONPOLY_GUARD, SETH_WENO6_EPS

#if defined _R_P_IS_R4P
integer, parameter :: KIND_CHECK = 1/merge(1, 0, R_P == R4P) !< Stops the build if PENF's R_P is not R4P.
! R4P: epsilon 1.2e-7, smallest normal 1.2e-38, largest 3.4e38.
real(R_P), parameter :: SETH_TENO_EPS = 1.e-30_R_P !< The papers' 1e-40 is below the smallest normal R4P number (and
                                                   !< flushed to zero on a GPU): 1e-30 is normal, and negligible
                                                   !< against any indicator R4P data can produce (a difference of one
                                                   !< ulp of unit data squares to 1.4e-14). The R4P kernels never
                                                   !< divide by it alone (their weights are ratios: seth_z_teno_body).
real(R_P), parameter :: SETH_POSITIVITY_BACKOFF = 5.e-7_R_P !< About 4 ulps of R4P, as 1e-15 is of R8P: 1 - 1e-15 is 1
                                                            !< in R4P, and the back-off would do nothing.
real(R_P), parameter :: SETH_NONPOLY_GUARD = 1.e-15_R_P !< The first-order guard of the THINC and DR overrides,
                                                       !< (f_{i+1}-f_i)(f_i-f_{i-1}) <= 1e-15 (Takagi 2022, Eq. 17;
                                                       !< Liang 2022, Eq. 33): the papers' number in both kinds. It
                                                       !< bounds a product of two differences of the data (dimensional,
                                                       !< scaled by rmu^2 in the kernels), not a step against 1, and is
                                                       !< a normal R4P number.
real(R_P), parameter :: SETH_WENO6_EPS = 1.e-30_R_P !< Guard of seth_weno6_point's WENO-Z weights against 0/0 on
                                                    !< constant data: ADAM's 1e-40 is below the smallest normal R4P
                                                    !< number, 1e-30 is normal (as SETH_TENO_EPS). No eps small
                                                    !< enough to be negligible keeps (tau/(beta + eps))**2 finite
                                                    !< in R4P next to a flat sub-stencil (beta = 0, eps = 1e-30:
                                                    !< it overflows once tau passes 2e-11, the ratio itself once
                                                    !< tau passes 3e8), so the R4P kernel scales the weights by
                                                    !< the smallest beta + eps when tau exceeds it.
#elif defined _R_P_IS_R16P
#error "SETH: R16P is not supported (no device code, no caller)"
#else
integer, parameter :: KIND_CHECK = 1/merge(1, 0, R_P == R8P) !< Stops the build if PENF's R_P is not R8P.
! R8P: the values the kernels were verified with.
real(R_P), parameter :: SETH_TENO_EPS = 1.e-40_R_P !< eps of the TENO papers (Li, Fu & Adams 2020, eq. 9): normal in
                                                   !< R8P (smallest normal 2.2e-308).
real(R_P), parameter :: SETH_POSITIVITY_BACKOFF = 1.e-15_R_P !< Relative step by which the positivity limiter backs
                                                             !< theta off (doubled per retry): a few ulps of R8P.
real(R_P), parameter :: SETH_NONPOLY_GUARD = 1.e-15_R_P !< The first-order guard of the THINC and DR overrides,
                                                       !< (f_{i+1}-f_i)(f_i-f_{i-1}) <= 1e-15 (Takagi 2022, Eq. 17;
                                                       !< Liang 2022, Eq. 33), on the data scaled by rmu.
real(R_P), parameter :: SETH_WENO6_EPS = 1.e-40_R_P !< Guard of seth_weno6_point's WENO-Z weights against 0/0 on
                                                    !< constant data: ADAM's SEAM_WENO6_EPS, the value its seam
                                                    !< transfer was verified with (normal in R8P).
#endif
endmodule seth_kind_parameters
