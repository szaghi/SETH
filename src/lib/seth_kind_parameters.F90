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
public :: SETH_TENO_EPS, SETH_POSITIVITY_BACKOFF

#if defined _R_P_IS_R4P
integer, parameter :: KIND_CHECK = 1/merge(1, 0, R_P == R4P) !< Stops the build if PENF's R_P is not R4P.
#error "SETH: the R4P parameter set is not defined yet (roadmap issue #4, 6.2, phase 2)"
#elif defined _R_P_IS_R16P
#error "SETH: R16P is not supported (no device code, no caller)"
#else
integer, parameter :: KIND_CHECK = 1/merge(1, 0, R_P == R8P) !< Stops the build if PENF's R_P is not R8P.
! R8P: the values the kernels were verified with.
real(R_P), parameter :: SETH_TENO_EPS = 1.e-40_R_P !< eps of the TENO papers (Li, Fu & Adams 2020, eq. 9): normal in
                                                   !< R8P (smallest normal 2.2e-308).
real(R_P), parameter :: SETH_POSITIVITY_BACKOFF = 1.e-15_R_P !< Relative step by which the positivity limiter backs
                                                             !< theta off (doubled per retry): a few ulps of R8P.
#endif
endmodule seth_kind_parameters
