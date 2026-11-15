!< SETH benchmark: WENO-JS 5 kernels with run-time array bounds and/or loop extents (de-optimised copies of SETH's).
module seth_bench_dyn_kernels
!< SETH benchmark: WENO-JS 5 kernels with run-time array bounds and/or loop extents (de-optimised copies of SETH's).
!<
!< SETH's kernels have named constants for the bounds of their local arrays and for their loop extents. These copies of
!< the same body (seth_bench_weno_js_dyn_body.INC) give up one or both, to measure each on the device:
!<
!< | Routine | Local arrays | Loop extents |
!< |---|---|---|
!< | `seth_weno5_js` (the library) | constant | constant |
!< | `weno5_js_smax` | constant, over-dimensioned to SMAX = 5 (ADAM's S_MAX) | the dummy S |
!< | `weno5_js_auto_cloop` | automatic (bounds from the dummy S) | constant |
!< | `weno5_js_auto` | automatic | the dummy S |
!<
!< The tables are WENO-5's constants in every variant; S = 3 is passed at run time.
use penf, only : I4P, R8P
use seth_weno_js_tables
implicit none
private
public :: weno5_js_smax, weno5_js_auto_cloop, weno5_js_auto, descaling_s3

integer(I4P), parameter :: SMAX = 5_I4P !< ADAM's S_MAX.

contains
   pure subroutine weno5_js_smax(S, vp, vm, eps, rmu, rp, rm)
   !< Over-dimensioned local arrays (constant bounds), run-time loop extents.
   integer(I4P), intent(in) :: S
   integer(I4P), parameter  :: WEXP = 3_I4P
   real(R8P),    parameter  :: WA(1:2,0:2) = WENO5_JS_A, WP(1:2,0:2,0:2) = WENO5_JS_P, WD(0:2,0:2,0:2) = WENO5_JS_D
#define SB SMAX
#define SL S
#include "seth_bench_weno_js_dyn_body.INC"
#undef SB
#undef SL
   endsubroutine weno5_js_smax

   pure subroutine weno5_js_auto_cloop(S, vp, vm, eps, rmu, rp, rm)
   !< Automatic local arrays, constant loop extents.
   integer(I4P), intent(in) :: S
   integer(I4P), parameter  :: WEXP = 3_I4P
   real(R8P),    parameter  :: WA(1:2,0:2) = WENO5_JS_A, WP(1:2,0:2,0:2) = WENO5_JS_P, WD(0:2,0:2,0:2) = WENO5_JS_D
#define SB S
#define SL 3
#include "seth_bench_weno_js_dyn_body.INC"
#undef SB
#undef SL
   endsubroutine weno5_js_auto_cloop

   pure subroutine weno5_js_auto(S, vp, vm, eps, rmu, rp, rm)
   !< Automatic local arrays, run-time loop extents.
   integer(I4P), intent(in) :: S
   integer(I4P), parameter  :: WEXP = 3_I4P
   real(R8P),    parameter  :: WA(1:2,0:2) = WENO5_JS_A, WP(1:2,0:2,0:2) = WENO5_JS_P, WD(0:2,0:2,0:2) = WENO5_JS_D
#define SB S
#define SL S
#include "seth_bench_weno_js_dyn_body.INC"
#undef SB
#undef SL
   endsubroutine weno5_js_auto

   pure function descaling_s3(sigma, vp, vm) result(rmu)
   !< seth_weno_descaling with S a constant (3).
   integer(I4P), parameter  :: S = 3_I4P
   real(R8P),    intent(in) :: sigma, vp(1-S:S), vm(1-S:S)
   real(R8P)                :: rmu, mu
   integer(I4P)             :: m
   !$acc routine seq
   !$omp declare target
   mu = 0._R8P
   do m=1-S, S-1
      mu = mu + abs(vm(m+1))
      mu = mu + abs(vp(m))
   enddo
   mu = max(mu/real(2*(2*S-1), R8P), tiny(1._R8P))
   rmu = sigma*(1._R8P/mu) + (1._R8P - sigma)
   endfunction descaling_s3
endmodule seth_bench_dyn_kernels
