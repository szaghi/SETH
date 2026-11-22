!< SETH benchmark: the three symbols ADAM's FNL WENO kernels import from adam_weno_object, without the rest of it.
module adam_weno_object
!< SETH benchmark: the three symbols ADAM's FNL WENO kernels import from adam_weno_object, without the rest of it.
!<
!< S_MAX and S_MAX_M1 as ADAM defines them; weno_weights_exponent copied verbatim from ADAM's adam_weno_object.F90
!< (commit 302f093e). The kernel files next to it (adam_fnl_weno_kernels.F90, adam_fnl_weno_width_kernels.F90,
!< adam_fnl_weno_width_body.INC) are ADAM's, byte for byte.
!< Frozen history: ADAM deleted those files and weno_weights_exponent in ADAM commit ba37c4cb (ADAM issue 139), after
!< its FLUME moved to SETH's kernels; the copies stay unchanged as the benchmark's baseline.
use penf, only : I4P
implicit none
private
public :: S_MAX, S_MAX_M1, weno_weights_exponent
integer(I4P), parameter :: S_MAX=5    !< Maximum number/dimensions of stencils.
integer(I4P), parameter :: S_MAX_M1=4 !< Maximum number/dimensions of stencils minus 1.

contains
   pure function weno_weights_exponent(S) result(wexp)
   !< Return the exponent of the smoothness indicators in the WENO nonlinear weights, `a = d / (zeps + IS)**wexp`.
   !<
   !< The one definition used by the host weights, the device weights and `weno_object%wexp`: `wexp = S`, reduced to
   !< `S - 1` for `S > 4`. Before, the host used `S` and the device `2`, so the two backends ran different schemes for
   !< `S > 2`. Being a function of `S`, the exponent also follows the stencil of each ROR iteration.
   integer(I4P), intent(in) :: S    !< Number of stencils used.
   integer(I4P)             :: wexp !< Exponent of the smoothness indicators.
   !$acc routine seq
   !$omp declare target

   wexp = S
   if (S > 4) wexp = S - 1
   endfunction weno_weights_exponent
endmodule adam_weno_object
