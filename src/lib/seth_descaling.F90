!< SETH, the generic descaler of the smoothness indicators (scale-invariant WENO weights).
module seth_descaling
!< SETH, the generic descaler of the smoothness indicators (scale-invariant WENO weights).
!<
!< Scale-invariant weights (Don, Li, Wang & Wang 2022, JCP 448, 110724) read `(IS rmu) rmu` in place of `IS`, so that
!< the regulariser acts relative to the magnitude of the data. `seth_weno_descaling` returns `rmu = 1/mu`, `mu` the
!< mean |value| over both stencils of a face (floored at `tiny`): the generic magnitude, for callers without a better
!< one (a field near zero by construction should pass a magnitude of its own). `sigma = 0` gives `rmu = 1` exactly,
!< the Jiang-Shu weights; `sigma = 1` the scale-invariant ones. ADAM's weno_descaling, in SETH's stencil convention.
use penf, only : I4P, R8P
implicit none
private
public :: seth_weno_descaling

contains
   pure function seth_weno_descaling(S, sigma, vp, vm) result(rmu)
   !< Inverse descaler `1/mu` of the smoothness indicators of one face, from its own stencils.
   integer(I4P), intent(in) :: S         !< Candidates of the scheme (WENO-(2S-1)).
   real(R8P),    intent(in) :: sigma     !< Descaler switch: 0 Jiang-Shu, 1 scale-invariant.
   real(R8P),    intent(in) :: vp(1-S:S) !< Values reconstructed with a left bias.
   real(R8P),    intent(in) :: vm(1-S:S) !< Values reconstructed with a right bias.
   real(R8P)                :: rmu       !< Inverse of the descaler.
   real(R8P)                :: mu        !< Descaler.
   integer(I4P)             :: m         !< Counter.
   !$acc routine seq
   !$omp declare target

   mu = 0._R8P
   do m=1-S, S-1
      mu = mu + abs(vm(m+1)) ! ADAM's row 1 (vm shifted one cell)
      mu = mu + abs(vp(m))   ! ADAM's row 2
   enddo
   mu = max(mu/real(2*(2*S-1), R8P), tiny(1._R8P))
   rmu = sigma*(1._R8P/mu) + (1._R8P - sigma)
   endfunction seth_weno_descaling
endmodule seth_descaling
