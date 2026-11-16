!< SETH, the generic descaler of the smoothness indicators (scale-invariant WENO weights).
module seth_descaling
!< SETH, the generic descaler of the smoothness indicators (scale-invariant WENO weights).
!<
!< Scale-invariant weights (Don, Li, Wang & Wang 2022, JCP 448, 110724) read `(IS rmu) rmu` in place of `IS`, so that
!< the regulariser acts relative to the magnitude of the data. `seth_weno_descaling` returns `rmu = 1/mu`, `mu` the
!< mean |value| over both stencils of a face (floored at `tiny`): the generic magnitude, for callers without a better
!< one (a field near zero by construction should pass a magnitude of its own). `sigma = 0` gives `rmu = 1` exactly,
!< the Jiang-Shu weights; `sigma = 1` the scale-invariant ones. ADAM's weno_descaling, in SETH's stencil convention.
!<
!< `seth_weno_descaling_ai` is the affine-invariant descaler of Wang & Don (2022, Appl. Numer. Math. 181:630, Eq. (13)):
!< the mean deviation from the stencil's mean, mu = mean |v_j - mean v|, instead of the mean |v_j|. The indicators are
!< blind to a constant, so with this mu the weights are blind to it too: small fluctuations riding on a large offset
!< (pressure, density, temperature) keep the ENO property, which the mean |v| loses (it is then the offset, and eps
!< swamps the indicators). Each bias over its own stencil (left: vp(1-S:S-1), right: vm(2-S:S)), the two averaged; mu
!< floored at 64 eps mean|v| and tiny, so that exactly constant data (where the weights do not matter) cannot send the
!< indicators' round-off times rmu^2 to overflow. Scaling by 2^k is exact, so the result stays bitwise scale invariant.
use penf, only : I4P, R8P
implicit none
private
public :: seth_weno_descaling
public :: seth_weno_descaling_ai

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

   pure function seth_weno_descaling_ai(S, vp, vm) result(rmu)
   !< Inverse affine-invariant descaler `1/mu`, mu = the mean deviation from the mean of each bias's stencil (Wang & Don).
   integer(I4P), intent(in) :: S         !< Candidates of the scheme (WENO-(2S-1)).
   real(R8P),    intent(in) :: vp(1-S:S) !< Values reconstructed with a left bias.
   real(R8P),    intent(in) :: vm(1-S:S) !< Values reconstructed with a right bias.
   real(R8P)                :: rmu       !< Inverse of the descaler.
   real(R8P)                :: avp, avm  !< Means of the two stencils.
   real(R8P)                :: mu, mabs  !< Descaler; mean |v|, for the floor.
   integer(I4P)             :: m         !< Counter.
   !$acc routine seq
   !$omp declare target

   avp = 0._R8P ; avm = 0._R8P ; mabs = 0._R8P
   do m=1-S, S-1
      avp = avp + vp(m)
      avm = avm + vm(m+1)
      mabs = mabs + abs(vp(m)) + abs(vm(m+1))
   enddo
   avp = avp/real(2*S-1, R8P) ; avm = avm/real(2*S-1, R8P)
   mu = 0._R8P
   do m=1-S, S-1
      mu = mu + abs(vp(m) - avp) + abs(vm(m+1) - avm)
   enddo
   mu = mu/real(2*(2*S-1), R8P)
   mabs = mabs/real(2*(2*S-1), R8P)
   mu = max(mu, 64._R8P*epsilon(1._R8P)*mabs, tiny(1._R8P))
   rmu = 1._R8P/mu
   endfunction seth_weno_descaling_ai
endmodule seth_descaling
