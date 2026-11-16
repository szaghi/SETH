!< SETH, the Ducros shock sensors: where a hybrid's central flux must give way to the shock-capturing one.
module seth_shock_sensors
!< SETH, the Ducros shock sensors: where a hybrid's central flux must give way to the shock-capturing one.
!<
!< Pointwise, from the velocity gradient g(a,b) = d u_a / d x_b at a cell (the caller computes it; its stencil, the
!< mapping to faces and the dilation of the flags are the caller's policy):
!<
!<    call seth_ducros(eps, g, phi)              phi = theta^2/(theta^2 + |omega|^2 + eps)          (Ducros et al. 1999)
!<    call seth_ducros_compression(eps, g, psi)  psi = max(-theta/sqrt(theta^2 + |omega|^2 + eps), 0) (De Vanna et al.
!<                                                                                                   2023a, Eq. (9))
!<
!< theta = div u, omega = curl u, eps a floor in the units of a gradient squared (De Vanna: (u_inf/delta_0)^2; ADAM:
!< a fraction of (c/dx)^2). phi is even in theta and flags expansions too; psi flags compressions only. In compression
!< psi = sqrt(phi) when both use the same eps: a threshold of one is not a threshold of the other (0.1, 0.2, 0.4 of psi
!< are 0.01, 0.04, 0.16 of phi). ADAM's compute_ducros is seth_ducros plus a density-jump override (its policy). Pure,
!< `acc routine seq`, `omp declare target`.
use penf, only : R_P
implicit none
private
public :: seth_ducros, seth_ducros_compression

contains
   pure subroutine seth_ducros(eps, g, phi)
   !< The Ducros sensor, theta^2/(theta^2 + |omega|^2 + eps).
   real(R_P), intent(in)  :: eps    !< Floor of the denominator.
   real(R_P), intent(in)  :: g(3,3) !< Velocity gradient, g(a,b) = d u_a / d x_b.
   real(R_P), intent(out) :: phi    !< Sensor in [0, 1].
   real(R_P)              :: th, wx, wy, wz
   !$acc routine seq
   !$omp declare target

   th = g(1,1) + g(2,2) + g(3,3)
   wx = g(3,2) - g(2,3) ; wy = g(1,3) - g(3,1) ; wz = g(2,1) - g(1,2)
   phi = th*th/(th*th + wx*wx + wy*wy + wz*wz + eps)
   endsubroutine seth_ducros

   pure subroutine seth_ducros_compression(eps, g, psi)
   !< The compression-only Ducros sensor of De Vanna et al., max(-theta/sqrt(theta^2 + |omega|^2 + eps), 0).
   real(R_P), intent(in)  :: eps    !< Floor of the denominator.
   real(R_P), intent(in)  :: g(3,3) !< Velocity gradient, g(a,b) = d u_a / d x_b.
   real(R_P), intent(out) :: psi    !< Sensor in [0, 1].
   real(R_P)              :: th, wx, wy, wz
   !$acc routine seq
   !$omp declare target

   th = g(1,1) + g(2,2) + g(3,3)
   wx = g(3,2) - g(2,3) ; wy = g(1,3) - g(3,1) ; wz = g(2,1) - g(1,2)
   psi = max(-th/sqrt(th*th + wx*wx + wy*wy + wz*wz + eps), 0._R_P)
   endsubroutine seth_ducros_compression
endmodule seth_shock_sensors
