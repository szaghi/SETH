!< SETH verification: scalar conservation laws on a periodic grid (linear advection, Burgers), around the kernels.
module seth_vv_scalar
!< SETH verification: scalar conservation laws on a periodic grid (linear advection, Burgers), around the kernels.
!<
!< u_t + f(u)_x = 0 on N periodic cells of width h, point values at the cell centres, conservative finite differences
!< with global Lax-Friedrichs splitting f+- = (f +- alpha u)/2, alpha = 1.1 max |f'(u)| for Burgers (alpha = a for
!< advection, so the splitting is the upwind one). The 10% margin matters: with alpha = max|u| exactly,
!< d f-/dx = (u - alpha) u_x/2 vanishes to second order at the maximum of u, a second-order critical point where TENO
!< drops to second order (measured: TENO6 on Burgers before the shock 5.7, 5.0, 4.7, 4.2, 4.1 from N = 80 to 1280
!< without it; 5.7, 5.8, 6.0, 6.0 with it). f+ is reconstructed with the left bias (vp), f- with the right bias (vm),
!< the flux F_{i+1/2} = rp + rm. SSP-RK3 in time.
use seth
use seth_test_schemes, only : reconstruct, VLO, VHI
implicit none
private
public :: FLUX_ADVECTION, FLUX_BURGERS, scalar_step

integer(I4P), parameter :: FLUX_ADVECTION = 1, FLUX_BURGERS = 2

contains
   subroutine residual(id, kind, u, n, h, res, par)
   integer(I4P), intent(in)           :: id, kind, n
   real(R8P),    intent(in)           :: u(0:n-1), h
   real(R8P),    intent(out)          :: res(0:n-1)
   real(R8P),    intent(in), optional :: par
   real(R8P)                          :: f(0:n-1), fh(0:n-1), alpha, vp(VLO:VHI), vm(VLO:VHI), rp, rm
   integer(I4P)                       :: i, m, j
   if (kind == FLUX_ADVECTION) then
      f = u ; alpha = 1._R8P
   else
      f = 0.5_R8P*u*u ; alpha = 1.1_R8P*maxval(abs(u))
   endif
   do i=0, n-1 ! face i+1/2
      do m=VLO, VHI
         j = modulo(i+m, n)
         vp(m) = 0.5_R8P*(f(j) + alpha*u(j))
         vm(m) = 0.5_R8P*(f(j) - alpha*u(j))
      enddo
      call reconstruct(id, vp, vm, rp, rm, par=par)
      fh(i) = rp + rm
   enddo
   do i=0, n-1
      res(i) = -(fh(i) - fh(modulo(i-1, n)))/h
   enddo
   endsubroutine residual

   subroutine scalar_step(id, kind, u, n, h, dt, par)
   !< One SSP-RK3 step.
   integer(I4P), intent(in)           :: id, kind, n
   real(R8P),    intent(inout)        :: u(0:n-1)
   real(R8P),    intent(in)           :: h, dt
   real(R8P),    intent(in), optional :: par
   real(R8P)                          :: u1(0:n-1), u2(0:n-1), r(0:n-1)
   call residual(id, kind, u, n, h, r, par)  ; u1 = u + dt*r
   call residual(id, kind, u1, n, h, r, par) ; u2 = 0.75_R8P*u + 0.25_R8P*(u1 + dt*r)
   call residual(id, kind, u2, n, h, r, par) ; u = u/3._R8P + 2._R8P/3._R8P*(u2 + dt*r)
   endsubroutine scalar_step
endmodule seth_vv_scalar
