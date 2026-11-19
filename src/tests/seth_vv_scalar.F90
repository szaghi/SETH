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
use selah, only : stepper, box, box_interior, SCHEME_SSP_33, FORM_SHU_OSHER_2R ! the time integrator (SELAH #2)
implicit none
private
public :: FLUX_ADVECTION, FLUX_BURGERS, scalar_step

integer(I4P), parameter :: FLUX_ADVECTION = 1, FLUX_BURGERS = 2

contains
   subroutine residual(id, kind, u, n, h, res, par)
   integer(I4P), intent(in)           :: id, kind, n
   real(R_P),    intent(in)           :: u(0:n-1), h
   real(R_P),    intent(out)          :: res(0:n-1)
   real(R_P),    intent(in), optional :: par
   real(R_P)                          :: f(0:n-1), fh(0:n-1), alpha, vp(VLO:VHI), vm(VLO:VHI), rp, rm
   integer(I4P)                       :: i, m, j
   if (kind == FLUX_ADVECTION) then
      f = u ; alpha = 1._R_P
   else
      f = 0.5_R_P*u*u ; alpha = 1.1_R_P*maxval(abs(u))
   endif
   do i=0, n-1 ! face i+1/2
      do m=VLO, VHI
         j = modulo(i+m, n)
         vp(m) = 0.5_R_P*(f(j) + alpha*u(j))
         vm(m) = 0.5_R_P*(f(j) - alpha*u(j))
      enddo
      call reconstruct(id, vp, vm, rp, rm, par=par)
      fh(i) = rp + rm
   enddo
   do i=0, n-1
      res(i) = -(fh(i) - fh(modulo(i-1, n)))/h
   enddo
   endsubroutine residual

   subroutine scalar_step(id, kind, u, n, h, dt, par)
   !< One SSP-RK3 step: SELAH's SSP(3,3) in its two-register Shu-Osher form (the stage state is u itself), on the
   !< periodic rank-1 field u(0:n-1) with no ghost cells.
   integer(I4P), intent(in)           :: id, kind, n
   real(R_P),    intent(inout)        :: u(0:n-1)
   real(R_P),    intent(in)           :: h, dt
   real(R_P),    intent(in), optional :: par
   type(stepper)                      :: st
   type(box)                          :: bx
   real(R_P), pointer                 :: r(:)
   integer(I4P)                       :: k
   bx = box_interior(lb=[0_I4P], ub=[n-1], ngc=0_I4P, components=.false.)
   call st%initialize(SCHEME_SSP_33, FORM_SHU_OSHER_2R, lb=[0_I4P], ub=[n-1], components=.false.)
   call st%begin_step(bx, u)
   do k=1, st%stages()
      call st%begin_stage(k, bx, dt, u)
      call st%rhs_register_view(r)
      call residual(id, kind, u, n, h, r, par)
      call st%end_stage(k, bx, dt, u)
   enddo
   call st%finish_step(bx, dt, u)
   call st%destroy
   endsubroutine scalar_step
endmodule seth_vv_scalar
