!< Tutorial 9: a hybrid. A sensor picks, face by face, the linear scheme where the data is smooth and TENO5 elsewhere.
program t09_hybrid
!run t09_hybrid t09_hybrid
use seth
use selah, only : stepper, box, box_interior, SCHEME_SSP_33, FORM_SHU_OSHER_2R   ! the time integrator: SELAH's
implicit none
integer,   parameter :: N = 200, NG = 3
real(R_P), parameter :: A = 1._R_P, CFL = 0.4_R_P
integer,   parameter :: TENO5 = 1, LINEAR5 = 2, HYBRID = 3
character(8)         :: label(3) = [character(8) :: 'TENO5', 'linear5', 'hybrid']
real(R_P)            :: u0(N), u(N), h, x
integer              :: s, i
integer(I8P)         :: nfaces, nnonlinear

h = 1._R_P/N
do i=1, N                                       ! the profile of tutorial 3
   x = (i - 0.5_R_P)*h
   u0(i) = merge(1._R_P, 0._R_P, x > 0.2_R_P .and. x < 0.4_R_P) + exp(-300*(x - 0.7_R_P)**2)
enddo
print '(A)', 'one period of linear advection, 200 cells; the hybrid with Li, Chen & Ren''s sensor'
print '(A)', 'scheme      L1 error  undershoot  overshoot   nonlinear faces'
do s=1, 3
   u = u0 ; nfaces = 0 ; nnonlinear = 0
   call advect(s, u)
   print '(A8,3ES11.2,F14.1,A)', label(s), sum(abs(u - u0))*h, max(-minval(u), 0._R_P), max(maxval(u) - 1, 0._R_P), &
                                 100*real(nnonlinear, R_P)/max(nfaces, 1_I8P), ' %'
enddo

contains
   subroutine residual(s, u, res)
   integer,   intent(in)  :: s
   real(R_P), intent(in)  :: u(N)
   real(R_P), intent(out) :: res(N)
   real(R_P)              :: ug(1-NG:N+NG), fp(1-NG:N+NG), fm(1-NG:N+NG), f(0:N), rp, rm, rmu, sp, sm
   integer                :: i
   ug(1:N) = u ; ug(1-NG:0) = u(N+1-NG:N) ; ug(N+1:N+NG) = u(1:NG)
   fp = 0.5_R_P*(A*ug + abs(A)*ug)
   fm = 0.5_R_P*(A*ug - abs(A)*ug)
   rmu = 1/max(maxval(u) - minval(u), tiny(1._R_P))  ! the field's range: every sensor scale and offset invariant
   do i=0, N
      select case(s)
      case(TENO5)
         call seth_teno5(fp(i-2:i+3), fm(i-2:i+3), SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)
         nnonlinear = nnonlinear + 1
      case(LINEAR5)
         call seth_linear5(fp(i-2:i+3), fm(i-2:i+3), rp, rm)
      case(HYBRID)
         !region hybrid
         ! the sensor on the data around the face, once per bias (the right bias reads the mirrored cells)
         call seth_sensor_lcr(ug(i-2:i+3), rmu, SETH_SENSOR_LCR_PSIC, sp)
         call seth_sensor_lcr(ug(i+3:i-2:-1), rmu, SETH_SENSOR_LCR_PSIC, sm)
         if (sp == 1 .and. sm == 1) then        ! both smooth: the linear scheme, no weights
            call seth_linear5(fp(i-2:i+3), fm(i-2:i+3), rp, rm)
         else
            call seth_teno5(fp(i-2:i+3), fm(i-2:i+3), SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)
            nnonlinear = nnonlinear + 1
         endif
         !endregion hybrid
      endselect
      nfaces = nfaces + 1
      f(i) = rp + rm
   enddo
   res = -(f(1:N) - f(0:N-1))/h
   endsubroutine residual

   subroutine advect(s, u)
   integer,   intent(in)    :: s
   real(R_P), intent(inout) :: u(N)
   real(R_P)                :: t, dt
   type(stepper)            :: st
   type(box)                :: bx
   real(R_P), pointer       :: r(:)
   integer                  :: k
   bx = box_interior(lb=[1], ub=[N], ngc=0, components=.false.)
   call st%initialize(SCHEME_SSP_33, FORM_SHU_OSHER_2R, lb=[1], ub=[N], components=.false.)
   t = 0
   do while (t < 1._R_P)
      dt = min(CFL*h/abs(A), 1._R_P - t)
      call st%begin_step(bx, u)
      do k=1, st%stages()                          ! SELAH's SSP(3,3): the stage state is u itself
         call st%begin_stage(k, bx, dt, u)
         call st%rhs_register_view(r)
         call residual(s, u, r)
         call st%end_stage(k, bx, dt, u)
      enddo
      call st%finish_step(bx, dt, u)
      t = t + dt
   enddo
   call st%destroy
   endsubroutine advect
endprogram t09_hybrid
