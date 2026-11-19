!< Tutorial 11: on the GPU. The advection solver of tutorial 3 with TENO5, its loops offloaded with OpenACC.
program t11_gpu
!run t11_gpu t11_gpu
use seth
use fundal, only : dev_alloc, dev_free, dev_memcpy_to_device, dev_memcpy_from_device
use selah, only : stepper, box, box_interior, SCHEME_SSP_33, FORM_SHU_OSHER_2R   ! the time integrator: SELAH's
implicit none
integer,   parameter :: N = 4096, NG = 3
real(R_P), parameter :: A = 1._R_P, CFL = 0.4_R_P
real(R_P), allocatable :: u0(:), uh(:), ug(:), fp(:), fm(:), f(:)
real(R_P), pointer   :: u(:) => null(), r(:)
type(stepper)        :: st
type(box)            :: bx
real(R_P)            :: h, t, dt, x
integer              :: i, k, steps, ierr

allocate(u0(N), uh(N), ug(1-NG:N+NG), fp(1-NG:N+NG), fm(1-NG:N+NG), f(0:N))
h = 1._R_P/N
do i=1, N
   x = (i - 0.5_R_P)*h
   u0(i) = merge(1._R_P, 0._R_P, x > 0.2_R_P .and. x < 0.4_R_P) + exp(-300*(x - 0.7_R_P)**2)
enddo
!region data
! the field in device memory (FUNDAL: on the GPU with OpenACC, on the host without), and SELAH's registers with it:
! one copy in, one copy out
call dev_alloc(fptr_dev=u, ubounds=[N], lbounds=[1], ierr=ierr)
call dev_memcpy_to_device(dst=u, src=u0)
bx = box_interior(lb=[1], ub=[N], ngc=0, components=.false.)
call st%initialize(SCHEME_SSP_33, FORM_SHU_OSHER_2R, lb=[1], ub=[N], components=.false.)
!$acc data create(ug, fp, fm, f)
t = 0 ; steps = 0
do while (t < 1._R_P)
   dt = min(CFL*h/abs(A), 1._R_P - t)
   call st%begin_step(bx, u)
   do k=1, st%stages()                           ! SELAH's SSP(3,3): the stage state is u itself
      call st%begin_stage(k, bx, dt, u)
      call st%rhs_register_view(r)
      call residual(u, r)
      call st%end_stage(k, bx, dt, u)
   enddo
   call st%finish_step(bx, dt, u)
   t = t + dt ; steps = steps + 1
enddo
!$acc end data
call dev_memcpy_from_device(dst=uh, src=u)
!endregion data
print '(A,I0,A,I0,A)', 'one period of linear advection, TENO5, ', N, ' cells, ', steps, ' steps'
print '(A,ES10.3)', '  L1 error  ', sum(abs(uh - u0))*h
print '(A,L2)', '  in [-1e-6, 1 + 1e-6]', minval(uh) > -1.e-6_R_P .and. maxval(uh) < 1 + 1.e-6_R_P
call st%destroy
call dev_free(u)

contains
   !region residual
   subroutine residual(u, res)
   !< res = -(F_{i+1/2} - F_{i-1/2})/h, every loop on the device; the kernel called inside the face loop. u and res are
   !< device memory (FUNDAL's and SELAH's: deviceptr), the work arrays live in the data region (present).
   real(R_P), intent(in)  :: u(N)
   real(R_P), intent(out) :: res(N)
   real(R_P)              :: vp(-2:3), vm(-2:3), rp, rm
   integer                :: i, m
   !$acc parallel loop deviceptr(u) present(ug)
   do i=1-NG, N+NG                              ! periodic ghost cells
      ug(i) = u(modulo(i - 1, N) + 1)
   enddo
   !$acc parallel loop present(ug, fp, fm)
   do i=1-NG, N+NG
      fp(i) = 0.5_R_P*(A*ug(i) + abs(A)*ug(i))
      fm(i) = 0.5_R_P*(A*ug(i) - abs(A)*ug(i))
   enddo
   !$acc parallel loop present(fp, fm, f) private(vp, vm, rp, rm)
   do i=0, N
      do m=-2, 3                                ! the face's stencils, in arrays of constant bounds
         vp(m) = fp(i+m) ; vm(m) = fm(i+m)
      enddo
      call seth_teno5(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)   ! acc routine seq: one thread per face
      f(i) = rp + rm
   enddo
   !$acc parallel loop present(f) deviceptr(res)
   do i=1, N
      res(i) = -(f(i) - f(i-1))/h
   enddo
   endsubroutine residual
   !endregion residual
endprogram t11_gpu
