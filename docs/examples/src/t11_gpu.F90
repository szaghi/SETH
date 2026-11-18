!< Tutorial 11: on the GPU. The advection solver of tutorial 3 with TENO5, its loops offloaded with OpenACC.
program t11_gpu
!run t11_gpu t11_gpu
use seth
implicit none
integer,   parameter :: N = 4096, NG = 3
real(R_P), parameter :: A = 1._R_P, CFL = 0.4_R_P
real(R_P), allocatable :: u(:), u0(:), u1(:), u2(:), res(:), ug(:), fp(:), fm(:), f(:)
real(R_P)            :: h, t, dt, x
integer              :: i, steps

allocate(u(N), u0(N), u1(N), u2(N), res(N), ug(1-NG:N+NG), fp(1-NG:N+NG), fm(1-NG:N+NG), f(0:N))
h = 1._R_P/N
do i=1, N
   x = (i - 0.5_R_P)*h
   u0(i) = merge(1._R_P, 0._R_P, x > 0.2_R_P .and. x < 0.4_R_P) + exp(-300*(x - 0.7_R_P)**2)
enddo
u = u0
!region data
! the fields live on the device for the whole run: one copy in, one copy out
!$acc data copy(u) create(u1, u2, res, ug, fp, fm, f)
t = 0 ; steps = 0
do while (t < 1._R_P)
   dt = min(CFL*h/abs(A), 1._R_P - t)
   call residual(u, res)
   !$acc parallel loop present(u, u1, res)
   do i=1, N
      u1(i) = u(i) + dt*res(i)
   enddo
   call residual(u1, res)
   !$acc parallel loop present(u, u1, u2, res)
   do i=1, N
      u2(i) = 0.75_R_P*u(i) + 0.25_R_P*(u1(i) + dt*res(i))
   enddo
   call residual(u2, res)
   !$acc parallel loop present(u, u2, res)
   do i=1, N
      u(i) = u(i)/3 + 2._R_P/3*(u2(i) + dt*res(i))
   enddo
   t = t + dt ; steps = steps + 1
enddo
!$acc end data
!endregion data
print '(A,I0,A,I0,A)', 'one period of linear advection, TENO5, ', N, ' cells, ', steps, ' steps'
print '(A,ES10.3)', '  L1 error  ', sum(abs(u - u0))*h
print '(A,L2)', '  in [-1e-6, 1 + 1e-6]', minval(u) > -1.e-6_R_P .and. maxval(u) < 1 + 1.e-6_R_P

contains
   !region residual
   subroutine residual(u, res)
   !< res = -(F_{i+1/2} - F_{i-1/2})/h, every loop on the device; the kernel called inside the face loop.
   real(R_P), intent(in)  :: u(N)
   real(R_P), intent(out) :: res(N)
   real(R_P)              :: vp(-2:3), vm(-2:3), rp, rm
   integer                :: i, m
   !$acc parallel loop present(u, ug)
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
   !$acc parallel loop present(f, res)
   do i=1, N
      res(i) = -(f(i) - f(i-1))/h
   enddo
   endsubroutine residual
   !endregion residual
endprogram t11_gpu
