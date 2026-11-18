!< Tutorial 3: a conservation law. Linear advection by finite volumes: flux splitting, face fluxes, SSP-RK3.
program t03_advection
!run t03_advection t03_advection
use seth
implicit none
integer,   parameter :: N = 200, NG = 3         ! cells, ghost cells
real(R_P), parameter :: A = 1._R_P              ! advection speed: u_t + (A u)_x = 0
real(R_P), parameter :: CFL = 0.4_R_P
integer,   parameter :: WENO5 = 1, TENO5 = 2, LINEAR5 = 3
character(8)         :: label(3) = [character(8) :: 'WENO5-JS', 'TENO5', 'linear5']
real(R_P)            :: u0(N), u(N), h
integer              :: s, i

h = 1._R_P/N
do i=1, N                                       ! a square wave and a smooth bump, sampled at the cell centres
   u0(i) = profile((i - 0.5_R_P)*h)
enddo
print '(A)', 'one period of linear advection, 200 cells, CFL 0.4'
print '(A)', 'scheme      L1 error  undershoot  overshoot'
do s=1, 3
   u = u0
   call advect(s, u)
   print '(A8,3ES11.2)', label(s), sum(abs(u - u0))*h, max(-minval(u), 0._R_P), max(maxval(u) - 1, 0._R_P)
enddo

contains
   pure function profile(x) result(v)
   real(R_P), intent(in) :: x
   real(R_P)             :: v
   v = merge(1._R_P, 0._R_P, x > 0.2_R_P .and. x < 0.4_R_P) + exp(-300*(x - 0.7_R_P)**2)
   endfunction profile

   !region residual
   subroutine residual(s, u, res)
   !< res = -(F_{i+1/2} - F_{i-1/2})/h, the face fluxes from the scheme s.
   integer,   intent(in)  :: s
   real(R_P), intent(in)  :: u(N)
   real(R_P), intent(out) :: res(N)
   real(R_P)              :: ug(1-NG:N+NG), fp(1-NG:N+NG), fm(1-NG:N+NG), f(0:N), rp, rm
   integer                :: i
   ug(1:N) = u ; ug(1-NG:0) = u(N+1-NG:N) ; ug(N+1:N+NG) = u(1:NG)   ! periodic ghost cells
   fp = 0.5_R_P*(A*ug + abs(A)*ug)              ! Lax-Friedrichs splitting, f = f+ + f-
   fm = 0.5_R_P*(A*ug - abs(A)*ug)
   do i=0, N                                    ! the face i+1/2: f+ reconstructed from the left, f- from the right
      select case(s)
      case(WENO5)
         call seth_weno5_js(fp(i-2:i+3), fm(i-2:i+3), SETH_WENO_JS_EPS, 1._R_P, rp, rm)
      case(TENO5)
         call seth_teno5(fp(i-2:i+3), fm(i-2:i+3), SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)
      case default
         call seth_linear5(fp(i-2:i+3), fm(i-2:i+3), rp, rm)
      endselect
      f(i) = rp + rm
   enddo
   res = -(f(1:N) - f(0:N-1))/h
   endsubroutine residual
   !endregion residual

   !region rk3
   subroutine advect(s, u)
   !< One period, by the three-stage SSP Runge-Kutta scheme of Shu and Osher.
   integer,   intent(in)    :: s
   real(R_P), intent(inout) :: u(N)
   real(R_P)                :: u1(N), u2(N), res(N), t, dt
   t = 0
   do while (t < 1._R_P)
      dt = min(CFL*h/abs(A), 1._R_P - t)
      call residual(s, u, res)  ; u1 = u + dt*res
      call residual(s, u1, res) ; u2 = 0.75_R_P*u + 0.25_R_P*(u1 + dt*res)
      call residual(s, u2, res) ; u  = u/3 + 2._R_P/3*(u2 + dt*res)
      t = t + dt
   enddo
   endsubroutine advect
   !endregion rk3
endprogram t03_advection
