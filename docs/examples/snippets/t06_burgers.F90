!< Tutorial 6: a nonlinear flux. Burgers' equation, a sine steepening into a shock.
program t06_burgers
use seth
implicit none
real(R_P), parameter :: PI = acos(-1._R_P)
integer,   parameter :: N = 100, NG = 4          ! WENO-Z7 reads cells -3..4: four ghost cells
integer,   parameter :: WENO5 = 1, WENOZ7 = 2, TENO5 = 3, LINEAR5 = 4
character(8)         :: label(4) = [character(8) :: 'WENO5-JS', 'WENO-Z7', 'TENO5', 'linear5']
real(R_P)            :: u(N), h, tv0
integer              :: s, i

h = 1._R_P/N
print '(A)', 'u_t + (u^2/2)_x = 0, u0 = 0.5 + sin(2 pi x), 100 cells, t = 0.5 (shock from t = 1/(2 pi))'
print '(A)', 'scheme     extrema  total variation  mass conserved'
do s=1, 4
   do i=1, N                                    ! the cell averages of 0.5 + sin(2 pi x)
      u(i) = 0.5_R_P + (cos(2*PI*(i - 1)*h) - cos(2*PI*i*h))/(2*PI*h)
   enddo
   tv0 = sum(abs(u - cshift(u, 1)))
   call evolve(s, u, 0.5_R_P)
   print '(A8,I10,F17.4,L16)', label(s), extrema(u), sum(abs(u - cshift(u, 1))), abs(sum(u)*h - 0.5_R_P) < 1.e-13_R_P
enddo
print '(A8,I10,F17.4)', 'initial', 2, tv0

contains
   subroutine residual(s, u, res)
   !< res = -(F_{i+1/2} - F_{i-1/2})/h for f(u) = u^2/2, Lax-Friedrichs splitting with alpha = 1.1 max|u|.
   integer,   intent(in)  :: s
   real(R_P), intent(in)  :: u(N)
   real(R_P), intent(out) :: res(N)
   real(R_P)              :: ug(1-NG:N+NG), fp(1-NG:N+NG), fm(1-NG:N+NG), f(0:N), alpha, rp, rm
   integer                :: i
   ug(1:N) = u ; ug(1-NG:0) = u(N+1-NG:N) ; ug(N+1:N+NG) = u(1:NG)
   alpha = 1.1_R_P*maxval(abs(u))               ! f'(u) = u; a margin above max|f'| (see the text)
   fp = 0.5_R_P*(0.5_R_P*ug**2 + alpha*ug)
   fm = 0.5_R_P*(0.5_R_P*ug**2 - alpha*ug)
   do i=0, N
      select case(s)
      case(WENO5)
         call seth_weno5_js(fp(i-2:i+3), fm(i-2:i+3), SETH_WENO_JS_EPS, 1._R_P, rp, rm)
      case(WENOZ7)
         call seth_wenoz7(fp(i-3:i+4), fm(i-3:i+4), SETH_TENO_EPS, 1, rp, rm)
      case(TENO5)
         call seth_teno5(fp(i-2:i+3), fm(i-2:i+3), SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)
      case default
         call seth_linear5(fp(i-2:i+3), fm(i-2:i+3), rp, rm)
      endselect
      f(i) = rp + rm
   enddo
   res = -(f(1:N) - f(0:N-1))/h
   endsubroutine residual

   pure integer function extrema(u)
   !< The local extrema of the periodic field u: the exact solution has one maximum and one minimum.
   real(R_P), intent(in) :: u(N)
   real(R_P)             :: d(N)
   d = cshift(u, 1) - u                         ! d(i) = u(i+1) - u(i)
   extrema = count(d*cshift(d, -1) < 0)
   endfunction extrema

   subroutine evolve(s, u, tend)
   !< SSP-RK3 to tend, CFL 0.4 on max|u|.
   integer,   intent(in)    :: s
   real(R_P), intent(inout) :: u(N)
   real(R_P), intent(in)    :: tend
   real(R_P)                :: u1(N), u2(N), res(N), t, dt
   t = 0
   do while (t < tend)
      dt = min(0.4_R_P*h/maxval(abs(u)), tend - t)
      call residual(s, u, res)  ; u1 = u + dt*res
      call residual(s, u1, res) ; u2 = 0.75_R_P*u + 0.25_R_P*(u1 + dt*res)
      call residual(s, u2, res) ; u  = u/3 + 2._R_P/3*(u2 + dt*res)
      t = t + dt
   enddo
   endsubroutine evolve
endprogram t06_burgers
