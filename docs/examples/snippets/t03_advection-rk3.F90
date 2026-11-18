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
