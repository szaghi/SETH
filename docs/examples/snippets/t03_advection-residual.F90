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
