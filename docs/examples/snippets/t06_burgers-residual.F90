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
