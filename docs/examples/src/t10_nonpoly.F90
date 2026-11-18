!< Tutorial 10: non-polynomial TENO. The DR and THINC overrides keep an advected jump sharp.
program t10_nonpoly
!run t10_nonpoly t10_nonpoly
use seth
implicit none
integer,   parameter :: N = 200, NG = 6         ! the composite kernels read cells -5..6: six ghost cells
real(R_P), parameter :: A = 1._R_P, CFL = 0.4_R_P
integer,   parameter :: NS = 4
character(12)        :: label(NS) = [character(12) :: 'TENO5', 'TENO5-DR', 'TENO5-THINC', 'TENO5-RBF+DR']
real(R_P)            :: u0(N), u(N), h, x
integer              :: s, i

h = 1._R_P/N
do i=1, N                                       ! a square wave on [0.25, 0.75]
   x = (i - 0.5_R_P)*h
   u0(i) = merge(1._R_P, 0._R_P, x > 0.25_R_P .and. x < 0.75_R_P)
enddo
print '(A)', 'one period of linear advection of a square wave, 200 cells'
print '(A)', 'scheme          L1 error   in [0, 1]  cells per jump'
do s=1, NS
   u = u0
   call advect(s, u)
   print '(A12,ES11.2,L10,I14)', label(s), sum(abs(u - u0))*h, &
                                 minval(u) > -1.e-14_R_P .and. maxval(u) < 1 + 1.e-14_R_P, width(u)
enddo

contains
   pure integer function width(u)
   !< The cells of the left jump between 5% and 95% of its height.
   real(R_P), intent(in) :: u(N)
   width = count(u(1:N/2) > 0.05_R_P .and. u(1:N/2) < 0.95_R_P)
   endfunction width

   subroutine residual(s, u, res)
   integer,   intent(in)  :: s
   real(R_P), intent(in)  :: u(N)
   real(R_P), intent(out) :: res(N)
   real(R_P)              :: ug(1-NG:N+NG), fp(1-NG:N+NG), fm(1-NG:N+NG), f(0:N), rp, rm
   integer                :: i
   ug(1:N) = u ; ug(1-NG:0) = u(N+1-NG:N) ; ug(N+1:N+NG) = u(1:NG)
   fp = 0.5_R_P*(A*ug + abs(A)*ug)
   fm = 0.5_R_P*(A*ug - abs(A)*ug)
   do i=0, N
      !region kernels
      select case(s)
      case(1)
         call seth_teno5(fp(i-2:i+3), fm(i-2:i+3), SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)
      case(2)                                   ! twelve cells, -5..6
         call seth_teno5_dr(fp(i-5:i+6), fm(i-5:i+6), SETH_TENO_EPS, SETH_TENO5_CT, SETH_DR_T1, SETH_DR_T2, &
                            1._R_P, rp, rm)
      case(3)
         call seth_teno5_thinc(fp(i-5:i+6), fm(i-5:i+6), SETH_TENO_EPS, SETH_TENO5_CT, SETH_THINC_BETA, 1._R_P, &
                               rp, rm)
      case(4)                                   ! the RBF shape parameters in the grid's units: s = 3 h
         call seth_teno5_rbf_dr(fp(i-5:i+6), fm(i-5:i+6), SETH_TENO_EPS, SETH_TENO5_CT, SETH_RBF_CR, &
                                SETH_RBF_ALPHA1*h, SETH_RBF_ALPHA2*h, SETH_DR_T1, SETH_DR_T2, 1._R_P, rp, rm)
      endselect
      !endregion kernels
      f(i) = rp + rm
   enddo
   res = -(f(1:N) - f(0:N-1))/h
   endsubroutine residual

   subroutine advect(s, u)
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
endprogram t10_nonpoly
