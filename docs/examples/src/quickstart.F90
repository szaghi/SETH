!< Quick start: the two states at every face of a step, from WENO5-JS and from its linear background.
program quickstart
!run quickstart quickstart
!region use
use seth
implicit none
!endregion use
integer,   parameter :: N = 10           ! cells 1..N on [0, 1]
real(R_P)            :: u(-2:N+3)        ! cell averages, 3 ghost cells per side
real(R_P)            :: rp, rm, lp, lm
integer              :: i
! a step: 1 left of x = 0.5, 0 right of it; the ghost cells copy the boundary cells
u = 0._R_P
u(-2:N/2) = 1._R_P
print '(A)', ' face x   WENO5 left  WENO5 right  linear5 left  linear5 right'
do i=0, N                                ! the face x_{i+1/2} between cells i and i+1
   !region face
   ! vp and vm are the six cells around the face, cell 0 the one left of it: u(i-2:i+3)
   call seth_weno5_js(u(i-2:i+3), u(i-2:i+3), SETH_WENO_JS_EPS, 1._R_P, rp, rm)
   call seth_linear5(u(i-2:i+3), u(i-2:i+3), lp, lm)
   !endregion face
   print '(F7.2,4F13.4)', i/real(N, R_P), shown([rp, rm, lp, lm])
enddo
contains
   pure function shown(v) result(r)
   !< v rounded to 4 decimals for printing, a negative zero printed as 0 (WENO's states differ from 0 and 1 by 1e-19).
   real(R_P), intent(in) :: v(:)
   real(R_P)             :: r(size(v))
   r = anint(v*1e4_R_P)/1e4_R_P + 0._R_P
   endfunction shown
endprogram quickstart
