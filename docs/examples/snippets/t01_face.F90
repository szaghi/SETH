!< Tutorial 1: one face by hand. The two states of WENO5-JS at one face, on smooth data and across a jump.
program t01_face
use seth
implicit none
real(R_P), parameter :: PI = acos(-1._R_P)
real(R_P), parameter :: H = 0.02_R_P            ! cell width
real(R_P), parameter :: XF = 0.1_R_P            ! the face
real(R_P)            :: u(-2:3), s(-2:3), rp, rm, lp, lm, wrp, wrm
integer              :: m

! cell averages of sin(2 pi x) on the six cells around the face: cell m is [XF + (m-1) H, XF + m H]
do m=-2, 3
   u(m) = (cos(2*PI*(XF + (m - 1)*H)) - cos(2*PI*(XF + m*H)))/(2*PI*H)
enddo
s = u                                           ! kept for the mirror, below
call seth_weno5_js_wr(u, u, SETH_WENO_JS_EPS, 1._R_P, rp, rm, wrp, wrm)
call seth_linear5(u, u, lp, lm)
print '(A)', 'smooth data, sin(2 pi x)'
print '(A,F10.6)', '  exact value at the face  ', sin(2*PI*XF)
print '(A,2F10.6)', '  WENO5-JS   rp, rm        ', rp, rm
print '(A,2F10.6)', '  linear5    rp, rm        ', lp, lm
print '(A,2F10.6)', '  weights ratios wrp, wrm  ', wrp, wrm

! a step between cells 1 and 2: the left state's stencil (cells -2..2) crosses it, the right state's (-1..3) too
u(-2:1) = 1._R_P
u(2:3)  = 0._R_P
call seth_weno5_js_wr(u, u, SETH_WENO_JS_EPS, 1._R_P, rp, rm, wrp, wrm)
call seth_linear5(u, u, lp, lm)
print '(A)', 'a step between cells 1 and 2'
print '(A,2F10.6)', '  WENO5-JS   rp, rm        ', rp, rm
print '(A,2F10.6)', '  linear5    rp, rm        ', lp, lm
print '(A,2F10.6)', '  weights ratios wrp, wrm  ', wrp, wrm

! vm(m) = vp(1-m): the right state from the mirrored data is the left state from the data
call seth_weno5_js(s, s(3:-2:-1), SETH_WENO_JS_EPS, 1._R_P, rp, rm)
print '(A)', 'the smooth data and its mirror'
print '(A,2F10.6)', '  rp, rm                   ', rp, rm
print '(A,L2)', '  equal to round-off       ', abs(rp - rm) <= 4*epsilon(1._R_P)
endprogram t01_face
