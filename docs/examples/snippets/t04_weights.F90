!< Tutorial 4: the weights. What eps and the descaler rmu do to WENO5-JS on a small jump.
program t04_weights
use seth
implicit none
real(R_P) :: step(-2:3), u(-2:3), scale, rp, rm, rmu, mean, o1
integer   :: k

! the step of tutorial 1, between cells 1 and 2: the right state's stencil crosses it
step(-2:1) = 1._R_P
step(2:3)  = 0._R_P

print '(A)', 'the step times a scale: overshoot of the right state, as a fraction of the jump'
print '(A)', '     scale     rmu = 1   rmu = descaling'
scale = 1._R_P
do k=1, 4
   u = scale*step
   call seth_weno5_js(u, u, SETH_WENO_JS_EPS, 1._R_P, rp, rm)        ! Jiang and Shu's weights
   o1 = over(rm, scale, scale)
   rmu = seth_weno_descaling(3, 1._R_P, u, u)                         ! 1/mean|v| of the stencils
   call seth_weno5_js(u, u, SETH_WENO_JS_EPS, rmu, rp, rm)
   print '(ES10.0,F12.4,F15.4)', scale, o1, over(rm, scale, scale)
   scale = scale*1.e-2_R_P
enddo

print '(A)', 'a step of 1e-7 on the offset 1: overshoot of the right state, as a fraction of the jump'
u = 1._R_P + 1.e-7_R_P*step
rmu = seth_weno_descaling(3, 1._R_P, u, u)
call seth_weno5_js(u, u, SETH_WENO_JS_EPS, rmu, rp, rm)
print '(A,F10.4)', '  descaling (mean |v|)           ', over(rm, 1 + 1.e-7_R_P, 1.e-7_R_P)
mean = sum(u)/size(u)                           ! subtract the stencil's mean, descale, add the mean back
rmu = seth_weno_descaling_ai(3, u - mean, u - mean)
call seth_weno5_js(u - mean, u - mean, SETH_WENO_JS_EPS, rmu, rp, rm)
rp = rp + mean ; rm = rm + mean
print '(A,F10.4)', '  mean subtracted, affine descaler', over(rm, 1 + 1.e-7_R_P, 1.e-7_R_P)
contains
   pure function over(r, top, jump) result(o)
   !< The overshoot of the state r above the top of the step, as a fraction of the jump (0 if r is not above it).
   real(R_P), intent(in) :: r, top, jump
   real(R_P)             :: o
   o = max(r - top, 0._R_P)/jump
   endfunction over
endprogram t04_weights
