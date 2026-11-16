!< SETH test: across a jump the WENO-JS states stay within the data's bounds (no Gibbs oscillation).
program seth_test_weno_js_eno
!< SETH test: across a jump the WENO-JS states stay within the data's bounds (no Gibbs oscillation).
!<
!< A unit step placed at every cell boundary of the stencil, both orientations; the states must lie in [0, 1] up to
!< 1e-6 (the Jiang-Shu weights of the candidates crossing the jump are O(eps^wexp), not zero).
use seth
implicit none
real(R_P), parameter   :: TOL = 1.e-6_R_P
real(R_P), allocatable :: vp(:)
real(R_P)              :: rp, rm, worst
integer(I4P)           :: S, jump, m, o
logical                :: passed

passed = .true.
do S=2, 5
   if (allocated(vp)) deallocate(vp)
   allocate(vp(1-S:S))
   worst = 0._R_P
   do o=1, 2
      do jump=1-S, S ! the value changes between cell jump-1 and cell jump
         do m=1-S, S
            vp(m) = merge(1._R_P, 0._R_P, (m >= jump) .eqv. (o == 1))
         enddo
         select case(S)
         case(2) ; call seth_weno3_js(vp, vp, SETH_WENO_JS_EPS, 1._R_P, rp, rm)
         case(3) ; call seth_weno5_js(vp, vp, SETH_WENO_JS_EPS, 1._R_P, rp, rm)
         case(4) ; call seth_weno7_js(vp, vp, SETH_WENO_JS_EPS, 1._R_P, rp, rm)
         case(5) ; call seth_weno9_js(vp, vp, SETH_WENO_JS_EPS, 1._R_P, rp, rm)
         endselect
         worst = max(worst, -rp, rp - 1._R_P, -rm, rm - 1._R_P)
      enddo
   enddo
   print '(A,I0,A,ES9.2)', 'S=', S, ' worst overshoot=', worst
   if (worst > TOL) passed = .false.
enddo
if (.not.passed) error stop 'a WENO-JS kernel overshoots across a jump'
print '(A)', 'seth_test_weno_js_eno: ok'
endprogram seth_test_weno_js_eno
