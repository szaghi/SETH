!< SETH test: the monotonicity-preserving bound of Suresh & Huynh (seth_mp_bound) on WENO states.
program seth_test_mp
!< SETH test: the monotonicity-preserving bound of Suresh & Huynh (seth_mp_bound) on WENO states.
!<
!< 1. Monotone data (random increasing and decreasing stencils): the bounded WENO5-JS states lie between v_0 and v_1,
!<    the values on the two sides of the face (to 4 eps of the data).
!< 2. A unit step at every cell boundary, both orientations, reconstructed by WENO9-JS, WENO-Z9 (p = 1, the scheme
!<    that overshoots) and TENO6: the bounded states lie in [0, 1] to 1e-12.
!< 3. Well-resolved smooth monotone data: the bound leaves the state bitwise unchanged.
!< 4. f = sin^4(pi x) on [0, 1], periodic, cell averages, N = 20 .. 640: the L1 and Linf errors of the left state with
!<    and without the bound, their orders printed; asserted: the bound changes no face (its curvature terms leave the
!<    accurate states at the extrema alone), every face it would change lies within 2 h of an
!<    extremum of f (x = 0 and 1/2), and the L1 order with the bound at the finest pair is at least 4.5. Balsara & Shu
!<    (2000, Table V) report order loss at the flat extrema of this function; this measures how much (WENO5-JS: the
!<    bound never acts on these grids, the orders are those of the plain scheme).
!<    The grids are the first MP_NG of N = 20 .. 1280 (seth_test_parameters: all of them in R8P; in R4P up to N = 80,
!<    where the cell averages, differences of the primitive, are still accurate).
!< 5. Scaling by 2^k and the mirror: bitwise.
use seth
use seth_test_parameters, only : EPS_ZT, MP_NG, TOL_BOUNDS
use seth_test_rng, only : seed_rng
use seth_test_functions, only : gavg
implicit none
real(R_P), parameter :: PI = acos(-1._R_P)
real(R_P)            :: vp(-4:5), vm(-4:5), rp, rm, rp0, rm0, x, worst, sc, s1, s2
integer(I4P)         :: n_, m, jump, side, k, j, nn, i, nchanged = 0
integer(I4P), parameter :: NG = 7_I4P
integer(I4P)         :: ngrid(NG) = [20, 40, 80, 160, 320, 640, 1280]
real(R_P)            :: e1(NG,2), ei(NG,2)
logical              :: passed

passed = .true.
call seed_rng(1009_I4P)
! 1. monotone data
worst = 0._R_P
do n_=1, 20000
   x = 0._R_P
   do m=-4, 5
      call random_number(sc) ; x = x + sc**3 ; vp(m) = x
   enddo
   if (mod(n_, 2) == 0) vp = -vp
   vm = vp
   call seth_weno5_js(vp(-2:3), vm(-2:3), SETH_WENO_JS_EPS, 1._R_P, rp, rm)
   call seth_mp_bound(vp(-2:3), vm(-2:3), SETH_MP_ALPHA, rp, rm)
   sc = 4._R_P*epsilon(1._R_P)*maxval(abs(vp))
   worst = max(worst, min(vp(0), vp(1)) - rp, rp - max(vp(0), vp(1)), min(vp(0), vp(1)) - rm, rm - max(vp(0), vp(1)))
   if (rp < min(vp(0), vp(1)) - sc .or. rp > max(vp(0), vp(1)) + sc .or. &
       rm < min(vp(0), vp(1)) - sc .or. rm > max(vp(0), vp(1)) + sc) then
      print '(A,4ES12.4)', 'monotone data: bounded state outside [v_0, v_1]: ', vp(0), vp(1), rp, rm ; passed = .false.
      exit
   endif
enddo
print '(A,ES10.2)', '1. monotone data, worst excursion beyond [v_0, v_1] (negative: inside): ', worst
! 2. steps
worst = 0._R_P
do side=1, 2
   do jump=-3, 5
      do m=-4, 5
         vp(m) = merge(1._R_P, 0._R_P, (m >= jump) .eqv. (side == 1))
      enddo
      vm = vp
      do k=1, 3
         select case(k)
         case(1) ; call seth_weno9_js(vp, vm, SETH_WENO_JS_EPS, 1._R_P, rp, rm)
         case(2) ; call seth_wenoz9(vp, vm, EPS_ZT, 1_I4P, rp, rm)
         case(3) ; call seth_teno6(vp(-2:3), vm(-2:3), SETH_TENO_EPS, SETH_TENO6_CT, rp, rm)
         endselect
         call seth_mp_bound(vp(-2:3), vm(-2:3), SETH_MP_ALPHA, rp, rm)
         worst = max(worst, -rp, rp - 1._R_P, -rm, rm - 1._R_P)
      enddo
   enddo
enddo
print '(A,ES10.2)', '2. steps (WENO9-JS, WENO-Z9 p = 1, TENO6), worst excursion beyond [0, 1]: ', worst
if (worst > TOL_BOUNDS) passed = .false.
! 3. well-resolved smooth monotone data: unchanged
k = 0
do n_=1, 2000
   call random_number(x)
   do m=-2, 3
      vp(m) = gavg(x + (m-1)*0.01_R_P, x + m*0.01_R_P)
   enddo
   if ((vp(1) - vp(0))*(vp(0) - vp(-1)) <= 0._R_P .or. (vp(2) - vp(1))*(vp(1) - vp(0)) <= 0._R_P) cycle
   if ((vp(0) - vp(-1))*(vp(-1) - vp(-2)) <= 0._R_P .or. (vp(3) - vp(2))*(vp(2) - vp(1)) <= 0._R_P) cycle
   vm = vp
   call seth_weno5_js(vp(-2:3), vm(-2:3), SETH_WENO_JS_EPS, 1._R_P, rp0, rm0)
   rp = rp0 ; rm = rm0
   call seth_mp_bound(vp(-2:3), vm(-2:3), SETH_MP_ALPHA, rp, rm)
   if (rp /= rp0 .or. rm /= rm0) k = k + 1
enddo
print '(A,I0)', '3. smooth monotone stencils changed by the bound: ', k
if (k > 0) passed = .false.
! 4. sin^4: orders with and without the bound, and where the bound acts
do i=1, MP_NG
   nn = ngrid(i)
   e1(i,:) = 0._R_P ; ei(i,:) = 0._R_P
   do j=0, nn-1 ! face (j+1)/nn, left state from the cells j-2 .. j+3 (cell j = [j/nn, (j+1)/nn])
      do m=-2, 3
         vp(m) = sin4avg(real(j+m, R_P)/nn, real(j+m+1, R_P)/nn)
      enddo
      vm = vp
      call seth_weno5_js(vp(-2:3), vm(-2:3), SETH_WENO_JS_EPS, 1._R_P, rp0, rm0)
      rp = rp0 ; rm = rm0
      call seth_mp_bound(vp(-2:3), vm(-2:3), SETH_MP_ALPHA, rp, rm)
      x = real(j+1, R_P)/nn
      e1(i,1) = e1(i,1) + abs(rp0 - sin(PI*x)**4)/nn ; ei(i,1) = max(ei(i,1), abs(rp0 - sin(PI*x)**4))
      e1(i,2) = e1(i,2) + abs(rp - sin(PI*x)**4)/nn  ; ei(i,2) = max(ei(i,2), abs(rp - sin(PI*x)**4))
      if (rp /= rp0) nchanged = nchanged + 1
      if (rp /= rp0 .and. min(abs(x - 0.5_R_P), x, 1._R_P - x) > 2._R_P/nn) then
         print '(A,I5,F9.5)', '   the bound changed a state away from the extrema: N, x = ', nn, x ; passed = .false.
      endif
   enddo
enddo
print '(A)', '4. sin^4(pi x), WENO5-JS left state:      N     L1 plain  order    L1 bound  order  Linf plain  Linf bound'
do i=1, MP_NG
   if (i == 1) then
      print '(A,I6,2(ES12.3,A7),2ES12.3)', '   ', ngrid(i), e1(i,1), '', e1(i,2), '', ei(i,1), ei(i,2)
   else
      print '(A,I6,2(ES12.3,F7.2),2ES12.3)', '   ', ngrid(i), e1(i,1), log(e1(i-1,1)/e1(i,1))/log(2._R_P), e1(i,2), &
            log(e1(i-1,2)/e1(i,2))/log(2._R_P), ei(i,1), ei(i,2)
   endif
enddo
print '(A,I0)', '   faces the bound changed, all grids: ', nchanged
if (nchanged > 0) then ! the curvature terms must leave accurate states near the extrema alone (Suresh & Huynh)
   print '(A)', '   ^ the bound changed an accurate state' ; passed = .false.
endif
if (log(e1(MP_NG-1,2)/e1(MP_NG,2))/log(2._R_P) < 4.5_R_P) then
   print '(A)', '   ^ L1 order with the bound below 4.5' ; passed = .false.
endif
! 5. scaling and mirror, bitwise
do n_=1, 2000
   do m=-2, 3
      call random_number(x) ; vp(m) = x
   enddo
   do m=-2, 3
      vm(m) = vp(1-m)
   enddo
   call seth_weno5_js(vp(-2:3), vm(-2:3), SETH_WENO_JS_EPS, 1._R_P, rp, rm)
   rm = rp ! WENO-JS mirrors to round-off only: give the bound the same state on both sides
   s1 = rp ; s2 = rm
   call seth_mp_bound(vp(-2:3), vm(-2:3), SETH_MP_ALPHA, rp, rm)
   if (rp /= rm) then ; print '(A)', '5. mirror mismatch' ; passed = .false. ; exit ; endif
   sc = 2._R_P**(mod(n_, 41) - 20)
   s1 = sc*s1 ; s2 = sc*s2
   call seth_mp_bound(sc*vp(-2:3), sc*vm(-2:3), SETH_MP_ALPHA, s1, s2)
   if (s1 /= sc*rp .or. s2 /= sc*rm) then ; print '(A)', '5. scaling mismatch' ; passed = .false. ; exit ; endif
enddo
if (.not.passed) error stop 'MP bound test failed'
print '(A)', 'seth_test_mp: ok'

contains
   pure function sin4avg(a, b) result(r)
   !< Average of sin^4(pi x) over [a, b]: the primitive of sin^4 is 3x/8 - sin(2 pi x)/(4 pi) + sin(4 pi x)/(32 pi).
   real(R_P), intent(in) :: a, b
   real(R_P)             :: r
   r = ((3._R_P*b/8._R_P - sin(2._R_P*PI*b)/(4._R_P*PI) + sin(4._R_P*PI*b)/(32._R_P*PI)) - &
        (3._R_P*a/8._R_P - sin(2._R_P*PI*a)/(4._R_P*PI) + sin(4._R_P*PI*a)/(32._R_P*PI)))/(b - a)
   endfunction sin4avg
endprogram seth_test_mp
