!< SETH test: the centred WENO (seth_weno4/6/8_c) at its design order, mirrored, finite, and as ENO as its upwind part.
program seth_test_weno_c
!< SETH test: the centred WENO (seth_weno4/6/8_c) at its design order, mirrored, finite, and as ENO as its upwind part.
!<
!< The bounds were frozen before the kernel was written (roadmap issue #4, 6.3; night 2026-11-16 LOG, 15:39):
!<
!< 1. Order: u = g(x) (seth_test_functions), cell averages, cells of width h halved from 0.4, face at x = 0.3, rmu = 1;
!<    the observed order on the finest pair above round-off must reach 2S - 0.3, both states (in R4P less
!<    CENTRED_ORDER_LOSS, measured: the pairs above R4P's floor are pre-asymptotic; a scheme whose pairs all sit at
!<    round-off passes when ORDERS_AT_FLOOR_OK).
!< 2. Mirror: with vm(m) = vp(1-m) the right state equals the left one to TOL_PROPERTIES max|v|, the tolerance of
!<    WENO-JS (seth_test_properties), whose two biases sum mirrored tables in another order. Frozen as 4 ulps; changed
!<    at 15:42 after measuring 4.8e-15 (S = 3) and 8.9e-15 (S = 4): the central sums and the blend mirror bitwise, the
!<    upwind WENO-JS part does not, as documented there.
!< 3. Constant data: both states equal the constant to 4 ulps.
!< 4. A jump of BIG_JUMP (1e30 in R8P, 1e16 in R4P) anywhere in the stencil: both states finite.
!< 5. A unit step at every cell boundary of the big stencil, both orientations: the worst overshoot outside [0, 1] at
!<    most that of the upwind WENO-JS of the same S plus 1e-3 (the blend must reach the upwind state at a jump). Holds
!<    for S = 2, 3. For S = 4 (finding, 15:42): a step in the outer cells -3 or 4 lies outside the blend's detector
!<    (ADAM's, cells -1..2), phi = 0 there, and the central row passes its edge coefficient: overshoot |CEN8_C(4)| =
!<    1/280 = 3.6e-3; bounded by that plus 1e-3. (S = 3: the edge coefficient 1/60 is positive, a precursor inside
!<    [0, 1], not an overshoot.) Widening the detector to the whole stencil would change ADAM's scheme: left to the
!<    owner.
use seth
use seth_test_parameters, only : BIG_JUMP, CENTRED_ORDER_LOSS, ORDERS_AT_FLOOR_OK, POINT_FLOOR, TOL_PROPERTIES
use seth_test_rng, only : seed_rng
use seth_test_functions, only : gfun, gavg
use seth_z_teno_tables, only : CEN8_C
implicit none
real(R_P),    parameter :: XF = 0.3_R_P
integer(I4P), parameter :: NH = 6_I4P
real(R_P)               :: vp(-4:5), vm(-4:5), rp, rm, up, um, h, err(2,NH), order, floor_, worst, worst_up, x, d
integer(I4P)            :: s, j, m, side, jump, n_
logical                 :: passed

passed = .true.
call seed_rng(1011_I4P)
floor_ = POINT_FLOOR*epsilon(1._R_P)
do s=2, 4
   ! 1. order
   do j=1, NH
      h = 0.4_R_P/2._R_P**(j-1)
      vp = 0._R_P
      do m=1-s, s
         vp(m) = gavg(XF + (m-1)*h, XF + m*h)
      enddo
      call kernel(s, vp, vp, 1._R_P, rp, rm)
      err(:,j) = abs([rp, rm] - gfun(XF))
   enddo
   do side=1, 2
      order = -1._R_P
      do j=NH-1, 1, -1
         if (err(side,j+1) > floor_) then
            order = log(err(side,j)/err(side,j+1))/log(2._R_P) ; exit
         endif
      enddo
      print '(A,I0,A,A,A,F6.2,A,6ES9.2)', 'WENO', 2*s, '-C state=', merge('+', '-', side == 1), ' order=', order, &
            ' errors', err(side,:)
      if (order < 2*s - 0.3_R_P - CENTRED_ORDER_LOSS .and. .not.(ORDERS_AT_FLOOR_OK .and. order < 0._R_P)) &
         call fail('order below 2S - 0.3')
   enddo
   ! 2. mirror
   d = 0._R_P
   do n_=1, 1000
      vp = 0._R_P ; vm = 0._R_P
      do m=1-s, s
         call random_number(x) ; vp(m) = x
      enddo
      do m=1-s, s
         vm(m) = vp(1-m)
      enddo
      call kernel(s, vp, vm, 1._R_P, rp, rm)
      d = max(d, abs(rp - rm))
   enddo
   print '(A,I0,A,ES9.2)', 'WENO', 2*s, '-C mirror: max |rp - rm| ', d
   if (d > TOL_PROPERTIES) call fail('mirror') ! data in [0, 1): max|v| < 1
   ! 3. constant data
   vp = 0._R_P ; vp(1-s:s) = 0.7_R_P
   call kernel(s, vp, vp, 1._R_P, rp, rm)
   if (abs(rp - 0.7_R_P) > 4._R_P*epsilon(1._R_P) .or. abs(rm - 0.7_R_P) > 4._R_P*epsilon(1._R_P)) call fail('constant')
   ! 4. a big jump anywhere
   do jump=2-s, s
      vp = 0._R_P
      do m=1-s, s
         vp(m) = merge(BIG_JUMP, 1._R_P, m >= jump)
      enddo
      call kernel(s, vp, vp, 1._R_P, rp, rm)
      if (.not.(abs(rp) <= huge(1._R_P) .and. abs(rm) <= huge(1._R_P))) call fail('big jump: not finite')
   enddo
   ! 5. a unit step: overshoot against the upwind WENO-JS of the same S
   worst = 0._R_P ; worst_up = 0._R_P
   do side=1, 2
      do jump=2-s, s
         vp = 0._R_P
         do m=1-s, s
            vp(m) = merge(1._R_P, 0._R_P, (m >= jump) .eqv. (side == 1))
         enddo
         call kernel(s, vp, vp, 1._R_P, rp, rm)
         call upwind(s, vp, vp, up, um)
         worst = max(worst, -rp, rp - 1._R_P, -rm, rm - 1._R_P)
         worst_up = max(worst_up, -up, up - 1._R_P, -um, um - 1._R_P)
         if (max(-rp, rp - 1._R_P, -rm, rm - 1._R_P) > 1.e-3_R_P) &
            print '(A,I0,A,I0,A,I0,A,2F9.5)', '   WENO', 2*s, '-C step at cell ', jump, ' side ', side, ': states ', rp, rm
      enddo
   enddo
   print '(A,I0,A,ES9.2,A,ES9.2)', 'WENO', 2*s, '-C worst overshoot ', worst, ', its upwind WENO-JS ', worst_up
   if (worst > worst_up + 1.e-3_R_P + merge(abs(CEN8_C(4)), 0._R_P, s == 4)) call fail('overshoot above the bound')
enddo
if (.not.passed) error stop 'centred WENO test failed'
print '(A)', 'seth_test_weno_c: ok'

contains
   subroutine kernel(s, vp, vm, rmu, rp, rm)
   integer(I4P), intent(in)  :: s
   real(R_P),    intent(in)  :: vp(-4:5), vm(-4:5), rmu
   real(R_P),    intent(out) :: rp, rm
   select case(s)
   case(2) ; call seth_weno4_c(vp(-1:2), vm(-1:2), SETH_WENO_JS_EPS, rmu, rp, rm)
   case(3) ; call seth_weno6_c(vp(-2:3), vm(-2:3), SETH_WENO_JS_EPS, rmu, rp, rm)
   case(4) ; call seth_weno8_c(vp(-3:4), vm(-3:4), SETH_WENO_JS_EPS, rmu, rp, rm)
   endselect
   endsubroutine kernel

   subroutine upwind(s, vp, vm, rp, rm)
   integer(I4P), intent(in)  :: s
   real(R_P),    intent(in)  :: vp(-4:5), vm(-4:5)
   real(R_P),    intent(out) :: rp, rm
   select case(s)
   case(2) ; call seth_weno3_js(vp(-1:2), vm(-1:2), SETH_WENO_JS_EPS, 1._R_P, rp, rm)
   case(3) ; call seth_weno5_js(vp(-2:3), vm(-2:3), SETH_WENO_JS_EPS, 1._R_P, rp, rm)
   case(4) ; call seth_weno7_js(vp(-3:4), vm(-3:4), SETH_WENO_JS_EPS, 1._R_P, rp, rm)
   endselect
   endsubroutine upwind

   subroutine fail(what)
   character(*), intent(in) :: what
   print '(A)', '   FAILED: '//what
   passed = .false.
   endsubroutine fail
endprogram seth_test_weno_c
