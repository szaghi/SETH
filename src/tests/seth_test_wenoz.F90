!< SETH test: WENO-Z3/5/7/9 converge at their design order, stay bounded across a jump, and mirror exactly.
program seth_test_wenoz
!< SETH test: WENO-Z3/5/7/9 converge at their design order, stay bounded across a jump, and mirror exactly.
!<
!< Every check runs on the reconstruction kernels and on the interpolation ones (`_interp`).
!<
!< 1. Order: u = g(x), the generic function of seth_test_functions (not exp: its averages are a constant multiple of
!<    its point values), cell averages (reconstruction) or point values at the cell centres (interpolation), cells of
!<    width h halved from 0.4, face at x = 0.3, p = 1 and p = 2; the observed order on the finest pair above round-off
!<    must reach 2S-1 - 0.5.
!< 2. ENO: a unit step at every cell boundary of the stencil, both orientations; the states lie in [0, 1] to 1e-6.
!< 3. Mirror: with vm(m) = vp(1-m) the right state equals the left one bitwise (the same arithmetic on the mirror).
use seth
use seth_test_parameters, only : EPS_ZT, ORDERS_AT_FLOOR_OK, POINT_FLOOR, POINT_ORDER_LOSS
use seth_test_rng, only : seed_rng
use seth_test_functions, only : gfun, gavg
use seth_vv_cases, only : vv_report, vv_open
implicit none
real(R_P), parameter :: EPS = EPS_ZT !< WENO-Z regulariser of the tests.
real(R_P), parameter :: XF = 0.3_R_P
integer(I4P), parameter :: NH = 6_I4P
real(R_P)            :: vp(-4:5), vm(-4:5), rp, rm, h, err(2,NH), order, floor_, worst, x
integer(I4P)         :: s, o, p, j, m, side, jump, n_, t
logical              :: passed, interp
character(7)         :: kind
integer              :: ur

passed = .true.
call seed_rng(1001_I4P)
if (vv_report()) call vv_open('orders-wenoz.dat', ur, 'kind order p state observed-order')
floor_ = POINT_FLOOR * epsilon(1._R_P) ! R8P: 100 eps |g|, |g| < 2 near XF
do t=1, 2
   interp = t == 2
   kind = merge('interp ', 'recon  ', interp)
   do s=2, 5 ! WENO-Z(2s-1)
      o = 2*s - 1
      do p=1, 2
         do j=1, NH
            h = merge(0.8_R_P, 0.4_R_P, s == 5) / 2._R_P**(j-1)
            vp = 0._R_P
            do m=1-s, s
               if (interp) then
                  vp(m) = gfun(XF + (m - 0.5_R_P)*h)
               else
                  vp(m) = gavg(XF + (m-1)*h, XF + m*h)
               endif
            enddo
            call kernel(s, p, interp, vp, vp, rp, rm)
            err(:,j) = abs([rp, rm] - gfun(XF))
         enddo
         do side=1, 2
            order = -1._R_P
            do j=NH-1, 1, -1
               if (err(side,j+1) > floor_) then
                  order = log(err(side,j)/err(side,j+1)) / log(2._R_P) ; exit
               endif
            enddo
            print '(A,A,I0,A,I0,A,A,A,F6.2)', kind, 'WENO-Z', o, ' p=', p, ' state=', merge('+','-',side==1), &
                  ' order=', order
            if (vv_report()) write(ur, '(A,1X,I2,1X,I1,1X,A1,F8.3)') trim(kind), o, p, merge('+','-',side==1), order
            if (order < o - 0.5_R_P - POINT_ORDER_LOSS .and. .not.(ORDERS_AT_FLOOR_OK .and. order < 0._R_P)) passed = .false.
         enddo
      enddo
      ! ENO
      worst = 0._R_P
      do side=1, 2
         do jump=1-s, s
            vp = 0._R_P
            do m=1-s, s
               vp(m) = merge(1._R_P, 0._R_P, (m >= jump) .eqv. (side == 1))
            enddo
            call kernel(s, 1, interp, vp, vp, rp, rm)
            worst = max(worst, -rp, rp - 1._R_P, -rm, rm - 1._R_P)
         enddo
      enddo
      print '(A,A,I0,A,ES9.2)', kind, 'WENO-Z', o, ' worst overshoot=', worst
      if (worst > 1.e-6_R_P) passed = .false.
      ! mirror
      do n_=1, 1000
         vp = 0._R_P ; vm = 0._R_P
         do m=1-s, s
            call random_number(x) ; vp(m) = x
         enddo
         do m=1-s, s
            vm(m) = vp(1-m)
         enddo
         call kernel(s, 1, interp, vp, vm, rp, rm)
         if (rp /= rm) then
            print '(A,A,I0,A,2ES24.16)', kind, 'WENO-Z', o, ' mirror mismatch: ', rp, rm ; passed = .false. ; exit
         endif
      enddo
   enddo
enddo
if (vv_report()) close(ur)
if (.not.passed) error stop 'WENO-Z test failed'
print '(A)', 'seth_test_wenoz: ok'

contains
   subroutine kernel(s, p, interp, vp, vm, rp, rm)
   integer(I4P), intent(in)  :: s, p
   logical,      intent(in)  :: interp
   real(R_P),    intent(in)  :: vp(-4:5), vm(-4:5)
   real(R_P),    intent(out) :: rp, rm
   if (interp) then
      select case(s)
      case(2) ; call seth_wenoz3_interp(vp(-1:2), vm(-1:2), EPS, p, rp, rm)
      case(3) ; call seth_wenoz5_interp(vp(-2:3), vm(-2:3), EPS, p, rp, rm)
      case(4) ; call seth_wenoz7_interp(vp(-3:4), vm(-3:4), EPS, p, rp, rm)
      case(5) ; call seth_wenoz9_interp(vp(-4:5), vm(-4:5), EPS, p, rp, rm)
      endselect
   else
      select case(s)
      case(2) ; call seth_wenoz3(vp(-1:2), vm(-1:2), EPS, p, rp, rm)
      case(3) ; call seth_wenoz5(vp(-2:3), vm(-2:3), EPS, p, rp, rm)
      case(4) ; call seth_wenoz7(vp(-3:4), vm(-3:4), EPS, p, rp, rm)
      case(5) ; call seth_wenoz9(vp(-4:5), vm(-4:5), EPS, p, rp, rm)
      endselect
   endif
   endsubroutine kernel
endprogram seth_test_wenoz
