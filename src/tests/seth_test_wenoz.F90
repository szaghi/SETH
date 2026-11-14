!< SETH test: WENO-Z3/5/7 converge at their design order, stay bounded across a jump, and mirror exactly.
program seth_test_wenoz
!< SETH test: WENO-Z3/5/7 converge at their design order, stay bounded across a jump, and mirror exactly.
!<
!< 1. Order: u = exp(x), cell averages of width h halved from 0.4, face at x = 0.3, p = 1 and p = 2; the observed
!<    order on the finest pair above round-off must reach 2S-1 - 0.5.
!< 2. ENO: a unit step at every cell boundary of the stencil, both orientations; the states lie in [0, 1] to 1e-6.
!< 3. Mirror: with vm(m) = vp(1-m) the right state equals the left one bitwise (the same arithmetic on the mirror).
use seth
implicit none
real(R8P), parameter :: EPS = 1.e-40_R8P !< WENO-Z regulariser of the tests.
real(R8P), parameter :: XF = 0.3_R8P
integer(I4P), parameter :: NH = 6_I4P
real(R8P)            :: vp(-3:4), vm(-3:4), rp, rm, h, err(2,NH), order, floor_, worst, x
integer(I4P)         :: s, o, p, j, m, side, jump, n_
logical              :: passed

passed = .true.
floor_ = 1.e3_R8P * epsilon(1._R8P) * exp(XF)
do s=2, 4 ! WENO-Z(2s-1)
   o = 2*s - 1
   do p=1, 2
      do j=1, NH
         h = 0.4_R8P / 2._R8P**(j-1)
         vp = 0._R8P
         do m=1-s, s
            vp(m) = (exp(XF + m*h) - exp(XF + (m-1)*h)) / h
         enddo
         call kernel(s, p, vp, vp, rp, rm)
         err(:,j) = abs([rp, rm] - exp(XF))
      enddo
      do side=1, 2
         order = -1._R8P
         do j=NH-1, 1, -1
            if (err(side,j+1) > floor_) then
               order = log(err(side,j)/err(side,j+1)) / log(2._R8P) ; exit
            endif
         enddo
         print '(A,I0,A,I0,A,A,A,F6.2)', 'WENO-Z', o, ' p=', p, ' state=', merge('+','-',side==1), ' order=', order
         if (order < o - 0.5_R8P) passed = .false.
      enddo
   enddo
   ! ENO
   worst = 0._R8P
   do side=1, 2
      do jump=1-s, s
         vp = 0._R8P
         do m=1-s, s
            vp(m) = merge(1._R8P, 0._R8P, (m >= jump) .eqv. (side == 1))
         enddo
         call kernel(s, 1, vp, vp, rp, rm)
         worst = max(worst, -rp, rp - 1._R8P, -rm, rm - 1._R8P)
      enddo
   enddo
   print '(A,I0,A,ES9.2)', 'WENO-Z', o, ' worst overshoot=', worst
   if (worst > 1.e-6_R8P) passed = .false.
   ! mirror
   do n_=1, 1000
      vp = 0._R8P ; vm = 0._R8P
      do m=1-s, s
         call random_number(x) ; vp(m) = x
      enddo
      do m=1-s, s
         vm(m) = vp(1-m)
      enddo
      call kernel(s, 1, vp, vm, rp, rm)
      if (rp /= rm) then
         print '(A,I0,A,2ES24.16)', 'WENO-Z', o, ' mirror mismatch: ', rp, rm ; passed = .false. ; exit
      endif
   enddo
enddo
if (.not.passed) error stop 'WENO-Z test failed'
print '(A)', 'seth_test_wenoz: ok'

contains
   subroutine kernel(s, p, vp, vm, rp, rm)
   integer(I4P), intent(in)  :: s, p
   real(R8P),    intent(in)  :: vp(-3:4), vm(-3:4)
   real(R8P),    intent(out) :: rp, rm
   select case(s)
   case(2) ; call seth_wenoz3(vp(-1:2), vm(-1:2), EPS, p, rp, rm)
   case(3) ; call seth_wenoz5(vp(-2:3), vm(-2:3), EPS, p, rp, rm)
   case(4) ; call seth_wenoz7(vp(-3:4), vm(-3:4), EPS, p, rp, rm)
   endselect
   endsubroutine kernel
endprogram seth_test_wenoz
