!< SETH test: the linear kernels (the backgrounds of the WENO and TENO schemes, the smooth branch of a hybrid).
program seth_test_linear
!< SETH test: the linear kernels (the backgrounds of the WENO and TENO schemes, the smooth branch of a hybrid).
!<
!< 1. Order on g(x) of seth_test_functions (cell averages, point values for the interpolations; face at 0.3), both
!<    states, >= design - 0.3: 3, 5, 7, 9, 6 and 6 (TENO8's background).
!< 2. Equal to the nonlinear scheme on smooth data, where that scheme is its linear background (TENO5 = WENO-5
!<    background, TENO6, TENO8, and the TENO5/6 interpolations): to 8 eps of the data (other summation order).
!< 3. vm(m) = vp(1-m) gives rm == rp bitwise.
use seth
use seth_test_rng, only : seed_rng
use seth_test_functions, only : gfun, gavg
implicit none
real(R_P),    parameter :: XF = 0.3_R_P
integer(I4P), parameter :: NH = 6_I4P, NK = 11_I4P
character(16), parameter :: NAMES(NK) = [character(16) :: 'linear3', 'linear5', 'linear7', 'linear9', &
                                         'linear3_interp', 'linear5_interp', 'linear7_interp', 'linear9_interp', &
                                         'linear6', 'linear6_interp', 'linear8']
real(R_P),    parameter :: DESIGN(NK) = [3._R_P, 5._R_P, 7._R_P, 9._R_P, 3._R_P, 5._R_P, 7._R_P, 9._R_P, 6._R_P, &
                                         6._R_P, 6._R_P]
logical,      parameter :: INTERP(NK) = [.false., .false., .false., .false., .true., .true., .true., .true., .false., &
                                         .true., .false.]
real(R_P)               :: vp(-4:5), vm(-4:5), rp, rm, sp, sm, h, err(2,NH), order, worst, x
integer(I4P)            :: k, j, m, side, n_
logical                 :: passed

passed = .true.
call seed_rng(1011_I4P)
do k=1, NK
   ! 1. order
   do j=1, NH
      h = 0.4_R_P/2._R_P**(j-1)
      do m=-4, 5
         if (INTERP(k)) then
            vp(m) = gfun(XF + (m - 0.5_R_P)*h)
         else
            vp(m) = gavg(XF + (m-1)*h, XF + m*h)
         endif
      enddo
      call kernel(k, vp, vp, rp, rm)
      err(:,j) = abs([rp, rm] - gfun(XF))
   enddo
   do side=1, 2
      order = -1._R_P
      do j=NH-1, 1, -1
         if (err(side,j+1) > 2.e2_R_P*epsilon(1._R_P)) then
            order = log(err(side,j)/err(side,j+1))/log(2._R_P) ; exit
         endif
      enddo
      print '(A16,A,A,A,F6.2)', NAMES(k), ' state=', merge('+','-',side==1), ' order=', order
      if (order < DESIGN(k) - 0.3_R_P) passed = .false.
   enddo
   ! 3. mirror
   do n_=1, 500
      do m=-4, 5
         call random_number(x) ; vp(m) = x
      enddo
      do m=-4, 5
         vm(m) = vp(1-m)
      enddo
      call kernel(k, vp, vm, rp, rm)
      if (rp /= rm) then ; print '(A)', NAMES(k)//': mirror mismatch' ; passed = .false. ; exit ; endif
   enddo
enddo
! 2. the nonlinear scheme on smooth data is its background
worst = 0._R_P
do n_=1, 2000
   call random_number(x)
   do m=-4, 5
      vp(m) = 1._R_P + x*m + 0.1_R_P*x*m*m + 0.01_R_P*x*m**3
   enddo
   call seth_linear5(vp(-2:3), vp(-2:3), rp, rm)
   call seth_teno5(vp(-2:3), vp(-2:3), SETH_TENO_EPS, SETH_TENO5_CT, sp, sm) ; call cmp(rp, rm, sp, sm)
   call seth_linear6(vp(-2:3), vp(-2:3), rp, rm)
   call seth_teno6(vp(-2:3), vp(-2:3), SETH_TENO_EPS, SETH_TENO6_CT, sp, sm) ; call cmp(rp, rm, sp, sm)
   call seth_linear8(vp(-3:4), vp(-3:4), rp, rm)
   call seth_teno8(vp(-3:4), vp(-3:4), SETH_TENO_EPS, SETH_TENO8_CT, sp, sm) ; call cmp(rp, rm, sp, sm)
   call seth_linear5_interp(vp(-2:3), vp(-2:3), rp, rm)
   call seth_teno5_interp(vp(-2:3), vp(-2:3), SETH_TENO_EPS, SETH_TENO5_CT, sp, sm) ; call cmp(rp, rm, sp, sm)
   call seth_linear6_interp(vp(-2:3), vp(-2:3), rp, rm)
   call seth_teno6_interp(vp(-2:3), vp(-2:3), SETH_TENO_EPS, SETH_TENO6_CT, sp, sm) ; call cmp(rp, rm, sp, sm)
enddo
print '(A,ES10.2)', 'linear kernels against TENO on smooth data, worst scaled difference: ', worst
if (worst > 8._R_P*epsilon(1._R_P)) passed = .false.
if (.not.passed) error stop 'linear kernel test failed'
print '(A)', 'seth_test_linear: ok'

contains
   subroutine kernel(k, vp, vm, rp, rm)
   integer(I4P), intent(in)  :: k
   real(R_P),    intent(in)  :: vp(-4:5), vm(-4:5)
   real(R_P),    intent(out) :: rp, rm
   select case(k)
   case(1)  ; call seth_linear3(vp(-1:2), vm(-1:2), rp, rm)
   case(2)  ; call seth_linear5(vp(-2:3), vm(-2:3), rp, rm)
   case(3)  ; call seth_linear7(vp(-3:4), vm(-3:4), rp, rm)
   case(4)  ; call seth_linear9(vp(-4:5), vm(-4:5), rp, rm)
   case(5)  ; call seth_linear3_interp(vp(-1:2), vm(-1:2), rp, rm)
   case(6)  ; call seth_linear5_interp(vp(-2:3), vm(-2:3), rp, rm)
   case(7)  ; call seth_linear7_interp(vp(-3:4), vm(-3:4), rp, rm)
   case(8)  ; call seth_linear9_interp(vp(-4:5), vm(-4:5), rp, rm)
   case(9)  ; call seth_linear6(vp(-2:3), vm(-2:3), rp, rm)
   case(10) ; call seth_linear6_interp(vp(-2:3), vm(-2:3), rp, rm)
   case(11) ; call seth_linear8(vp(-3:4), vm(-3:4), rp, rm)
   endselect
   endsubroutine kernel

   subroutine cmp(rp, rm, sp, sm)
   real(R_P), intent(in) :: rp, rm, sp, sm
   worst = max(worst, abs(rp - sp)/maxval(abs(vp)), abs(rm - sm)/maxval(abs(vp)))
   endsubroutine cmp
endprogram seth_test_linear
