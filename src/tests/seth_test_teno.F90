!< SETH test: TENO5, TENO5-opt, TENO6, TENO6-opt, TENO6-NENO6 and the TENO5/6 interpolations: order, linear scheme on
!< smooth data, cutoff across a jump, no overflow on a large jump, exact mirroring.
program seth_test_teno
!< SETH test: TENO5, TENO5-opt, TENO6, TENO6-opt, TENO6-NENO6 and the TENO5/6 interpolations: order, linear scheme on
!< smooth data, cutoff across a jump, no overflow on a large jump, exact mirroring.
!<
!< 1. Order on u = g(x), the generic function of seth_test_functions (cell averages, or point values at the cell
!<    centres for the interpolations; face at 0.3), >= order - 0.5: TENO5 5, TENO6 6, TENO6-opt 5 and TENO5-opt 4 (Fu
!<    et al.'s optimised weights lose one order), TENO6-NENO6 4 (the NENO6 paper's weights lose two), the interpolations 5 and
!<    6. C_T: the sourced values.
!< 2. Smooth data keeps every candidate (delta = 1): the state equals the linear scheme sum d_k q_k / sum d to 4 eps of
!<    the data magnitude.
!< 3. A unit step at every cell boundary, both orientations: the states lie in [0, 1] to 1e-12 (the candidates
!<    crossing the jump are cut, the kept ones reproduce constants).
!< 4. A step of height 1e30: the states are finite and in [0, 1e30] to round-off (gamma = g**6 alone would
!<    overflow).
!< 5. vm(m) = vp(1-m) gives rm == rp bitwise.
!< 6. TENO5-A and TENO8-A with a2 = 0 equal TENO5 at C_T = 1e-5 and TENO8 at 1e-7 bitwise (the cutoff 10^-floor(a1)).
use seth
use seth_test_parameters, only : BIG_JUMP, ORDERS_AT_FLOOR_OK, POINT_FLOOR, POINT_ORDER_LOSS, TOL_BOUNDS
use seth_test_rng, only : seed_rng
use seth_z_teno_tables
use seth_test_functions, only : gfun, gavg
use seth_vv_cases, only : vv_report, vv_open
implicit none
real(R_P),    parameter :: XF = 0.3_R_P
integer(I4P), parameter :: NH = 6_I4P
real(R_P)               :: vp(-2:3), vm(-2:3), rp, rm, h, err(2,NH), order, floor_, worst, x, lin, ct, big
integer(I4P)            :: sc, j, m, side, jump, n_
logical                 :: passed
integer                 :: ur
integer(I4P), parameter :: NSC = 8_I4P
character(12), parameter :: NAMES(NSC) = ['TENO5       ', 'TENO6       ', 'TENO6-opt   ', 'TENO6-NENO6 ', &
                                          'TENO5-interp', 'TENO6-interp', 'TENO5-opt   ', &
                                          'TENO5-AZ    ']
real(R_P),    parameter :: DESIGN(NSC) = [5._R_P, 6._R_P, 5._R_P, 4._R_P, 5._R_P, 6._R_P, 4._R_P, 5._R_P]
real(R_P),    parameter :: CTS(NSC) = [SETH_TENO5_CT, SETH_TENO6_CT, SETH_TENO6_OPT_CT, SETH_TENO6_NENO6_CT, &
                                       SETH_TENO5_CT, SETH_TENO6_CT, SETH_TENO5_OPT_CT, &
                                       SETH_TENO5_AZ_CT]

passed = .true.
call seed_rng(1005_I4P)
if (vv_report()) call vv_open('orders-teno.dat', ur, 'scheme state observed-order design')
floor_ = POINT_FLOOR * epsilon(1._R_P) ! R8P: 100 eps |g|, |g| < 2 near XF
do sc=1, NSC
   ct = CTS(sc)
   ! 1. order
   do j=1, NH
      h = 0.4_R_P / 2._R_P**(j-1)
      do m=-2, 3
         if (sc == 5 .or. sc == 6) then
            vp(m) = gfun(XF + (m - 0.5_R_P)*h)
         else
            vp(m) = gavg(XF + (m-1)*h, XF + m*h)
         endif
      enddo
      call kernel(sc, vp, vp, ct, rp, rm)
      err(:,j) = abs([rp, rm] - gfun(XF))
   enddo
   do side=1, 2
      order = -1._R_P
      do j=NH-1, 1, -1
         if (err(side,j+1) > floor_) then
            order = log(err(side,j)/err(side,j+1)) / log(2._R_P) ; exit
         endif
      enddo
      print '(A,A,A,A,F6.2)', trim(NAMES(sc)), ' state=', merge('+','-',side==1), ' order=', order
      if (vv_report()) write(ur, '(A,1X,A1,F8.3,F5.1)') trim(NAMES(sc)), merge('+','-',side==1), order, DESIGN(sc)
      if (order < DESIGN(sc) - 0.5_R_P - POINT_ORDER_LOSS .and. .not.(ORDERS_AT_FLOOR_OK .and. order < 0._R_P)) passed = .false.
   enddo
   ! 2. linear scheme on smooth data
   worst = 0._R_P
   do n_=1, 1000
      call random_number(x)
      do m=-2, 3
         vp(m) = 1._R_P + x*m + 0.1_R_P*x*m*m + 0.01_R_P*x*m**3
      enddo
      call kernel(sc, vp, vp, ct, rp, rm)
      lin = linear(sc, vp)
      worst = max(worst, abs(rp - lin)/maxval(abs(vp)))
   enddo
   print '(A,A,ES9.2)', trim(NAMES(sc)), ' smooth data vs linear scheme: ', worst
   if (worst > 4._R_P*epsilon(1._R_P)) passed = .false.
   ! 3. unit step
   worst = 0._R_P
   do side=1, 2
      do jump=-1, 3
         do m=-2, 3
            vp(m) = merge(1._R_P, 0._R_P, (m >= jump) .eqv. (side == 1))
         enddo
         call kernel(sc, vp, vp, ct, rp, rm)
         worst = max(worst, -rp, rp - 1._R_P, -rm, rm - 1._R_P)
      enddo
   enddo
   print '(A,A,ES9.2)', trim(NAMES(sc)), ' worst overshoot: ', worst
   if (worst > TOL_BOUNDS) passed = .false.
   ! 4. large step
   big = BIG_JUMP
   do jump=-1, 3
      do m=-2, 3
         vp(m) = merge(big, 0._R_P, m >= jump)
      enddo
      call kernel(sc, vp, vp, ct, rp, rm)
      x = 4._R_P*epsilon(1._R_P)*big ! round-off of a convex combination of big
      if (.not.(rp >= -x .and. rp <= big + x .and. rm >= -x .and. rm <= big + x)) then
         print '(A,A,2ES12.4)', trim(NAMES(sc)), ' large step out of bounds or NaN: ', rp, rm ; passed = .false.
      endif
   enddo
   ! 5. mirror
   do n_=1, 1000
      do m=-2, 3
         call random_number(x) ; vp(m) = x
      enddo
      do m=-2, 3
         vm(m) = vp(1-m)
      enddo
      call kernel(sc, vp, vm, ct, rp, rm)
      if (rp /= rm) then
         print '(A,A)', trim(NAMES(sc)), ' mirror mismatch' ; passed = .false. ; exit
      endif
   enddo
enddo
if (vv_report()) close(ur)
! 6. the adaptive cutoff with a2 = 0 is the fixed cutoff 10^-floor(a1), bitwise (TENO5-A against TENO5 at 1e-5,
!    TENO8-A against TENO8 at 1e-7), on smooth, rough and stepped stencils
block
real(R_P)    :: w8(-3:4), z8(-3:4), r1, r2, s1, s2
integer(I4P) :: nd
nd = 0
do n_=1, 3000
   do m=-3, 4
      call random_number(x)
      w8(m) = merge(x, merge(1._R_P, 0._R_P, m >= mod(n_, 7) - 3), mod(n_, 3) /= 0)
      call random_number(x) ; z8(m) = x
   enddo
   call seth_teno5_a(w8(-2:3), z8(-2:3), SETH_TENO_EPS, SETH_TENO5_A_CR, SETH_TENO5_A_XI, 5._R_P, 0._R_P, 1._R_P, r1, r2)
   call seth_teno5(w8(-2:3), z8(-2:3), SETH_TENO_EPS, 1.e-5_R_P, s1, s2)
   if (r1 /= s1 .or. r2 /= s2) nd = nd + 1
   call seth_teno8_a(w8, z8, SETH_TENO_EPS, SETH_TENO8_A_CR, SETH_TENO8_A_XI, 7._R_P, 0._R_P, 1._R_P, r1, r2)
   call seth_teno8(w8, z8, SETH_TENO_EPS, 1.e-7_R_P, s1, s2)
   if (r1 /= s1 .or. r2 /= s2) nd = nd + 1
enddo
print '(A,I0)', 'adaptive cutoff with a2 = 0 against the fixed cutoff, differences: ', nd
if (nd > 0) passed = .false.
endblock
if (.not.passed) error stop 'TENO test failed'
print '(A)', 'seth_test_teno: ok'

contains
   subroutine kernel(sc, vp, vm, ct, rp, rm)
   integer(I4P), intent(in)  :: sc
   real(R_P),    intent(in)  :: vp(-2:3), vm(-2:3), ct
   real(R_P),    intent(out) :: rp, rm
   select case(sc)
   case(1) ; call seth_teno5(vp, vm, SETH_TENO_EPS, ct, rp, rm)
   case(2) ; call seth_teno6(vp, vm, SETH_TENO_EPS, ct, rp, rm)
   case(3) ; call seth_teno6_opt(vp, vm, SETH_TENO_EPS, ct, rp, rm)
   case(4) ; call seth_teno6_neno6(vp, vm, SETH_TENO_EPS, ct, rp, rm)
   case(5) ; call seth_teno5_interp(vp, vm, SETH_TENO_EPS, ct, rp, rm)
   case(6) ; call seth_teno6_interp(vp, vm, SETH_TENO_EPS, ct, rp, rm)
   case(7) ; call seth_teno5_opt(vp, vm, SETH_TENO_EPS, ct, rp, rm)
   case(8) ; call seth_teno5_az(vp, vm, SETH_TENO_EPS, ct, 1._R_P, rp, rm) ! Tang's form
   endselect
   endsubroutine kernel

   function linear(sc, v) result(r)
   !< The linear scheme of the descriptor: sum_k d_k q_k / sum d.
   integer(I4P), intent(in) :: sc
   real(R_P),    intent(in) :: v(-2:3)
   real(R_P)                :: r
   select case(sc)
   case(1) ; r = lin_of(TENO5_NC, TENO5_MAXW, TENO5_CLO, TENO5_CW, TENO5_P, TENO5_D, v)
   case(2) ; r = lin_of(TENO6_NC, TENO6_MAXW, TENO6_CLO, TENO6_CW, TENO6_P, TENO6_D, v)
   case(3) ; r = lin_of(TENO6_OPT_NC, TENO6_OPT_MAXW, TENO6_OPT_CLO, TENO6_OPT_CW, TENO6_OPT_P, TENO6_OPT_D, v)
   case(4) ; r = lin_of(TENO6_NENO6_NC, TENO6_NENO6_MAXW, TENO6_NENO6_CLO, TENO6_NENO6_CW, TENO6_NENO6_P, &
                        TENO6_NENO6_D, v)
   case(5) ; r = lin_of(TENO5_NC, TENO5_MAXW, TENO5_CLO, TENO5_CW, TENO5_INTERP_P, TENO5_INTERP_D, v)
   case(6) ; r = lin_of(TENO6_NC, TENO6_MAXW, TENO6_CLO, TENO6_CW, TENO6_INTERP_P, TENO6_INTERP_D, v)
   case(8) ; r = lin_of(TENO5_AZ_NC, TENO5_AZ_MAXW, TENO5_AZ_CLO, TENO5_AZ_CW, TENO5_AZ_P, TENO5_AZ_D, v)
   case(7) ; r = lin_of(TENO5_OPT_NC, TENO5_OPT_MAXW, TENO5_OPT_CLO, TENO5_OPT_CW, TENO5_OPT_P, TENO5_OPT_D, v)
   endselect
   endfunction linear

   function lin_of(nc, maxw, clo, cw, cp, d, v) result(r)
   integer(I4P), intent(in) :: nc, maxw, clo(0:nc-1), cw(0:nc-1)
   real(R_P),    intent(in) :: cp(0:maxw-1,0:nc-1), d(0:nc-1), v(-2:3)
   real(R_P)                :: r, q
   integer(I4P)             :: k, a
   r = 0._R_P
   do k=0, nc-1
      q = 0._R_P
      do a=0, cw(k)-1
         q = q + cp(a,k)*v(clo(k)+a)
      enddo
      r = r + (d(k)/sum(d))*q
   enddo
   endfunction lin_of
endprogram seth_test_teno
