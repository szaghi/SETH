!< SETH test: TENO5, TENO6, TENO6-opt, TENO6-NENO6 and the TENO5/6 interpolations: order, linear scheme on smooth data,
!< cutoff across a jump, no overflow on a large jump, exact mirroring.
program seth_test_teno
!< SETH test: TENO5, TENO6, TENO6-opt, TENO6-NENO6 and the TENO5/6 interpolations: order, linear scheme on smooth data,
!< cutoff across a jump, no overflow on a large jump, exact mirroring.
!<
!< 1. Order on u = g(x), the generic function of seth_test_functions (cell averages, or point values at the cell
!<    centres for the interpolations; face at 0.3), >= order - 0.5: TENO5 5, TENO6 6, TENO6-opt 5 (Fu et al.'s
!<    optimised weights lose one order), TENO6-NENO6 4 (the NENO6 paper's weights lose two), the interpolations 5 and
!<    6. C_T: the sourced values.
!< 2. Smooth data keeps every candidate (delta = 1): the state equals the linear scheme sum d_k q_k / sum d to 4 eps of
!<    the data magnitude.
!< 3. A unit step at every cell boundary, both orientations: the states lie in [0, 1] to 1e-12 (the candidates
!<    crossing the jump are cut, the kept ones reproduce constants).
!< 4. A step of height 1e30: the states are finite and in [0, 1e30] to round-off (gamma = g**6 alone would
!<    overflow).
!< 5. vm(m) = vp(1-m) gives rm == rp bitwise.
use seth
use seth_test_rng, only : seed_rng
use seth_z_teno_tables
use seth_test_functions, only : gfun, gavg
implicit none
real(R8P),    parameter :: XF = 0.3_R8P
integer(I4P), parameter :: NH = 6_I4P
real(R8P)               :: vp(-2:3), vm(-2:3), rp, rm, h, err(2,NH), order, floor_, worst, x, lin, ct, big
integer(I4P)            :: sc, j, m, side, jump, n_
logical                 :: passed
integer(I4P), parameter :: NSC = 6_I4P
character(12), parameter :: NAMES(NSC) = ['TENO5       ', 'TENO6       ', 'TENO6-opt   ', 'TENO6-NENO6 ', &
                                          'TENO5-interp', 'TENO6-interp']
real(R8P),    parameter :: DESIGN(NSC) = [5._R8P, 6._R8P, 5._R8P, 4._R8P, 5._R8P, 6._R8P]
real(R8P),    parameter :: CTS(NSC) = [SETH_TENO5_CT, SETH_TENO6_CT, SETH_TENO6_OPT_CT, SETH_TENO6_NENO6_CT, &
                                       SETH_TENO5_CT, SETH_TENO6_CT]

passed = .true.
call seed_rng(1005_I4P)
floor_ = 2.e2_R8P * epsilon(1._R8P) ! 100 eps |g|, |g| < 2 near XF
do sc=1, NSC
   ct = CTS(sc)
   ! 1. order
   do j=1, NH
      h = 0.4_R8P / 2._R8P**(j-1)
      do m=-2, 3
         if (sc > 4) then
            vp(m) = gfun(XF + (m - 0.5_R8P)*h)
         else
            vp(m) = gavg(XF + (m-1)*h, XF + m*h)
         endif
      enddo
      call kernel(sc, vp, vp, ct, rp, rm)
      err(:,j) = abs([rp, rm] - gfun(XF))
   enddo
   do side=1, 2
      order = -1._R8P
      do j=NH-1, 1, -1
         if (err(side,j+1) > floor_) then
            order = log(err(side,j)/err(side,j+1)) / log(2._R8P) ; exit
         endif
      enddo
      print '(A,A,A,A,F6.2)', trim(NAMES(sc)), ' state=', merge('+','-',side==1), ' order=', order
      if (order < DESIGN(sc) - 0.5_R8P) passed = .false.
   enddo
   ! 2. linear scheme on smooth data
   worst = 0._R8P
   do n_=1, 1000
      call random_number(x)
      do m=-2, 3
         vp(m) = 1._R8P + x*m + 0.1_R8P*x*m*m + 0.01_R8P*x*m**3
      enddo
      call kernel(sc, vp, vp, ct, rp, rm)
      lin = linear(sc, vp)
      worst = max(worst, abs(rp - lin)/maxval(abs(vp)))
   enddo
   print '(A,A,ES9.2)', trim(NAMES(sc)), ' smooth data vs linear scheme: ', worst
   if (worst > 4._R8P*epsilon(1._R8P)) passed = .false.
   ! 3. unit step
   worst = 0._R8P
   do side=1, 2
      do jump=-1, 3
         do m=-2, 3
            vp(m) = merge(1._R8P, 0._R8P, (m >= jump) .eqv. (side == 1))
         enddo
         call kernel(sc, vp, vp, ct, rp, rm)
         worst = max(worst, -rp, rp - 1._R8P, -rm, rm - 1._R8P)
      enddo
   enddo
   print '(A,A,ES9.2)', trim(NAMES(sc)), ' worst overshoot: ', worst
   if (worst > 1.e-12_R8P) passed = .false.
   ! 4. large step
   big = 1.e30_R8P
   do jump=-1, 3
      do m=-2, 3
         vp(m) = merge(big, 0._R8P, m >= jump)
      enddo
      call kernel(sc, vp, vp, ct, rp, rm)
      x = 4._R8P*epsilon(1._R8P)*big ! round-off of a convex combination of 1e30
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
if (.not.passed) error stop 'TENO test failed'
print '(A)', 'seth_test_teno: ok'

contains
   subroutine kernel(sc, vp, vm, ct, rp, rm)
   integer(I4P), intent(in)  :: sc
   real(R8P),    intent(in)  :: vp(-2:3), vm(-2:3), ct
   real(R8P),    intent(out) :: rp, rm
   select case(sc)
   case(1) ; call seth_teno5(vp, vm, SETH_TENO_EPS, ct, rp, rm)
   case(2) ; call seth_teno6(vp, vm, SETH_TENO_EPS, ct, rp, rm)
   case(3) ; call seth_teno6_opt(vp, vm, SETH_TENO_EPS, ct, rp, rm)
   case(4) ; call seth_teno6_neno6(vp, vm, SETH_TENO_EPS, ct, rp, rm)
   case(5) ; call seth_teno5_interp(vp, vm, SETH_TENO_EPS, ct, rp, rm)
   case(6) ; call seth_teno6_interp(vp, vm, SETH_TENO_EPS, ct, rp, rm)
   endselect
   endsubroutine kernel

   function linear(sc, v) result(r)
   !< The linear scheme of the descriptor: sum_k d_k q_k / sum d.
   integer(I4P), intent(in) :: sc
   real(R8P),    intent(in) :: v(-2:3)
   real(R8P)                :: r
   select case(sc)
   case(1) ; r = lin_of(TENO5_NC, TENO5_MAXW, TENO5_CLO, TENO5_CW, TENO5_P, TENO5_D, v)
   case(2) ; r = lin_of(TENO6_NC, TENO6_MAXW, TENO6_CLO, TENO6_CW, TENO6_P, TENO6_D, v)
   case(3) ; r = lin_of(TENO6_OPT_NC, TENO6_OPT_MAXW, TENO6_OPT_CLO, TENO6_OPT_CW, TENO6_OPT_P, TENO6_OPT_D, v)
   case(4) ; r = lin_of(TENO6_NENO6_NC, TENO6_NENO6_MAXW, TENO6_NENO6_CLO, TENO6_NENO6_CW, TENO6_NENO6_P, &
                        TENO6_NENO6_D, v)
   case(5) ; r = lin_of(TENO5_NC, TENO5_MAXW, TENO5_CLO, TENO5_CW, TENO5_INTERP_P, TENO5_INTERP_D, v)
   case(6) ; r = lin_of(TENO6_NC, TENO6_MAXW, TENO6_CLO, TENO6_CW, TENO6_INTERP_P, TENO6_INTERP_D, v)
   endselect
   endfunction linear

   function lin_of(nc, maxw, clo, cw, cp, d, v) result(r)
   integer(I4P), intent(in) :: nc, maxw, clo(0:nc-1), cw(0:nc-1)
   real(R8P),    intent(in) :: cp(0:maxw-1,0:nc-1), d(0:nc-1), v(-2:3)
   real(R8P)                :: r, q
   integer(I4P)             :: k, a
   r = 0._R8P
   do k=0, nc-1
      q = 0._R8P
      do a=0, cw(k)-1
         q = q + cp(a,k)*v(clo(k)+a)
      enddo
      r = r + (d(k)/sum(d))*q
   enddo
   endfunction lin_of
endprogram seth_test_teno
