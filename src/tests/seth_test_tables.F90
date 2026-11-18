!< SETH test: the compiled coefficient tables have the properties the generators checked (what ships, re-verified).
program seth_test_tables
!< SETH test: the compiled coefficient tables have the properties the generators checked (what ships, re-verified).
!<
!< The generators check the tables in exact arithmetic before writing them; this test re-checks the compiled constants,
!< so that a hand edit of a generated file, an emission or rounding bug, or a stale table cannot ship unnoticed. Two
!< TENO mutants (a tau coefficient, an indicator pivot) passed every other test: they change only TENO's cut
!< decisions, which no test stencil sits near. Checked, in double precision:
!<
!< - WENO-Z and TENO descriptors: every candidate exact for polynomials below its width (cell averages, face at 0);
!<   sum d = 1 and the linear combination exact below the scheme's order; every indicator, global one included, 0 on
!<   constants and 1 on unit-slope data; tau of its order on smooth data (cell averages of the generic function g of
!<   seth_test_functions, h halved; not exp, whose averages are a constant multiple of its point values).
!< - The interpolation descriptors (_INTERP_P, _INTERP_D) of WENO-Z3/5/7/9, TENO5, TENO6: every candidate exact on
!<   point values at the cell centres below its width; the weights positive, summing to 1, the combination exact
!<   below the scheme's order; tau (the reconstruction indicators on point values) of its order on point values of g.
!< - WENO-JS tables: each interface's weights sum to 1; the candidates and the optimal combination exact below 2S-1
!<   (cell averages; point values for the interpolation tables); the indicators 0 on constants and WENO<N>_JS_D_SCALE
!<   on unit slope.
!< - The non-polynomial tables of issue #6 (seth_nonpoly_tables, written by scripts/gen_nonpoly_tables.py): every
!<   polynomial row NP_RBF_P(:,mask) zero outside its stencil and exact on cell averages below its width (the masks
!<   1, 2, 4 three points, 3 and 5 four points, 0 and 7 five points; 6, TENO5's renormalised weights, below 3); the pair
!<   rows equal to the pair weights of TENO5's candidates (issue #6, section 1.2: (1,1,0) q0/2 + q1/2, (1,0,1)
!<   3 q0/4 + q2/4) and the row of the mask 6 to 3 q1/4 + q2/4 (decision D4); every s^2 row NP_RBF_C(:,mask) zero
!<   outside its stencil and summing to zero (Liang & Fu 2024, Eq. 25: the correction vanishes on constants); and the
!<   DR steepness literals SETH_DR_T1, SETH_DR_T2 within one ulp of the intrinsic tanh(0.91), tanh(1.23) (step 2.2).
use seth
use seth_test_parameters, only : ORDERS_AT_FLOOR_OK, TAU_FLOOR, TAU_ORDER_LOSS, TOL_TABLES
use seth_weno_js_tables
use seth_z_teno_tables
use seth_nonpoly_tables, only : NP_RBF_P, NP_RBF_C, SETH_DR_T1, SETH_DR_T2
use seth_test_functions, only : gfun, gavg
implicit none
real(R_P), parameter :: TOL = TOL_TABLES
logical              :: passed

passed = .true.
call descriptor('wenoz3', 3, 3._R_P, WENOZ3_NC, WENOZ3_MAXW, WENOZ3_CLO, WENOZ3_CW, WENOZ3_P, WENOZ3_D, WENOZ3_L, &
                WENOZ3_DQ, WENOZ3_TC, WENOZ3_TG, WENOZ3_GLO, WENOZ3_GW, WENOZ3_LG, WENOZ3_DG)
call descriptor('wenoz5', 5, 5._R_P, WENOZ5_NC, WENOZ5_MAXW, WENOZ5_CLO, WENOZ5_CW, WENOZ5_P, WENOZ5_D, WENOZ5_L, &
                WENOZ5_DQ, WENOZ5_TC, WENOZ5_TG, WENOZ5_GLO, WENOZ5_GW, WENOZ5_LG, WENOZ5_DG)
call descriptor('wenoz7', 7, 7._R_P, WENOZ7_NC, WENOZ7_MAXW, WENOZ7_CLO, WENOZ7_CW, WENOZ7_P, WENOZ7_D, WENOZ7_L, &
                WENOZ7_DQ, WENOZ7_TC, WENOZ7_TG, WENOZ7_GLO, WENOZ7_GW, WENOZ7_LG, WENOZ7_DG)
call descriptor('wenoz9', 9, 7._R_P, WENOZ9_NC, WENOZ9_MAXW, WENOZ9_CLO, WENOZ9_CW, WENOZ9_P, WENOZ9_D, WENOZ9_L, &
                WENOZ9_DQ, WENOZ9_TC, WENOZ9_TG, WENOZ9_GLO, WENOZ9_GW, WENOZ9_LG, WENOZ9_DG)
call descriptor('teno5', 5, 5._R_P, TENO5_NC, TENO5_MAXW, TENO5_CLO, TENO5_CW, TENO5_P, TENO5_D, TENO5_L, &
                TENO5_DQ, TENO5_TC, TENO5_TG, TENO5_GLO, TENO5_GW, TENO5_LG, TENO5_DG)
call descriptor('teno6', 6, 6._R_P, TENO6_NC, TENO6_MAXW, TENO6_CLO, TENO6_CW, TENO6_P, TENO6_D, TENO6_L, &
                TENO6_DQ, TENO6_TC, TENO6_TG, TENO6_GLO, TENO6_GW, TENO6_LG, TENO6_DG)
call descriptor('teno6_opt', 5, 6._R_P, TENO6_OPT_NC, TENO6_OPT_MAXW, TENO6_OPT_CLO, TENO6_OPT_CW, TENO6_OPT_P, &
                TENO6_OPT_D, TENO6_OPT_L, TENO6_OPT_DQ, TENO6_OPT_TC, TENO6_OPT_TG, TENO6_OPT_GLO, TENO6_OPT_GW, &
                TENO6_OPT_LG, TENO6_OPT_DG)
call descriptor('teno6_neno6', 4, 6._R_P, TENO6_NENO6_NC, TENO6_NENO6_MAXW, TENO6_NENO6_CLO, TENO6_NENO6_CW, &
                TENO6_NENO6_P, TENO6_NENO6_D, TENO6_NENO6_L, TENO6_NENO6_DQ, TENO6_NENO6_TC, TENO6_NENO6_TG, &
                TENO6_NENO6_GLO, TENO6_NENO6_GW, TENO6_NENO6_LG, TENO6_NENO6_DG)
call descriptor('teno5_opt', 4, 5._R_P, TENO5_OPT_NC, TENO5_OPT_MAXW, TENO5_OPT_CLO, TENO5_OPT_CW, TENO5_OPT_P, &
                TENO5_OPT_D, TENO5_OPT_L, TENO5_OPT_DQ, TENO5_OPT_TC, TENO5_OPT_TG, TENO5_OPT_GLO, TENO5_OPT_GW, &
                TENO5_OPT_LG, TENO5_OPT_DG)
call descriptor('teno5_az', 5, 4._R_P, TENO5_AZ_NC, TENO5_AZ_MAXW, TENO5_AZ_CLO, TENO5_AZ_CW, TENO5_AZ_P, &
                TENO5_AZ_D, TENO5_AZ_L, TENO5_AZ_DQ, TENO5_AZ_TC, TENO5_AZ_TG, TENO5_AZ_GLO, TENO5_AZ_GW, &
                TENO5_AZ_LG, TENO5_AZ_DG)
call descriptor('teno8', 6, 6._R_P, TENO8_NC, TENO8_MAXW, TENO8_CLO, TENO8_CW, TENO8_P, &
                TENO8_D, TENO8_L, TENO8_DQ, TENO8_TC, TENO8_TG, TENO8_GLO, TENO8_GW, &
                TENO8_LG, TENO8_DG)
call interp_descriptor('wenoz3_interp', 3, 3._R_P, WENOZ3_NC, WENOZ3_MAXW, WENOZ3_CLO, WENOZ3_CW, WENOZ3_INTERP_P, &
                       WENOZ3_INTERP_D, WENOZ3_L, WENOZ3_DQ, WENOZ3_TC, WENOZ3_TG, WENOZ3_GLO, WENOZ3_GW, WENOZ3_LG, &
                       WENOZ3_DG)
call interp_descriptor('wenoz5_interp', 5, 5._R_P, WENOZ5_NC, WENOZ5_MAXW, WENOZ5_CLO, WENOZ5_CW, WENOZ5_INTERP_P, &
                       WENOZ5_INTERP_D, WENOZ5_L, WENOZ5_DQ, WENOZ5_TC, WENOZ5_TG, WENOZ5_GLO, WENOZ5_GW, WENOZ5_LG, &
                       WENOZ5_DG)
call interp_descriptor('wenoz7_interp', 7, 7._R_P, WENOZ7_NC, WENOZ7_MAXW, WENOZ7_CLO, WENOZ7_CW, WENOZ7_INTERP_P, &
                       WENOZ7_INTERP_D, WENOZ7_L, WENOZ7_DQ, WENOZ7_TC, WENOZ7_TG, WENOZ7_GLO, WENOZ7_GW, WENOZ7_LG, &
                       WENOZ7_DG)
call interp_descriptor('wenoz9_interp', 9, 7._R_P, WENOZ9_NC, WENOZ9_MAXW, WENOZ9_CLO, WENOZ9_CW, WENOZ9_INTERP_P, &
                       WENOZ9_INTERP_D, WENOZ9_L, WENOZ9_DQ, WENOZ9_TC, WENOZ9_TG, WENOZ9_GLO, WENOZ9_GW, WENOZ9_LG, &
                       WENOZ9_DG)
call interp_descriptor('teno5_interp', 5, 5._R_P, TENO5_NC, TENO5_MAXW, TENO5_CLO, TENO5_CW, TENO5_INTERP_P, &
                       TENO5_INTERP_D, TENO5_L, TENO5_DQ, TENO5_TC, TENO5_TG, TENO5_GLO, TENO5_GW, TENO5_LG, TENO5_DG)
call interp_descriptor('teno6_interp', 6, 6._R_P, TENO6_NC, TENO6_MAXW, TENO6_CLO, TENO6_CW, TENO6_INTERP_P, &
                       TENO6_INTERP_D, TENO6_L, TENO6_DQ, TENO6_TC, TENO6_TG, TENO6_GLO, TENO6_GW, TENO6_LG, TENO6_DG)
call js('weno3_js', 2, WENO3_JS_A, WENO3_JS_P, WENO3_JS_D, WENO3_JS_D_SCALE, WENO3_JS_INTERP_A, WENO3_JS_INTERP_P)
call js('weno5_js', 3, WENO5_JS_A, WENO5_JS_P, WENO5_JS_D, WENO5_JS_D_SCALE, WENO5_JS_INTERP_A, WENO5_JS_INTERP_P)
call js('weno7_js', 4, WENO7_JS_A, WENO7_JS_P, WENO7_JS_D, WENO7_JS_D_SCALE, WENO7_JS_INTERP_A, WENO7_JS_INTERP_P)
call js('weno9_js', 5, WENO9_JS_A, WENO9_JS_P, WENO9_JS_D, WENO9_JS_D_SCALE, WENO9_JS_INTERP_A, WENO9_JS_INTERP_P)
call nonpoly
if (.not.passed) error stop 'a compiled table lost a property'
print '(A)', 'seth_test_tables: ok'

contains
   pure function avg(m, n) result(a)
   !< Average of x**n over cell m = [m-1, m].
   integer(I4P), intent(in) :: m, n
   real(R_P)                :: a
   a = (real(m, R_P)**(n+1) - real(m-1, R_P)**(n+1))/(n+1)
   endfunction avg

   subroutine nonpoly()
   !< The tables of the non-polynomial kernels (issue #6): rows by mask, the pair weights, the s^2 rows, the DR literals.
   integer(I4P), parameter :: LO(0:7) = [-2_I4P, -1_I4P, 0_I4P, -1_I4P, -2_I4P, -2_I4P, -2_I4P, -2_I4P] !< stencils
   integer(I4P), parameter :: HI(0:7) = [ 2_I4P,  1_I4P, 2_I4P,  2_I4P,  0_I4P,  1_I4P,  2_I4P,  2_I4P]
   integer(I4P), parameter :: EXACT(0:7) = [5_I4P, 3_I4P, 3_I4P, 4_I4P, 3_I4P, 4_I4P, 3_I4P, 5_I4P] !< exact below
   real(R_P)               :: q(-2:2,0:2), pair(-2:2), e
   integer(I4P)            :: k, m, n, j
   character(1)            :: c
   ! TENO5's candidates on -2:2 (cells by offset from the face's upwind cell)
   q = 0._R_P
   do k=0, 2
      do j=0, 2
         q(TENO5_CLO(k)+j,k) = TENO5_P(j,k)
      enddo
   enddo
   do k=0, 7
      write(c, '(I1)') k
      do m=-2, 2
         if ((m < LO(k) .or. m > HI(k)) .and. (NP_RBF_P(m,k) /= 0._R_P .or. NP_RBF_C(m,k) /= 0._R_P)) &
            call fail('nonpoly row of mask '//c//': a cell outside its stencil')
      enddo
      do n=0, EXACT(k)-1 ! cell m = [m-1, m] with the face at 0
         e = -merge(1._R_P, 0._R_P, n == 0)
         do m=-2, 2
            e = e + NP_RBF_P(m,k)*avg(m, n)
         enddo
         if (abs(e) > TOL) call fail('nonpoly row of mask '//c//': not exact below its width')
      enddo
      if (abs(sum(NP_RBF_C(:,k))) > TOL) call fail('nonpoly s^2 row of mask '//c//': does not sum to zero')
   enddo
   do k=0, 2
      if (any(abs(NP_RBF_P(:,2**k) - q(:,k)) > TOL)) call fail('nonpoly single-candidate row differs from TENO5''s')
   enddo
   pair = 0.5_R_P*q(:,0) + 0.5_R_P*q(:,1)
   if (any(abs(NP_RBF_P(:,3) - pair) > TOL)) call fail('nonpoly pair row (1,1,0) /= q0/2 + q1/2')
   pair = 0.75_R_P*q(:,0) + 0.25_R_P*q(:,2)
   if (any(abs(NP_RBF_P(:,5) - pair) > TOL)) call fail('nonpoly pair row (1,0,1) /= 3 q0/4 + q2/4')
   pair = 0.75_R_P*q(:,1) + 0.25_R_P*q(:,2)
   if (any(abs(NP_RBF_P(:,6) - pair) > TOL)) call fail('nonpoly row (0,1,1) /= 3 q1/4 + q2/4 (D4)')
   if (any(abs(NP_RBF_P(:,7) - NP_RBF_P(:,0)) > 0._R_P)) call fail('nonpoly rows of the masks 0 and 7 differ')
   if (abs(SETH_DR_T1 - tanh(0.91_R_P)) > spacing(SETH_DR_T1)) call fail('SETH_DR_T1 is not tanh(0.91) to one ulp')
   if (abs(SETH_DR_T2 - tanh(1.23_R_P)) > spacing(SETH_DR_T2)) call fail('SETH_DR_T2 is not tanh(1.23) to one ulp')
   print '(A,2ES10.2)', '   nonpoly: 8 rows, 8 s^2 rows, the pairs, D4, DR literals (ulp offsets) ', &
         (SETH_DR_T1 - tanh(0.91_R_P))/spacing(SETH_DR_T1), (SETH_DR_T2 - tanh(1.23_R_P))/spacing(SETH_DR_T2)
   endsubroutine nonpoly

   subroutine fail(what)
   character(*), intent(in) :: what
   print '(A)', '   FAILED: '//what
   passed = .false.
   endsubroutine fail

   subroutine descriptor(name, order, tau_order, nc, maxw, clo, cw, cp, d, cl, dq, tc, tg, glo, gw, lg, dg)
   character(*), intent(in) :: name
   integer,      intent(in) :: order
   real(R_P),    intent(in) :: tau_order
   integer(I4P), intent(in) :: nc, maxw, clo(0:nc-1), cw(0:nc-1), glo, gw
   real(R_P),    intent(in) :: cp(0:maxw-1,0:nc-1), d(0:nc-1), cl(0:maxw-1,0:maxw-1,0:nc-1), dq(0:maxw-1,0:nc-1)
   real(R_P),    intent(in) :: tc(0:nc-1), tg, lg(0:gw-1,0:gw-1), dg(0:gw-1)
   real(R_P)                :: s, lin, u(-8:8), b, tau(4), h, o
   integer(I4P)             :: k, j, n, m, lev
   print '(A)', name
   ! candidates exact below their width
   do k=0, nc-1
      do n=0, cw(k)-1
         s = 0._R_P
         do j=0, cw(k)-1
            s = s + cp(j,k)*avg(clo(k)+j, n)
         enddo
         if (abs(s - merge(1._R_P, 0._R_P, n == 0)) > TOL*10._R_P**n) call fail('candidate not exact')
      enddo
   enddo
   ! linear weights
   ! 16 eps: the NENO6 paper's published decimals sum to 1 + 1.3e-15 (the kernels divide by sum d delta anyway)
   if (abs(sum(d) - 1._R_P) > 16._R_P*epsilon(1._R_P)) call fail('linear weights do not sum to 1')
   do n=0, order-1
      lin = 0._R_P
      do k=0, nc-1
         do j=0, cw(k)-1
            lin = lin + d(k)*cp(j,k)*avg(clo(k)+j, n)
         enddo
      enddo
      if (abs(lin - merge(1._R_P, 0._R_P, n == 0)) > TOL*10._R_P**n) call fail('linear scheme not of its order')
   enddo
   ! indicators: 0 on constants, 1 on unit slope
   do k=0, nc-1
      do m=-8, 8 ; u(m) = 1._R_P ; enddo
      if (abs(beta(k, u, nc, maxw, clo, cw, cl, dq)) > TOL) call fail('indicator not blind to constants')
      do m=-8, 8 ; u(m) = real(m, R_P) ; enddo
      if (abs(beta(k, u, nc, maxw, clo, cw, cl, dq) - 1._R_P) > TOL) call fail('indicator not 1 on unit slope')
   enddo
   if (tg /= 0._R_P) then
      do m=-8, 8 ; u(m) = 1._R_P ; enddo
      if (abs(betag(u, gw, glo, lg, dg)) > TOL) call fail('global indicator not blind to constants')
      do m=-8, 8 ; u(m) = real(m, R_P) ; enddo
      if (abs(betag(u, gw, glo, lg, dg) - 1._R_P) > TOL) call fail('global indicator not 1 on unit slope')
   endif
   ! tau of its order on smooth data
   do lev=1, 4
      h = 0.2_R_P/2._R_P**(lev-1)
      do m=-8, 8 ! cell averages of g on cells of width h
         u(m) = gavg((m-1)*h, m*h)
      enddo
      tau(lev) = 0._R_P
      do k=0, nc-1
         tau(lev) = tau(lev) + tc(k)*beta(k, u, nc, maxw, clo, cw, cl, dq)
      enddo
      if (tg /= 0._R_P) tau(lev) = tau(lev) + tg*betag(u, gw, glo, lg, dg)
      tau(lev) = abs(tau(lev))
   enddo
   o = tau_observed(tau)
   print '(A,F6.2,A,F4.1,A)', '   tau order ', o, ' (required ', tau_order, ')'
   if (o < tau_order - 0.3_R_P - TAU_ORDER_LOSS .and. .not.(ORDERS_AT_FLOOR_OK .and. o < 0._R_P)) &
      call fail('tau below its order')
   endsubroutine descriptor

   subroutine interp_descriptor(name, order, tau_order, nc, maxw, clo, cw, cp, d, cl, dq, tc, tg, glo, gw, lg, dg)
   !< The interpolation tables of a descriptor (the indicators and tau are the reconstruction ones, on point values).
   character(*), intent(in) :: name
   integer,      intent(in) :: order
   real(R_P),    intent(in) :: tau_order
   integer(I4P), intent(in) :: nc, maxw, clo(0:nc-1), cw(0:nc-1), glo, gw
   real(R_P),    intent(in) :: cp(0:maxw-1,0:nc-1), d(0:nc-1), cl(0:maxw-1,0:maxw-1,0:nc-1), dq(0:maxw-1,0:nc-1)
   real(R_P),    intent(in) :: tc(0:nc-1), tg, lg(0:gw-1,0:gw-1), dg(0:gw-1)
   real(R_P)                :: s, lin, u(-8:8), tau(4), h, o
   integer(I4P)             :: k, j, n, m, lev
   print '(A)', name
   do k=0, nc-1 ! candidates exact on point values (cell centres m - 1/2) below their width
      do n=0, cw(k)-1
         s = 0._R_P
         do j=0, cw(k)-1
            s = s + cp(j,k)*(clo(k) + j - 0.5_R_P)**n
         enddo
         if (abs(s - merge(1._R_P, 0._R_P, n == 0)) > TOL*10._R_P**n) call fail('interpolation candidate not exact')
      enddo
   enddo
   if (any(d <= 0._R_P)) call fail('an interpolation weight is not positive')
   if (abs(sum(d) - 1._R_P) > 16._R_P*epsilon(1._R_P)) call fail('interpolation weights do not sum to 1')
   do n=0, order-1
      lin = 0._R_P
      do k=0, nc-1
         do j=0, cw(k)-1
            lin = lin + d(k)*cp(j,k)*(clo(k) + j - 0.5_R_P)**n
         enddo
      enddo
      if (abs(lin - merge(1._R_P, 0._R_P, n == 0)) > TOL*10._R_P**n) call fail('interpolation not of its order')
   enddo
   do lev=1, 4 ! tau on point values of g
      h = 0.2_R_P/2._R_P**(lev-1)
      do m=-8, 8
         u(m) = gfun((m - 0.5_R_P)*h)
      enddo
      tau(lev) = 0._R_P
      do k=0, nc-1
         tau(lev) = tau(lev) + tc(k)*beta(k, u, nc, maxw, clo, cw, cl, dq)
      enddo
      if (tg /= 0._R_P) tau(lev) = tau(lev) + tg*betag(u, gw, glo, lg, dg)
      tau(lev) = abs(tau(lev))
   enddo
   o = tau_observed(tau)
   print '(A,F6.2,A,F4.1,A)', '   tau order on point values ', o, ' (required ', tau_order, ')'
   if (o < tau_order - 0.3_R_P - TAU_ORDER_LOSS .and. .not.(ORDERS_AT_FLOOR_OK .and. o < 0._R_P)) &
      call fail('tau below its order on point values')
   endsubroutine interp_descriptor

   pure function beta(k, u, nc, maxw, clo, cw, cl, dq) result(b)
   !< Indicator of candidate k from its factors: sum_i DQ(i) (sum_{a>=i} L(a,i) u_a)^2.
   integer(I4P), intent(in) :: k, nc, maxw, clo(0:nc-1), cw(0:nc-1)
   real(R_P),    intent(in) :: u(-8:8), cl(0:maxw-1,0:maxw-1,0:nc-1), dq(0:maxw-1,0:nc-1)
   real(R_P)                :: b, ls
   integer(I4P)             :: i, a
   b = 0._R_P
   do i=0, cw(k)-1
      ls = 0._R_P
      do a=i, cw(k)-1
         ls = ls + cl(a,i,k)*u(clo(k)+a)
      enddo
      b = b + dq(i,k)*ls*ls
   enddo
   endfunction beta

   pure function betag(u, gw, glo, lg, dg) result(b)
   !< The global indicator from its factors.
   integer(I4P), intent(in) :: gw, glo
   real(R_P),    intent(in) :: u(-8:8), lg(0:gw-1,0:gw-1), dg(0:gw-1)
   real(R_P)                :: b, ls
   integer(I4P)             :: i, a
   b = 0._R_P
   do i=0, gw-1
      ls = 0._R_P
      do a=i, gw-1
         ls = ls + lg(a,i)*u(glo+a)
      enddo
      b = b + dg(i)*ls*ls
   enddo
   endfunction betag

   subroutine js(name, s, a, p, d, scale, ai, pi_)
   character(*), intent(in) :: name
   integer,      intent(in) :: s
   real(R_P),    intent(in) :: a(1:2,0:s-1), p(1:2,0:s-1,0:s-1), d(0:s-1,0:s-1,0:s-1), scale
   real(R_P),    intent(in) :: ai(1:2,0:s-1), pi_(1:2,0:s-1,0:s-1)
   real(R_P)                :: face, lin, val, b, u(-8:8)
   integer(I4P)             :: f, k, c, n, m, s2, s3, t
   print '(A)', name
   do t=1, 2 ! 1: reconstruction (cell averages), 2: interpolation (point values)
      do f=1, 2
         face = merge(-0.5_R_P, 0.5_R_P, f == 1) ! row 1: face i-1/2 of the centre cell, row 2: i+1/2
         if (abs(sum(merge(a(f,:), ai(f,:), t == 1)) - 1._R_P) > 4._R_P*epsilon(1._R_P)) call fail('weights sum')
         do n=0, 2*s-2
            lin = 0._R_P
            do k=0, s-1
               do c=0, s-1
                  m = k - c
                  if (t == 1) then
                     val = ((real(m, R_P) + 0.5_R_P)**(n+1) - (real(m, R_P) - 0.5_R_P)**(n+1))/(n+1)
                     lin = lin + a(f,k)*p(f,c,k)*val
                  else
                     lin = lin + ai(f,k)*pi_(f,c,k)*real(m, R_P)**n
                  endif
               enddo
            enddo
            if (abs(lin - face**n) > TOL*10._R_P**n) call fail('optimal combination not exact')
         enddo
      enddo
   enddo
   do k=0, s-1 ! indicators: d(s3,s2,k) v(k-s3) v(k-s2)
      do m=-8, 8 ; u(m) = 1._R_P ; enddo
      b = 0._R_P
      do s2=0, s-1 ; do s3=0, s-1 ; b = b + d(s3,s2,k)*u(k-s3)*u(k-s2) ; enddo ; enddo
      if (abs(b) > TOL*scale) call fail('indicator not blind to constants')
      do m=-8, 8 ; u(m) = real(m, R_P) ; enddo
      b = 0._R_P
      do s2=0, s-1 ; do s3=0, s-1 ; b = b + d(s3,s2,k)*u(k-s3)*u(k-s2) ; enddo ; enddo
      if (abs(b - scale) > TOL*scale*100._R_P) call fail('indicator not D_SCALE on unit slope')
   enddo
   endsubroutine js

   pure function tau_observed(tau) result(o)
   !< The order of tau on the finest pair of h (0.2, 0.1, 0.05, 0.025) whose finer tau is above TAU_FLOOR (R8P 0: the
   !< pair 0.05, 0.025); -1 if none. tau cancels to h^order out of indicators of order h^2: in R4P its finest values
   !< are round-off.
   real(R_P), intent(in) :: tau(4)
   real(R_P)             :: o
   integer(I4P)          :: j
   o = -1._R_P
   do j=3, 1, -1
      if (tau(j+1) > TAU_FLOOR) then
         o = log(tau(j)/tau(j+1))/log(2._R_P) ; return
      endif
   enddo
   endfunction tau_observed
endprogram seth_test_tables
