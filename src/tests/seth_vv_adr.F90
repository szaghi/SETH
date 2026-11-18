!< SETH verification: the approximate dispersion relation (ADR) of every nonlinear reconstruction.
program seth_vv_adr
!< SETH verification: the approximate dispersion relation (ADR) of every nonlinear reconstruction.
!<
!< The procedure of Pirozzoli (2006), as the TENO and WENO papers use it (llm-wiki weno-teno
!< `definitions/approximate-dispersion-relation.md`): for each reduced wavenumber theta_n = 2 pi n/N, n = 1 .. N/2 - 1,
!< N = 200, the nonlinear finite-difference operator D_j = F_{j+1/2} - F_{j-1/2} of u_t + u_x = 0 (upwind: F = the
!< left-biased state) is applied to the single mode u_j = cos(theta_n j) on a periodic grid (h = 1), and its Fourier
!< coefficient at theta_n, G = (2/N) sum_j D_j exp(-i theta_n j), gives the modified wavenumber
!< theta~ = -i G: Re theta~ the dispersion (ideal theta), Im theta~ the dissipation (ideal 0, negative dissipative,
!< positive anti-dissipative; the sign convention of the wiki page and of seth_test_spectral). The linear background
!< of each scheme (its optimal weights) is computed analytically from its tables, as in seth_test_spectral.
!<
!< The ADR of a nonlinear scheme depends on the test signal: unit amplitude here, every scheme with its default
!< parameters (WENO-JS eps = 1e-6, which is not scale invariant, so its ADR changes with the amplitude; WENO-Z and TENO
!< eps = 1e-40, scale invariant).
!<
!< Asserted: (1) TENO is its linear background for the well-resolved waves: |theta~ - theta~_lin| <= 1e-12 for every
!< theta_n <= pi/4 (TENO5, TENO6: the cutoff keeps every candidate; Fu, Hu & Adams 2016, Fig. 3); (2) every scheme of
!< order >= 5 approaches its linear background at low wavenumbers, |theta~ - theta~_lin| <= 1e-3 theta for theta <= 0.3
!< (the third-order weights do not at unit amplitude: printed); (3) dissipation is never positive beyond round-off,
!< Im theta~ <= 1e-12, except for TENO6, whose ADR is measured anti-dissipative (Im up to 2e-4) in a narrow band
!< theta/pi = 0.46 .. 0.49 where its cutoff starts to act (it departs from its linear background at theta/pi = 0.46):
!< for TENO6 the anti-dissipative band must lie beyond theta = 0.4 pi and stay below 1e-3. Grid-scale packets of that
!< wavelength (about 4 h) are what TENO6 leaves undamped behind the Burgers shock (seth_vv_burgers). TENO8 (Liang & Fu's
!< background, whose linear scheme is dissipative everywhere) is measured anti-dissipative too, Im up to 3.2e-3 in
!< theta/pi = 0.36 .. 0.45 where its cutoff starts to act, more for smaller C_T (2.2e-3 at 1e-5, 9.6e-3 at 1e-10):
!< its band must lie beyond 0.3 pi and stay below 5e-3. TENO8-A, whose adaptive cutoff reaches 1e-10 on these waves
!< (Liang & Fu 2024b), measures Im up to 9.6e-3 in theta/pi = 0.38 .. 0.51: beyond 0.3 pi, below 1.2e-2. TENO5-A has the
!< ADR of TENO5 exactly: Ren's sensor reads the crests of a unit wave as discontinuities (eta -> 0) and keeps
!< C_T = 1e-5 there. (4) The ADR is one evaluation of the operator on one mode; for
!< every scheme with an anti-dissipative band, each mode of the band is advected for 4 periods (u_t + u_x = 0,
!< SSP-RK3, CFL 0.4) and its amplitude must never exceed the initial one: the band is not a growth of that mode.
!< G6 of issue #6 (steps 3.1, 3.2): TENO5-ST equals linear5 for theta <= pi/4 (TOL_ADR) and is nowhere anti-dissipative.
!< TENO5-RBF (s = 3 rbf_dx = 0.03, 0.06, 0.15, 0.3, frozen; its sensor raises s above s1 where it reads a crest as
!< non-smooth, so the nonlinear ADR is not exactly its row's): inside the validated range of D7 (s <= 0.12: the first
!< two), the anti-dissipative band's edge within 10 % of the predicted sqrt(30) s and its peak within twice the
!< predicted (200/3) s^6 (seth_test_nonpoly, spec 18, for the linear row). Outside it reported: measured 2026-11-18, at
!< s = 0.15 the peak is 1.6e-2, 22 times the prediction, at s = 0.3 0.198 (4 times), the band reaching theta = 1.8.
!< Report mode: adr-<scheme>.dat (theta, Re and Im of theta~, of theta~_lin); adr-teno5_st.dat, adr-teno5_rbf-s<k>.dat.
use seth
use seth_test_parameters, only : TOL_ADR
use seth_test_schemes
use seth_weno_js_tables
use seth_z_teno_tables
use seth_nonpoly_tables, only : NP_RBF_P, NP_RBF_C
use seth_vv_cases, only : vv_report, vv_open, vv_schemes, vv_label
use seth_vv_euler, only : PATH_SPLIT
implicit none
integer(I4P), parameter   :: N = 200_I4P
real(R_P),    parameter   :: PI = acos(-1._R_P)
integer(I4P), allocatable :: ids(:)
real(R_P)                 :: c(-4:5), th, worst_teno, worst_low, maxim, thanti, gmax
complex(R_P)              :: z, zl
integer(I4P)              :: k, nn, u
logical                   :: passed, lin_dissipative

passed = .true.
call vv_schemes(PATH_SPLIT, .not.vv_report(), ids)
print '(A16,3A14)', 'scheme', 'TENO vs lin', 'low-k vs lin', 'max Im'
do k=1, size(ids)
   call background(vv_label(ids(k)), c, lin_dissipative)
   worst_teno = 0._R_P ; worst_low = 0._R_P ; maxim = -huge(1._R_P) ; thanti = huge(1._R_P)
   if (vv_report()) call vv_open('adr-'//vv_label(ids(k))//'.dat', u, &
                                 'theta Re(theta~) Im(theta~) Re(theta~_lin) Im(theta~_lin) (N = 200, unit amplitude)')
   do nn=1, N/2-1
      th = 2._R_P*PI*nn/N
      z = adr(ids(k), th)
      zl = xit(c, th)
      if (th <= PI/4._R_P) worst_teno = max(worst_teno, abs(z - zl))
      if (th <= 0.3_R_P) worst_low = max(worst_low, abs(z - zl)/th)
      maxim = max(maxim, aimag(z))
      if (aimag(z) > TOL_ADR) thanti = min(thanti, th)
      if (vv_report()) write(u, '(5ES24.15)') th, real(z, R_P), aimag(z), real(zl, R_P), aimag(zl)
   enddo
   if (vv_report()) close(u)
   print '(A16,3ES14.4)', vv_label(ids(k)), worst_teno, worst_low, maxim
   if (any(vv_label(ids(k)) == [character(16) :: 'teno5', 'teno6', 'teno5_opt', 'teno5_az', 'teno8']) .and. &
       worst_teno > TOL_ADR) &
      call fail(vv_label(ids(k))//': not its linear background for theta <= pi/4')
   if (SCHEMES(ids(k))%s >= 3 .and. worst_low > 1.e-3_R_P) &
      call fail(vv_label(ids(k))//': far from its linear background at low wavenumbers')
   if (vv_label(ids(k)) == 'teno6') then
      if (thanti < 0.4_R_P*PI .or. maxim > 1.e-3_R_P) call fail('teno6: anti-dissipation below 0.4 pi or above 1e-3')
   elseif (vv_label(ids(k)) == 'teno8') then
      if (thanti < 0.3_R_P*PI .or. maxim > 5.e-3_R_P) call fail('teno8: anti-dissipation below 0.3 pi or above 5e-3')
   elseif (vv_label(ids(k)) == 'teno8_a') then
      if (thanti < 0.3_R_P*PI .or. maxim > 1.2e-2_R_P) call fail('teno8_a: anti-dissipation below 0.3 pi or above 1.2e-2')
   elseif (lin_dissipative .and. maxim > TOL_ADR) then
      call fail(vv_label(ids(k))//': anti-dissipative')
   endif
   if (maxim > TOL_ADR) then ! (4) the band is not a growth
      do nn=1, N/2-1
         th = 2._R_P*PI*nn/N
         if (aimag(adr(ids(k), th)) <= TOL_ADR) cycle
         gmax = growth(ids(k), nn)
         print '(A16,A,F6.3,A,F12.9)', vv_label(ids(k)), ': theta/pi ', th/PI, ', 4 periods, max amplitude ', gmax
         if (gmax > 1._R_P + TOL_ADR) call fail(vv_label(ids(k))//': a mode of the anti-dissipative band grows')
      enddo
   endif
enddo
call nonpoly_adr()
if (.not.passed) error stop 'ADR verification failed'
print '(A)', 'seth_vv_adr: ok'

contains
   subroutine nonpoly_adr()
   !< G6 of issue #6: TENO5-ST against linear5, TENO5-RBF at s^2 = 9e-4, 3.6e-3, 2.25e-2, 9e-2 against its linear
   !< five-point row P + s^2 C (s = 3 rbf_dx: the sensor keeps s = s1 on a unit wave).
   real(R_P), parameter :: SS(4) = [0.03_R_P, 0.06_R_P, 0.15_R_P, 0.3_R_P]
   real(R_P)            :: cl(-4:5), dx0, wl, wim, imlin, thb
   integer(I4P)         :: ist, irb, j, n_, ua
   character(1)         :: cj
   ist = by_name('teno5_st') ; irb = by_name('teno5_rbf')
   dx0 = rbf_dx
   cl = 0._R_P ; cl(-2:2) = NP_RBF_P(:,7)
   wl = 0._R_P ; wim = -huge(1._R_P)
   if (vv_report()) call vv_open('adr-teno5_st.dat', ua, &
                                 'theta Re(theta~) Im(theta~) Re(theta~_lin) Im(theta~_lin) (N = 200, lin: linear5)')
   do n_=1, N/2-1
      th = 2._R_P*PI*n_/N
      z = adr(ist, th) ; zl = xit(cl, th)
      if (th <= PI/4._R_P) wl = max(wl, abs(z - zl))
      wim = max(wim, aimag(z))
      if (vv_report()) write(ua, '(5ES24.15)') th, real(z, R_P), aimag(z), real(zl, R_P), aimag(zl)
   enddo
   if (vv_report()) close(ua)
   print '(A,2ES12.4)', 'teno5_st: |theta~ - linear5| for theta <= pi/4, max Im ', wl, wim
   if (wl > TOL_ADR .or. wim > TOL_ADR) call fail('teno5_st: not linear5 at theta <= pi/4, or anti-dissipative')
   do j=1, 4
      rbf_dx = SS(j)/SETH_RBF_ALPHA1
      cl = 0._R_P ; cl(-2:2) = NP_RBF_P(:,7) + SS(j)**2*NP_RBF_C(:,7)
      wl = 0._R_P ; wim = -huge(1._R_P) ; imlin = -huge(1._R_P) ; thb = 0._R_P
      write(cj, '(I1)') j
      if (vv_report()) call vv_open('adr-teno5_rbf-s'//cj//'.dat', ua, &
         'theta Re(theta~) Im(theta~) Re(row) Im(row) (N = 200; lin: the row P + s^2 C, s in the first line)')
      if (vv_report()) write(ua, '(A,ES24.15)') '# s', SS(j)
      do n_=1, N/2-1
         th = 2._R_P*PI*n_/N
         z = adr(irb, th) ; zl = xit(cl, th)
         if (th <= PI/4._R_P) wl = max(wl, abs(z - zl))
         wim = max(wim, aimag(z)) ; imlin = max(imlin, aimag(zl))
         if (aimag(z) > TOL_ADR) thb = th
         if (vv_report()) write(ua, '(5ES24.15)') th, real(z, R_P), aimag(z), real(zl, R_P), aimag(zl)
      enddo
      if (vv_report()) close(ua)
      print '(A,F6.3,A,ES11.3,A,2ES11.3,A,ES11.3,A,F7.4)', 'teno5_rbf s =', SS(j), ': |theta~ - row| (theta <= pi/4)', &
            wl, '; max Im nonlinear, row', wim, imlin, '; predicted peak', 200._R_P/3._R_P*SS(j)**6, &
            '; band edge theta', thb
      if (j <= 2 .and. 200._R_P/3._R_P*SS(j)**6 <= 10._R_P*TOL_ADR) then
         print '(A)', '   (the predicted band is below the round-off of this precision: not checked)'
      elseif (j <= 2) then ! inside the validated range (D7)
         if (abs(thb/(sqrt(30._R_P)*SS(j)) - 1._R_P) > 0.1_R_P) call fail('teno5_rbf: band edge not sqrt(30) s')
         if (wim > 2._R_P*200._R_P/3._R_P*SS(j)**6) call fail('teno5_rbf: band peak above twice the prediction')
      endif
   enddo
   rbf_dx = dx0
   endsubroutine nonpoly_adr

   function by_name(name) result(id)
   character(*), intent(in) :: name
   integer(I4P)             :: id, j
   id = 0
   do j=1, NALL
      if (trim(SCHEMES(j)%name) == name) id = j
   enddo
   if (id == 0) error stop 'seth_vv_adr: scheme not in the registry'
   endfunction by_name

   function growth(id, nm) result(amax)
   !< Largest amplitude of mode nm advected by u_t + u_x = 0 for 4 periods (periodic, SSP-RK3, CFL 0.4, left state).
   integer(I4P), intent(in) :: id, nm
   real(R_P)                :: amax, u0(0:N-1), u1(0:N-1), u2(0:N-1), h, dt
   integer(I4P)             :: j, step
   h = 1._R_P/N ; dt = 0.4_R_P*h
   do j=0, N-1
      u0(j) = cos(2._R_P*PI*nm*j/N)
   enddo
   amax = maxval(abs(u0))
   do step=1, nint(4._R_P/dt)
      u1 = u0 + dt*rhs(id, u0)
      u2 = 0.75_R_P*u0 + 0.25_R_P*(u1 + dt*rhs(id, u1))
      u0 = u0/3._R_P + (2._R_P/3._R_P)*(u2 + dt*rhs(id, u2))
      amax = max(amax, maxval(abs(u0)))
   enddo
   endfunction growth

   function rhs(id, w) result(r)
   !< -d/dx of the left-biased flux of u_t + u_x = 0 (periodic), scheme id.
   integer(I4P), intent(in) :: id
   real(R_P),    intent(in) :: w(0:N-1)
   real(R_P)                :: r(0:N-1), f(0:N-1), v(VLO:VHI), zero(VLO:VHI), rp, rm
   integer(I4P)             :: jj, m
   zero = 0._R_P
   do jj=0, N-1 ! face jj+1/2
      do m=VLO, VHI
         v(m) = w(modulo(jj+m, N))
      enddo
      call reconstruct(id, v, zero, rp, rm)
      f(jj) = rp
   enddo
   do jj=0, N-1
      r(jj) = -(f(jj) - f(modulo(jj-1, N)))*N
   enddo
   endfunction rhs

   function adr(id, th) result(z)
   !< theta~ of scheme id at theta by the single-mode procedure.
   integer(I4P), intent(in) :: id
   real(R_P),    intent(in) :: th
   complex(R_P)             :: z, g
   real(R_P)                :: uu(0:N-1), f(0:N-1), v(VLO:VHI), zero(VLO:VHI), rp, rm
   integer(I4P)             :: j, m
   do j=0, N-1
      uu(j) = cos(th*j)
   enddo
   zero = 0._R_P
   do j=0, N-1 ! face j+1/2
      do m=VLO, VHI
         v(m) = uu(modulo(j+m, N))
      enddo
      call reconstruct(id, v, zero, rp, rm)
      f(j) = rp
   enddo
   g = (0._R_P, 0._R_P)
   do j=0, N-1
      g = g + (f(j) - f(modulo(j-1, N)))*exp(cmplx(0._R_P, -th*j, R_P))
   enddo
   g = g*2._R_P/N
   z = cmplx(0._R_P, -1._R_P, R_P)*g
   endfunction adr

   subroutine background(name, c, dissipative)
   !< Linear face coefficients of the scheme's background (optimal weights), and whether it is dissipative.
   character(*), intent(in)  :: name
   real(R_P),    intent(out) :: c(-4:5)
   logical,      intent(out) :: dissipative
   dissipative = .true.
   select case(name)
   case('weno3_js') ; call js(2, WENO3_JS_A, WENO3_JS_P, c)
   case('weno5_js') ; call js(3, WENO5_JS_A, WENO5_JS_P, c)
   case('weno7_js') ; call js(4, WENO7_JS_A, WENO7_JS_P, c)
   case('weno9_js') ; call js(5, WENO9_JS_A, WENO9_JS_P, c)
   case('wenoz3') ; call dsc(WENOZ3_NC, WENOZ3_MAXW, WENOZ3_CLO, WENOZ3_CW, WENOZ3_P, WENOZ3_D, c)
   case('wenoz5') ; call dsc(WENOZ5_NC, WENOZ5_MAXW, WENOZ5_CLO, WENOZ5_CW, WENOZ5_P, WENOZ5_D, c)
   case('wenoz7') ; call dsc(WENOZ7_NC, WENOZ7_MAXW, WENOZ7_CLO, WENOZ7_CW, WENOZ7_P, WENOZ7_D, c)
   case('wenoz9') ; call dsc(WENOZ9_NC, WENOZ9_MAXW, WENOZ9_CLO, WENOZ9_CW, WENOZ9_P, WENOZ9_D, c)
   case('teno5') ; call dsc(TENO5_NC, TENO5_MAXW, TENO5_CLO, TENO5_CW, TENO5_P, TENO5_D, c)
   case('teno5_dr', 'teno5_thinc') ; call dsc(TENO5_NC, TENO5_MAXW, TENO5_CLO, TENO5_CW, TENO5_P, TENO5_D, c) ! TENO5's
   case('teno6') ; call dsc(TENO6_NC, TENO6_MAXW, TENO6_CLO, TENO6_CW, TENO6_P, TENO6_D, c)
   case('teno6_opt')
      call dsc(TENO6_OPT_NC, TENO6_OPT_MAXW, TENO6_OPT_CLO, TENO6_OPT_CW, TENO6_OPT_P, TENO6_OPT_D, c)
      dissipative = .false.
   case('teno6_neno6')
      call dsc(TENO6_NENO6_NC, TENO6_NENO6_MAXW, TENO6_NENO6_CLO, TENO6_NENO6_CW, TENO6_NENO6_P, TENO6_NENO6_D, c)
      dissipative = .false.
   case('teno5_opt')
      call dsc(TENO5_OPT_NC, TENO5_OPT_MAXW, TENO5_OPT_CLO, TENO5_OPT_CW, TENO5_OPT_P, TENO5_OPT_D, c)
   case('teno5_az')
      call dsc(TENO5_AZ_NC, TENO5_AZ_MAXW, TENO5_AZ_CLO, TENO5_AZ_CW, TENO5_AZ_P, TENO5_AZ_D, c)
   case('teno8')
      call dsc(TENO8_NC, TENO8_MAXW, TENO8_CLO, TENO8_CW, TENO8_P, TENO8_D, c)
   case('teno5_a') ; call dsc(TENO5_NC, TENO5_MAXW, TENO5_CLO, TENO5_CW, TENO5_P, TENO5_D, c)
   case('teno8_a') ; call dsc(TENO8_NC, TENO8_MAXW, TENO8_CLO, TENO8_CW, TENO8_P, TENO8_D, c)
   case('weno4_c') ; c = 0._R_P ; c(CEN4_LO:CEN4_HI) = CEN4_C ; dissipative = .false. ! the central row (phi -> 0)
   case('weno6_c') ; c = 0._R_P ; c(CEN6_LO:CEN6_HI) = CEN6_C ; dissipative = .false.
   case('weno8_c') ; c = 0._R_P ; c(CEN8_LO:CEN8_HI) = CEN8_C ; dissipative = .false.
   case default ; error stop 'background: unknown scheme'
   end select
   endsubroutine background

   subroutine js(s, a, p, c)
   !< Linear face coefficients of WENO-JS, the left-biased state (ADAM's row 2: P(2, cell j, stencil k) on cell k-j).
   integer(I4P), intent(in)  :: s
   real(R_P),    intent(in)  :: a(1:2,0:s-1), p(1:2,0:s-1,0:s-1)
   real(R_P),    intent(out) :: c(-4:5)
   integer(I4P)              :: kk, j
   c = 0._R_P
   do kk=0, s-1
      do j=0, s-1
         c(kk-j) = c(kk-j) + a(2,kk)*p(2,j,kk)
      enddo
   enddo
   endsubroutine js

   subroutine dsc(nc, maxw, clo, cw, cp, d, c)
   !< Linear face coefficients of a WENO-Z/TENO descriptor: sum_k d_k/sum d times candidate k.
   integer(I4P), intent(in)  :: nc, maxw, clo(0:nc-1), cw(0:nc-1)
   real(R_P),    intent(in)  :: cp(0:maxw-1,0:nc-1), d(0:nc-1)
   real(R_P),    intent(out) :: c(-4:5)
   integer(I4P)              :: kk, j
   c = 0._R_P
   do kk=0, nc-1
      do j=0, cw(kk)-1
         c(clo(kk)+j) = c(clo(kk)+j) + d(kk)/sum(d)*cp(j,kk)
      enddo
   enddo
   endsubroutine dsc

   pure function xit(c, xi) result(z)
   !< Modified wavenumber of the finite-difference derivative with face coefficients c (seth_test_spectral).
   real(R_P), intent(in) :: c(-4:5), xi
   complex(R_P)          :: z, s
   integer(I4P)          :: m
   s = (0._R_P, 0._R_P)
   do m=-4, 5
      s = s + c(m)*exp(cmplx(0._R_P, m*xi, R_P))
   enddo
   z = cmplx(0._R_P, -1._R_P, R_P)*s*(1._R_P - exp(cmplx(0._R_P, -xi, R_P)))
   endfunction xit

   subroutine fail(what)
   character(*), intent(in) :: what
   print '(A)', '   FAILED: '//what
   passed = .false.
   endsubroutine fail
endprogram seth_vv_adr
