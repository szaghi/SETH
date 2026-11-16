!< SETH test: the spectral properties of every scheme's linear background, from its compile-time constants.
program seth_test_spectral
!< SETH test: the spectral properties of every scheme's linear background, from its compile-time constants.
!<
!< The finite-difference derivative of the linear scheme, D u_i = (F_{i+1/2} - F_{i-1/2})/h with F_{i+1/2} =
!< sum_m c_m u_{i+m}, applied to u_j = exp(i j xi) gives (i/h) xi~(xi), xi~ = -i sum_m c_m e^{i m xi} (1 - e^{-i xi})
!< (Martin et al. 2006; Pirozzoli 2002; llm-wiki high-order-capturing-schemes `definitions/modified-wavenumber.md`).
!< Re xi~ is the dispersion, Im xi~ the dissipation (negative: dissipative, for right-running waves). Asserted:
!<
!< 1. No anti-dissipation: Im xi~ <= 1e-14 on (0, pi], for the max-order backgrounds (WENO-JS = WENO-Z, TENO5,
!<    TENO6); the optimised TENO6 weights are printed (bandwidth-optimised schemes can be anti-dissipative at low
!<    wavenumbers, Fu et al. 2016).
!< 2. Low-wavenumber accuracy: |xi~ - xi| ~ xi^(o+1) for a scheme of order o: the slope between xi = 0.2 and 0.1 is at
!<    least o + 1 - 0.3.
!< 3. Resolution: the largest xi/pi with |Re xi~ - xi|/xi below 1% is 0.35 for the 5th-order and 0.42 for the
!<    7th-order WENO background (Martin et al. 2006, via llm-wiki `methods/bandwidth-optimized-weno.md`), within 0.01.
!<
!< The interpolation tables are not a derivative operator and are left out.
use seth
use seth_test_parameters, only : TOL_ANTIDISSIPATION
use seth_weno_js_tables
use seth_z_teno_tables
use seth_vv_cases, only : vv_report, vv_open
implicit none
integer(I4P), parameter :: NX = 20000_I4P
real(R_P),    parameter :: PI = acos(-1._R_P)
real(R_P)               :: c(-4:5), imax, slope, eff
integer(I4P)            :: k
logical                 :: passed

passed = .true.
print '(A16,A6,3A14)', 'background', 'order', 'max Im xi~', 'slope at 0', 'resolution'
do k=1, 12
   select case(k)
   case(1) ; call js(1, WENO1_JS_A, WENO1_JS_P, c) ; call check('weno1_js', 1, .true., -1._R_P)
   case(2) ; call js(2, WENO3_JS_A, WENO3_JS_P, c) ; call check('weno3_js', 3, .true., -1._R_P)
   case(3) ; call js(3, WENO5_JS_A, WENO5_JS_P, c) ; call check('weno5_js', 5, .true., 0.35_R_P)
   case(4) ; call js(4, WENO7_JS_A, WENO7_JS_P, c) ; call check('weno7_js', 7, .true., 0.42_R_P)
   case(5) ; call js(5, WENO9_JS_A, WENO9_JS_P, c) ; call check('weno9_js', 9, .true., -1._R_P)
   case(6) ; call dsc(TENO5_NC, TENO5_MAXW, TENO5_CLO, TENO5_CW, TENO5_P, TENO5_D, c)
             call check('teno5', 5, .true., 0.35_R_P)
   case(7) ; call dsc(TENO6_NC, TENO6_MAXW, TENO6_CLO, TENO6_CW, TENO6_P, TENO6_D, c)
             call check('teno6', 6, .true., -1._R_P)
   case(8) ; call dsc(TENO6_OPT_NC, TENO6_OPT_MAXW, TENO6_OPT_CLO, TENO6_OPT_CW, TENO6_OPT_P, TENO6_OPT_D, c)
             call check('teno6_opt', 5, .false., -1._R_P)
   case(9) ; call dsc(TENO6_NENO6_NC, TENO6_NENO6_MAXW, TENO6_NENO6_CLO, TENO6_NENO6_CW, TENO6_NENO6_P, &
                      TENO6_NENO6_D, c)
             call check('teno6_neno6', 4, .false., -1._R_P)
   case(10); call dsc(WENOZ7_NC, WENOZ7_MAXW, WENOZ7_CLO, WENOZ7_CW, WENOZ7_P, WENOZ7_D, c)
             call check('wenoz7 (= js7)', 7, .true., 0.42_R_P)
   case(11); call dsc(TENO5_OPT_NC, TENO5_OPT_MAXW, TENO5_OPT_CLO, TENO5_OPT_CW, TENO5_OPT_P, TENO5_OPT_D, c)
             call check('teno5_opt', 4, .true., -1._R_P)
   case(12); call dsc(TENO8_NC, TENO8_MAXW, TENO8_CLO, TENO8_CW, TENO8_P, TENO8_D, c)
             call check('teno8', 6, .true., -1._R_P)
   endselect
enddo
if (.not.passed) error stop 'spectral test failed'
print '(A)', 'seth_test_spectral: ok'

contains
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
   !< Modified wavenumber of the finite-difference derivative with face coefficients c.
   real(R_P), intent(in) :: c(-4:5), xi
   complex(R_P)          :: z, s
   integer(I4P)          :: m
   s = (0._R_P, 0._R_P)
   do m=-4, 5
      s = s + c(m)*exp(cmplx(0._R_P, m*xi, R_P))
   enddo
   z = cmplx(0._R_P, -1._R_P, R_P)*s*(1._R_P - exp(cmplx(0._R_P, -xi, R_P)))
   endfunction xit

   subroutine check(name, o, dissipative, eff_ref)
   character(*), intent(in) :: name
   integer(I4P), intent(in) :: o
   logical,      intent(in) :: dissipative
   real(R_P),    intent(in) :: eff_ref
   complex(R_P)             :: z
   real(R_P)                :: xi
   integer(I4P)             :: i
   integer                  :: u
   imax = -huge(1._R_P) ; eff = -1._R_P
   if (vv_report()) then
      call vv_open('spectral-'//name(1:index(name//' ', ' ')-1)//'.dat', u, 'xi Re(xi~) Im(xi~) (linear background)')
      do i=1, NX, 100
         xi = PI*i/NX ; z = xit(c, xi) ; write(u, '(3ES24.15)') xi, real(z, R_P), aimag(z)
      enddo
      close(u)
   endif
   do i=1, NX
      xi = PI*i/NX
      z = xit(c, xi)
      imax = max(imax, aimag(z))
      if (eff < 0._R_P .and. abs(real(z, R_P) - xi)/xi >= 0.01_R_P) eff = real(i-1, R_P)/NX
   enddo
   slope = log(abs(xit(c, 0.2_R_P) - 0.2_R_P)/abs(xit(c, 0.1_R_P) - 0.1_R_P))/log(2._R_P)
   print '(A16,I6,ES14.2,2F14.3)', name, o, imax, slope, eff
   if (dissipative .and. imax > TOL_ANTIDISSIPATION) then ; print '(A)', '   ^ anti-dissipative' ; passed = .false. ; endif
   if (slope < o + 1 - 0.3_R_P) then ; print '(A)', '   ^ low-wavenumber error slope below o + 1' ; passed = .false. ; endif
   if (eff_ref > 0._R_P .and. abs(eff - eff_ref) > 0.01_R_P) then
      print '(A,F6.3)', '   ^ resolution differs from Martin et al. 2006: ', eff_ref ; passed = .false.
   endif
   endsubroutine check
endprogram seth_test_spectral
