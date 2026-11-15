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
!< wavelength (about 4 h) are what TENO6 leaves undamped behind the Burgers shock (seth_vv_burgers).
!< Report mode: adr-<scheme>.dat (theta, Re and Im of theta~, of theta~_lin).
use seth
use seth_test_schemes
use seth_weno_js_tables
use seth_z_teno_tables
use seth_vv_cases, only : vv_report, vv_open, vv_schemes, vv_label
use seth_vv_euler, only : PATH_SPLIT
implicit none
integer(I4P), parameter   :: N = 200_I4P
real(R8P),    parameter   :: PI = acos(-1._R8P)
integer(I4P), allocatable :: ids(:)
real(R8P)                 :: c(-4:5), th, worst_teno, worst_low, maxim, thanti
complex(R8P)              :: z, zl
integer(I4P)              :: k, nn, u
logical                   :: passed, lin_dissipative

passed = .true.
call vv_schemes(PATH_SPLIT, .not.vv_report(), ids)
print '(A16,3A14)', 'scheme', 'TENO vs lin', 'low-k vs lin', 'max Im'
do k=1, size(ids)
   call background(vv_label(ids(k)), c, lin_dissipative)
   worst_teno = 0._R8P ; worst_low = 0._R8P ; maxim = -huge(1._R8P) ; thanti = huge(1._R8P)
   if (vv_report()) call vv_open('adr-'//vv_label(ids(k))//'.dat', u, &
                                 'theta Re(theta~) Im(theta~) Re(theta~_lin) Im(theta~_lin) (N = 200, unit amplitude)')
   do nn=1, N/2-1
      th = 2._R8P*PI*nn/N
      z = adr(ids(k), th)
      zl = xit(c, th)
      if (th <= PI/4._R8P) worst_teno = max(worst_teno, abs(z - zl))
      if (th <= 0.3_R8P) worst_low = max(worst_low, abs(z - zl)/th)
      maxim = max(maxim, aimag(z))
      if (aimag(z) > 1.e-12_R8P) thanti = min(thanti, th)
      if (vv_report()) write(u, '(5ES24.15)') th, real(z, R8P), aimag(z), real(zl, R8P), aimag(zl)
   enddo
   if (vv_report()) close(u)
   print '(A16,3ES14.4)', vv_label(ids(k)), worst_teno, worst_low, maxim
   if (any(vv_label(ids(k)) == [character(16) :: 'teno5', 'teno6']) .and. worst_teno > 1.e-12_R8P) &
      call fail(vv_label(ids(k))//': not its linear background for theta <= pi/4')
   if (SCHEMES(ids(k))%s >= 3 .and. worst_low > 1.e-3_R8P) &
      call fail(vv_label(ids(k))//': far from its linear background at low wavenumbers')
   if (vv_label(ids(k)) == 'teno6') then
      if (thanti < 0.4_R8P*PI .or. maxim > 1.e-3_R8P) call fail('teno6: anti-dissipation below 0.4 pi or above 1e-3')
   elseif (lin_dissipative .and. maxim > 1.e-12_R8P) then
      call fail(vv_label(ids(k))//': anti-dissipative')
   endif
enddo
if (.not.passed) error stop 'ADR verification failed'
print '(A)', 'seth_vv_adr: ok'

contains
   function adr(id, th) result(z)
   !< theta~ of scheme id at theta by the single-mode procedure.
   integer(I4P), intent(in) :: id
   real(R8P),    intent(in) :: th
   complex(R8P)             :: z, g
   real(R8P)                :: uu(0:N-1), f(0:N-1), v(VLO:VHI), zero(VLO:VHI), rp, rm
   integer(I4P)             :: j, m
   do j=0, N-1
      uu(j) = cos(th*j)
   enddo
   zero = 0._R8P
   do j=0, N-1 ! face j+1/2
      do m=VLO, VHI
         v(m) = uu(modulo(j+m, N))
      enddo
      call reconstruct(id, v, zero, rp, rm)
      f(j) = rp
   enddo
   g = (0._R8P, 0._R8P)
   do j=0, N-1
      g = g + (f(j) - f(modulo(j-1, N)))*exp(cmplx(0._R8P, -th*j, R8P))
   enddo
   g = g*2._R8P/N
   z = cmplx(0._R8P, -1._R8P, R8P)*g
   endfunction adr

   subroutine background(name, c, dissipative)
   !< Linear face coefficients of the scheme's background (optimal weights), and whether it is dissipative.
   character(*), intent(in)  :: name
   real(R8P),    intent(out) :: c(-4:5)
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
   case('teno6') ; call dsc(TENO6_NC, TENO6_MAXW, TENO6_CLO, TENO6_CW, TENO6_P, TENO6_D, c)
   case('teno6_opt')
      call dsc(TENO6_OPT_NC, TENO6_OPT_MAXW, TENO6_OPT_CLO, TENO6_OPT_CW, TENO6_OPT_P, TENO6_OPT_D, c)
      dissipative = .false.
   case('teno6_neno6')
      call dsc(TENO6_NENO6_NC, TENO6_NENO6_MAXW, TENO6_NENO6_CLO, TENO6_NENO6_CW, TENO6_NENO6_P, TENO6_NENO6_D, c)
      dissipative = .false.
   case default ; error stop 'background: unknown scheme'
   end select
   endsubroutine background

   subroutine js(s, a, p, c)
   !< Linear face coefficients of WENO-JS, the left-biased state (ADAM's row 2: P(2, cell j, stencil k) on cell k-j).
   integer(I4P), intent(in)  :: s
   real(R8P),    intent(in)  :: a(1:2,0:s-1), p(1:2,0:s-1,0:s-1)
   real(R8P),    intent(out) :: c(-4:5)
   integer(I4P)              :: kk, j
   c = 0._R8P
   do kk=0, s-1
      do j=0, s-1
         c(kk-j) = c(kk-j) + a(2,kk)*p(2,j,kk)
      enddo
   enddo
   endsubroutine js

   subroutine dsc(nc, maxw, clo, cw, cp, d, c)
   !< Linear face coefficients of a WENO-Z/TENO descriptor: sum_k d_k/sum d times candidate k.
   integer(I4P), intent(in)  :: nc, maxw, clo(0:nc-1), cw(0:nc-1)
   real(R8P),    intent(in)  :: cp(0:maxw-1,0:nc-1), d(0:nc-1)
   real(R8P),    intent(out) :: c(-4:5)
   integer(I4P)              :: kk, j
   c = 0._R8P
   do kk=0, nc-1
      do j=0, cw(kk)-1
         c(clo(kk)+j) = c(clo(kk)+j) + d(kk)/sum(d)*cp(j,kk)
      enddo
   enddo
   endsubroutine dsc

   pure function xit(c, xi) result(z)
   !< Modified wavenumber of the finite-difference derivative with face coefficients c (seth_test_spectral).
   real(R8P), intent(in) :: c(-4:5), xi
   complex(R8P)          :: z, s
   integer(I4P)          :: m
   s = (0._R8P, 0._R8P)
   do m=-4, 5
      s = s + c(m)*exp(cmplx(0._R8P, m*xi, R8P))
   enddo
   z = cmplx(0._R8P, -1._R8P, R8P)*s*(1._R8P - exp(cmplx(0._R8P, -xi, R8P)))
   endfunction xit

   subroutine fail(what)
   character(*), intent(in) :: what
   print '(A)', '   FAILED: '//what
   passed = .false.
   endsubroutine fail
endprogram seth_vv_adr
