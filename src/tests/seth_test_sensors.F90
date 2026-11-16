!< SETH test: the smoothness sensors (seth_sensor_teno, _ren, _lcr, _tau) against Li, Chen & Ren's static tests.
program seth_test_sensors
!< SETH test: the smoothness sensors (seth_sensor_teno, _ren, _lcr, _tau) against Li, Chen & Ren's static tests.
!<
!< 1. Monochromatic waves (Li, Chen & Ren 2022, Fig. 9): f_j = sin(k j), j = 0 .. 191, unit amplitude (their
!<    normalisation by max(|f|, 1) leaves it unchanged); the fraction of interfaces judged smooth (s = 1) for k = 0.01
!<    .. 3.00. Printed: the largest k below which every interface is smooth for every smaller k. Asserted (the paper,
!<    within 0.05): Fu's TENO detector 0.64, Li, Chen & Ren's 2.06; Ren's fails "at very low wavenumbers" (below 0.3).
!< 2. Test (b) of their section 2.2: sin(12 pi x) - 2 for x < 0, sin(24 pi x + pi/2) + 2 for x >= 0, 192 cells on
!<    [-1, 1]: Li, Chen & Ren's detector flags the jump and nothing farther than 3 h from it; Fu's flags the faster wave
!<    (k = pi/4 > 0.64) too; Ren's flags the extrema. Asserted for the first; the others printed.
!< 3. Properties: s in [0, 1]; constant data gives s = 1; a unit step anywhere in the stencil gives s < 1 (TENO, Ren,
!<    Li-Chen-Ren on cells -2..3; Wan & Xia's on its cells -2..2).
use seth
use seth_vv_cases, only : vv_report, vv_open
implicit none
integer(I4P), parameter :: N = 192_I4P, NS = 4_I4P
real(R8P),    parameter :: PI = acos(-1._R8P)
character(8), parameter :: NAMES(NS) = ['teno    ', 'ren     ', 'lcr     ', 'tau     ']
real(R8P)               :: f(0:N-1), v(-2:3), s, k, frac(NS), kok(NS), x, sb(NS)
integer                 :: uk, ub
integer(I4P)            :: is, j, ik, nsm, nf, m, jump
logical                 :: passed, ok(NS), flag

passed = .true.
! 1. monochromatic waves
kok = 0._R8P ; ok = .true.
if (vv_report()) call vv_open('sensors-k.dat', uk, 'k fraction-smooth: teno ren lcr tau (N = 192, unit sine)')
do ik=1, 300
   k = 0.01_R8P*ik
   do j=0, N-1
      f(j) = sin(k*j)
   enddo
   do is=1, NS
      nsm = 0 ; nf = 0
      do j=2, N-4
         v = f(j-2:j+3)
         call sensor(is, v, s)
         nf = nf + 1
         if (s >= 1._R8P) nsm = nsm + 1
      enddo
      frac(is) = real(nsm, R8P)/nf
      if (ok(is) .and. nsm == nf) then
         kok(is) = k
      else
         ok(is) = .false.
      endif
   enddo
   if (vv_report()) write(uk, '(5ES16.7)') k, frac
enddo
if (vv_report()) close(uk)
print '(A)', '1. largest wavenumber below which every interface of a sine is judged smooth:'
do is=1, NS
   print '(A10,F6.2)', NAMES(is), kok(is)
enddo
if (abs(kok(1) - 0.64_R8P) > 0.05_R8P) then ; print '(A)', '   ^ Fu: not 0.64 +- 0.05' ; passed = .false. ; endif
if (abs(kok(3) - 2.06_R8P) > 0.05_R8P) then ; print '(A)', '   ^ Li, Chen & Ren: not 2.06 +- 0.05' ; passed = .false. ; endif
if (kok(2) > 0.3_R8P) then ; print '(A)', '   ^ Ren: smooth beyond 0.3' ; passed = .false. ; endif
! 2. test (b)
do j=0, N-1
   x = -1._R8P + (j + 0.5_R8P)*2._R8P/N
   if (x < 0._R8P) then
      f(j) = sin(12._R8P*PI*x) - 2._R8P
   else
      f(j) = sin(24._R8P*PI*x + PI/2._R8P) + 2._R8P
   endif
enddo
print '(A)', '2. test (b): interfaces judged non-smooth (s < 1), total and farther than 3 h from the jump:'
if (vv_report()) then
   call vv_open('sensors-b.dat', ub, 'x(face) f(j) s: teno ren lcr tau (test (b) of Li, Chen & Ren 2022)')
   do j=2, N-4
      do is=1, NS
         v = f(j-2:j+3) ; call sensor(is, v, sb(is))
      enddo
      write(ub, '(6ES16.7)') -1._R8P + (j + 1._R8P)*2._R8P/N, f(j), sb
   enddo
   close(ub)
endif
do is=1, NS
   nf = 0 ; nsm = 0 ; flag = .false.
   do j=2, N-4
      v = f(j-2:j+3)
      call sensor(is, v, s)
      x = -1._R8P + (j + 1._R8P)*2._R8P/N ! the face j+1/2
      if (s < 1._R8P) then
         nf = nf + 1
         if (abs(x) > 3._R8P*2._R8P/N + 1.e-12_R8P) nsm = nsm + 1
         if (abs(x) < 1.e-12_R8P) flag = .true.
      endif
   enddo
   print '(A10,2I6,A,L2)', NAMES(is), nf, nsm, '   the jump flagged:', flag
   if (is == 3 .and. (nsm > 0 .or. .not.flag)) then
      print '(A)', '   ^ Li, Chen & Ren: a smooth face flagged, or the jump missed' ; passed = .false.
   endif
enddo
! 3. properties
do is=1, NS
   v = 3.7_R8P
   call sensor(is, v, s)
   if (s /= 1._R8P) then ; print '(A)', NAMES(is)//': constant data not smooth' ; passed = .false. ; endif
   do jump=merge(-1, -1, is /= 4), merge(3, 2, is /= 4)
      do m=-2, 3
         v(m) = merge(1._R8P, 0._R8P, m >= jump)
      enddo
      call sensor(is, v, s)
      if (.not.(s >= 0._R8P .and. s < 1._R8P)) then
         print '(A,I3,ES12.4)', NAMES(is)//': a unit step not flagged, jump at ', jump, s ; passed = .false.
      endif
   enddo
enddo
if (.not.passed) error stop 'sensor test failed'
print '(A)', 'seth_test_sensors: ok'

contains
   subroutine sensor(is, v, s)
   integer(I4P), intent(in)  :: is
   real(R8P),    intent(in)  :: v(-2:3)
   real(R8P),    intent(out) :: s
   real(R8P)                 :: rmu
   rmu = 1._R8P/max(maxval(abs(v)), 1._R8P) ! Li, Chen & Ren's normalisation (Eq. (61))
   select case(is)
   case(1) ; call seth_sensor_teno(v, rmu, SETH_SENSOR_TENO_CT, s)
   case(2) ; call seth_sensor_ren(v, rmu, SETH_SENSOR_REN_PSIC, s)
   case(3) ; call seth_sensor_lcr(v, rmu, SETH_SENSOR_LCR_PSIC, s)
   case(4) ; call seth_sensor_tau(v, rmu, SETH_SENSOR_TAU_ETA, s)
   endselect
   if (.not.(s >= 0._R8P .and. s <= 1._R8P)) error stop 'sensor outside [0, 1]'
   endsubroutine sensor
endprogram seth_test_sensors
