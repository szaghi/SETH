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
!< 4. Critical points (roadmap Phase 5, item (c)): cos(pi x) (order 1) and sin(pi x)^3 (order 2) around x = 0, 2 .. 1024
!<    cells per unit length: is a face of the five nearest the point flagged? Asserted: Li, Chen & Ren's never; no
!<    sensor at 256 cells per unit or more. Printed: the resolution from which each sensor leaves the point alone.
!< 5. Scaled and offset data (item (e)): test (b) as a f + b, a = 1e-6 .. 1e6, b = 0, 100, with three normalisations
!<    (rmu = 1; the paper's per-stencil 1/max(|v|, 1); the field's range, 1/(max f - min f)): the fraction of faces
!<    whose decision changes. Asserted: the range normalisation changes at most one face (invariant up to round-off);
!<    without normalisation the absolute regularisers change decisions at a = 1e-6 (the control).
!< Report mode: sensors-cp.dat and sensors-scaling.dat too.
use seth
use seth_test_parameters, only : SENSOR_SCALES, TOL_X
use seth_vv_cases, only : vv_report, vv_open
implicit none
integer(I4P), parameter :: N = 192_I4P, NS = 4_I4P
real(R_P),    parameter :: PI = acos(-1._R_P)
character(8), parameter :: NAMES(NS) = ['teno    ', 'ren     ', 'lcr     ', 'tau     ']
real(R_P)               :: f(0:N-1), v(-2:3), s, k, frac(NS), kok(NS), x, sb(NS)
integer                 :: uk, ub
integer(I4P)            :: is, j, ik, nsm, nf, m, jump
logical                 :: passed, ok(NS), flag

passed = .true.
! 1. monochromatic waves
kok = 0._R_P ; ok = .true.
if (vv_report()) call vv_open('sensors-k.dat', uk, 'k fraction-smooth: teno ren lcr tau (N = 192, unit sine)')
do ik=1, 300
   k = 0.01_R_P*ik
   do j=0, N-1
      f(j) = sin(k*j)
   enddo
   do is=1, NS
      nsm = 0 ; nf = 0
      do j=2, N-4
         v = f(j-2:j+3)
         call sensor(is, v, s)
         nf = nf + 1
         if (s >= 1._R_P) nsm = nsm + 1
      enddo
      frac(is) = real(nsm, R_P)/nf
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
if (abs(kok(1) - 0.64_R_P) > 0.05_R_P) then ; print '(A)', '   ^ Fu: not 0.64 +- 0.05' ; passed = .false. ; endif
if (abs(kok(3) - 2.06_R_P) > 0.05_R_P) then ; print '(A)', '   ^ Li, Chen & Ren: not 2.06 +- 0.05' ; passed = .false. ; endif
if (kok(2) > 0.3_R_P) then ; print '(A)', '   ^ Ren: smooth beyond 0.3' ; passed = .false. ; endif
! 2. test (b)
do j=0, N-1
   x = -1._R_P + (j + 0.5_R_P)*2._R_P/N
   if (x < 0._R_P) then
      f(j) = sin(12._R_P*PI*x) - 2._R_P
   else
      f(j) = sin(24._R_P*PI*x + PI/2._R_P) + 2._R_P
   endif
enddo
print '(A)', '2. test (b): interfaces judged non-smooth (s < 1), total and farther than 3 h from the jump:'
if (vv_report()) then
   call vv_open('sensors-b.dat', ub, 'x(face) f(j) s: teno ren lcr tau (test (b) of Li, Chen & Ren 2022)')
   do j=2, N-4
      do is=1, NS
         v = f(j-2:j+3) ; call sensor(is, v, sb(is))
      enddo
      write(ub, '(6ES16.7)') -1._R_P + (j + 1._R_P)*2._R_P/N, f(j), sb
   enddo
   close(ub)
endif
do is=1, NS
   nf = 0 ; nsm = 0 ; flag = .false.
   do j=2, N-4
      v = f(j-2:j+3)
      call sensor(is, v, s)
      x = -1._R_P + (j + 1._R_P)*2._R_P/N ! the face j+1/2
      if (s < 1._R_P) then
         nf = nf + 1
         if (abs(x) > 3._R_P*2._R_P/N + TOL_X) nsm = nsm + 1
         if (abs(x) < TOL_X) flag = .true.
      endif
   enddo
   print '(A10,2I6,A,L2)', NAMES(is), nf, nsm, '   the jump flagged:', flag
   if (is == 3 .and. (nsm > 0 .or. .not.flag)) then
      print '(A)', '   ^ Li, Chen & Ren: a smooth face flagged, or the jump missed' ; passed = .false.
   endif
enddo
! 3. properties
do is=1, NS
   v = 3.7_R_P
   call sensor(is, v, s)
   if (s /= 1._R_P) then ; print '(A)', NAMES(is)//': constant data not smooth' ; passed = .false. ; endif
   do jump=merge(-1, -1, is /= 4), merge(3, 2, is /= 4)
      do m=-2, 3
         v(m) = merge(1._R_P, 0._R_P, m >= jump)
      enddo
      call sensor(is, v, s)
      if (.not.(s >= 0._R_P .and. s < 1._R_P)) then
         print '(A,I3,ES12.4)', NAMES(is)//': a unit step not flagged, jump at ', jump, s ; passed = .false.
      endif
   enddo
enddo
! 4. critical points under refinement
call critical_points()
! 5. scaled and offset data
call scaling()
if (.not.passed) error stop 'sensor test failed'
print '(A)', 'seth_test_sensors: ok'

contains
   subroutine critical_points()
   !< f = cos(pi x) (critical point of order 1 at x = 0) and f = sin(pi x)^3 (order 2), cells of width h = 1/m centred
   !< at (j + 1/2) h, the five faces j h, j = -2 .. 2, around x = 0: flagged (some s < 1) or not, for m = 2 .. 1024
   !< cells per unit length (a period of cos(pi x) is 2 m cells). Printed per sensor: the flags, and the smallest m
   !< from which on no face near the point is flagged.
   integer(I4P), parameter :: NM = 10
   integer(I4P)            :: ncp, im, m, jf, is2, mres(NS)
   real(R_P)               :: hh, vv(-2:3), ss
   logical                 :: fl(NS,NM)
   integer                 :: uc
   if (vv_report()) call vv_open('sensors-cp.dat', uc, 'n_cp m(cells per unit) flagged: teno ren lcr tau '// &
                                 '(1 = some face of the five around the critical point has s < 1)')
   do ncp=1, 2
      do im=1, NM
         m = 2**im ; hh = 1._R_P/m
         do is2=1, NS
            fl(is2,im) = .false.
            do jf=-2, 2 ! face jf h: stencil cells jf-3 .. jf+2, i.e. v(-2:3) = f(cells jf-3+1 .. ) with cell 0 = jf-1
               do j=-2, 3
                  x = (jf - 1 + j + 0.5_R_P)*hh
                  vv(j) = merge(cos(PI*x), sin(PI*x)**3, ncp == 1)
               enddo
               call sensor(is2, vv, ss)
               if (ss < 1._R_P) fl(is2,im) = .true.
            enddo
         enddo
         if (vv_report()) write(uc, '(2I6,4I3)') ncp, m, merge(1, 0, fl(:,im))
      enddo
      do is2=1, NS
         mres(is2) = -1
         do im=NM, 1, -1
            if (fl(is2,im)) exit
            mres(is2) = 2**im
         enddo
      enddo
      print '(A,I1,A)', '4. critical point of order ', ncp, ': flagged at m = 2, 4, .. 1024 cells per unit; '// &
                        'smooth from m ='
      do is2=1, NS
         print '(A10,1X,10L2,I7)', NAMES(is2), fl(is2,:), mres(is2)
      enddo
      if (any(fl(3,:))) then ; print '(A)', '   ^ Li, Chen & Ren: a critical point flagged' ; passed = .false. ; endif
      if (any(mres < 0) .or. any(fl(:,NM-2:NM))) then
         print '(A)', '   ^ a sensor still flags the critical point at 256 cells per unit or more' ; passed = .false.
      endif
   enddo
   if (vv_report()) close(uc)
   endsubroutine critical_points

   subroutine scaling()
   !< Test (b) transformed to a f + b, a = 1e-6 .. 1e6, b = 0, 100: the fraction of faces whose decision (s = 1 or
   !< not) differs from that on the untransformed data, with three normalisations of the stencil: none (rmu = 1), the
   !< paper's per-stencil 1/max(|v|_max, 1), and the field's range 1/(max f - min f) (a global descaler).
   real(R_P), parameter :: AS(5) = SENSOR_SCALES, BS(2) = [0._R_P, 100._R_P] ! R8P: 1e-6, 1e-3, 1, 1e3, 1e6
   character(8), parameter :: NORMS(3) = ['none    ', 'paper   ', 'range   ']
   real(R_P)            :: g(0:N-1), vv(-2:3), r, ss
   logical              :: base(NS,3,2:N-4), dec
   integer(I4P)         :: ia, ib, in, is2, jj, nmis(NS,3)
   integer              :: ul
   do jj=0, N-1
      x = -1._R_P + (jj + 0.5_R_P)*2._R_P/N
      f(jj) = merge(sin(12._R_P*PI*x) - 2._R_P, sin(24._R_P*PI*x + PI/2._R_P) + 2._R_P, x < 0._R_P)
   enddo
   if (vv_report()) call vv_open('sensors-scaling.dat', ul, 'a b normalisation mismatched-fraction: teno ren lcr tau')
   print '(A)', '5. test (b) as a f + b: faces whose decision changes, fraction (none / paper / range normalisation)'
   r = maxval(f) - minval(f)
   do in=1, 3
      do jj=2, N-4
         vv = f(jj-2:jj+3)
         do is2=1, NS
            call sensor_n(is2, in, vv, r, ss)
            base(is2,in,jj) = ss >= 1._R_P
         enddo
      enddo
   enddo
   do ib=1, 2
      do ia=1, 5
         g = AS(ia)*f + BS(ib)
         r = maxval(g) - minval(g)
         nmis = 0
         do in=1, 3
            do jj=2, N-4
               vv = g(jj-2:jj+3)
               do is2=1, NS
                  call sensor_n(is2, in, vv, r, ss)
                  dec = ss >= 1._R_P
                  if (dec .neqv. base(is2,in,jj)) nmis(is2,in) = nmis(is2,in) + 1
               enddo
            enddo
         enddo
         print '(A,ES8.1,A,F6.1,3(A,4F6.3))', '   a ', AS(ia), ' b ', BS(ib), ' | none', real(nmis(:,1))/(N-5), &
               ' | paper', real(nmis(:,2))/(N-5), ' | range', real(nmis(:,3))/(N-5)
         if (any(nmis(:,3) > 1)) then
            print '(A)', '   ^ the range normalisation is not scale and offset invariant' ; passed = .false.
         endif
         if (ia == 1 .and. ib == 1 .and. .not.(nmis(1,1) > 0 .and. nmis(2,1) > 0)) then
            print '(A)', '   ^ the absolute regularisers did not act at a = 1e-6 (the control)' ; passed = .false.
         endif
         if (vv_report()) then
            do in=1, 3
               write(ul, '(2ES12.4,1X,A,4F9.5)') AS(ia), BS(ib), trim(NORMS(in)), real(nmis(:,in))/(N-5)
            enddo
         endif
      enddo
   enddo
   if (vv_report()) close(ul)
   endsubroutine scaling

   subroutine sensor_n(is2, in, v, r, s)
   !< Sensor is2 with normalisation in (1: none, 2: the paper's per-stencil max(|v|, 1), 3: 1/r, r the field's range).
   integer(I4P), intent(in)  :: is2, in
   real(R_P),    intent(in)  :: v(-2:3), r
   real(R_P),    intent(out) :: s
   real(R_P)                 :: rmu
   select case(in)
   case(1) ; rmu = 1._R_P
   case(2) ; rmu = 1._R_P/max(maxval(abs(v)), 1._R_P)
   case default ; rmu = 1._R_P/r
   endselect
   select case(is2)
   case(1) ; call seth_sensor_teno(v, rmu, SETH_SENSOR_TENO_CT, s)
   case(2) ; call seth_sensor_ren(v, rmu, SETH_SENSOR_REN_PSIC, s)
   case(3) ; call seth_sensor_lcr(v, rmu, SETH_SENSOR_LCR_PSIC, s)
   case(4) ; call seth_sensor_tau(v, rmu, SETH_SENSOR_TAU_ETA, s)
   endselect
   endsubroutine sensor_n

   subroutine sensor(is, v, s)
   integer(I4P), intent(in)  :: is
   real(R_P),    intent(in)  :: v(-2:3)
   real(R_P),    intent(out) :: s
   real(R_P)                 :: rmu
   rmu = 1._R_P/max(maxval(abs(v)), 1._R_P) ! Li, Chen & Ren's normalisation (Eq. (61))
   select case(is)
   case(1) ; call seth_sensor_teno(v, rmu, SETH_SENSOR_TENO_CT, s)
   case(2) ; call seth_sensor_ren(v, rmu, SETH_SENSOR_REN_PSIC, s)
   case(3) ; call seth_sensor_lcr(v, rmu, SETH_SENSOR_LCR_PSIC, s)
   case(4) ; call seth_sensor_tau(v, rmu, SETH_SENSOR_TAU_ETA, s)
   endselect
   if (.not.(s >= 0._R_P .and. s <= 1._R_P)) error stop 'sensor outside [0, 1]'
   endsubroutine sensor
endprogram seth_test_sensors
