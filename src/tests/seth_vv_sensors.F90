!< SETH verification: the four sensors in flows with shocks, as diagnostics and as the switch of a cost hybrid.
program seth_vv_sensors
!< SETH verification: the four sensors in flows with shocks, as diagnostics and as the switch of a cost hybrid.
!<
!< Roadmap issue #4, Phase 5, item (d). Two flows of the 1-D Euler harness (seth_vv_euler, splitting path):
!< Shu-Osher (400 cells on [-5, 5], t = 1.8, transmissive) and the Woodward-Colella blast waves (400 cells, t = 0.038,
!< reflecting walls, the positivity limiter on), against WENO5-JS on 4 times the cells (16 in the report mode).
!<
!< 1. Diagnostic: the pure TENO5 run, 60 frames; at each frame every sensor reads the density of every face (as the
!<    harness's hybrid does: cells i-2..i+3, normalised by max(|rho|, 1), the published thresholds, no dilation).
!<    Printed: the fraction of faces flagged, averaged over the frames. Report: the flag field in x and t per sensor
!<    (sensors-<case>-flags.dat).
!< 2. The cost hybrid (HYB_COST: the linear background where the sensor, dilated by two faces, says smooth; the
!<    characteristic nonlinear scheme elsewhere) of TENO5 and TENO6 driven by each sensor: L1(rho) against the
!<    reference, its ratio to the pure scheme, the fraction of faces flagged (the nonlinear share), positivity.
!< Asserted: every run positive (the blast with the limiter); every hybrid within 3 times the pure scheme's L1 (a
!< sensor that misses a shock blows up, it does not lose a few per cent); the diagnostic fractions ordered as the
!< static tests predict on Shu-Osher (Ren's flags the most, then Fu's TENO detector, then Li, Chen & Ren's).
!< Report: sensors-flows.dat (case scheme sensor L1 ratio flagged-fraction positive).
use seth
use seth_test_schemes
use seth_vv_euler
use seth_vv_cases
implicit none
integer(I4P),  parameter :: NSN = 4, NF = 60
character(8),  parameter :: SNAMES(NSN) = ['teno    ', 'ren     ', 'lcr     ', 'tau     ']
integer(I4P),  parameter :: SIDS(NSN) = [SENS_TENO, SENS_REN, SENS_LCR, SENS_TAU]
character(16), parameter :: BASES(2) = [character(16) :: 'teno5', 'teno6']
integer(I4P)             :: us
logical                  :: passed
passed = .true.
if (vv_report()) call vv_open('sensors-flows.dat', us, 'case scheme sensor L1(rho) ratio-to-pure flagged-fraction positive')
call study('shu-osher', init_shu_osher, -5._R_P, 5._R_P, BC_TRANSMISSIVE, 1.8_R_P, .false.)
call study('blast', init_blast, 0._R_P, 1._R_P, BC_REFLECTIVE, 0.038_R_P, .true.)
hybrid = HYB_NONE ; hyb_sensor = SENS_LCR ; limiter = .false.
if (vv_report()) close(us)
if (.not.passed) error stop 'sensor flow verification failed'
print '(A)', 'seth_vv_sensors: ok'

contains
   subroutine study(name, init, xl, xr, bc, tend, lim)
   character(*), intent(in) :: name
   interface
      pure function init(x) result(w)
      import :: R_P
      real(R_P), intent(in) :: x
      real(R_P)             :: w(3)
      endfunction init
   endinterface
   real(R_P),    intent(in) :: xl, xr, tend
   integer(I4P), intent(in) :: bc
   logical,      intent(in) :: lim
   integer(I4P), parameter  :: N = 400
   real(R_P), allocatable   :: qref(:,:)
   real(R_P)                :: tf(NF), fr(3,N,NF), w0(3,N), q(3,1-NG:N+NG), h, v(-2:3), s, mean(NSN), l1(1), l1p, frac
   logical                  :: pos(1), flag(N-5,NF,NSN)
   integer(I4P)             :: ids(1), k, i, f, is, ib, uf, nref
   h = (xr - xl)/N
   nref = merge(16, 4, vv_report())*N
   limiter = lim
   ! 1. diagnostic flags on the pure TENO5 run
   hybrid = HYB_NONE
   ids(1) = scheme_index('teno5')
   do i=1, N
      w0(:,i) = init(xl + (i - 0.5_R_P)*h)
   enddo
   tf = [(tend*k/real(NF, R_P), k=1, NF)]
   call euler_run(ids(1), PATH_SPLIT, w0, N, h, bc, tend, 0.5_R_P, q, pos(1), fr, tf, NF)
   if (.not.pos(1)) call fail(name//': the pure TENO5 run lost positivity')
   mean = 0._R_P
   do f=1, NF
      do i=3, N-3 ! face i+1/2, cells i-2 .. i+3
         v = fr(1,i-2:i+3,f)
         do is=1, NSN
            call sensor(is, v, s)
            flag(i-2,f,is) = s < 1._R_P
         enddo
      enddo
      do is=1, NSN
         mean(is) = mean(is) + count(flag(:,f,is))/real((N-5)*NF, R_P)
      enddo
   enddo
   print '(A)', '1. '//name//', pure TENO5: fraction of faces flagged, mean over 60 frames (teno ren lcr tau)'
   print '(10X,4F8.4)', mean
   if (name == 'shu-osher' .and. .not.(mean(2) > mean(1) .and. mean(1) > mean(3))) &
      call fail(name//': the flagged fractions are not ordered ren > teno > lcr')
   if (vv_report()) then
      call vv_open('sensors-'//name//'-flags.dat', uf, 'flags of '//name//' (pure TENO5): line 1 the face x; then '// &
                   'per frame and sensor (teno ren lcr tau): t sensor flags(1 = s < 1)')
      write(uf, '(*(ES14.6))') (xl + i*h, i=3, N-3)
      do f=1, NF
         do is=1, NSN
            write(uf, '(ES14.6,I3,*(I2))') tf(f), is, merge(1, 0, flag(:,f,is))
         enddo
      enddo
      close(uf)
   endif
   ! 2. the cost hybrid driven by each sensor
   print '(A)', '2. '//name//', cost hybrid: L1(rho), ratio to the pure scheme, flagged fraction, positive'
   do ib=1, 2
      ids(1) = scheme_index(BASES(ib))
      hybrid = HYB_NONE
      call reference_problem('sens-'//name, init, xl, xr, bc, tend, 0.5_R_P, N, nref, 3_I4P, ids, PATH_SPLIT, l1, pos, &
                             files=.false., qref_in=qref)
      l1p = l1(1)
      print '(A10,A6,ES12.4,8X,8X,L3)', trim(BASES(ib)), 'pure', l1p, pos(1)
      if (.not.pos(1)) call fail(name//' '//trim(BASES(ib))//': the pure scheme lost positivity')
      if (vv_report()) write(us, '(A,1X,A,1X,A,ES24.15,2F10.5,1X,L1)') name, trim(BASES(ib)), 'pure', l1p, 1._R_P, &
                                                                       1._R_P, pos(1)
      do is=1, NSN
         hybrid = HYB_COST ; hyb_sensor = SIDS(is) ; nflagged = 0 ; nhybfaces = 0
         call reference_problem('sens-'//name, init, xl, xr, bc, tend, 0.5_R_P, N, nref, 3_I4P, ids, PATH_SPLIT, l1, &
                                pos, files=.false., qref_in=qref)
         frac = real(nflagged, R_P)/max(real(nhybfaces, R_P), 1._R_P)
         print '(A10,A6,ES12.4,F8.3,F8.3,L3)', trim(BASES(ib)), trim(SNAMES(is)), l1(1), l1(1)/l1p, frac, pos(1)
         if (vv_report()) write(us, '(A,1X,A,1X,A,ES24.15,2F10.5,1X,L1)') name, trim(BASES(ib)), trim(SNAMES(is)), &
                                                                          l1(1), l1(1)/l1p, frac, pos(1)
         if (.not.pos(1)) call fail(name//' '//trim(BASES(ib))//' '//trim(SNAMES(is))//': lost positivity')
         if (pos(1) .and. l1(1) > 3._R_P*l1p) call fail(name//' '//trim(BASES(ib))//' '//trim(SNAMES(is))// &
                                                        ': L1 above 3 times the pure scheme')
      enddo
   enddo
   hybrid = HYB_NONE ; hyb_sensor = SENS_LCR
   endsubroutine study

   subroutine sensor(is, v, s)
   !< Sensor is on the stencil v, normalised as the harness's hybrid does, at its published threshold.
   integer(I4P), intent(in)  :: is
   real(R_P),    intent(in)  :: v(-2:3)
   real(R_P),    intent(out) :: s
   real(R_P)                 :: rmu
   rmu = 1._R_P/max(maxval(abs(v)), 1._R_P)
   select case(is)
   case(1) ; call seth_sensor_teno(v, rmu, SETH_SENSOR_TENO_CT, s)
   case(2) ; call seth_sensor_ren(v, rmu, SETH_SENSOR_REN_PSIC, s)
   case(3) ; call seth_sensor_lcr(v, rmu, SETH_SENSOR_LCR_PSIC, s)
   case(4) ; call seth_sensor_tau(v, rmu, SETH_SENSOR_TAU_ETA, s)
   endselect
   endsubroutine sensor

   pure function scheme_index(name) result(id)
   character(*), intent(in) :: name
   integer(I4P)             :: id
   do id=1, NSCHEMES
      if (trim(SCHEMES(id)%name) == trim(name)) return
   enddo
   id = -1
   endfunction scheme_index

   subroutine fail(what)
   character(*), intent(in) :: what
   print '(A)', '   FAILED: '//what
   passed = .false.
   endsubroutine fail
endprogram seth_vv_sensors
