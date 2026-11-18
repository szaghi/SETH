!< SETH test: every kernel executed on the device (OpenACC) agrees with the same kernel on the host.
program seth_test_device
!< SETH test: every kernel executed on the device (OpenACC) agrees with the same kernel on the host.
!<
!< Each of the 63 cases runs one kernel inside `!$acc parallel loop` over N deterministic stencils (smooth cubics,
!< jumps, an integer LCG), then the same loop on the host, and compares: states within 1e-13 of the stencil's
!< magnitude, weights ratios within 1e-9 (the device and the host compile the kernel separately and may contract
!< multiply-adds differently; the tolerances are those of seth_test_golden, for the same reasons). One case calls
!< seth_weno_descaling on the device. Built without OpenACC the directives are ignored and both loops run on the host
!< (the comparison is then trivial); built with tests-nvf-oac it is the device test. Run it on an idle GPU
!< (CUDA_VISIBLE_DEVICES). Generated case by case: one block per kernel, the same body. The positivity limiter
!< (seth_positivity_euler, nv = 3 and 5) is checked after them on random cell pairs (limiter_case).
use seth
use seth_test_parameters, only : EPS_ZT, POS_RANGE_EXP, POSITIVITY_FLOOR, TOL_DEVICE_FLUX, TOL_DEVICE_STATE, &
                                 TOL_DEVICE_WR
implicit none
integer(I4P), parameter :: N = 4096_I4P
integer(I4P), parameter :: NK = 63_I4P
character(26), parameter :: NAMES(NK) = [ &
   'seth_weno1_js             ', &
   'seth_weno1_js_wr          ', &
   'seth_weno1_js_interp      ', &
   'seth_weno1_js_interp_wr   ', &
   'seth_weno3_js             ', &
   'seth_weno3_js_wr          ', &
   'seth_weno3_js_interp      ', &
   'seth_weno3_js_interp_wr   ', &
   'seth_weno5_js             ', &
   'seth_weno5_js_wr          ', &
   'seth_weno5_js_interp      ', &
   'seth_weno5_js_interp_wr   ', &
   'seth_weno7_js             ', &
   'seth_weno7_js_wr          ', &
   'seth_weno7_js_interp      ', &
   'seth_weno7_js_interp_wr   ', &
   'seth_weno9_js             ', &
   'seth_weno9_js_wr          ', &
   'seth_weno9_js_interp      ', &
   'seth_weno9_js_interp_wr   ', &
   'seth_weno5_js+descaling   ', &
   'seth_wenoz3               ', &
   'seth_wenoz5               ', &
   'seth_wenoz7               ', &
   'seth_wenoz9               ', &
   'seth_teno5                ', &
   'seth_teno6                ', &
   'seth_teno6_opt            ', &
   'seth_teno6_neno6          ', &
   'seth_wenoz3_interp        ', &
   'seth_wenoz5_interp        ', &
   'seth_wenoz7_interp        ', &
   'seth_wenoz9_interp        ', &
   'seth_teno5_interp         ', &
   'seth_teno6_interp         ', &
   'seth_teno5_opt            ', &
   'seth_teno5_az             ', &
   'seth_teno8                ', &
   'seth_teno5_a              ', &
   'seth_teno8_a              ', &
   'seth_weno5_js+descaling_ai', &
   'seth_weno9_js+mp_bound    ', &
   'seth_linear5              ', &
   'seth_linear8              ', &
   'seth_sensor_teno          ', &
   'seth_sensor_lcr           ', &
   'seth_weno4_c              ', &
   'seth_weno6_c              ', &
   'seth_weno8_c              ', &
   'seth_sensor_teno_eps      ', &
   'seth_teno5_mask           ', &
   'seth_teno5_dr             ', &
   'seth_teno5_thinc          ', &
   'seth_teno5_st             ', &
   'seth_teno5_rbf            ', &
   'seth_teno5_rbf_dr         ', &
   'seth_dr/thinc_detect+zeta ', &
   'seth_dr_override          ', &
   'seth_thinc_override       ', &
   'seth_teno5_face_masks     ', &
   'seth_teno5_dr_stored      ', &
   'seth_teno5_thinc_stored   ', &
   'seth_teno5_rbf_dr_stored  ']
real(R_P)               :: v(-5:6,N), w(-5:6,N), rp(N), rm(N), wrp(N), wrm(N), hp(N), hm(N), hwp(N), hwm(N)
real(R_P)               :: mag, ds, dw, worst_s, worst_w
integer(I4P)            :: k, i, nbit, nfail
logical                 :: passed

call stencils(v, w)
passed = .true.
#ifdef _OPENACC
print '(A)', 'OpenACC build: kernels on the device'
#else
print '(A)', 'no OpenACC: both loops on the host'
#endif
do k=1, NK
   select case(k)
   case(1) ! seth_weno1_js
      block
      real(R_P) :: a(0:1), b(0:1), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(0:1,i) ; b = w(0:1,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno1_js(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(0:1,i) ; b = w(0:1,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno1_js(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      endblock
   case(2) ! seth_weno1_js_wr
      block
      real(R_P) :: a(0:1), b(0:1), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(0:1,i) ; b = w(0:1,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno1_js_wr(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(0:1,i) ; b = w(0:1,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno1_js_wr(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      endblock
   case(3) ! seth_weno1_js_interp
      block
      real(R_P) :: a(0:1), b(0:1), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(0:1,i) ; b = w(0:1,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno1_js_interp(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(0:1,i) ; b = w(0:1,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno1_js_interp(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      endblock
   case(4) ! seth_weno1_js_interp_wr
      block
      real(R_P) :: a(0:1), b(0:1), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(0:1,i) ; b = w(0:1,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno1_js_interp_wr(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(0:1,i) ; b = w(0:1,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno1_js_interp_wr(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      endblock
   case(5) ! seth_weno3_js
      block
      real(R_P) :: a(-1:2), b(-1:2), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno3_js(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno3_js(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      endblock
   case(6) ! seth_weno3_js_wr
      block
      real(R_P) :: a(-1:2), b(-1:2), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno3_js_wr(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno3_js_wr(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      endblock
   case(7) ! seth_weno3_js_interp
      block
      real(R_P) :: a(-1:2), b(-1:2), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno3_js_interp(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno3_js_interp(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      endblock
   case(8) ! seth_weno3_js_interp_wr
      block
      real(R_P) :: a(-1:2), b(-1:2), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno3_js_interp_wr(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno3_js_interp_wr(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      endblock
   case(9) ! seth_weno5_js
      block
      real(R_P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno5_js(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno5_js(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      endblock
   case(10) ! seth_weno5_js_wr
      block
      real(R_P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno5_js_wr(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno5_js_wr(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      endblock
   case(11) ! seth_weno5_js_interp
      block
      real(R_P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno5_js_interp(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno5_js_interp(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      endblock
   case(12) ! seth_weno5_js_interp_wr
      block
      real(R_P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno5_js_interp_wr(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno5_js_interp_wr(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      endblock
   case(13) ! seth_weno7_js
      block
      real(R_P) :: a(-3:4), b(-3:4), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno7_js(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno7_js(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      endblock
   case(14) ! seth_weno7_js_wr
      block
      real(R_P) :: a(-3:4), b(-3:4), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno7_js_wr(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno7_js_wr(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      endblock
   case(15) ! seth_weno7_js_interp
      block
      real(R_P) :: a(-3:4), b(-3:4), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno7_js_interp(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno7_js_interp(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      endblock
   case(16) ! seth_weno7_js_interp_wr
      block
      real(R_P) :: a(-3:4), b(-3:4), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno7_js_interp_wr(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno7_js_interp_wr(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      endblock
   case(17) ! seth_weno9_js
      block
      real(R_P) :: a(-4:5), b(-4:5), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno9_js(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno9_js(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      endblock
   case(18) ! seth_weno9_js_wr
      block
      real(R_P) :: a(-4:5), b(-4:5), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno9_js_wr(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno9_js_wr(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      endblock
   case(19) ! seth_weno9_js_interp
      block
      real(R_P) :: a(-4:5), b(-4:5), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno9_js_interp(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno9_js_interp(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      endblock
   case(20) ! seth_weno9_js_interp_wr
      block
      real(R_P) :: a(-4:5), b(-4:5), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno9_js_interp_wr(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno9_js_interp_wr(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      endblock
   case(21) ! seth_weno5_js+descaling
      block
      real(R_P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         rmu = seth_weno_descaling(3_I4P, 1._R_P, a, b)
         call seth_weno5_js(a, b, SETH_WENO_JS_EPS, RMU, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         rmu = seth_weno_descaling(3_I4P, 1._R_P, a, b)
         call seth_weno5_js(a, b, SETH_WENO_JS_EPS, RMU, rp(i), rm(i))
      enddo
      endblock
   case(22) ! seth_wenoz3
      block
      real(R_P) :: a(-1:2), b(-1:2), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_wenoz3(a, b, EPS_ZT, 2_I4P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_wenoz3(a, b, EPS_ZT, 2_I4P, rp(i), rm(i))
      enddo
      endblock
   case(23) ! seth_wenoz5
      block
      real(R_P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_wenoz5(a, b, EPS_ZT, 2_I4P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_wenoz5(a, b, EPS_ZT, 2_I4P, rp(i), rm(i))
      enddo
      endblock
   case(24) ! seth_wenoz7
      block
      real(R_P) :: a(-3:4), b(-3:4), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_wenoz7(a, b, EPS_ZT, 2_I4P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_wenoz7(a, b, EPS_ZT, 2_I4P, rp(i), rm(i))
      enddo
      endblock
   case(25) ! seth_wenoz9
      block
      real(R_P) :: a(-4:5), b(-4:5), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_wenoz9(a, b, EPS_ZT, 2_I4P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_wenoz9(a, b, EPS_ZT, 2_I4P, rp(i), rm(i))
      enddo
      endblock
   case(26) ! seth_teno5
      block
      real(R_P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno5(a, b, SETH_TENO_EPS, SETH_TENO5_CT, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno5(a, b, SETH_TENO_EPS, SETH_TENO5_CT, rp(i), rm(i))
      enddo
      endblock
   case(27) ! seth_teno6
      block
      real(R_P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno6(a, b, SETH_TENO_EPS, SETH_TENO6_CT, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno6(a, b, SETH_TENO_EPS, SETH_TENO6_CT, rp(i), rm(i))
      enddo
      endblock
   case(28) ! seth_teno6_opt
      block
      real(R_P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno6_opt(a, b, SETH_TENO_EPS, SETH_TENO6_OPT_CT, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno6_opt(a, b, SETH_TENO_EPS, SETH_TENO6_OPT_CT, rp(i), rm(i))
      enddo
      endblock
   case(29) ! seth_teno6_neno6
      block
      real(R_P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno6_neno6(a, b, SETH_TENO_EPS, SETH_TENO6_NENO6_CT, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno6_neno6(a, b, SETH_TENO_EPS, SETH_TENO6_NENO6_CT, rp(i), rm(i))
      enddo
      endblock
   case(30) ! seth_wenoz3_interp
      block
      real(R_P) :: a(-1:2), b(-1:2), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_wenoz3_interp(a, b, EPS_ZT, 2_I4P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_wenoz3_interp(a, b, EPS_ZT, 2_I4P, rp(i), rm(i))
      enddo
      endblock
   case(31) ! seth_wenoz5_interp
      block
      real(R_P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_wenoz5_interp(a, b, EPS_ZT, 2_I4P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_wenoz5_interp(a, b, EPS_ZT, 2_I4P, rp(i), rm(i))
      enddo
      endblock
   case(32) ! seth_wenoz7_interp
      block
      real(R_P) :: a(-3:4), b(-3:4), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_wenoz7_interp(a, b, EPS_ZT, 2_I4P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_wenoz7_interp(a, b, EPS_ZT, 2_I4P, rp(i), rm(i))
      enddo
      endblock
   case(33) ! seth_wenoz9_interp
      block
      real(R_P) :: a(-4:5), b(-4:5), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_wenoz9_interp(a, b, EPS_ZT, 2_I4P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_wenoz9_interp(a, b, EPS_ZT, 2_I4P, rp(i), rm(i))
      enddo
      endblock
   case(34) ! seth_teno5_interp
      block
      real(R_P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno5_interp(a, b, SETH_TENO_EPS, SETH_TENO5_CT, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno5_interp(a, b, SETH_TENO_EPS, SETH_TENO5_CT, rp(i), rm(i))
      enddo
      endblock
   case(35) ! seth_teno6_interp
      block
      real(R_P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno6_interp(a, b, SETH_TENO_EPS, SETH_TENO6_CT, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno6_interp(a, b, SETH_TENO_EPS, SETH_TENO6_CT, rp(i), rm(i))
      enddo
      endblock
   case(36) ! seth_teno5_opt
      block
      real(R_P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno5_opt(a, b, SETH_TENO_EPS, SETH_TENO5_OPT_CT, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno5_opt(a, b, SETH_TENO_EPS, SETH_TENO5_OPT_CT, rp(i), rm(i))
      enddo
      endblock
   case(37) ! seth_teno5_az
      block
      real(R_P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno5_az(a, b, SETH_TENO_EPS, SETH_TENO5_AZ_CT, seth_weno_descaling(3_I4P, 1._R_P, a, b), rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno5_az(a, b, SETH_TENO_EPS, SETH_TENO5_AZ_CT, seth_weno_descaling(3_I4P, 1._R_P, a, b), rp(i), rm(i))
      enddo
      endblock
   case(38) ! seth_teno8
      block
      real(R_P) :: a(-3:4), b(-3:4), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno8(a, b, SETH_TENO_EPS, SETH_TENO8_CT, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno8(a, b, SETH_TENO_EPS, SETH_TENO8_CT, rp(i), rm(i))
      enddo
      endblock
   case(39) ! seth_teno5_a
      block
      real(R_P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno5_a(a, b, SETH_TENO_EPS, SETH_TENO5_A_CR, SETH_TENO5_A_XI, SETH_TENO5_A_A1, SETH_TENO5_A_A2, &
                           1._R_P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno5_a(a, b, SETH_TENO_EPS, SETH_TENO5_A_CR, SETH_TENO5_A_XI, SETH_TENO5_A_A1, SETH_TENO5_A_A2, &
                           1._R_P, rp(i), rm(i))
      enddo
      endblock
   case(40) ! seth_teno8_a
      block
      real(R_P) :: a(-3:4), b(-3:4), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno8_a(a, b, SETH_TENO_EPS, SETH_TENO8_A_CR, SETH_TENO8_A_XI, SETH_TENO8_A_A1, SETH_TENO8_A_A2, &
                           1._R_P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno8_a(a, b, SETH_TENO_EPS, SETH_TENO8_A_CR, SETH_TENO8_A_XI, SETH_TENO8_A_A1, SETH_TENO8_A_A2, &
                           1._R_P, rp(i), rm(i))
      enddo
      endblock
   case(41) ! seth_weno5_js+descaling_ai
      block
      real(R_P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         rmu = seth_weno_descaling_ai(3_I4P, a, b)
         call seth_weno5_js(a, b, SETH_WENO_JS_EPS, rmu, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         rmu = seth_weno_descaling_ai(3_I4P, a, b)
         call seth_weno5_js(a, b, SETH_WENO_JS_EPS, rmu, rp(i), rm(i))
      enddo
      endblock
   case(42) ! seth_weno9_js+mp_bound
      block
      real(R_P) :: a(-4:5), b(-4:5), c(-2:3), e(-2:3), rmu
      rmu = 1._R_P
      ! rmu is private and so undefined on the device: the kernels take the literal 1 (passing rmu put the WENO9-JS
      ! states up to 0.51 off the host on the GPU)
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, c, e, rmu)
      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno9_js(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
         c = a(-2:3) ; e = b(-2:3)
         call seth_mp_bound(c, e, SETH_MP_ALPHA, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno9_js(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
         c = a(-2:3) ; e = b(-2:3)
         call seth_mp_bound(c, e, SETH_MP_ALPHA, rp(i), rm(i))
      enddo
      endblock
   case(43) ! seth_linear5
      block
      real(R_P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_linear5(a, b, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_linear5(a, b, rp(i), rm(i))
      enddo
      endblock
   case(44) ! seth_linear8
      block
      real(R_P) :: a(-3:4), b(-3:4), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_linear8(a, b, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_linear8(a, b, rp(i), rm(i))
      enddo
      endblock
   case(45) ! seth_sensor_teno: the sensor value in rp, rm = rp
      block
      real(R_P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_sensor_teno(a, 1._R_P/max(maxval(abs(a)), 1._R_P), SETH_SENSOR_TENO_CT, rp(i))
         rm(i) = rp(i)
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_sensor_teno(a, 1._R_P/max(maxval(abs(a)), 1._R_P), SETH_SENSOR_TENO_CT, rp(i))
         rm(i) = rp(i)
      enddo
      endblock
   case(46) ! seth_sensor_lcr: the sensor value in rp, rm = rp
      block
      real(R_P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_sensor_lcr(a, 1._R_P/max(maxval(abs(a)), 1._R_P), SETH_SENSOR_LCR_PSIC, rp(i))
         rm(i) = rp(i)
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_sensor_lcr(a, 1._R_P/max(maxval(abs(a)), 1._R_P), SETH_SENSOR_LCR_PSIC, rp(i))
         rm(i) = rp(i)
      enddo
      endblock
   case(47) ! seth_weno4_c
      block
      real(R_P) :: a(-1:2), b(-1:2), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno4_c(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno4_c(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      endblock
   case(48) ! seth_weno6_c
      block
      real(R_P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno6_c(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno6_c(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      endblock
   case(49) ! seth_weno8_c
      block
      real(R_P) :: a(-3:4), b(-3:4), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno8_c(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_weno8_c(a, b, SETH_WENO_JS_EPS, 1._R_P, rp(i), rm(i))
      enddo
      endblock
   case(50) ! seth_sensor_teno_eps with Fu's multi-dimensional eps, 1e-3: the sensor value in rp, rm = rp
      block
      real(R_P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R_P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_sensor_teno_eps(a, 1._R_P/max(maxval(abs(a)), 1._R_P), SETH_SENSOR_TENO_CT, SETH_SENSOR_TENO_EPS_MD, rp(i))
         rm(i) = rp(i)
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_sensor_teno_eps(a, 1._R_P/max(maxval(abs(a)), 1._R_P), SETH_SENSOR_TENO_CT, SETH_SENSOR_TENO_EPS_MD, rp(i))
         rm(i) = rp(i)
      enddo
      endblock
   case(51) ! seth_teno5_mask: the masks of the left-biased stencil and of the mirrored right-biased one in rp, rm
      block
      real(R_P)    :: a(-2:2), b(-2:2)
      integer(I4P) :: mk, mk2, m
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, mk, mk2)
      do i=1, N
         a = v(-2:2,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         do m=-2, 2
            b(m) = w(1-m,i)
         enddo
         call seth_teno5_mask(a, SETH_TENO_EPS, SETH_TENO5_CT, mk)
         call seth_teno5_mask(b, SETH_TENO_EPS, SETH_TENO5_CT, mk2)
         rp(i) = real(mk, R_P) ; rm(i) = real(mk2, R_P)
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:2,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         do m=-2, 2
            b(m) = w(1-m,i)
         enddo
         call seth_teno5_mask(a, SETH_TENO_EPS, SETH_TENO5_CT, mk)
         call seth_teno5_mask(b, SETH_TENO_EPS, SETH_TENO5_CT, mk2)
         rp(i) = real(mk, R_P) ; rm(i) = real(mk2, R_P)
      enddo
      endblock
   case(52) ! seth_teno5_dr, the papers' constants
      block
      real(R_P) :: a(-5:6), b(-5:6)
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b)
      do i=1, N
         a = v(:,i) ; b = w(:,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno5_dr(a, b, SETH_TENO_EPS, SETH_TENO5_CT, SETH_DR_T1, SETH_DR_T2, 1._R_P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(:,i) ; b = w(:,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno5_dr(a, b, SETH_TENO_EPS, SETH_TENO5_CT, SETH_DR_T1, SETH_DR_T2, 1._R_P, rp(i), rm(i))
      enddo
      endblock
   case(53) ! seth_teno5_thinc, the paper's constants
      block
      real(R_P) :: a(-5:6), b(-5:6)
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b)
      do i=1, N
         a = v(:,i) ; b = w(:,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno5_thinc(a, b, SETH_TENO_EPS, SETH_TENO5_CT, SETH_THINC_BETA, 1._R_P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(:,i) ; b = w(:,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno5_thinc(a, b, SETH_TENO_EPS, SETH_TENO5_CT, SETH_THINC_BETA, 1._R_P, rp(i), rm(i))
      enddo
      endblock
   case(54, 55) ! seth_teno5_st; seth_teno5_rbf with s1 = s2 = 0.09 (the paper's alpha = 3 at dx = 0.03)
      block
      real(R_P) :: a(-2:3), b(-2:3), sv
      sv = merge(0.09_R_P, 0._R_P, k == 55)
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         if (k == 54) then
            call seth_teno5_st(a, b, SETH_TENO_EPS, SETH_TENO5_CT, rp(i), rm(i))
         else
            call seth_teno5_rbf(a, b, SETH_TENO_EPS, SETH_TENO5_CT, SETH_RBF_CR, sv, sv, rp(i), rm(i))
         endif
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         if (k == 54) then
            call seth_teno5_st(a, b, SETH_TENO_EPS, SETH_TENO5_CT, rp(i), rm(i))
         else
            call seth_teno5_rbf(a, b, SETH_TENO_EPS, SETH_TENO5_CT, SETH_RBF_CR, sv, sv, rp(i), rm(i))
         endif
      enddo
      endblock
   case(56) ! seth_teno5_rbf_dr, s1 = s2 = 0.09, the papers' DR constants
      block
      real(R_P) :: a(-5:6), b(-5:6)
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b)
      do i=1, N
         a = v(:,i) ; b = w(:,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno5_rbf_dr(a, b, SETH_TENO_EPS, SETH_TENO5_CT, SETH_RBF_CR, 0.09_R_P, 0.09_R_P, SETH_DR_T1, &
                                SETH_DR_T2, 1._R_P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(:,i) ; b = w(:,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno5_rbf_dr(a, b, SETH_TENO_EPS, SETH_TENO5_CT, SETH_RBF_CR, 0.09_R_P, 0.09_R_P, SETH_DR_T1, &
                                SETH_DR_T2, 1._R_P, rp(i), rm(i))
      enddo
      endblock
   case(57) ! seth_dr_detect, seth_thinc_detect, seth_teno5_zeta on the masks of seth_teno5_mask: the levels coded in
      block ! rp (DR: lp + 3 lm), rm (THINC: lp + 2 lm, + 4 (zeta + 1) of the cell 0); a detector must fire somewhere
      real(R_P)    :: a(-2:2), b(-2:2)
      integer(I4P) :: mp(-3:3), mm(-3:3), lp, lm, tp, tm, c, m
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, mp, mm, lp, lm, tp, tm, c, m)
      do i=1, N
         wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         do c=-3, 3
            do m=-2, 2
               a(m) = v(c+m,i) ; b(m) = w(1-c-m,i)
            enddo
            call seth_teno5_mask(a, SETH_TENO_EPS, SETH_TENO5_CT, mp(c))
            call seth_teno5_mask(b, SETH_TENO_EPS, SETH_TENO5_CT, mm(c))
         enddo
         call seth_dr_detect(mp, mm, lp, lm)
         call seth_thinc_detect(mp, mm, tp, tm)
         rp(i) = real(lp + 3*lm, R_P) ; rm(i) = real(tp + 2*tm + 4*(seth_teno5_zeta(mp(0)) + 1), R_P)
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         do c=-3, 3
            do m=-2, 2
               a(m) = v(c+m,i) ; b(m) = w(1-c-m,i)
            enddo
            call seth_teno5_mask(a, SETH_TENO_EPS, SETH_TENO5_CT, mp(c))
            call seth_teno5_mask(b, SETH_TENO_EPS, SETH_TENO5_CT, mm(c))
         enddo
         call seth_dr_detect(mp, mm, lp, lm)
         call seth_thinc_detect(mp, mm, tp, tm)
         rp(i) = real(lp + 3*lm, R_P) ; rm(i) = real(tp + 2*tm + 4*(seth_teno5_zeta(mp(0)) + 1), R_P)
      enddo
      print '(A,I5,A,I5,A)', '   detectors on the device: DR fires on ', count(mod(nint(hp), 3) > 0), &
            ' stencils, THINC on ', count(mod(nint(hm), 2) > 0), ' (left states)'
      if (count(mod(nint(hp), 3) > 0) == 0 .or. count(mod(nint(hm), 2) > 0) == 0) then
         print '(A)', '   FAILED: a detector never fires on the test stencils' ; passed = .false.
      endif
      endblock
   case(58, 59) ! seth_dr_override (levels 0, 1, 2 by stencil), seth_thinc_override (0, 1), states preset to the mean
      block      ! of f_i and f_{i+1}; at least one state changed by each override
      real(R_P)    :: a(-1:2), b(-1:2)
      integer(I4P) :: lv
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, lv)
      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         rp(i) = 0.5_R_P*(a(0) + a(1)) ; rm(i) = 0.5_R_P*(b(0) + b(1))
         if (k == 58) then
            lv = mod(i, 3_I4P)
            call seth_dr_override(a, b, lv, lv, SETH_DR_T1, SETH_DR_T2, 1._R_P, rp(i), rm(i))
         else
            lv = mod(i, 2_I4P)
            call seth_thinc_override(a, b, lv, lv, SETH_THINC_BETA, 1._R_P, rp(i), rm(i))
         endif
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         rp(i) = 0.5_R_P*(a(0) + a(1)) ; rm(i) = 0.5_R_P*(b(0) + b(1))
         if (k == 58) then
            lv = mod(i, 3_I4P)
            call seth_dr_override(a, b, lv, lv, SETH_DR_T1, SETH_DR_T2, 1._R_P, rp(i), rm(i))
         else
            lv = mod(i, 2_I4P)
            call seth_thinc_override(a, b, lv, lv, SETH_THINC_BETA, 1._R_P, rp(i), rm(i))
         endif
      enddo
      if (count(hp /= 0.5_R_P*(v(0,:) + v(1,:))) == 0) then
         print '(A)', '   FAILED: the override never acted on the device' ; passed = .false.
      endif
      endblock
   case(60) ! seth_teno5_face_masks: the two masks of the face coded in rp (left, cell 0) and rm (right, cell 1)
      block
      real(R_P)    :: a(-2:3), b(-2:3)
      integer(I4P) :: ml, mr
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, ml, mr)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno5_face_masks(a, b, SETH_TENO_EPS, SETH_TENO5_CT, ml, mr)
         rp(i) = real(ml, R_P) ; rm(i) = real(mr, R_P)
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         call seth_teno5_face_masks(a, b, SETH_TENO_EPS, SETH_TENO5_CT, ml, mr)
         rp(i) = real(ml, R_P) ; rm(i) = real(mr, R_P)
      enddo
      endblock
   case(61:63) ! the stored-mask pair: the masks of the faces -3..3 by seth_teno5_face_masks, then seth_teno5_dr_stored,
      block    ! seth_teno5_thinc_stored, seth_teno5_rbf_dr_stored (s1 = s2 = 0.09)
      real(R_P)    :: a(-2:3), b(-2:3)
      integer(I4P) :: msl(-3:3), msr(-2:4), mp(-3:3), mm(-3:3), c, m
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, msl, msr, mp, mm, c, m)
      do i=1, N
         wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         do c=-3, 3 ! pass 1 at the faces c+1/2 of the stencil
            do m=-2, 3
               a(m) = v(c+m,i) ; b(m) = w(c+m,i)
            enddo
            call seth_teno5_face_masks(a, b, SETH_TENO_EPS, SETH_TENO5_CT, msl(c), msr(c+1))
         enddo
         do c=-3, 3
            mp(c) = msl(c) ; mm(c) = msr(1-c)
         enddo
         a = v(-2:3,i) ; b = w(-2:3,i)
         if (k == 61) then
            call seth_teno5_dr_stored(a, b, mp, mm, SETH_TENO_EPS, SETH_TENO5_CT, SETH_DR_T1, SETH_DR_T2, 1._R_P, &
                                      rp(i), rm(i))
         elseif (k == 62) then
            call seth_teno5_thinc_stored(a, b, mp, mm, SETH_TENO_EPS, SETH_TENO5_CT, SETH_THINC_BETA, 1._R_P, &
                                         rp(i), rm(i))
         else
            call seth_teno5_rbf_dr_stored(a, b, mp, mm, SETH_TENO_EPS, SETH_TENO5_CT, SETH_RBF_CR, 0.09_R_P, 0.09_R_P, &
                                          SETH_DR_T1, SETH_DR_T2, 1._R_P, rp(i), rm(i))
         endif
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         wrp(i) = -1._R_P ; wrm(i) = -1._R_P
         do c=-3, 3 ! pass 1 at the faces c+1/2 of the stencil
            do m=-2, 3
               a(m) = v(c+m,i) ; b(m) = w(c+m,i)
            enddo
            call seth_teno5_face_masks(a, b, SETH_TENO_EPS, SETH_TENO5_CT, msl(c), msr(c+1))
         enddo
         do c=-3, 3
            mp(c) = msl(c) ; mm(c) = msr(1-c)
         enddo
         a = v(-2:3,i) ; b = w(-2:3,i)
         if (k == 61) then
            call seth_teno5_dr_stored(a, b, mp, mm, SETH_TENO_EPS, SETH_TENO5_CT, SETH_DR_T1, SETH_DR_T2, 1._R_P, &
                                      rp(i), rm(i))
         elseif (k == 62) then
            call seth_teno5_thinc_stored(a, b, mp, mm, SETH_TENO_EPS, SETH_TENO5_CT, SETH_THINC_BETA, 1._R_P, &
                                         rp(i), rm(i))
         else
            call seth_teno5_rbf_dr_stored(a, b, mp, mm, SETH_TENO_EPS, SETH_TENO5_CT, SETH_RBF_CR, 0.09_R_P, 0.09_R_P, &
                                          SETH_DR_T1, SETH_DR_T2, 1._R_P, rp(i), rm(i))
         endif
      enddo
      endblock
   endselect
   nbit = 0 ; nfail = 0 ; worst_s = 0._R_P ; worst_w = 0._R_P
   do i=1, N
      mag = max(maxval(abs(v(-4:5,i))), maxval(abs(w(-4:5,i))), tiny(1._R_P)) ! -4:5: as before the widening
      ds = max(abs(rp(i) - hp(i)), abs(rm(i) - hm(i)))/mag
      dw = max(abs(wrp(i) - hwp(i)), abs(wrm(i) - hwm(i)))
      if (rp(i) == hp(i) .and. rm(i) == hm(i) .and. wrp(i) == hwp(i) .and. wrm(i) == hwm(i)) nbit = nbit + 1
      worst_s = max(worst_s, ds) ; worst_w = max(worst_w, dw)
      if (.not.(ds <= TOL_DEVICE_STATE .and. dw <= TOL_DEVICE_WR)) nfail = nfail + 1
   enddo
   print '(A26,A,I5,A,I5,A,ES9.2,A,ES9.2)', NAMES(k), ': bitwise ', nbit, '/', N, ', worst state ', worst_s, &
         ', worst wr ', worst_w
   if (nfail > 0) passed = .false.
enddo
call limiter_case(3_I4P)
call limiter_case(5_I4P)
if (.not.passed) error stop 'device results differ from the host beyond round-off'
print '(A)', 'seth_test_device: ok'

contains
   subroutine limiter_case(nv)
   !< seth_positivity_euler on the device against the host: NL deterministic cell pairs (an integer LCG: rho and p
   !< log-uniform in [1e-8, 1], velocities in [-10, 10]), Lax-Friedrichs low-order flux at CFL 1/2, a perturbed
   !< high-order flux. theta within 1e-10 and the flux within 1e-10 of its magnitude (the closed form's square root and
   !< the multiply-adds may round differently), and every device flux admissible (checked on the host).
   integer(I4P), intent(in) :: nv
   integer(I4P), parameter  :: NL = 4096_I4P
   real(R_P),    parameter  :: GM = 1.4_R_P, EP = SETH_POSITIVITY_EPS
   real(R_P)                :: ql(nv,NL), qr(nv,NL), fh(nv,NL), fl(nv,NL), lam(NL), fd(nv,NL), td(NL)
   real(R_P)                :: fc(nv,NL), tc(NL)
   real(R_P)                :: a, b, worst, mag
   integer(I8P)             :: st
   integer(I4P)             :: i, k, nbit, nlim, nfail
   st = 20261115_I8P
   do i=1, NL
      call lim_state(nv, st, ql(:,i)) ; call lim_state(nv, st, qr(:,i))
      a = max(lim_speed(nv, ql(:,i)), lim_speed(nv, qr(:,i))) ; lam(i) = 1._R_P/a
      fl(:,i) = 0.5_R_P*(lim_flux(nv, ql(:,i)) + lim_flux(nv, qr(:,i))) - 0.5_R_P*a*(qr(:,i) - ql(:,i))
      do k=1, nv
         b = lim_rnd(st)
         fh(k,i) = fl(k,i) + (2._R_P*b - 1._R_P)*(abs(fl(k,i)) + abs(ql(k,i)) + abs(qr(k,i))) &
                   *10._R_P**(3._R_P*lim_rnd(st) - 2._R_P)
      enddo
   enddo
   !$acc parallel loop gang vector copyin(ql, qr, fh, fl, lam) copyout(fd, td)
   do i=1, NL
      call seth_positivity_euler(nv, GM, lam(i), EP, EP, ql(:,i), qr(:,i), fh(:,i), fl(:,i), fd(:,i), td(i))
   enddo
   do i=1, NL
      call seth_positivity_euler(nv, GM, lam(i), EP, EP, ql(:,i), qr(:,i), fh(:,i), fl(:,i), fc(:,i), tc(i))
   enddo
   nbit = 0 ; nlim = 0 ; nfail = 0 ; worst = 0._R_P
   do i=1, NL
      if (tc(i) < 1._R_P) nlim = nlim + 1
      if (td(i) == tc(i) .and. all(fd(:,i) == fc(:,i))) nbit = nbit + 1
      mag = max(maxval(abs(fh(:,i))), maxval(abs(fl(:,i))), tiny(1._R_P))
      worst = max(worst, abs(td(i) - tc(i)), maxval(abs(fd(:,i) - fc(:,i)))/mag)
      if (.not.(abs(td(i) - tc(i)) <= TOL_DEVICE_FLUX .and. maxval(abs(fd(:,i) - fc(:,i))) <= TOL_DEVICE_FLUX*mag)) &
         nfail = nfail + 1
      ! to the round-off of the states as evaluated here (seth_test_positivity, adm_ro)
      if (.not.(lim_adm_ro(nv, ql(:,i) - lam(i)*fd(:,i), abs(ql(:,i)) + abs(lam(i)*fd(:,i))) .and. &
                lim_adm_ro(nv, qr(:,i) + lam(i)*fd(:,i), abs(qr(:,i)) + abs(lam(i)*fd(:,i))))) then
         if (seth_euler_admissible(nv, GM, ql(:,i) - lam(i)*fl(:,i), EP, EP) .and. &
             seth_euler_admissible(nv, GM, qr(:,i) + lam(i)*fl(:,i), EP, EP)) nfail = nfail + 1
      endif
   enddo
   print '(A,I1,A,I5,A,I5,A,I5,A,ES9.2)', 'seth_positivity_euler nv=', nv, ': limited ', nlim, ', bitwise ', nbit, &
         '/', &
         NL, ', worst ', worst
   if (nfail > 0) passed = .false.
   endsubroutine limiter_case

   function lim_rnd(st) result(r)
   !< The next number of the integer LCG (Park-Miller), in (0, 1).
   integer(I8P), intent(inout) :: st
   real(R_P)                   :: r
   st = mod(16807_I8P*st, 2147483647_I8P) ; r = real(st, R_P)/2147483647._R_P
   endfunction lim_rnd

   subroutine lim_state(nv, st, q)
   !< rho and p log-uniform in [1e-8, 1], velocities in [-10, 10] (gamma = 1.4).
   integer(I4P), intent(in)    :: nv
   integer(I8P), intent(inout) :: st
   real(R_P),    intent(out)   :: q(nv)
   real(R_P)                   :: rho, p
   integer(I4P)                :: kk
   rho = 10._R_P**(-POS_RANGE_EXP*lim_rnd(st)) ; p = 10._R_P**(-POS_RANGE_EXP*lim_rnd(st)) ! R8P: from 1e-8
   q(1) = rho
   do kk=2, nv-1
      q(kk) = rho*(20._R_P*lim_rnd(st) - 10._R_P)
   enddo
   q(nv) = p/0.4_R_P + 0.5_R_P*dot_product(q(2:nv-1), q(2:nv-1))/rho
   endsubroutine lim_state

   pure function lim_flux(nv, q) result(fx)
   !< The Euler flux normal to direction 1 (gamma = 1.4).
   integer(I4P), intent(in) :: nv
   real(R_P),    intent(in) :: q(nv)
   real(R_P)                :: fx(nv), u, p
   u = q(2)/q(1) ; p = 0.4_R_P*(q(nv) - 0.5_R_P*dot_product(q(2:nv-1), q(2:nv-1))/q(1))
   fx = u*q ; fx(2) = fx(2) + p ; fx(nv) = fx(nv) + p*u
   endfunction lim_flux

   pure logical function lim_adm_ro(nv, q, mag) result(ok)
   !< rho and p above the floor 1e-13 to the round-off of the one-sided state q (mag = |q_cell| + |lam f|): the
   !< allowance of seth_test_positivity's adm_ro, gamma = 1.4.
   integer(I4P), intent(in) :: nv
   real(R_P),    intent(in) :: q(nv), mag(nv)
   real(R_P)                :: ke, tol, d(nv)
   integer(I4P)             :: k
   ok = .false.
   d = 4._R_P*epsilon(1._R_P)*mag
   if (.not.(q(1) >= POSITIVITY_FLOOR - d(1))) return
   ke = 0._R_P ; tol = d(nv)
   do k=2, nv-1
      ke = ke + q(k)*q(k)
      tol = tol + abs(q(k)/q(1))*d(k)
   enddo
   ke = 0.5_R_P*ke/q(1)
   tol = 0.4_R_P*(tol + ke/q(1)*d(1) + 4._R_P*epsilon(1._R_P)*(abs(q(nv)) + ke))
   ok = 0.4_R_P*(q(nv) - ke) >= POSITIVITY_FLOOR - tol
   endfunction lim_adm_ro

   pure function lim_speed(nv, q) result(c)
   !< |u_1| + c (gamma = 1.4).
   integer(I4P), intent(in) :: nv
   real(R_P),    intent(in) :: q(nv)
   real(R_P)                :: c, p
   p = 0.4_R_P*(q(nv) - 0.5_R_P*dot_product(q(2:nv-1), q(2:nv-1))/q(1))
   c = abs(q(2)/q(1)) + sqrt(1.4_R_P*p/q(1))
   endfunction lim_speed

   subroutine stencils(v, w)
   !< N stencil pairs, exact in binary: smooth cubics, jumps, an integer LCG, cycling.
   real(R_P), intent(out) :: v(-5:6,N), w(-5:6,N)
   integer(I8P)           :: state
   integer(I4P)           :: i, m, j
   state = 20261115_I8P
   do i=1, N
      j = mod(i, 16)
      do m=-4, 5
         select case(mod(i, 3))
         case(0)
            v(m,i) = 1._R_P + (j/8._R_P)*m - (j/64._R_P)*m*m + (1._R_P/512._R_P)*m**3
            w(m,i) = -0.5_R_P + (j/16._R_P)*m + (j/128._R_P)*m*m
         case(1)
            v(m,i) = merge(1._R_P, 0._R_P, m >= mod(i, 10) - 4)
            w(m,i) = merge(-2._R_P, 0.5_R_P, m >= mod(i, 9) - 4)
         case default
            state = mod(16807_I8P*state, 2147483647_I8P) ; v(m,i) = real(state, R_P)/2147483647._R_P - 0.5_R_P
            state = mod(16807_I8P*state, 2147483647_I8P) ; w(m,i) = real(state, R_P)/2147483647._R_P - 0.3_R_P
         endselect
      enddo
      ! the cells -5 and 6 of the composite kernels: the same formula, or (the LCG) constant continuation, so that the
      ! cells -4..5 and the LCG's stream are those written before the widening
      do m=-5, 6, 11
         select case(mod(i, 3))
         case(0)
            v(m,i) = 1._R_P + (j/8._R_P)*m - (j/64._R_P)*m*m + (1._R_P/512._R_P)*m**3
            w(m,i) = -0.5_R_P + (j/16._R_P)*m + (j/128._R_P)*m*m
         case(1)
            v(m,i) = merge(1._R_P, 0._R_P, m >= mod(i, 10) - 4)
            w(m,i) = merge(-2._R_P, 0.5_R_P, m >= mod(i, 9) - 4)
         case default
            v(m,i) = v(merge(-4, 5, m < 0),i) ; w(m,i) = w(merge(-4, 5, m < 0),i)
         endselect
      enddo
   enddo
   endsubroutine stencils
endprogram seth_test_device
