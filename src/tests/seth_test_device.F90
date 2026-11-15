!< SETH test: every kernel executed on the device (OpenACC) agrees with the same kernel on the host.
program seth_test_device
!< SETH test: every kernel executed on the device (OpenACC) agrees with the same kernel on the host.
!<
!< Each of the 35 cases runs one kernel inside `!$acc parallel loop` over N deterministic stencils (smooth cubics,
!< jumps, an integer LCG), then the same loop on the host, and compares: states within 1e-13 of the stencil's
!< magnitude, weights ratios within 1e-9 (the device and the host compile the kernel separately and may contract
!< multiply-adds differently; the tolerances are those of seth_test_golden, for the same reasons). One case calls
!< seth_weno_descaling on the device. Built without OpenACC the directives are ignored and both loops run on the host
!< (the comparison is then trivial); built with tests-nvf-oac it is the device test. Run it on an idle GPU
!< (CUDA_VISIBLE_DEVICES). Generated case by case: one block per kernel, the same body. The positivity limiter
!< (seth_positivity_euler, nv = 3 and 5) is checked after them on random cell pairs (limiter_case).
use seth
implicit none
integer(I4P), parameter :: N = 4096_I4P
integer(I4P), parameter :: NK = 35_I4P
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
   'seth_teno6_interp         ']
real(R8P)               :: v(-4:5,N), w(-4:5,N), rp(N), rm(N), wrp(N), wrm(N), hp(N), hm(N), hwp(N), hwm(N)
real(R8P)               :: mag, ds, dw, worst_s, worst_w
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
      real(R8P) :: a(0:1), b(0:1), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(0:1,i) ; b = w(0:1,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno1_js(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(0:1,i) ; b = w(0:1,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno1_js(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
      endblock
   case(2) ! seth_weno1_js_wr
      block
      real(R8P) :: a(0:1), b(0:1), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(0:1,i) ; b = w(0:1,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno1_js_wr(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(0:1,i) ; b = w(0:1,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno1_js_wr(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      endblock
   case(3) ! seth_weno1_js_interp
      block
      real(R8P) :: a(0:1), b(0:1), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(0:1,i) ; b = w(0:1,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno1_js_interp(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(0:1,i) ; b = w(0:1,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno1_js_interp(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
      endblock
   case(4) ! seth_weno1_js_interp_wr
      block
      real(R8P) :: a(0:1), b(0:1), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(0:1,i) ; b = w(0:1,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno1_js_interp_wr(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(0:1,i) ; b = w(0:1,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno1_js_interp_wr(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      endblock
   case(5) ! seth_weno3_js
      block
      real(R8P) :: a(-1:2), b(-1:2), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno3_js(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno3_js(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
      endblock
   case(6) ! seth_weno3_js_wr
      block
      real(R8P) :: a(-1:2), b(-1:2), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno3_js_wr(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno3_js_wr(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      endblock
   case(7) ! seth_weno3_js_interp
      block
      real(R8P) :: a(-1:2), b(-1:2), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno3_js_interp(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno3_js_interp(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
      endblock
   case(8) ! seth_weno3_js_interp_wr
      block
      real(R8P) :: a(-1:2), b(-1:2), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno3_js_interp_wr(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno3_js_interp_wr(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      endblock
   case(9) ! seth_weno5_js
      block
      real(R8P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno5_js(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno5_js(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
      endblock
   case(10) ! seth_weno5_js_wr
      block
      real(R8P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno5_js_wr(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno5_js_wr(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      endblock
   case(11) ! seth_weno5_js_interp
      block
      real(R8P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno5_js_interp(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno5_js_interp(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
      endblock
   case(12) ! seth_weno5_js_interp_wr
      block
      real(R8P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno5_js_interp_wr(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno5_js_interp_wr(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      endblock
   case(13) ! seth_weno7_js
      block
      real(R8P) :: a(-3:4), b(-3:4), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno7_js(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno7_js(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
      endblock
   case(14) ! seth_weno7_js_wr
      block
      real(R8P) :: a(-3:4), b(-3:4), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno7_js_wr(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno7_js_wr(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      endblock
   case(15) ! seth_weno7_js_interp
      block
      real(R8P) :: a(-3:4), b(-3:4), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno7_js_interp(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno7_js_interp(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
      endblock
   case(16) ! seth_weno7_js_interp_wr
      block
      real(R8P) :: a(-3:4), b(-3:4), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno7_js_interp_wr(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno7_js_interp_wr(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      endblock
   case(17) ! seth_weno9_js
      block
      real(R8P) :: a(-4:5), b(-4:5), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno9_js(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno9_js(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
      endblock
   case(18) ! seth_weno9_js_wr
      block
      real(R8P) :: a(-4:5), b(-4:5), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno9_js_wr(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno9_js_wr(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      endblock
   case(19) ! seth_weno9_js_interp
      block
      real(R8P) :: a(-4:5), b(-4:5), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno9_js_interp(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno9_js_interp(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
      endblock
   case(20) ! seth_weno9_js_interp_wr
      block
      real(R8P) :: a(-4:5), b(-4:5), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno9_js_interp_wr(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_weno9_js_interp_wr(a, b, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i), wrp(i), wrm(i))
      enddo
      endblock
   case(21) ! seth_weno5_js+descaling
      block
      real(R8P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         rmu = seth_weno_descaling(3_I4P, 1._R8P, a, b)
         call seth_weno5_js(a, b, SETH_WENO_JS_EPS, RMU, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         rmu = seth_weno_descaling(3_I4P, 1._R8P, a, b)
         call seth_weno5_js(a, b, SETH_WENO_JS_EPS, RMU, rp(i), rm(i))
      enddo
      endblock
   case(22) ! seth_wenoz3
      block
      real(R8P) :: a(-1:2), b(-1:2), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_wenoz3(a, b, 1.e-40_R8P, 2_I4P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_wenoz3(a, b, 1.e-40_R8P, 2_I4P, rp(i), rm(i))
      enddo
      endblock
   case(23) ! seth_wenoz5
      block
      real(R8P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_wenoz5(a, b, 1.e-40_R8P, 2_I4P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_wenoz5(a, b, 1.e-40_R8P, 2_I4P, rp(i), rm(i))
      enddo
      endblock
   case(24) ! seth_wenoz7
      block
      real(R8P) :: a(-3:4), b(-3:4), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_wenoz7(a, b, 1.e-40_R8P, 2_I4P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_wenoz7(a, b, 1.e-40_R8P, 2_I4P, rp(i), rm(i))
      enddo
      endblock
   case(25) ! seth_wenoz9
      block
      real(R8P) :: a(-4:5), b(-4:5), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_wenoz9(a, b, 1.e-40_R8P, 2_I4P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_wenoz9(a, b, 1.e-40_R8P, 2_I4P, rp(i), rm(i))
      enddo
      endblock
   case(26) ! seth_teno5
      block
      real(R8P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_teno5(a, b, SETH_TENO_EPS, SETH_TENO5_CT, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_teno5(a, b, SETH_TENO_EPS, SETH_TENO5_CT, rp(i), rm(i))
      enddo
      endblock
   case(27) ! seth_teno6
      block
      real(R8P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_teno6(a, b, SETH_TENO_EPS, SETH_TENO6_CT, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_teno6(a, b, SETH_TENO_EPS, SETH_TENO6_CT, rp(i), rm(i))
      enddo
      endblock
   case(28) ! seth_teno6_opt
      block
      real(R8P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_teno6_opt(a, b, SETH_TENO_EPS, SETH_TENO6_OPT_CT, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_teno6_opt(a, b, SETH_TENO_EPS, SETH_TENO6_OPT_CT, rp(i), rm(i))
      enddo
      endblock
   case(29) ! seth_teno6_neno6
      block
      real(R8P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_teno6_neno6(a, b, SETH_TENO_EPS, SETH_TENO6_NENO6_CT, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_teno6_neno6(a, b, SETH_TENO_EPS, SETH_TENO6_NENO6_CT, rp(i), rm(i))
      enddo
      endblock
   case(30) ! seth_wenoz3_interp
      block
      real(R8P) :: a(-1:2), b(-1:2), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_wenoz3_interp(a, b, 1.e-40_R8P, 2_I4P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-1:2,i) ; b = w(-1:2,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_wenoz3_interp(a, b, 1.e-40_R8P, 2_I4P, rp(i), rm(i))
      enddo
      endblock
   case(31) ! seth_wenoz5_interp
      block
      real(R8P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_wenoz5_interp(a, b, 1.e-40_R8P, 2_I4P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_wenoz5_interp(a, b, 1.e-40_R8P, 2_I4P, rp(i), rm(i))
      enddo
      endblock
   case(32) ! seth_wenoz7_interp
      block
      real(R8P) :: a(-3:4), b(-3:4), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_wenoz7_interp(a, b, 1.e-40_R8P, 2_I4P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-3:4,i) ; b = w(-3:4,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_wenoz7_interp(a, b, 1.e-40_R8P, 2_I4P, rp(i), rm(i))
      enddo
      endblock
   case(33) ! seth_wenoz9_interp
      block
      real(R8P) :: a(-4:5), b(-4:5), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_wenoz9_interp(a, b, 1.e-40_R8P, 2_I4P, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-4:5,i) ; b = w(-4:5,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_wenoz9_interp(a, b, 1.e-40_R8P, 2_I4P, rp(i), rm(i))
      enddo
      endblock
   case(34) ! seth_teno5_interp
      block
      real(R8P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_teno5_interp(a, b, SETH_TENO_EPS, SETH_TENO5_CT, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_teno5_interp(a, b, SETH_TENO_EPS, SETH_TENO5_CT, rp(i), rm(i))
      enddo
      endblock
   case(35) ! seth_teno6_interp
      block
      real(R8P) :: a(-2:3), b(-2:3), rmu
      rmu = 1._R8P
      !$acc parallel loop gang vector copyin(v, w) copyout(rp, rm, wrp, wrm) private(a, b, rmu)
      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_teno6_interp(a, b, SETH_TENO_EPS, SETH_TENO6_CT, rp(i), rm(i))
      enddo
      hp = rp ; hm = rm ; hwp = wrp ; hwm = wrm

      do i=1, N
         a = v(-2:3,i) ; b = w(-2:3,i) ; wrp(i) = -1._R8P ; wrm(i) = -1._R8P
         call seth_teno6_interp(a, b, SETH_TENO_EPS, SETH_TENO6_CT, rp(i), rm(i))
      enddo
      endblock
   endselect
   nbit = 0 ; nfail = 0 ; worst_s = 0._R8P ; worst_w = 0._R8P
   do i=1, N
      mag = max(maxval(abs(v(:,i))), maxval(abs(w(:,i))), tiny(1._R8P))
      ds = max(abs(rp(i) - hp(i)), abs(rm(i) - hm(i)))/mag
      dw = max(abs(wrp(i) - hwp(i)), abs(wrm(i) - hwm(i)))
      if (rp(i) == hp(i) .and. rm(i) == hm(i) .and. wrp(i) == hwp(i) .and. wrm(i) == hwm(i)) nbit = nbit + 1
      worst_s = max(worst_s, ds) ; worst_w = max(worst_w, dw)
      if (.not.(ds <= 1.e-13_R8P .and. dw <= 1.e-9_R8P)) nfail = nfail + 1
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
   real(R8P),    parameter  :: GM = 1.4_R8P, EP = SETH_POSITIVITY_EPS
   real(R8P)                :: ql(nv,NL), qr(nv,NL), fh(nv,NL), fl(nv,NL), lam(NL), fd(nv,NL), td(NL)
   real(R8P)                :: fc(nv,NL), tc(NL)
   real(R8P)                :: a, b, worst, mag
   integer(I8P)             :: st
   integer(I4P)             :: i, k, nbit, nlim, nfail
   st = 20261115_I8P
   do i=1, NL
      call lim_state(nv, st, ql(:,i)) ; call lim_state(nv, st, qr(:,i))
      a = max(lim_speed(nv, ql(:,i)), lim_speed(nv, qr(:,i))) ; lam(i) = 1._R8P/a
      fl(:,i) = 0.5_R8P*(lim_flux(nv, ql(:,i)) + lim_flux(nv, qr(:,i))) - 0.5_R8P*a*(qr(:,i) - ql(:,i))
      do k=1, nv
         b = lim_rnd(st)
         fh(k,i) = fl(k,i) + (2._R8P*b - 1._R8P)*(abs(fl(k,i)) + abs(ql(k,i)) + abs(qr(k,i))) &
                   *10._R8P**(3._R8P*lim_rnd(st) - 2._R8P)
      enddo
   enddo
   !$acc parallel loop gang vector copyin(ql, qr, fh, fl, lam) copyout(fd, td)
   do i=1, NL
      call seth_positivity_euler(nv, GM, lam(i), EP, EP, ql(:,i), qr(:,i), fh(:,i), fl(:,i), fd(:,i), td(i))
   enddo
   do i=1, NL
      call seth_positivity_euler(nv, GM, lam(i), EP, EP, ql(:,i), qr(:,i), fh(:,i), fl(:,i), fc(:,i), tc(i))
   enddo
   nbit = 0 ; nlim = 0 ; nfail = 0 ; worst = 0._R8P
   do i=1, NL
      if (tc(i) < 1._R8P) nlim = nlim + 1
      if (td(i) == tc(i) .and. all(fd(:,i) == fc(:,i))) nbit = nbit + 1
      mag = max(maxval(abs(fh(:,i))), maxval(abs(fl(:,i))), tiny(1._R8P))
      worst = max(worst, abs(td(i) - tc(i)), maxval(abs(fd(:,i) - fc(:,i)))/mag)
      if (.not.(abs(td(i) - tc(i)) <= 1.e-10_R8P .and. maxval(abs(fd(:,i) - fc(:,i))) <= 1.e-10_R8P*mag)) &
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
   real(R8P)                   :: r
   st = mod(16807_I8P*st, 2147483647_I8P) ; r = real(st, R8P)/2147483647._R8P
   endfunction lim_rnd

   subroutine lim_state(nv, st, q)
   !< rho and p log-uniform in [1e-8, 1], velocities in [-10, 10] (gamma = 1.4).
   integer(I4P), intent(in)    :: nv
   integer(I8P), intent(inout) :: st
   real(R8P),    intent(out)   :: q(nv)
   real(R8P)                   :: rho, p
   integer(I4P)                :: kk
   rho = 10._R8P**(-8._R8P*lim_rnd(st)) ; p = 10._R8P**(-8._R8P*lim_rnd(st))
   q(1) = rho
   do kk=2, nv-1
      q(kk) = rho*(20._R8P*lim_rnd(st) - 10._R8P)
   enddo
   q(nv) = p/0.4_R8P + 0.5_R8P*dot_product(q(2:nv-1), q(2:nv-1))/rho
   endsubroutine lim_state

   pure function lim_flux(nv, q) result(fx)
   !< The Euler flux normal to direction 1 (gamma = 1.4).
   integer(I4P), intent(in) :: nv
   real(R8P),    intent(in) :: q(nv)
   real(R8P)                :: fx(nv), u, p
   u = q(2)/q(1) ; p = 0.4_R8P*(q(nv) - 0.5_R8P*dot_product(q(2:nv-1), q(2:nv-1))/q(1))
   fx = u*q ; fx(2) = fx(2) + p ; fx(nv) = fx(nv) + p*u
   endfunction lim_flux

   pure logical function lim_adm_ro(nv, q, mag) result(ok)
   !< rho and p above the floor 1e-13 to the round-off of the one-sided state q (mag = |q_cell| + |lam f|): the
   !< allowance of seth_test_positivity's adm_ro, gamma = 1.4.
   integer(I4P), intent(in) :: nv
   real(R8P),    intent(in) :: q(nv), mag(nv)
   real(R8P)                :: ke, tol, d(nv)
   integer(I4P)             :: k
   ok = .false.
   d = 4._R8P*epsilon(1._R8P)*mag
   if (.not.(q(1) >= 1.e-13_R8P - d(1))) return
   ke = 0._R8P ; tol = d(nv)
   do k=2, nv-1
      ke = ke + q(k)*q(k)
      tol = tol + abs(q(k)/q(1))*d(k)
   enddo
   ke = 0.5_R8P*ke/q(1)
   tol = 0.4_R8P*(tol + ke/q(1)*d(1) + 4._R8P*epsilon(1._R8P)*(abs(q(nv)) + ke))
   ok = 0.4_R8P*(q(nv) - ke) >= 1.e-13_R8P - tol
   endfunction lim_adm_ro

   pure function lim_speed(nv, q) result(c)
   !< |u_1| + c (gamma = 1.4).
   integer(I4P), intent(in) :: nv
   real(R8P),    intent(in) :: q(nv)
   real(R8P)                :: c, p
   p = 0.4_R8P*(q(nv) - 0.5_R8P*dot_product(q(2:nv-1), q(2:nv-1))/q(1))
   c = abs(q(2)/q(1)) + sqrt(1.4_R8P*p/q(1))
   endfunction lim_speed

   subroutine stencils(v, w)
   !< N stencil pairs, exact in binary: smooth cubics, jumps, an integer LCG, cycling.
   real(R8P), intent(out) :: v(-4:5,N), w(-4:5,N)
   integer(I8P)           :: state
   integer(I4P)           :: i, m, j
   state = 20261115_I8P
   do i=1, N
      j = mod(i, 16)
      do m=-4, 5
         select case(mod(i, 3))
         case(0)
            v(m,i) = 1._R8P + (j/8._R8P)*m - (j/64._R8P)*m*m + (1._R8P/512._R8P)*m**3
            w(m,i) = -0.5_R8P + (j/16._R8P)*m + (j/128._R8P)*m*m
         case(1)
            v(m,i) = merge(1._R8P, 0._R8P, m >= mod(i, 10) - 4)
            w(m,i) = merge(-2._R8P, 0.5_R8P, m >= mod(i, 9) - 4)
         case default
            state = mod(16807_I8P*state, 2147483647_I8P) ; v(m,i) = real(state, R8P)/2147483647._R8P - 0.5_R8P
            state = mod(16807_I8P*state, 2147483647_I8P) ; w(m,i) = real(state, R8P)/2147483647._R8P - 0.3_R8P
         endselect
      enddo
   enddo
   endsubroutine stencils
endprogram seth_test_device
