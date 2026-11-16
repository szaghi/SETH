!< SETH benchmark: SETH's compile-time-constant WENO kernels against ADAM's table-driven FNL kernels.
program seth_bench_weno
!< SETH benchmark: SETH's compile-time-constant WENO kernels against ADAM's table-driven FNL kernels.
!<
!< n faces of a 1-D field (cells 1-NG..n+NG; face i between cells i and i+1), two biased fields fp, fm (a smooth wave
!< with a jump every 2048 cells). Every variant reads the same stencils and writes the two states of every face:
!<
!<  1 SETH seth_weno5_js            (coefficients compile-time constants)
!<  2 ADAM weno_reconstruct_upwind_dev, S=3   (run-time S, tables passed as assumed-shape arrays)
!<  3 ADAM weno_reconstruct_upwind_s3_dev     (S a constant, tables still arrays)
!<  4-6 the same at WENO-7 (S=4)
!<  7 SETH seth_wenoz5 (p=1), 8 SETH seth_teno6 (C_T=1e-7): for reference
!<  9, 10 ADAM's host weno_reconstruct_upwind (what FLUME's CPU backend calls), S=3 and S=4: CPU builds only (no
!<        acc routine); in an OpenACC build they are skipped
!< 11.. SETH kernels added after v0.1.0, each against its cost budget (roadmap issue #4, gate PERF1): 11 seth_teno5,
!<      12 seth_teno5_opt, 13 seth_teno5_az (rmu = 1), 14 seth_teno8, 15 seth_teno5_a, 16 seth_teno8_a,
!<      17 seth_weno5_js with seth_weno_descaling, 18 with seth_weno_descaling_ai, 19 seth_weno5_js + seth_mp_bound,
!<      20 seth_linear5, 21 seth_linear7
!<
!< ADAM's variants pack ADAM's rows as FLUME does (v(2,m) = fp(i+m), v(1,m) = fm(i+m+1)). Built with -acc the loops run
!< on the device (timed with !$acc wait + system_clock); without, on one CPU core. Each variant: 2 warm-up passes, then
!< `reps` timed passes; the median and the minimum are reported. Checks: the SETH results equal the host SETH results
!< on a sample of faces, and every ADAM variant equals the SETH one of its order, to round-off (64 eps of the data).
!<
!< Usage: seth_bench_weno [n_faces [reps]]   (defaults 16777216 and 20)
use, intrinsic :: iso_fortran_env, only : int64
use seth
use seth_weno_js_tables
use adam_fnl_weno_kernels,       only : weno_reconstruct_upwind_dev
use adam_fnl_weno_width_kernels, only : weno_reconstruct_upwind_s3_dev, weno_reconstruct_upwind_s4_dev
use seth_adam_reference,         only : weno_reconstruct_upwind
implicit none
integer(I4P), parameter :: NG = 4_I4P         !< Ghost cells (WENO-7 reads i-3..i+4).
integer(I4P), parameter :: NV = 21_I4P        !< Variants.
integer(I4P), parameter :: SAMPLE = 65536_I4P !< Faces checked against the host.
real(R8P),    parameter :: TOL = 64._R8P*epsilon(1._R8P)
character(31), parameter :: NAMES(NV) = [character(31) :: 'SETH  weno5_js', 'ADAM  generic S=3', 'ADAM  width s3', &
                                         'SETH  weno7_js', 'ADAM  generic S=4', 'ADAM  width s4', &
                                         'SETH  wenoz5 (p=1)', 'SETH  teno6 (C_T=1e-7)', &
                                         'ADAM  host S=3 (FLUME CPU)', 'ADAM  host S=4 (FLUME CPU)', &
                                         'SETH  teno5 (C_T=1e-5)', 'SETH  teno5_opt (C_T=1e-5)', &
                                         'SETH  teno5_az (C_T=1e-5)', 'SETH  teno8 (C_T=1e-7)', &
                                         'SETH  teno5_a', 'SETH  teno8_a', &
                                         'SETH  weno5_js + descaling', 'SETH  weno5_js + descaling_ai', &
                                         'SETH  weno5_js + mp_bound', 'SETH  linear5', 'SETH  linear7']
real(R8P), allocatable :: fp(:), fm(:), rp(:), rm(:), rp5(:), rm5(:), rp7(:), rm7(:)
real(R8P), allocatable :: a3(:,:,:), p3(:,:,:,:), d3(:,:,:,:), a4(:,:,:), p4(:,:,:,:), d4(:,:,:,:)
real(R8P), allocatable :: times(:)
real(R8P)              :: med(NV), tmin(NV), diff, vmax
integer(I4P)           :: n, reps, v, r, i, c
integer(int64)         :: t0, t1, rate
character(32)          :: arg
logical                :: passed

n = 16777216_I4P ; reps = 20_I4P
if (command_argument_count() >= 1) then ; call get_command_argument(1, arg) ; read(arg,*) n ; endif
if (command_argument_count() >= 2) then ; call get_command_argument(2, arg) ; read(arg,*) reps ; endif
allocate(fp(1-NG:n+NG), fm(1-NG:n+NG), rp(n), rm(n), rp5(n), rm5(n), rp7(n), rm7(n), times(reps))
do c=1-NG, n+NG
   fp(c) = 1._R8P + 0.5_R8P*sin(0.01_R8P*c) + merge(1._R8P, 0._R8P, mod(c+1000000, 4096) < 2048)
   fm(c) = 0.5_R8P*cos(0.013_R8P*c) - merge(0.7_R8P, 0._R8P, mod(c+1000000, 4096) < 2048)
enddo
allocate(a3(1:2,0:2,1:3), p3(1:2,0:2,0:2,1:3), d3(0:2,0:2,0:2,1:3))
allocate(a4(1:2,0:3,1:4), p4(1:2,0:3,0:3,1:4), d4(0:3,0:3,0:3,1:4))
a3 = 0._R8P ; p3 = 0._R8P ; d3 = 0._R8P ; a4 = 0._R8P ; p4 = 0._R8P ; d4 = 0._R8P
a3(:,:,3) = WENO5_JS_A ; p3(:,:,:,3) = WENO5_JS_P ; d3(:,:,:,3) = WENO5_JS_D
a4(:,:,4) = WENO7_JS_A ; p4(:,:,:,4) = WENO7_JS_P ; d4(:,:,:,4) = WENO7_JS_D
!$acc enter data copyin(fp, fm, a3, p3, d3, a4, p4, d4) create(rp, rm)

call system_clock(count_rate=rate)
passed = .true.
print '(A,I0,A,I0,A)', 'faces: ', n, ', timed passes: ', reps, ' per variant'
med = -1._R8P ; tmin = -1._R8P
do v=1, NV
#ifdef _OPENACC
   if (v == 9 .or. v == 10) cycle ! host routines: no acc routine
#endif
   do r=1, 2
      call run(v)
   enddo
   do r=1, reps
      call system_clock(t0)
      call run(v)
      !$acc wait
      call system_clock(t1)
      times(r) = real(t1 - t0, R8P)/real(rate, R8P)
   enddo
   call sort(times)
   med(v) = times((reps+1)/2) ; tmin(v) = times(1)
   !$acc update self(rp, rm)
   call check(v)
enddo
!$acc exit data delete(fp, fm, a3, p3, d3, a4, p4, d4, rp, rm)
print '(A)', ''
print '(A31,4A13)', 'variant', 'median ms', 'min ms', 'ns/face', 'vs ADAM gen.'
do v=1, NV
   if (med(v) < 0._R8P) cycle
   print '(A31,3F13.3,F13.2)', NAMES(v), 1.e3_R8P*med(v), 1.e3_R8P*tmin(v), 1.e9_R8P*med(v)/n, &
                               med(merge(2, 5, v <= 3 .or. v == 7 .or. v == 9))/med(v)
enddo
if (.not.passed) error stop 'a benchmark check failed'
print '(A)', 'seth_bench_weno: checks ok'

contains
   subroutine run(v)
   !< One pass of variant v over all faces.
   integer(I4P), intent(in) :: v
   real(R8P)                :: vp5(-2:3), vm5(-2:3), vp7(-3:4), vm7(-3:4), w3(1:2,-2:2), w4(1:2,-3:3), vr(1:2)
   integer(I4P)             :: i, m
   select case(v)
   case(1)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp5, vm5)
      do i=1, n
         do m=-2, 3
            vp5(m) = fp(i+m) ; vm5(m) = fm(i+m)
         enddo
         call seth_weno5_js(vp5, vm5, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
   case(2)
      !$acc parallel loop gang vector present(fp, fm, rp, rm, a3, p3, d3) private(w3, vr)
      do i=1, n
         do m=-2, 2
            w3(1,m) = fm(i+m+1) ; w3(2,m) = fp(i+m)
         enddo
         call weno_reconstruct_upwind_dev(S=3_I4P, weno_a=a3, weno_p=p3, weno_d=d3, weno_zeps=SETH_WENO_JS_EPS, &
                                          weno_rmu=1._R8P, V=w3, VR=vr)
         rp(i) = vr(2) ; rm(i) = vr(1)
      enddo
   case(3)
      !$acc parallel loop gang vector present(fp, fm, rp, rm, a3, p3, d3) private(w3, vr)
      do i=1, n
         do m=-2, 2
            w3(1,m) = fm(i+m+1) ; w3(2,m) = fp(i+m)
         enddo
         call weno_reconstruct_upwind_s3_dev(weno_a=a3, weno_p=p3, weno_d=d3, weno_zeps=SETH_WENO_JS_EPS, &
                                             weno_rmu=1._R8P, V=w3, VR=vr)
         rp(i) = vr(2) ; rm(i) = vr(1)
      enddo
   case(4)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp7, vm7)
      do i=1, n
         do m=-3, 4
            vp7(m) = fp(i+m) ; vm7(m) = fm(i+m)
         enddo
         call seth_weno7_js(vp7, vm7, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
   case(5)
      !$acc parallel loop gang vector present(fp, fm, rp, rm, a4, p4, d4) private(w4, vr)
      do i=1, n
         do m=-3, 3
            w4(1,m) = fm(i+m+1) ; w4(2,m) = fp(i+m)
         enddo
         call weno_reconstruct_upwind_dev(S=4_I4P, weno_a=a4, weno_p=p4, weno_d=d4, weno_zeps=SETH_WENO_JS_EPS, &
                                          weno_rmu=1._R8P, V=w4, VR=vr)
         rp(i) = vr(2) ; rm(i) = vr(1)
      enddo
   case(6)
      !$acc parallel loop gang vector present(fp, fm, rp, rm, a4, p4, d4) private(w4, vr)
      do i=1, n
         do m=-3, 3
            w4(1,m) = fm(i+m+1) ; w4(2,m) = fp(i+m)
         enddo
         call weno_reconstruct_upwind_s4_dev(weno_a=a4, weno_p=p4, weno_d=d4, weno_zeps=SETH_WENO_JS_EPS, &
                                             weno_rmu=1._R8P, V=w4, VR=vr)
         rp(i) = vr(2) ; rm(i) = vr(1)
      enddo
   case(7)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp5, vm5)
      do i=1, n
         do m=-2, 3
            vp5(m) = fp(i+m) ; vm5(m) = fm(i+m)
         enddo
         call seth_wenoz5(vp5, vm5, 1.e-40_R8P, 1_I4P, rp(i), rm(i))
      enddo
   case(8)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp5, vm5)
      do i=1, n
         do m=-2, 3
            vp5(m) = fp(i+m) ; vm5(m) = fm(i+m)
         enddo
         call seth_teno6(vp5, vm5, SETH_TENO_EPS, 1.e-7_R8P, rp(i), rm(i))
      enddo
   case(11)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp5, vm5)
      do i=1, n
         do m=-2, 3
            vp5(m) = fp(i+m) ; vm5(m) = fm(i+m)
         enddo
         call seth_teno5(vp5, vm5, SETH_TENO_EPS, SETH_TENO5_CT, rp(i), rm(i))
      enddo
   case(12)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp5, vm5)
      do i=1, n
         do m=-2, 3
            vp5(m) = fp(i+m) ; vm5(m) = fm(i+m)
         enddo
         call seth_teno5_opt(vp5, vm5, SETH_TENO_EPS, SETH_TENO5_OPT_CT, rp(i), rm(i))
      enddo
   case(13)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp5, vm5)
      do i=1, n
         do m=-2, 3
            vp5(m) = fp(i+m) ; vm5(m) = fm(i+m)
         enddo
         call seth_teno5_az(vp5, vm5, SETH_TENO_EPS, SETH_TENO5_AZ_CT, 1._R8P, rp(i), rm(i))
      enddo
   case(14)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp7, vm7)
      do i=1, n
         do m=-3, 4
            vp7(m) = fp(i+m) ; vm7(m) = fm(i+m)
         enddo
         call seth_teno8(vp7, vm7, SETH_TENO_EPS, SETH_TENO8_CT, rp(i), rm(i))
      enddo
   case(15)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp5, vm5)
      do i=1, n
         do m=-2, 3
            vp5(m) = fp(i+m) ; vm5(m) = fm(i+m)
         enddo
         call seth_teno5_a(vp5, vm5, SETH_TENO_EPS, SETH_TENO5_A_CR, SETH_TENO5_A_XI, SETH_TENO5_A_A1, SETH_TENO5_A_A2, &
                           1._R8P, rp(i), rm(i))
      enddo
   case(16)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp7, vm7)
      do i=1, n
         do m=-3, 4
            vp7(m) = fp(i+m) ; vm7(m) = fm(i+m)
         enddo
         call seth_teno8_a(vp7, vm7, SETH_TENO_EPS, SETH_TENO8_A_CR, SETH_TENO8_A_XI, SETH_TENO8_A_A1, SETH_TENO8_A_A2, &
                           1._R8P, rp(i), rm(i))
      enddo
   case(17)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp5, vm5)
      do i=1, n
         do m=-2, 3
            vp5(m) = fp(i+m) ; vm5(m) = fm(i+m)
         enddo
         call seth_weno5_js(vp5, vm5, SETH_WENO_JS_EPS, seth_weno_descaling(3_I4P, 1._R8P, vp5, vm5), rp(i), rm(i))
      enddo
   case(18)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp5, vm5)
      do i=1, n
         do m=-2, 3
            vp5(m) = fp(i+m) ; vm5(m) = fm(i+m)
         enddo
         call seth_weno5_js(vp5, vm5, SETH_WENO_JS_EPS, seth_weno_descaling_ai(3_I4P, vp5, vm5), rp(i), rm(i))
      enddo
   case(19)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp5, vm5)
      do i=1, n
         do m=-2, 3
            vp5(m) = fp(i+m) ; vm5(m) = fm(i+m)
         enddo
         call seth_weno5_js(vp5, vm5, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
         call seth_mp_bound(vp5, vm5, SETH_MP_ALPHA, rp(i), rm(i))
      enddo
   case(20)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp5, vm5)
      do i=1, n
         do m=-2, 3
            vp5(m) = fp(i+m) ; vm5(m) = fm(i+m)
         enddo
         call seth_linear5(vp5, vm5, rp(i), rm(i))
      enddo
   case(21)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp7, vm7)
      do i=1, n
         do m=-3, 4
            vp7(m) = fp(i+m) ; vm7(m) = fm(i+m)
         enddo
         call seth_linear7(vp7, vm7, rp(i), rm(i))
      enddo
   case(9)
      do i=1, n
         do m=-2, 2
            w3(1,m) = fm(i+m+1) ; w3(2,m) = fp(i+m)
         enddo
         call weno_reconstruct_upwind(3_I4P, a3, p3, d3, SETH_WENO_JS_EPS, 1._R8P, w3, vr)
         rp(i) = vr(2) ; rm(i) = vr(1)
      enddo
   case(10)
      do i=1, n
         do m=-3, 3
            w4(1,m) = fm(i+m+1) ; w4(2,m) = fp(i+m)
         enddo
         call weno_reconstruct_upwind(4_I4P, a4, p4, d4, SETH_WENO_JS_EPS, 1._R8P, w4, vr)
         rp(i) = vr(2) ; rm(i) = vr(1)
      enddo
   endselect
   endsubroutine run

   subroutine check(v)
   !< SETH variants against host SETH on a sample; ADAM variants against the SETH variant of their order.
   integer(I4P), intent(in) :: v
   real(R8P)                :: vp5(-2:3), vm5(-2:3), vp7(-3:4), vm7(-3:4), hp, hm
   integer(I4P)             :: i, m, k, stride
   diff = 0._R8P
   select case(v)
   case(1, 4, 7, 8, 11:)
      stride = max(1_I4P, n/SAMPLE)
      do k=1, n, stride
         i = k
         do m=-2, 3
            vp5(m) = fp(i+m) ; vm5(m) = fm(i+m)
         enddo
         do m=-3, 4
            vp7(m) = fp(i+m) ; vm7(m) = fm(i+m)
         enddo
         select case(v)
         case(1) ; call seth_weno5_js(vp5, vm5, SETH_WENO_JS_EPS, 1._R8P, hp, hm)
         case(4) ; call seth_weno7_js(vp7, vm7, SETH_WENO_JS_EPS, 1._R8P, hp, hm)
         case(7) ; call seth_wenoz5(vp5, vm5, 1.e-40_R8P, 1_I4P, hp, hm)
         case(8) ; call seth_teno6(vp5, vm5, SETH_TENO_EPS, 1.e-7_R8P, hp, hm)
         case(11) ; call seth_teno5(vp5, vm5, SETH_TENO_EPS, SETH_TENO5_CT, hp, hm)
         case(12) ; call seth_teno5_opt(vp5, vm5, SETH_TENO_EPS, SETH_TENO5_OPT_CT, hp, hm)
         case(13) ; call seth_teno5_az(vp5, vm5, SETH_TENO_EPS, SETH_TENO5_AZ_CT, 1._R8P, hp, hm)
         case(14) ; call seth_teno8(vp7, vm7, SETH_TENO_EPS, SETH_TENO8_CT, hp, hm)
         case(15) ; call seth_teno5_a(vp5, vm5, SETH_TENO_EPS, SETH_TENO5_A_CR, SETH_TENO5_A_XI, SETH_TENO5_A_A1, &
                                      SETH_TENO5_A_A2, 1._R8P, hp, hm)
         case(17) ; call seth_weno5_js(vp5, vm5, SETH_WENO_JS_EPS, seth_weno_descaling(3_I4P, 1._R8P, vp5, vm5), hp, hm)
         case(18) ; call seth_weno5_js(vp5, vm5, SETH_WENO_JS_EPS, seth_weno_descaling_ai(3_I4P, vp5, vm5), hp, hm)
         case(19) ; call seth_weno5_js(vp5, vm5, SETH_WENO_JS_EPS, 1._R8P, hp, hm)
                     call seth_mp_bound(vp5, vm5, SETH_MP_ALPHA, hp, hm)
         case(20) ; call seth_linear5(vp5, vm5, hp, hm)
         case(21) ; call seth_linear7(vp7, vm7, hp, hm)
         case(16) ; call seth_teno8_a(vp7, vm7, SETH_TENO_EPS, SETH_TENO8_A_CR, SETH_TENO8_A_XI, SETH_TENO8_A_A1, &
                                      SETH_TENO8_A_A2, 1._R8P, hp, hm)
         endselect
         vmax = max(maxval(abs(vp7)), maxval(abs(vm7)))
         diff = max(diff, abs(rp(i) - hp)/vmax, abs(rm(i) - hm)/vmax)
      enddo
      if (v == 1) then ; rp5 = rp ; rm5 = rm ; endif
      if (v == 4) then ; rp7 = rp ; rm7 = rm ; endif
      print '(A31,A,ES9.2)', NAMES(v), ': device vs host SETH, max scaled difference ', diff
   case(2, 3, 9)
      do i=1, n
         diff = max(diff, abs(rp(i) - rp5(i))/3._R8P, abs(rm(i) - rm5(i))/3._R8P)
      enddo
      print '(A31,A,ES9.2)', NAMES(v), ': vs SETH weno5_js, max scaled difference ', diff
   case(5, 6, 10)
      do i=1, n
         diff = max(diff, abs(rp(i) - rp7(i))/3._R8P, abs(rm(i) - rm7(i))/3._R8P)
      enddo
      print '(A31,A,ES9.2)', NAMES(v), ': vs SETH weno7_js, max scaled difference ', diff
   endselect
   if (diff > TOL) passed = .false.
   endsubroutine check

   subroutine sort(x)
   !< Insertion sort (a few tens of values).
   real(R8P), intent(inout) :: x(:)
   real(R8P)                :: t
   integer(I4P)             :: i, j
   do i=2, size(x)
      t = x(i) ; j = i - 1
      do while (j >= 1)
         if (x(j) <= t) exit
         x(j+1) = x(j) ; j = j - 1
      enddo
      x(j+1) = t
   enddo
   endsubroutine sort
endprogram seth_bench_weno
