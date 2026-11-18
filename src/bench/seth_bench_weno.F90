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
!<      20 seth_linear5, 21 seth_linear7, 22 seth_sensor_lcr alone (s stored in rp, rm), 23 the cost hybrid of WENO-Z5:
!<      seth_sensor_lcr, then seth_linear5 where s = 1 and seth_wenoz5 elsewhere (the flagged fraction of this field),
!<      24 seth_keep_pe6 and 25 ADAM's compute_kep_flux on the state rho = 1 + fp^2, u = fm, v = fp, w = 0.1, p = 1 + fm^2
!<      (the flux of direction 1 in rp, its energy component in rm); 26, 27, 28 seth_weno4_c, seth_weno6_c,
!<      seth_weno8_c (the centred WENO, rmu = 1), 29 seth_weno3_js (the upwind part of seth_weno4_c, its budget's
!<      reference)
!< 30.. issue #6 (non-polynomial reconstructions): 30 seth_teno5_mask of both biases (the mask of vp(-2:2) and of the
!<      mirrored vm(1-m), m = -2..2, in rp and rm), against its budget of 0.70 seth_teno5; 31 seth_teno5_dr and
!<      32 seth_teno5_thinc (self-contained, -5:6, the papers' constants), against 3.0 seth_teno5 (owner decision D3);
!<      33, 34 the same without the early exit (every mask of cells -3..2 or -3..3 of both states computed, then the
!<      library's detector and override: the work the early exit saves); 35, 36 the stored-mask variant (notes/
!<      kernel-api-design.md, section 11.4), the library's pair since 2026-11-18 (issue #6, step 2.6): a first pass,
!<      seth_teno5_face_masks at every face, stores the left mask of cell i and the right mask of cell i+1, a second
!<      reconstructs each face with seth_teno5_dr_stored / seth_teno5_thinc_stored; both passes timed. On this scalar
!<      field its states equal the self-contained kernels' (checked);
!<      37 seth_teno5_st (budget 1.05 seth_teno5), 38 seth_teno5_rbf with s1 = s2 = 0.09 (budget 1.25), 39 the
!<      complete scheme seth_teno5_rbf_dr (s1 = s2 = 0.09; budget 3.0, the self-contained kernels')
!< 40.. issue #6, the pieces alone (2026-11-18): 40 seth_teno5_mask of one bias (the budget's unit: 0.70 seth_teno5,
!<      which reconstructs both states), 41 seth_dr_detect and 42 seth_thinc_detect on the stored masks of variants 35
!<      and 36 (budget 0.10), 43 seth_dr_override at level 1 on every face (0.15) and 44 seth_thinc_override at level 1
!<      on every face (0.30), the states preset to the mean of f_i and f_{i+1}; 45 seth_teno5_rbf_dr without the early
!<      exit and 46 with stored masks (seth_teno5_face_masks, seth_teno5_rbf_dr_stored; step 4.1)
!<
!< ADAM's variants pack ADAM's rows as FLUME does (v(2,m) = fp(i+m), v(1,m) = fm(i+m+1)). Built with -acc the loops run
!< on the device (timed with !$acc wait + system_clock); without, on one CPU core. Each variant: 2 warm-up passes, then
!< `reps` timed passes; the median and the minimum are reported. Checks: the SETH results equal the host SETH results
!< on a sample of faces, and every ADAM variant equals the SETH one of its order, to round-off (64 eps of the data).
!<
!< Usage: seth_bench_weno [n_faces [reps [jump_spacing]]]   (defaults 16777216, 20 and 2048: the field jumps every
!< jump_spacing cells; 0: no jump, a smooth field; the composites' cost against the fraction of faces near a jump, PERF7)
use, intrinsic :: iso_fortran_env, only : int64
use seth
use seth_weno_js_tables
use adam_fnl_weno_kernels,       only : weno_reconstruct_upwind_dev
use adam_fnl_weno_width_kernels, only : weno_reconstruct_upwind_s3_dev, weno_reconstruct_upwind_s4_dev
use seth_adam_reference,         only : weno_reconstruct_upwind, compute_kep_flux
implicit none
integer(I4P), parameter :: NG = 6_I4P         !< Ghost cells (the composites of issue #6 read i-5..i+6).
integer(I4P), parameter :: NV = 46_I4P        !< Variants.
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
                                         'SETH  weno5_js + mp_bound', 'SETH  linear5', 'SETH  linear7', &
                                         'SETH  sensor_lcr', 'SETH  hybrid lcr+linear5/wenoz5', &
                                         'SETH  keep_pe6', 'ADAM  compute_kep_flux', &
                                         'SETH  weno4_c', 'SETH  weno6_c', 'SETH  weno8_c', 'SETH  weno3_js', &
                                         'SETH  teno5_mask (both biases)', 'SETH  teno5_dr', 'SETH  teno5_thinc', &
                                         'SETH  teno5_dr, no early exit', 'SETH  teno5_thinc, no early exit', &
                                         'SETH  teno5_dr, stored masks', 'SETH  teno5_thinc, stored masks', &
                                         'SETH  teno5_st', 'SETH  teno5_rbf (s1 = s2 = 0.09)', 'SETH  teno5_rbf_dr', &
                                         'SETH  teno5_mask (one bias)', 'SETH  dr_detect (stored masks)', &
                                         'SETH  thinc_detect (stored masks)', 'SETH  dr_override (level 1)', &
                                         'SETH  thinc_override (level 1)', 'SETH  teno5_rbf_dr, no early exit', &
                                         'SETH  teno5_rbf_dr, stored masks']
real(R8P), allocatable :: fp(:), fm(:), rp(:), rm(:), rp5(:), rm5(:), rp7(:), rm7(:)
real(R8P), allocatable :: a3(:,:,:), p3(:,:,:,:), d3(:,:,:,:), a4(:,:,:), p4(:,:,:,:), d4(:,:,:,:)
real(R8P), allocatable :: times(:)
integer(I4P), allocatable :: ml(:), mr(:) !< Stored masks of the stored-mask variants (cells -3..n+4).
real(R8P)              :: med(NV), tmin(NV), diff, vmax
integer(I4P)           :: n, reps, v, r, i, c, js
integer(int64)         :: t0, t1, rate
character(32)          :: arg
logical                :: passed

n = 16777216_I4P ; reps = 20_I4P
if (command_argument_count() >= 1) then ; call get_command_argument(1, arg) ; read(arg,*) n ; endif
if (command_argument_count() >= 2) then ; call get_command_argument(2, arg) ; read(arg,*) reps ; endif
js = 2048_I4P
if (command_argument_count() >= 3) then ; call get_command_argument(3, arg) ; read(arg,*) js ; endif
allocate(fp(1-NG:n+NG), fm(1-NG:n+NG), rp(n), rm(n), rp5(n), rm5(n), rp7(n), rm7(n), times(reps))
allocate(ml(-3:n+4), mr(-3:n+4))
do c=1-NG, n+NG
   fp(c) = 1._R8P + 0.5_R8P*sin(0.01_R8P*c)
   fm(c) = 0.5_R8P*cos(0.013_R8P*c)
   if (js > 0) then
      fp(c) = fp(c) + merge(1._R8P, 0._R8P, mod(c+1000000, 2*js) < js)
      fm(c) = fm(c) - merge(0.7_R8P, 0._R8P, mod(c+1000000, 2*js) < js)
   endif
enddo
allocate(a3(1:2,0:2,1:3), p3(1:2,0:2,0:2,1:3), d3(0:2,0:2,0:2,1:3))
allocate(a4(1:2,0:3,1:4), p4(1:2,0:3,0:3,1:4), d4(0:3,0:3,0:3,1:4))
a3 = 0._R8P ; p3 = 0._R8P ; d3 = 0._R8P ; a4 = 0._R8P ; p4 = 0._R8P ; d4 = 0._R8P
a3(:,:,3) = WENO5_JS_A ; p3(:,:,:,3) = WENO5_JS_P ; d3(:,:,:,3) = WENO5_JS_D
a4(:,:,4) = WENO7_JS_A ; p4(:,:,:,4) = WENO7_JS_P ; d4(:,:,:,4) = WENO7_JS_D
!$acc enter data copyin(fp, fm, a3, p3, d3, a4, p4, d4) create(rp, rm, ml, mr)

call system_clock(count_rate=rate)
passed = .true.
print '(A,I0,A,I0,A,I0,A)', 'faces: ', n, ', timed passes: ', reps, ' per variant, a jump every ', js, ' cells (0: none)'
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
!$acc exit data delete(fp, fm, a3, p3, d3, a4, p4, d4, rp, rm, ml, mr)
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
   real(R8P)                :: ws(5,-2:3), fk(5), vp4(-1:2), vm4(-1:2), wa(-2:2), wb(-2:2), vpw(-5:6), vmw(-5:6)
   integer(I4P)             :: i, m, mka, mkb, mpw(-3:3), mmw(-3:3), lpw, lmw
   select case(v)
   case(31, 32)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vpw, vmw)
      do i=1, n
         do m=-5, 6
            vpw(m) = fp(i+m) ; vmw(m) = fm(i+m)
         enddo
         if (v == 31) then
            call seth_teno5_dr(vpw, vmw, SETH_TENO_EPS, SETH_TENO5_CT, SETH_DR_T1, SETH_DR_T2, 1._R8P, rp(i), rm(i))
         else
            call seth_teno5_thinc(vpw, vmw, SETH_TENO_EPS, SETH_TENO5_CT, SETH_THINC_BETA, 1._R8P, rp(i), rm(i))
         endif
      enddo
   case(39)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vpw, vmw)
      do i=1, n
         do m=-5, 6
            vpw(m) = fp(i+m) ; vmw(m) = fm(i+m)
         enddo
         call seth_teno5_rbf_dr(vpw, vmw, SETH_TENO_EPS, SETH_TENO5_CT, SETH_RBF_CR, 0.09_R8P, 0.09_R8P, SETH_DR_T1, &
                                SETH_DR_T2, 1._R8P, rp(i), rm(i))
      enddo
   case(37, 38)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp5, vm5)
      do i=1, n
         do m=-2, 3
            vp5(m) = fp(i+m) ; vm5(m) = fm(i+m)
         enddo
         if (v == 37) then
            call seth_teno5_st(vp5, vm5, SETH_TENO_EPS, SETH_TENO5_CT, rp(i), rm(i))
         else
            call seth_teno5_rbf(vp5, vm5, SETH_TENO_EPS, SETH_TENO5_CT, SETH_RBF_CR, 0.09_R8P, 0.09_R8P, rp(i), rm(i))
         endif
      enddo
   case(35, 36)
      !$acc parallel loop gang vector present(fp, fm, ml, mr) private(vp5, vm5, mka, mkb)
      do i=-3, n+3 ! pass 1, seth_teno5_face_masks: at each face i+1/2, the left mask of cell i, the right of cell i+1
         do m=-2, 3
            vp5(m) = fp(i+m) ; vm5(m) = fm(i+m)
         enddo
         call seth_teno5_face_masks(vp5, vm5, SETH_TENO_EPS, SETH_TENO5_CT, mka, mkb)
         ml(i) = mka ; mr(i+1) = mkb
      enddo
      !$acc parallel loop gang vector present(fp, fm, rp, rm, ml, mr) private(vp5, vm5, mpw, mmw)
      do i=1, n ! pass 2, seth_teno5_dr_stored / seth_teno5_thinc_stored: the faces, from the stored masks
         do m=-2, 3
            vp5(m) = fp(i+m) ; vm5(m) = fm(i+m)
         enddo
         do m=-3, 3
            mpw(m) = ml(i+m) ; mmw(m) = mr(i+1-m)
         enddo
         if (v == 35) then
            call seth_teno5_dr_stored(vp5, vm5, mpw, mmw, SETH_TENO_EPS, SETH_TENO5_CT, SETH_DR_T1, SETH_DR_T2, &
                                      1._R8P, rp(i), rm(i))
         else
            call seth_teno5_thinc_stored(vp5, vm5, mpw, mmw, SETH_TENO_EPS, SETH_TENO5_CT, SETH_THINC_BETA, 1._R8P, &
                                         rp(i), rm(i))
         endif
      enddo
   case(33, 34)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vpw, vmw)
      do i=1, n
         do m=-5, 6
            vpw(m) = fp(i+m) ; vmw(m) = fm(i+m)
         enddo
         call composite_full(v == 34, vpw, vmw, rp(i), rm(i))
      enddo
   case(40)
      !$acc parallel loop gang vector present(fp, rp, rm) private(wa, mka)
      do i=1, n
         do m=-2, 2
            wa(m) = fp(i+m)
         enddo
         call seth_teno5_mask(wa, SETH_TENO_EPS, SETH_TENO5_CT, mka)
         rp(i) = real(mka, R8P) ; rm(i) = rp(i)
      enddo
   case(41, 42)
      !$acc parallel loop gang vector present(rp, rm, ml, mr) private(mpw, mmw, lpw, lmw)
      do i=1, n
         do m=-3, 3
            mpw(m) = ml(i+m) ; mmw(m) = mr(i+1-m)
         enddo
         if (v == 41) then
            call seth_dr_detect(mpw, mmw, lpw, lmw)
         else
            call seth_thinc_detect(mpw, mmw, lpw, lmw)
         endif
         rp(i) = real(lpw, R8P) ; rm(i) = real(lmw, R8P)
      enddo
   case(43, 44)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp4, vm4)
      do i=1, n
         do m=-1, 2
            vp4(m) = fp(i+m) ; vm4(m) = fm(i+m)
         enddo
         rp(i) = 0.5_R8P*(vp4(0) + vp4(1)) ; rm(i) = 0.5_R8P*(vm4(0) + vm4(1))
         if (v == 43) then
            call seth_dr_override(vp4, vm4, 1_I4P, 1_I4P, SETH_DR_T1, SETH_DR_T2, 1._R8P, rp(i), rm(i))
         else
            call seth_thinc_override(vp4, vm4, 1_I4P, 1_I4P, SETH_THINC_BETA, 1._R8P, rp(i), rm(i))
         endif
      enddo
   case(45)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vpw, vmw)
      do i=1, n
         do m=-5, 6
            vpw(m) = fp(i+m) ; vmw(m) = fm(i+m)
         enddo
         call rbf_dr_full(vpw, vmw, rp(i), rm(i))
      enddo
   case(46)
      !$acc parallel loop gang vector present(fp, fm, ml, mr) private(vp5, vm5, mka, mkb)
      do i=-3, n+3 ! pass 1, seth_teno5_face_masks: at each face i+1/2, the left mask of cell i, the right of cell i+1
         do m=-2, 3
            vp5(m) = fp(i+m) ; vm5(m) = fm(i+m)
         enddo
         call seth_teno5_face_masks(vp5, vm5, SETH_TENO_EPS, SETH_TENO5_CT, mka, mkb)
         ml(i) = mka ; mr(i+1) = mkb
      enddo
      !$acc parallel loop gang vector present(fp, fm, rp, rm, ml, mr) private(vp5, vm5, mpw, mmw)
      do i=1, n ! pass 2, seth_teno5_rbf_dr_stored: the faces, from the stored masks
         do m=-2, 3
            vp5(m) = fp(i+m) ; vm5(m) = fm(i+m)
         enddo
         do m=-3, 3
            mpw(m) = ml(i+m) ; mmw(m) = mr(i+1-m)
         enddo
         call seth_teno5_rbf_dr_stored(vp5, vm5, mpw, mmw, SETH_TENO_EPS, SETH_TENO5_CT, SETH_RBF_CR, 0.09_R8P, &
                                       0.09_R8P, SETH_DR_T1, SETH_DR_T2, 1._R8P, rp(i), rm(i))
      enddo
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
   case(22)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp5, vm5)
      do i=1, n
         do m=-2, 3
            vp5(m) = fp(i+m)
         enddo
         call seth_sensor_lcr(vp5, 1._R8P, SETH_SENSOR_LCR_PSIC, rp(i))
         rm(i) = rp(i)
      enddo
   case(23)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp5, vm5)
      do i=1, n
         do m=-2, 3
            vp5(m) = fp(i+m) ; vm5(m) = fm(i+m)
         enddo
         call hybrid5(vp5, vm5, rp(i), rm(i))
      enddo
   case(24, 25)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(ws, fk)
      do i=1, n
         do m=-2, 3
            ws(1,m) = 1._R8P + fp(i+m)**2 ; ws(2,m) = fm(i+m) ; ws(3,m) = fp(i+m) ; ws(4,m) = 0.1_R8P
            ws(5,m) = 1._R8P + fm(i+m)**2
         enddo
         if (v == 24) then
            call seth_keep_pe6(1_I4P, 1.4_R8P, ws, fk)
         else
            call compute_kep_flux(1_I4P, 1.4_R8P, ws, fk)
         endif
         rp(i) = fk(1) ; rm(i) = fk(5)
      enddo
   case(26)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp4, vm4)
      do i=1, n
         do m=-1, 2
            vp4(m) = fp(i+m) ; vm4(m) = fm(i+m)
         enddo
         call seth_weno4_c(vp4, vm4, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
   case(29)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp4, vm4)
      do i=1, n
         do m=-1, 2
            vp4(m) = fp(i+m) ; vm4(m) = fm(i+m)
         enddo
         call seth_weno3_js(vp4, vm4, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
   case(27)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp5, vm5)
      do i=1, n
         do m=-2, 3
            vp5(m) = fp(i+m) ; vm5(m) = fm(i+m)
         enddo
         call seth_weno6_c(vp5, vm5, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
   case(28)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp7, vm7)
      do i=1, n
         do m=-3, 4
            vp7(m) = fp(i+m) ; vm7(m) = fm(i+m)
         enddo
         call seth_weno8_c(vp7, vm7, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
   case(30)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(wa, wb, mka, mkb)
      do i=1, n
         do m=-2, 2
            wa(m) = fp(i+m) ; wb(m) = fm(i+1-m)
         enddo
         call seth_teno5_mask(wa, SETH_TENO_EPS, SETH_TENO5_CT, mka)
         call seth_teno5_mask(wb, SETH_TENO_EPS, SETH_TENO5_CT, mkb)
         rp(i) = real(mka, R8P) ; rm(i) = real(mkb, R8P)
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
   real(R8P)                :: vp5(-2:3), vm5(-2:3), vp7(-3:4), vm7(-3:4), hp, hm, ws(5,-2:3), fk(5)
   real(R8P)                :: vp4(-1:2), vm4(-1:2), wa(-2:2), wb(-2:2), vpw(-5:6), vmw(-5:6)
   integer(I4P)             :: i, m, k, stride, mka, mkb
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
         case(22) ; call seth_sensor_lcr(vp5, 1._R8P, SETH_SENSOR_LCR_PSIC, hp) ; hm = hp
         case(24, 25)
            do m=-2, 3
               ws(1,m) = 1._R8P + vp5(m)**2 ; ws(2,m) = vm5(m) ; ws(3,m) = vp5(m) ; ws(4,m) = 0.1_R8P
               ws(5,m) = 1._R8P + vm5(m)**2
            enddo
            if (v == 24) then
               call seth_keep_pe6(1_I4P, 1.4_R8P, ws, fk)
            else
               call compute_kep_flux(1_I4P, 1.4_R8P, ws, fk)
            endif
            hp = fk(1) ; hm = fk(5)
         case(23) ; call hybrid5(vp5, vm5, hp, hm)
         case(26)
            vp4 = vp5(-1:2) ; vm4 = vm5(-1:2)
            call seth_weno4_c(vp4, vm4, SETH_WENO_JS_EPS, 1._R8P, hp, hm)
         case(29)
            vp4 = vp5(-1:2) ; vm4 = vm5(-1:2)
            call seth_weno3_js(vp4, vm4, SETH_WENO_JS_EPS, 1._R8P, hp, hm)
         case(27) ; call seth_weno6_c(vp5, vm5, SETH_WENO_JS_EPS, 1._R8P, hp, hm)
         case(28) ; call seth_weno8_c(vp7, vm7, SETH_WENO_JS_EPS, 1._R8P, hp, hm)
         case(39)
            do m=-5, 6
               vpw(m) = fp(i+m) ; vmw(m) = fm(i+m)
            enddo
            call seth_teno5_rbf_dr(vpw, vmw, SETH_TENO_EPS, SETH_TENO5_CT, SETH_RBF_CR, 0.09_R8P, 0.09_R8P, &
                                   SETH_DR_T1, SETH_DR_T2, 1._R8P, hp, hm)
         case(37) ; call seth_teno5_st(vp5, vm5, SETH_TENO_EPS, SETH_TENO5_CT, hp, hm)
         case(38) ; call seth_teno5_rbf(vp5, vm5, SETH_TENO_EPS, SETH_TENO5_CT, SETH_RBF_CR, 0.09_R8P, 0.09_R8P, hp, hm)
         case(31:36)
            do m=-5, 6
               vpw(m) = fp(i+m) ; vmw(m) = fm(i+m)
            enddo
            if (v == 31 .or. v == 33 .or. v == 35) then
               call seth_teno5_dr(vpw, vmw, SETH_TENO_EPS, SETH_TENO5_CT, SETH_DR_T1, SETH_DR_T2, 1._R8P, hp, hm)
            else
               call seth_teno5_thinc(vpw, vmw, SETH_TENO_EPS, SETH_TENO5_CT, SETH_THINC_BETA, 1._R8P, hp, hm)
            endif
         case(40)
            call seth_teno5_mask(vp5(-2:2), SETH_TENO_EPS, SETH_TENO5_CT, mka)
            hp = real(mka, R8P) ; hm = hp
         case(41, 42)
            call host_levels(v == 42, i, mka, mkb)
            hp = real(mka, R8P) ; hm = real(mkb, R8P)
         case(43, 44)
            vp4 = vp5(-1:2) ; vm4 = vm5(-1:2)
            hp = 0.5_R8P*(vp4(0) + vp4(1)) ; hm = 0.5_R8P*(vm4(0) + vm4(1))
            if (v == 43) then
               call seth_dr_override(vp4, vm4, 1_I4P, 1_I4P, SETH_DR_T1, SETH_DR_T2, 1._R8P, hp, hm)
            else
               call seth_thinc_override(vp4, vm4, 1_I4P, 1_I4P, SETH_THINC_BETA, 1._R8P, hp, hm)
            endif
         case(45, 46)
            do m=-5, 6
               vpw(m) = fp(i+m) ; vmw(m) = fm(i+m)
            enddo
            call seth_teno5_rbf_dr(vpw, vmw, SETH_TENO_EPS, SETH_TENO5_CT, SETH_RBF_CR, 0.09_R8P, 0.09_R8P, &
                                   SETH_DR_T1, SETH_DR_T2, 1._R8P, hp, hm)
         case(30)
            do m=-2, 2
               wa(m) = vp5(m) ; wb(m) = vm5(1-m)
            enddo
            call seth_teno5_mask(wa, SETH_TENO_EPS, SETH_TENO5_CT, mka)
            call seth_teno5_mask(wb, SETH_TENO_EPS, SETH_TENO5_CT, mkb)
            hp = real(mka, R8P) ; hm = real(mkb, R8P)
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

   pure subroutine hybrid5(vp5, vm5, rp, rm)
   !< The cost hybrid of WENO-Z5 (p = 1): the linear background where seth_sensor_lcr says smooth on both biases.
   real(R8P), intent(in)  :: vp5(-2:3), vm5(-2:3)
   real(R8P), intent(out) :: rp, rm
   real(R8P)              :: sp, sm
   !$acc routine seq
   call seth_sensor_lcr(vp5, 1._R8P, SETH_SENSOR_LCR_PSIC, sp)
   call seth_sensor_lcr(vm5, 1._R8P, SETH_SENSOR_LCR_PSIC, sm)
   if (sp >= 1._R8P .and. sm >= 1._R8P) then
      call seth_linear5(vp5, vm5, rp, rm)
   else
      call seth_wenoz5(vp5, vm5, 1.e-40_R8P, 1_I4P, rp, rm)
   endif
   endsubroutine hybrid5

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
   subroutine host_levels(thinc, i, lp, lm)
   !< The detector levels of face i from masks computed here on the host (the reference of variants 41, 42).
   logical,      intent(in)  :: thinc
   integer(I4P), intent(in)  :: i
   integer(I4P), intent(out) :: lp, lm
   real(R8P)                 :: u(-2:2)
   integer(I4P)              :: mp(-3:3), mm(-3:3), k, m
   do k=-3, 3
      do m=-2, 2
         u(m) = fp(i+k+m)
      enddo
      call seth_teno5_mask(u, SETH_TENO_EPS, SETH_TENO5_CT, mp(k))
      do m=-2, 2
         u(m) = fm(i+1-k-m)
      enddo
      call seth_teno5_mask(u, SETH_TENO_EPS, SETH_TENO5_CT, mm(k))
   enddo
   if (thinc) then
      call seth_thinc_detect(mp, mm, lp, lm)
   else
      call seth_dr_detect(mp, mm, lp, lm)
   endif
   endsubroutine host_levels

   pure subroutine rbf_dr_full(vp, vm, rp, rm)
   !< seth_teno5_rbf_dr (s1 = s2 = 0.09) without the early exit: the RBF states, every mask of cells -3..2 of both
   !< states, the library's DR detector and override (variant 45); equal to seth_teno5_rbf_dr (checked).
   !$acc routine seq
   real(R8P), intent(in)  :: vp(-5:6), vm(-5:6)
   real(R8P), intent(out) :: rp, rm
   real(R8P)              :: a(-2:3), b(-2:3), u(-2:2), c(-1:2), d(-1:2)
   integer(I4P)           :: mp(-3:3), mm(-3:3), k, m, lp, lm
   do m=-2, 3
      a(m) = vp(m) ; b(m) = vm(m)
   enddo
   call seth_teno5_rbf(a, b, SETH_TENO_EPS, SETH_TENO5_CT, SETH_RBF_CR, 0.09_R8P, 0.09_R8P, rp, rm)
   mp = 7 ; mm = 7
   do k=-3, 2
      do m=-2, 2
         u(m) = vp(k+m)
      enddo
      call seth_teno5_mask(u, SETH_TENO_EPS, SETH_TENO5_CT, mp(k))
      do m=-2, 2
         u(m) = vm(1-k-m)
      enddo
      call seth_teno5_mask(u, SETH_TENO_EPS, SETH_TENO5_CT, mm(k))
   enddo
   do m=-1, 2
      c(m) = vp(m) ; d(m) = vm(m)
   enddo
   call seth_dr_detect(mp, mm, lp, lm)
   call seth_dr_override(c, d, lp, lm, SETH_DR_T1, SETH_DR_T2, 1._R8P, rp, rm)
   endsubroutine rbf_dr_full

   pure subroutine composite_full(thinc, vp, vm, rp, rm)
   !< TENO5-DR or TENO5-THINC without the early exit: every mask of cells -3..3 of both states, the library's detector
   !< and override (variants 33, 34); equal to seth_teno5_dr / seth_teno5_thinc (checked).
   !$acc routine seq
   logical,   intent(in)  :: thinc
   real(R8P), intent(in)  :: vp(-5:6), vm(-5:6)
   real(R8P), intent(out) :: rp, rm
   real(R8P)              :: a(-2:3), b(-2:3), u(-2:2), c(-1:2), d(-1:2)
   integer(I4P)           :: mp(-3:3), mm(-3:3), k, m, lp, lm
   do m=-2, 3
      a(m) = vp(m) ; b(m) = vm(m)
   enddo
   call seth_teno5(a, b, SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)
   mp = 7 ; mm = 7 ! DR does not read cell 3
   do k=-3, merge(3, 2, thinc)
      do m=-2, 2
         u(m) = vp(k+m)
      enddo
      call seth_teno5_mask(u, SETH_TENO_EPS, SETH_TENO5_CT, mp(k))
      do m=-2, 2
         u(m) = vm(1-k-m)
      enddo
      call seth_teno5_mask(u, SETH_TENO_EPS, SETH_TENO5_CT, mm(k))
   enddo
   do m=-1, 2
      c(m) = vp(m) ; d(m) = vm(m)
   enddo
   if (thinc) then
      call seth_thinc_detect(mp, mm, lp, lm)
      call seth_thinc_override(c, d, lp, lm, SETH_THINC_BETA, 1._R8P, rp, rm)
   else
      call seth_dr_detect(mp, mm, lp, lm)
      call seth_dr_override(c, d, lp, lm, SETH_DR_T1, SETH_DR_T2, 1._R8P, rp, rm)
   endif
   endsubroutine composite_full
endprogram seth_bench_weno
