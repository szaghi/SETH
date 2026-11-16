!< SETH benchmark: SETH's kernels in the working kind R_P, the same source built in R8P and in R4P.
program seth_bench_kinds
!< SETH benchmark: SETH's kernels in the working kind R_P, the same source built in R8P and in R4P.
!<
!< Roadmap issue #4, 6.2 (phase 5): what single precision buys, kernel by kernel. The SETH variants of
!< seth_bench_weno (whose ADAM comparisons are R8P by construction), on the same field and with the same timing, in
!< R_P: build it in both kinds (bench-kinds-* and bench-*-r4 modes) and compare the ns per face.
!<
!< n faces of a 1-D field (cells 1-NG..n+NG; face i between cells i and i+1), two biased fields fp, fm (a smooth wave
!< with a jump every 2048 cells, values within [-1.2, 2.5]). Every variant reads the same stencils and writes the two
!< states of every face:
!<
!<  1 seth_weno5_js, 2 seth_weno7_js, 3 seth_wenoz5 (p = 1), 4 seth_teno6, 5 seth_teno5, 6 seth_teno5_opt,
!<  7 seth_teno5_az (rmu = 1), 8 seth_teno8, 9 seth_teno5_a, 10 seth_teno8_a, 11 seth_weno5_js with
!<  seth_weno_descaling, 12 with seth_weno_descaling_ai, 13 seth_weno5_js + seth_mp_bound, 14 seth_linear5,
!<  15 seth_linear7, 16 seth_sensor_lcr alone (s stored in rp, rm), 17 the cost hybrid of WENO-Z5 (seth_sensor_lcr,
!<  then seth_linear5 where s = 1 on both biases and seth_wenoz5 elsewhere), 18 seth_keep_pe6 on the state
!<  rho = 1 + fp^2, u = fm, v = fp, w = 0.1, p = 1 + fm^2 (the flux of direction 1 in rp, its energy component in rm)
!<
!< Built with -acc the loops run on the device (timed with !$acc wait + system_clock); without, on one CPU core. Each
!< variant: 2 warm-up passes, then `reps` timed passes; the median and the minimum are reported, in the rows that
!< scripts/bench_compare.py reads. Check: the device results equal the host results on a sample of faces, to 64 ulps
!< of the data, but at faces where a host/device rounding difference flips a branch (a TENO cut-off, the hybrid's
!< switch): those are counted and printed, not failed, when they stay below 1 per 1000 sampled faces.
!<
!< Usage: seth_bench_kinds [n_faces [reps]]   (defaults 16777216 and 20)
use, intrinsic :: iso_fortran_env, only : int64
use seth
implicit none
integer(I4P), parameter :: NG = 4_I4P         !< Ghost cells (the 8-point stencils read i-3..i+4).
integer(I4P), parameter :: NV = 18_I4P        !< Variants.
integer(I4P), parameter :: SAMPLE = 65536_I4P !< Faces checked against the host.
character(31), parameter :: NAMES(NV) = [character(31) :: 'SETH  weno5_js', 'SETH  weno7_js', 'SETH  wenoz5 (p=1)', &
                                         'SETH  teno6', 'SETH  teno5', 'SETH  teno5_opt', 'SETH  teno5_az', &
                                         'SETH  teno8', 'SETH  teno5_a', 'SETH  teno8_a', &
                                         'SETH  weno5_js + descaling', 'SETH  weno5_js + descaling_ai', &
                                         'SETH  weno5_js + mp_bound', 'SETH  linear5', 'SETH  linear7', &
                                         'SETH  sensor_lcr', 'SETH  hybrid lcr+linear5/wenoz5', 'SETH  keep_pe6']
real(R_P), allocatable :: fp(:), fm(:), rp(:), rm(:)
real(R8P), allocatable :: times(:)
real(R8P)              :: med(NV), tmin(NV)
real(R_P)              :: tol
integer(I4P)           :: n, reps, v, r, c, nflip
integer(int64)         :: t0, t1, rate
character(32)          :: arg
logical                :: passed

n = 16777216_I4P ; reps = 20_I4P
if (command_argument_count() >= 1) then ; call get_command_argument(1, arg) ; read(arg,*) n ; endif
if (command_argument_count() >= 2) then ; call get_command_argument(2, arg) ; read(arg,*) reps ; endif
allocate(fp(1-NG:n+NG), fm(1-NG:n+NG), rp(n), rm(n), times(reps))
do c=1-NG, n+NG
   fp(c) = real(1._R8P + 0.5_R8P*sin(0.01_R8P*c) + merge(1._R8P, 0._R8P, mod(c+1000000, 4096) < 2048), R_P)
   fm(c) = real(0.5_R8P*cos(0.013_R8P*c) - merge(0.7_R8P, 0._R8P, mod(c+1000000, 4096) < 2048), R_P)
enddo
tol = 64._R_P*epsilon(1._R_P)
!$acc enter data copyin(fp, fm) create(rp, rm)

call system_clock(count_rate=rate)
passed = .true.
print '(A,I0,A,I0,A,I0,A)', 'faces: ', n, ', timed passes: ', reps, ' per variant, R_P of ', storage_size(1._R_P), &
                             ' bits'
med = -1._R8P ; tmin = -1._R8P
do v=1, NV
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
!$acc exit data delete(fp, fm, rp, rm)
print '(A)', ''
print '(A31,3A13)', 'variant', 'median ms', 'min ms', 'ns/face'
do v=1, NV
   print '(A31,3F13.3)', NAMES(v), 1.e3_R8P*med(v), 1.e3_R8P*tmin(v), 1.e9_R8P*med(v)/n
enddo
if (.not.passed) error stop 'a benchmark check failed'
print '(A)', 'seth_bench_kinds: checks ok'

contains
   pure subroutine kernel(v, vp7, vm7, hp, hm)
   !< Variant v on the 8-point stencils vp7, vm7 of one face (the 6-point kernels read their middle -2..3).
   !$acc routine seq
   integer(I4P), intent(in)  :: v
   real(R_P),    intent(in)  :: vp7(-3:4), vm7(-3:4)
   real(R_P),    intent(out) :: hp, hm
   real(R_P)                 :: vp5(-2:3), vm5(-2:3), sp, sm, ws(5,-2:3), fk(5)
   integer(I4P)              :: m
   vp5 = vp7(-2:3) ; vm5 = vm7(-2:3)
   select case(v)
   case(1)  ; call seth_weno5_js(vp5, vm5, SETH_WENO_JS_EPS, 1._R_P, hp, hm)
   case(2)  ; call seth_weno7_js(vp7, vm7, SETH_WENO_JS_EPS, 1._R_P, hp, hm)
   case(3)  ; call seth_wenoz5(vp5, vm5, SETH_TENO_EPS, 1_I4P, hp, hm)
   case(4)  ; call seth_teno6(vp5, vm5, SETH_TENO_EPS, SETH_TENO6_CT, hp, hm)
   case(5)  ; call seth_teno5(vp5, vm5, SETH_TENO_EPS, SETH_TENO5_CT, hp, hm)
   case(6)  ; call seth_teno5_opt(vp5, vm5, SETH_TENO_EPS, SETH_TENO5_OPT_CT, hp, hm)
   case(7)  ; call seth_teno5_az(vp5, vm5, SETH_TENO_EPS, SETH_TENO5_AZ_CT, 1._R_P, hp, hm)
   case(8)  ; call seth_teno8(vp7, vm7, SETH_TENO_EPS, SETH_TENO8_CT, hp, hm)
   case(9)  ; call seth_teno5_a(vp5, vm5, SETH_TENO_EPS, SETH_TENO5_A_CR, SETH_TENO5_A_XI, SETH_TENO5_A_A1, &
                                SETH_TENO5_A_A2, 1._R_P, hp, hm)
   case(10) ; call seth_teno8_a(vp7, vm7, SETH_TENO_EPS, SETH_TENO8_A_CR, SETH_TENO8_A_XI, SETH_TENO8_A_A1, &
                                SETH_TENO8_A_A2, 1._R_P, hp, hm)
   case(11) ; call seth_weno5_js(vp5, vm5, SETH_WENO_JS_EPS, seth_weno_descaling(3_I4P, 1._R_P, vp5, vm5), hp, hm)
   case(12) ; call seth_weno5_js(vp5, vm5, SETH_WENO_JS_EPS, seth_weno_descaling_ai(3_I4P, vp5, vm5), hp, hm)
   case(13) ; call seth_weno5_js(vp5, vm5, SETH_WENO_JS_EPS, 1._R_P, hp, hm)
              call seth_mp_bound(vp5, vm5, SETH_MP_ALPHA, hp, hm)
   case(14) ; call seth_linear5(vp5, vm5, hp, hm)
   case(15) ; call seth_linear7(vp7, vm7, hp, hm)
   case(16) ; call seth_sensor_lcr(vp5, 1._R_P, SETH_SENSOR_LCR_PSIC, hp) ; hm = hp
   case(17)
      call seth_sensor_lcr(vp5, 1._R_P, SETH_SENSOR_LCR_PSIC, sp)
      call seth_sensor_lcr(vm5, 1._R_P, SETH_SENSOR_LCR_PSIC, sm)
      if (sp >= 1._R_P .and. sm >= 1._R_P) then
         call seth_linear5(vp5, vm5, hp, hm)
      else
         call seth_wenoz5(vp5, vm5, SETH_TENO_EPS, 1_I4P, hp, hm)
      endif
   case default
      do m=-2, 3
         ws(1,m) = 1._R_P + vp5(m)**2 ; ws(2,m) = vm5(m) ; ws(3,m) = vp5(m) ; ws(4,m) = 0.1_R_P
         ws(5,m) = 1._R_P + vm5(m)**2
      enddo
      call seth_keep_pe6(1_I4P, 1.4_R_P, ws, fk)
      hp = fk(1) ; hm = fk(5)
   endselect
   endsubroutine kernel

   subroutine run(v)
   !< One pass of variant v over all faces. One loop per stencil width, so that each kernel reads only its stencil.
   integer(I4P), intent(in) :: v
   real(R_P)                :: vp7(-3:4), vm7(-3:4)
   integer(I4P)             :: i, m
   if (v == 2 .or. v == 8 .or. v == 10 .or. v == 15) then
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp7, vm7)
      do i=1, n
         do m=-3, 4
            vp7(m) = fp(i+m) ; vm7(m) = fm(i+m)
         enddo
         call kernel(v, vp7, vm7, rp(i), rm(i))
      enddo
   else
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp7, vm7)
      do i=1, n
         do m=-2, 3
            vp7(m) = fp(i+m) ; vm7(m) = fm(i+m)
         enddo
         vp7(-3) = 0._R_P ; vm7(-3) = 0._R_P ; vp7(4) = 0._R_P ; vm7(4) = 0._R_P
         call kernel(v, vp7, vm7, rp(i), rm(i))
      enddo
   endif
   endsubroutine run

   subroutine check(v)
   !< The results of variant v against the host on a sample of faces, scaled by the stencil's magnitude.
   integer(I4P), intent(in) :: v
   real(R_P)                :: vp7(-3:4), vm7(-3:4), hp, hm, vmax, d, diff
   integer(I4P)             :: i, m, nsamp
   diff = 0._R_P ; nflip = 0 ; nsamp = 0
   do i=1, n, max(1_I4P, n/SAMPLE)
      do m=-3, 4
         vp7(m) = fp(i+m) ; vm7(m) = fm(i+m)
      enddo
      call kernel(v, vp7, vm7, hp, hm)
      vmax = max(maxval(abs(vp7)), maxval(abs(vm7)))
      d = max(abs(rp(i) - hp), abs(rm(i) - hm))/vmax
      nsamp = nsamp + 1
      if (d > tol) then
         nflip = nflip + 1
      else
         diff = max(diff, d)
      endif
   enddo
   print '(A31,A,ES9.2,A,I0,A,I0)', NAMES(v), ': device vs host, max scaled difference ', diff, &
                                    ', faces beyond 64 ulps ', nflip, ' of ', nsamp
   if (nflip*1000 > nsamp) passed = .false.
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
endprogram seth_bench_kinds
