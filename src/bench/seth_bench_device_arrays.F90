!< SETH benchmark: what automatic arrays and run-time loop extents cost in a device kernel.
program seth_bench_device_arrays
!< SETH benchmark: what automatic arrays and run-time loop extents cost in a device kernel.
!<
!< Two rules of ADAM's device kernels: (1) no automatic local arrays (bounds from a dummy argument) in a routine called
!< on the device: over-dimension them with a constant; (2) loop extents known at compile time. SETH's WENO-JS, WENO-Z and
!< TENO kernels follow both by construction (one kernel per scheme, named constants); seth_positivity_euler and
!< seth_weno_descaling take their size as a dummy argument (the limiter's first version had automatic arrays). This program measures each rule:
!<
!< - WENO-JS 5: the library kernel against copies of its body that give up rule 1, rule 2 or both
!<   (seth_bench_dyn_kernels), S = 3 passed at run time;
!< - the descaler: seth_weno_descaling(S, ..) against a copy with S constant, each followed by seth_weno5_js;
!< - the positivity limiter, nv = 3 and 5: the library kernel (locals over-dimensioned to 5, loops to nv) against its
!<   first version (automatic arrays, array expressions over nv: seth_bench_positivity_auto), a bench copy of the same
!<   design and one with nv constant (both rules).
!<
!< Built with -acc the loops run on the device (timed with !$acc wait + system_clock); without, on one CPU core. Each
!< variant: 2 warm-up passes, then `reps` timed passes; median and minimum. S and nv are read from variables the compiler
!< cannot fold (volatile). Checks: every WENO variant equals the library kernel; every limiter variant gives the
!< library's theta and flux to 1e-8.
!<
!< Usage: seth_bench_device_arrays [n_faces [reps [auto]]]   (defaults 16777216 and 20; the limiter runs on n_faces/4)
!<
!< The variants with automatic arrays run only with the third argument `auto`: on the device they need a larger heap
!< than the default (NVCOMPILER_ACC_CUDA_HEAPSIZE=268435456 for the default size; without it the launch fails with
!< "DEVICE FORTRAN AUTO ALLOCATION FAILED") and take minutes.
use, intrinsic :: iso_fortran_env, only : int64
use seth
use seth_bench_dyn_kernels
use seth_bench_positivity_auto, only : positivity_euler_auto
use seth_bench_positivity_max, only : pos_max => positivity_euler
use seth_bench_positivity_nv3, only : pos_nv3 => positivity_euler
use seth_bench_positivity_nv5, only : pos_nv5 => positivity_euler
implicit none
integer(I4P), parameter :: NG = 4_I4P, NVAR = 14_I4P
real(R8P),    parameter :: G = 1.4_R8P, TOL = 64._R8P*epsilon(1._R8P)
character(46), parameter :: NAMES(NVAR) = [character(46) :: &
   'weno5_js  library (const arrays, const loops)', 'weno5_js  S_MAX arrays, run-time loops', &
   'weno5_js  automatic arrays, const loops', 'weno5_js  automatic arrays, run-time loops', &
   'weno5_js + descaling(S dummy)', 'weno5_js + descaling(S const)', &
   'positivity nv=3 library (NV_MAX, loops to nv)', 'positivity nv=3 automatic arrays (ea4959b)', &
   'positivity nv=3 bench NV_MAX, loops to nv', 'positivity nv=3 const arrays, const loops', &
   'positivity nv=5 library (NV_MAX, loops to nv)', 'positivity nv=5 automatic arrays (ea4959b)', &
   'positivity nv=5 bench NV_MAX, loops to nv', 'positivity nv=5 const arrays, const loops']
integer(I4P), parameter :: REF(NVAR) = [1, 1, 1, 1, 5, 5, 7, 7, 7, 7, 11, 11, 11, 11]
real(R8P), allocatable :: fp(:), fm(:), rp(:), rm(:), rp0(:), rm0(:)
real(R8P), allocatable :: ql(:,:), fh(:,:), fl(:,:), lam(:), fo(:,:), th(:), fo0(:,:), th0(:)
real(R8P), allocatable :: times(:)
real(R8P)              :: med(NVAR), tmin(NVAR), diff
integer(I4P)           :: n, np, reps, v, r, c, nlim
integer(I4P), volatile :: s_run, nv_run
integer(int64)         :: t0, t1, rate
character(32)          :: arg
logical                :: passed, with_auto

n = 16777216_I4P ; reps = 20_I4P
if (command_argument_count() >= 1) then ; call get_command_argument(1, arg) ; read(arg,*) n ; endif
if (command_argument_count() >= 2) then ; call get_command_argument(2, arg) ; read(arg,*) reps ; endif
with_auto = .false.
if (command_argument_count() >= 3) then ; call get_command_argument(3, arg) ; with_auto = trim(arg) == 'auto' ; endif
np = n/4_I4P
med = -1._R8P ; tmin = -1._R8P
s_run = 3_I4P ; nv_run = 3_I4P
allocate(fp(1-NG:n+NG), fm(1-NG:n+NG), rp(n), rm(n), rp0(n), rm0(n), times(reps))
allocate(ql(5,0:np), fh(5,np), fl(5,np), lam(np), fo(5,np), th(np), fo0(5,np), th0(np))
do c=1-NG, n+NG
   fp(c) = 1._R8P + 0.5_R8P*sin(0.01_R8P*c) + merge(1._R8P, 0._R8P, mod(c+1000000, 4096) < 2048)
   fm(c) = 0.5_R8P*cos(0.013_R8P*c) - merge(0.7_R8P, 0._R8P, mod(c+1000000, 4096) < 2048)
enddo
ql = 0._R8P ; fh = 0._R8P ; fl = 0._R8P ; fo = 0._R8P
!$acc enter data copyin(fp, fm) create(rp, rm, ql, fh, fl, lam, fo, th)

call system_clock(count_rate=rate)
passed = .true.
print '(A,I0,A,I0,A,I0,A)', 'faces: ', n, ' (limiter: ', np, '), timed passes: ', reps, ' per variant'
do v=1, NVAR
   if (v == 7 .or. v == 11) then
      nv_run = merge(3_I4P, 5_I4P, v == 7)
      call limiter_data(nv_run)
      !$acc update device(ql, fh, fl, lam)
   endif
   if (.not.with_auto .and. (v == 3 .or. v == 4 .or. v == 8 .or. v == 12)) cycle
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
   if (v <= 6) then
      !$acc update self(rp, rm)
      if (v == REF(v)) then
         rp0 = rp ; rm0 = rm
      else
         diff = max(maxval(abs(rp - rp0)), maxval(abs(rm - rm0)))
         if (diff > 4._R8P*TOL) then
            passed = .false. ; print '(A,ES9.2)', '   FAILED: '//trim(NAMES(v))//' differs from its reference by ', diff
         endif
      endif
   else
      !$acc update self(fo, th)
      if (v == REF(v)) then
         fo0 = fo ; th0 = th
         nlim = count(th < 1._R8P)
         print '(A,I0,A,I0,A,F5.1,A)', '   nv=', nv_run, ': limited faces ', nlim, ' (', 100._R8P*nlim/np, ' %)'
      else
         diff = max(maxval(abs(th - th0)), maxval(abs(fo - fo0))/maxval(abs(fo0)))
         if (diff > 1.e-8_R8P) then ! measured: 9e-10 on one face of 2^20 (nvfortran, CPU), below 1e-10 elsewhere
            passed = .false. ; print '(A,ES9.2)', '   FAILED: '//trim(NAMES(v))//' differs from the library by ', diff
         endif
      endif
   endif
enddo
!$acc exit data delete(fp, fm, rp, rm, ql, fh, fl, lam, fo, th)
print '(A)', ''
print '(A46,4A12)', 'variant', 'median ms', 'min ms', 'ns/face', 'vs first'
do v=1, NVAR
   if (med(v) < 0._R8P) cycle
   print '(A46,3F12.3,F12.2)', NAMES(v), 1.e3_R8P*med(v), 1.e3_R8P*tmin(v), &
                               1.e9_R8P*med(v)/merge(n, np, v <= 6), med(v)/med(REF(v))
enddo
if (.not.passed) error stop 'a benchmark check failed'
print '(A)', 'seth_bench_device_arrays: checks ok'

contains
   subroutine limiter_data(nv)
   !< Cell states (rho, p from 1e-6 to 1, velocities to 3), Lax-Friedrichs fluxes at CFL 1/2 and a high-order flux
   !< perturbed so that part of the faces is limited (host).
   integer(I4P), intent(in) :: nv
   real(R8P)                :: rho, p, a, x, fa(5), fb(5)
   integer(I4P)             :: i, k
   do i=0, np
      x = real(i, R8P)
      rho = 10._R8P**(-3._R8P + 3._R8P*sin(0.37_R8P*x)) ; p = 10._R8P**(-3._R8P + 3._R8P*cos(0.23_R8P*x))
      ql(1,i) = rho
      do k=2, nv-1
         ql(k,i) = rho*3._R8P*sin(0.11_R8P*x + k)
      enddo
      ql(nv,i) = p/(G - 1._R8P) + 0.5_R8P*sum(ql(2:nv-1,i)**2)/rho
   enddo
   do i=1, np
      a = max(speed(nv, ql(:,i-1)), speed(nv, ql(:,i)))
      lam(i) = 1._R8P/a
      fa = eflux(nv, ql(:,i-1)) ; fb = eflux(nv, ql(:,i))
      fl(1:nv,i) = 0.5_R8P*(fa(1:nv) + fb(1:nv)) - 0.5_R8P*a*(ql(1:nv,i) - ql(1:nv,i-1))
      fh(1:nv,i) = fl(1:nv,i) + 0.3_R8P*sin(1.7_R8P*i)*(abs(fa(1:nv)) + abs(fb(1:nv)))
   enddo
   endsubroutine limiter_data

   pure function speed(nv, q) result(s)
   integer(I4P), intent(in) :: nv
   real(R8P),    intent(in) :: q(5)
   real(R8P)                :: s, p
   p = (G - 1._R8P)*(q(nv) - 0.5_R8P*sum(q(2:nv-1)**2)/q(1))
   s = abs(q(2)/q(1)) + sqrt(G*p/q(1))
   endfunction speed

   pure function eflux(nv, q) result(fx)
   integer(I4P), intent(in) :: nv
   real(R8P),    intent(in) :: q(5)
   real(R8P)                :: fx(5), u, p
   u = q(2)/q(1)
   p = (G - 1._R8P)*(q(nv) - 0.5_R8P*sum(q(2:nv-1)**2)/q(1))
   fx = u*q
   fx(2) = fx(2) + p
   fx(nv) = fx(nv) + p*u
   endfunction eflux

   subroutine run(v)
   !< One pass of variant v.
   integer(I4P), intent(in) :: v
   real(R8P)                :: vp(-2:3), vm(-2:3), l(5), rr(5), h(5), lo(5), f(5), rmu
   integer(I4P)             :: i, m, k, s_, nv_
   s_ = s_run ; nv_ = nv_run
   select case(v)
   case(1)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp, vm)
      do i=1, n
         do m=-2, 3
            vp(m) = fp(i+m) ; vm(m) = fm(i+m)
         enddo
         call seth_weno5_js(vp, vm, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
   case(2)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp, vm)
      do i=1, n
         do m=-2, 3
            vp(m) = fp(i+m) ; vm(m) = fm(i+m)
         enddo
         call weno5_js_smax(s_, vp, vm, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
   case(3)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp, vm)
      do i=1, n
         do m=-2, 3
            vp(m) = fp(i+m) ; vm(m) = fm(i+m)
         enddo
         call weno5_js_auto_cloop(s_, vp, vm, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
   case(4)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp, vm)
      do i=1, n
         do m=-2, 3
            vp(m) = fp(i+m) ; vm(m) = fm(i+m)
         enddo
         call weno5_js_auto(s_, vp, vm, SETH_WENO_JS_EPS, 1._R8P, rp(i), rm(i))
      enddo
   case(5)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp, vm, rmu)
      do i=1, n
         do m=-2, 3
            vp(m) = fp(i+m) ; vm(m) = fm(i+m)
         enddo
         rmu = seth_weno_descaling(s_, 1._R8P, vp, vm)
         call seth_weno5_js(vp, vm, SETH_WENO_JS_EPS, rmu, rp(i), rm(i))
      enddo
   case(6)
      !$acc parallel loop gang vector present(fp, fm, rp, rm) private(vp, vm, rmu)
      do i=1, n
         do m=-2, 3
            vp(m) = fp(i+m) ; vm(m) = fm(i+m)
         enddo
         rmu = descaling_s3(1._R8P, vp, vm)
         call seth_weno5_js(vp, vm, SETH_WENO_JS_EPS, rmu, rp(i), rm(i))
      enddo
   case(7, 11)
      !$acc parallel loop gang vector present(ql, fh, fl, lam, fo, th) private(l, rr, h, lo, f)
      do i=1, np
         do k=1, 5
            l(k) = ql(k,i-1) ; rr(k) = ql(k,i) ; h(k) = fh(k,i) ; lo(k) = fl(k,i) ; f(k) = 0._R8P
         enddo
         call seth_positivity_euler(nv_, G, lam(i), SETH_POSITIVITY_EPS, SETH_POSITIVITY_EPS, l, rr, h, lo, f, th(i))
         do k=1, 5
            fo(k,i) = f(k)
         enddo
      enddo
   case(8, 12)
      !$acc parallel loop gang vector present(ql, fh, fl, lam, fo, th) private(l, rr, h, lo, f)
      do i=1, np
         do k=1, 5
            l(k) = ql(k,i-1) ; rr(k) = ql(k,i) ; h(k) = fh(k,i) ; lo(k) = fl(k,i) ; f(k) = 0._R8P
         enddo
         call positivity_euler_auto(nv_, G, lam(i), SETH_POSITIVITY_EPS, SETH_POSITIVITY_EPS, l, rr, h, lo, f, th(i))
         do k=1, 5
            fo(k,i) = f(k)
         enddo
      enddo
   case(9, 13)
      !$acc parallel loop gang vector present(ql, fh, fl, lam, fo, th) private(l, rr, h, lo, f)
      do i=1, np
         do k=1, 5
            l(k) = ql(k,i-1) ; rr(k) = ql(k,i) ; h(k) = fh(k,i) ; lo(k) = fl(k,i) ; f(k) = 0._R8P
         enddo
         call pos_max(nv_, G, lam(i), SETH_POSITIVITY_EPS, SETH_POSITIVITY_EPS, l, rr, h, lo, f, th(i))
         do k=1, 5
            fo(k,i) = f(k)
         enddo
      enddo
   case(10)
      !$acc parallel loop gang vector present(ql, fh, fl, lam, fo, th) private(l, rr, h, lo, f)
      do i=1, np
         do k=1, 5
            l(k) = ql(k,i-1) ; rr(k) = ql(k,i) ; h(k) = fh(k,i) ; lo(k) = fl(k,i) ; f(k) = 0._R8P
         enddo
         call pos_nv3(nv_, G, lam(i), SETH_POSITIVITY_EPS, SETH_POSITIVITY_EPS, l, rr, h, lo, f, th(i))
         do k=1, 5
            fo(k,i) = f(k)
         enddo
      enddo
   case(14)
      !$acc parallel loop gang vector present(ql, fh, fl, lam, fo, th) private(l, rr, h, lo, f)
      do i=1, np
         do k=1, 5
            l(k) = ql(k,i-1) ; rr(k) = ql(k,i) ; h(k) = fh(k,i) ; lo(k) = fl(k,i) ; f(k) = 0._R8P
         enddo
         call pos_nv5(nv_, G, lam(i), SETH_POSITIVITY_EPS, SETH_POSITIVITY_EPS, l, rr, h, lo, f, th(i))
         do k=1, 5
            fo(k,i) = f(k)
         enddo
      enddo
   endselect
   endsubroutine run

   subroutine sort(x)
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
endprogram seth_bench_device_arrays
