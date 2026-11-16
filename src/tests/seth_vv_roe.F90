!< SETH verification: the Roe flux splitting with entropy fix of the 1-D Euler harness (issue #5), on its own.
program seth_vv_roe
!< SETH verification: the Roe flux splitting with entropy fix of the 1-D Euler harness (issue #5), on its own.
!<
!< PATH_ROE (seth_vv_euler; notes/issue5-roe-sources.md) is the "Roe scheme with entropy-fix" of Fu, Hu & Adams (2016,
!< section 5.3) read as Jiang & Shu's (1996) RF: Shu & Osher's (1989) Algorithm 2.3 in characteristic fields. Checked:
!<
!< 1. Order: the density wave of seth_vv_euler_order (rho = 1 + 0.2 sin(pi x), u = p = 1, periodic, t = 2). No field
!<    changes sign there (u - c < 0 < u < u + c), so the fix must never act (nsonic = 0) and the path is the unsplit
!<    upwind reconstruction; asserted as in seth_vv_euler_order (TENO and WENO-Z at their order - 0.5, every scheme at
!<    least 2), and no face limited.
!< 2. Sod and Lax against the exact solution, N = 200, CFL 0.5, without the limiter: positive, and L1 of the density at
!<    most TUBE_RATIO times that of the splitting path (PATH_SPLIT).
!< 3. The transonic rarefaction of Toro's test 1, (1, 0.75, 1) | (0.125, 0, 0.1) at x = 0.3, t = 0.2: the left
!<    rarefaction contains the sonic point x/t = 0, where an upwind Roe flux without fix admits a stationary expansion
!<    shock (Shu & Osher, after Algorithm 2.1). Linf of the density error in the fan region x in [0.2, 0.35], with
!<    the fix and without it (roe_fix = .false.): asserted, the error with the fix below FAN_BOUND and the error
!<    without it at least FAN_FACTOR times larger (the negative control: the fix must matter).
!< 4. The 123 problem (two transonic rarefactions, near vacuum), with the fix: positive with the positivity limiter,
!<    L1 printed with and without it.
!< Quick schemes (WENO5-JS, WENO-Z5, TENO5, TENO6) in the test run; every reconstruction in the report mode, which
!< writes roe-summary.dat (case scheme L1-roe L1-split positive-without-limiter).
use seth
use seth_test_schemes
use seth_vv_euler
use seth_vv_cases
use seth_vv_riemann, only : riemann_sample
implicit none
real(R8P), parameter      :: PI = acos(-1._R8P)
! Recorded bounds (gfortran -O3, quick schemes, 2026-11-16): Sod/Lax L1 ratio to PATH_SPLIT at most 1.18 (TENO6 on Sod),
! the others 0.76-0.93; Linf in the transonic fan 1.32e-2 .. 1.66e-2 with the fix, 5.3 .. 10.6 times that without it.
real(R8P), parameter      :: TUBE_RATIO = 1.25_R8P, FAN_BOUND = 1.75e-2_R8P, FAN_FACTOR = 4._R8P
integer(I4P), allocatable :: ids(:), nfb(:)
real(R8P), allocatable    :: l1r(:), l1s(:), e1(:), ei(:)
logical, allocatable      :: posr(:), poss(:)
integer(I4P)              :: k, j, nl, n, us
real(R8P)                 :: o, need, efix, enofix
logical                   :: passed
passed = .true.
if (vv_report()) call vv_open('roe-summary.dat', us, 'case scheme L1-roe L1-split positive-without-limiter')
call vv_schemes(PATH_ROE, .not.vv_report(), ids)
allocate(l1r(size(ids)), l1s(size(ids)), posr(size(ids)), poss(size(ids)), nfb(size(ids)))

! 1. order on the density wave
print '(A)', '1. density wave, PATH_ROE: L1, Linf (finest), observed order'
nl = merge(5, 4, vv_report())
allocate(e1(nl), ei(nl))
limiter = .true.
do k=1, size(ids)
   nsonic = 0
   do j=1, nl
      n = 20*2**(j-1)
      nlimited = 0
      call wave(ids(k), n, e1(j), ei(j))
      if (nlimited > 0) call fail(vv_label(ids(k))//': the positivity limiter acted on a smooth flow')
   enddo
   o = order(e1)
   need = SCHEMES(ids(k))%order
   print '(A18,2ES12.4,F8.2,A,F4.1)', vv_label(ids(k)), e1(nl), ei(nl), o, '  expected ', need
   if (nsonic /= 0) call fail(vv_label(ids(k))//': the entropy fix acted on a flow without sonic points')
   if (o < 2._R8P) call fail(vv_label(ids(k))//': below second order')
   if (SCHEMES(ids(k))%family == FAM_TENO .or. (SCHEMES(ids(k))%family == FAM_Z .and. SCHEMES(ids(k))%s >= 3)) then
      if (o < need - 0.5_R8P) call fail(vv_label(ids(k))//': below its expected order')
   endif
enddo
limiter = .false.

! 2. Sod and Lax
call tube('sod', [1._R8P, 0._R8P, 1._R8P], [0.125_R8P, 0._R8P, 0.1_R8P], 0.5_R8P, 0.2_R8P)
call tube('lax', [0.445_R8P, 0.698_R8P, 3.528_R8P], [0.5_R8P, 0._R8P, 0.571_R8P], 0.5_R8P, 0.13_R8P)

! 3. transonic rarefaction
print '(A)', '3. transonic rarefaction (Toro test 1), N = 200: Linf(rho) in x = 0.2 .. 0.35, with and without the fix'
do k=1, size(ids)
   roe_fix = .true.  ; call toro1(ids(k), efix)
   roe_fix = .false. ; call toro1(ids(k), enofix)
   roe_fix = .true.
   print '(A18,2ES12.4,F8.2)', vv_label(ids(k)), efix, enofix, enofix/efix
   if (vv_report()) write(us, '(A,1X,A,2ES24.15,1X,L1)') 'toro1', vv_label(ids(k)), efix, enofix, .true.
   if (efix > FAN_BOUND) call fail(vv_label(ids(k))//': transonic rarefaction, error above the bound with the fix')
   if (enofix < FAN_FACTOR*efix) call fail(vv_label(ids(k))//': transonic rarefaction, the fix does not matter')
enddo

! 4. 123
print '(A)', '4. 123 problem, N = 200: L1(rho) PATH_ROE with the limiter, without it, positive without it'
limiter = .false.
call riemann_problem('roe-123', [1._R8P, -2._R8P, 0.4_R8P], [1._R8P, 2._R8P, 0.4_R8P], 0.5_R8P, 0.15_R8P, 0.5_R8P, &
                     200, ids, PATH_ROE, l1s, poss, nfb, files=.false.)
limiter = .true.
call riemann_problem('roe-123', [1._R8P, -2._R8P, 0.4_R8P], [1._R8P, 2._R8P, 0.4_R8P], 0.5_R8P, 0.15_R8P, 0.5_R8P, &
                     200, ids, PATH_ROE, l1r, posr, nfb, files=.false.)
limiter = .false.
do k=1, size(ids)
   print '(A18,2ES12.4,L3)', vv_label(ids(k)), l1r(k), l1s(k), poss(k)
   if (vv_report()) write(us, '(A,1X,A,2ES24.15,1X,L1)') '123', vv_label(ids(k)), l1r(k), l1s(k), poss(k)
   if (.not.posr(k)) call fail(vv_label(ids(k))//': 123 not positive with the limiter')
enddo

if (vv_report()) close(us)
if (.not.passed) error stop 'Roe path verification failed'
print '(A)', 'seth_vv_roe: ok'

contains
   subroutine tube(name, wl, wr, x0, tend)
   !< A shock tube on both PATH_ROE and PATH_SPLIT, without the limiter.
   character(*), intent(in) :: name
   real(R8P),    intent(in) :: wl(3), wr(3), x0, tend
   integer(I4P)             :: kk
   call riemann_problem(name, wl, wr, x0, tend, 0.5_R8P, 200, ids, PATH_ROE, l1r, posr, nfb, files=.false.)
   call riemann_problem(name, wl, wr, x0, tend, 0.5_R8P, 200, ids, PATH_SPLIT, l1s, poss, nfb, files=.false.)
   print '(A)', '2. '//name//', N = 200, no limiter: L1(rho) PATH_ROE, PATH_SPLIT, ratio'
   do kk=1, size(ids)
      print '(A18,2ES12.4,F8.3)', vv_label(ids(kk)), l1r(kk), l1s(kk), l1r(kk)/l1s(kk)
      if (vv_report()) write(us, '(A,1X,A,2ES24.15,1X,L1)') name, vv_label(ids(kk)), l1r(kk), l1s(kk), posr(kk)
      if (.not.posr(kk)) call fail(vv_label(ids(kk))//': '//name//' not positive')
      if (l1r(kk) > TUBE_RATIO*l1s(kk)) call fail(vv_label(ids(kk))//': '//name//' L1 above the bound')
   enddo
   endsubroutine tube

   subroutine toro1(id, einf)
   !< Toro's test 1 on [0, 1], N = 200, CFL 0.5, PATH_ROE: Linf of the density error in the fan region.
   integer(I4P), intent(in)  :: id
   real(R8P),    intent(out) :: einf
   integer(I4P), parameter   :: NN = 200
   real(R8P), parameter      :: WL(3) = [1._R8P, 0.75_R8P, 1._R8P], WR(3) = [0.125_R8P, 0._R8P, 0.1_R8P]
   real(R8P)                 :: w0(3,NN), q(3,1-NG:NN+NG), h, x, ex(3)
   logical                   :: pos
   integer(I4P)              :: i
   h = 1._R8P/NN
   do i=1, NN
      w0(:,i) = merge(WL, WR, (i - 0.5_R8P)*h < 0.3_R8P)
   enddo
   call euler_run(id, PATH_ROE, w0, NN, h, BC_TRANSMISSIVE, 0.2_R8P, 0.5_R8P, q, pos)
   einf = huge(1._R8P)
   if (.not.pos) return
   einf = 0._R8P
   do i=1, NN
      x = (i - 0.5_R8P)*h
      if (x < 0.2_R8P .or. x > 0.35_R8P) cycle
      ex = riemann_sample(WL, WR, GAMMA, (x - 0.3_R8P)/0.2_R8P)
      einf = max(einf, abs(q(1,i) - ex(1)))
   enddo
   endsubroutine toro1

   subroutine wave(id, n, l1, linf)
   integer(I4P), intent(in)  :: id, n
   real(R8P),    intent(out) :: l1, linf
   real(R8P)                 :: q(3,1-NG:n+NG), h, t, dt, sc, e(n), x
   logical                   :: ok
   integer(I4P)              :: i
   h = 2._R8P/n
   do i=1, n
      x = -1._R8P + (i - 0.5_R8P)*h
      q(:,i) = prim2cons([1._R8P + 0.2_R8P*sin(PI*x), 1._R8P, 1._R8P])
   enddo
   sc = (h/0.1_R8P)**((SCHEMES(id)%order - 3._R8P)/3._R8P)
   t = 0._R8P
   do while (t < 2._R8P*(1._R8P - 1.e-14_R8P))
      dt = min(0.4_R8P*sc*euler_dt(q, n, h, 1._R8P), 2._R8P - t)
      call euler_step(id, PATH_ROE, q, n, h, dt, BC_PERIODIC, ok)
      if (.not.ok) error stop 'seth_vv_roe: a smooth wave lost positivity'
      t = t + dt
   enddo
   do i=1, n
      x = -1._R8P + (i - 0.5_R8P)*h
      e(i) = q(1,i) - (1._R8P + 0.2_R8P*sin(PI*x))
   enddo
   l1 = l1_norm(e, h) ; linf = maxval(abs(e))
   endsubroutine wave

   pure function order(e) result(o)
   !< Observed order on the finest pair whose finer error is above 1e-10 (as seth_vv_euler_order).
   real(R8P), intent(in) :: e(:)
   real(R8P)             :: o
   integer(I4P)          :: j
   o = -1._R8P
   do j=size(e)-1, 1, -1
      if (e(j+1) > 1.e-10_R8P) then
         o = log(e(j)/e(j+1))/log(2._R8P) ; return
      endif
   enddo
   endfunction order

   subroutine fail(what)
   character(*), intent(in) :: what
   print '(A)', '   FAILED: '//what
   passed = .false.
   endsubroutine fail
endprogram seth_vv_roe
