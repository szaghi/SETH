!< SETH verification: the order of accuracy of the 1-D Euler solver on a smooth exact solution, both paths.
program seth_vv_euler_order
!< SETH verification: the order of accuracy of the 1-D Euler solver on a smooth exact solution, both paths.
!<
!< A density wave advected at constant velocity and pressure: rho = 1 + 0.2 sin(pi x), u = 1, p = 1 on [-1, 1],
!< periodic, one period (t = 2); the exact solution is the initial one. Every scheme of both paths, N = 20 .. 320
!< (20 .. 160 in the test run), dt = 0.4 h (h/h0)^((o-3)/3), h0 = 2/20, so that the O(dt^3) error of SSP-RK3 falls with
!< the spatial one (o: the scheme's order, capped at 6 on the interpolation path). L1 and Linf of the density.
!<
!< Expected orders: the order of the reconstruction on the splitting path, min(order, 6) on the interpolation path
!< (the flux corrections are O(h^6), seth_vv_euler). The wave has first-order critical points (rho' = 0 at x = +-1/2),
!< where the WENO-JS weights lose accuracy (Henrick et al. 2005) and WENO-Z3 drops: those are printed, not asserted.
!< The positivity limiter is on, and must never act (a smooth flow). Asserted: every WENO-Z (but Z3) and TENO scheme
!< at its expected order - 0.5 in L1 on the finest pair above 1e-10 (round-off accumulates below); every scheme at
!< least order 2; no limited face.
!< Report mode: euler-order.dat (scheme N L1 Linf).
use seth
use seth_test_schemes
use seth_vv_euler
use seth_vv_cases
implicit none
real(R8P), parameter      :: PI = acos(-1._R8P)
integer(I4P), allocatable :: ids(:)
integer(I4P)              :: path, k, j, nl, n, u
real(R8P), allocatable    :: e1(:), ei(:)
real(R8P)                 :: o, need
logical                   :: passed
passed = .true.
nl = merge(5, 4, vv_report())
allocate(e1(nl), ei(nl))
if (vv_report()) call vv_open('euler-order.dat', u, 'scheme N L1(rho) Linf(rho) (density wave, t = 2)')
limiter = .true.
do path=PATH_SPLIT, PATH_INTERP
   call vv_schemes(path, .not.vv_report(), ids)
   do k=1, size(ids)
      do j=1, nl
         n = 20*2**(j-1)
         nlimited = 0
         call wave(ids(k), path, n, e1(j), ei(j))
         if (nlimited > 0) call fail(vv_label(ids(k))//': the positivity limiter acted on a smooth flow')
         if (vv_report()) write(u, '(A,1X,I6,2ES24.15)') vv_label(ids(k)), n, e1(j), ei(j)
      enddo
      o = order(e1)
      need = SCHEMES(ids(k))%order ; if (path == PATH_INTERP) need = min(need, 6._R8P)
      print '(A18,2ES12.4,F8.2,A,F4.1)', vv_label(ids(k)), e1(nl), ei(nl), o, '  expected ', need
      if (o < 2._R8P) call fail(vv_label(ids(k))//': below second order')
      if (SCHEMES(ids(k))%family == FAM_TENO .or. (SCHEMES(ids(k))%family == FAM_Z .and. SCHEMES(ids(k))%s >= 3)) then
         if (o < need - 0.5_R8P) call fail(vv_label(ids(k))//': below its expected order')
      endif
   enddo
enddo
if (vv_report()) close(u)
if (.not.passed) error stop 'Euler order verification failed'
print '(A)', 'seth_vv_euler_order: ok'

contains
   subroutine wave(id, path, n, l1, linf)
   integer(I4P), intent(in)  :: id, path, n
   real(R8P),    intent(out) :: l1, linf
   real(R8P)                 :: q(3,1-NG:n+NG), h, t, dt, sc, e(n), x, o
   logical                   :: ok
   integer(I4P)              :: i
   h = 2._R8P/n
   do i=1, n
      x = -1._R8P + (i - 0.5_R8P)*h
      q(:,i) = prim2cons([1._R8P + 0.2_R8P*sin(PI*x), 1._R8P, 1._R8P])
   enddo
   o = SCHEMES(id)%order ; if (path == PATH_INTERP) o = min(o, 6._R8P)
   sc = (h/0.1_R8P)**((o - 3._R8P)/3._R8P)
   t = 0._R8P
   do while (t < 2._R8P*(1._R8P - 1.e-14_R8P))
      dt = min(0.4_R8P*sc*euler_dt(q, n, h, 1._R8P), 2._R8P - t)
      call euler_step(id, path, q, n, h, dt, BC_PERIODIC, ok)
      if (.not.ok) error stop 'seth_vv_euler_order: a smooth wave lost positivity'
      t = t + dt
   enddo
   do i=1, n
      x = -1._R8P + (i - 0.5_R8P)*h
      e(i) = q(1,i) - (1._R8P + 0.2_R8P*sin(PI*x))
   enddo
   l1 = l1_norm(e, h) ; linf = maxval(abs(e))
   endsubroutine wave

   pure function order(e) result(o)
   !< Observed order on the finest pair whose finer error is above 1e-10: below it the round-off accumulated over the
   !< tens of thousands of RK3 steps of the time-accurate dt (a plateau at 1e-12 .. 1e-11, measured) takes over.
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
endprogram seth_vv_euler_order
