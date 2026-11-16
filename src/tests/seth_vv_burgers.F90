!< SETH verification: the inviscid Burgers equation before and after the shock forms, against the exact solution.
program seth_vv_burgers
!< SETH verification: the inviscid Burgers equation before and after the shock forms, against the exact solution.
!<
!< u_t + (u^2/2)_x = 0 on [-1, 1], periodic, u0 = 0.5 + sin(pi x) (Shu 1998, section 4; the classic WENO accuracy
!< test). In the frame moving at 0.5, v = u - 0.5 solves Burgers from sin(pi y): the characteristics y = xi + t
!< sin(pi xi) cross first at y = +-1 when t = 1/pi, and by the odd symmetry of v the shock then stays at y = +-1. The
!< exact v at y in (0, 1) is sin(pi xi) with xi the root of xi + t sin(pi xi) = y on the rising branch [0, xi_m] of the
!< map (xi_m = acos(-1/(pi t))/pi after the shock, 1 before): the characteristics beyond xi_m have entered the shock.
!< Bisection to 1e-15; v(-y) = -v(y).
!<
!< 1. Smooth, t = 0.5/pi: L1 and Linf errors over all cells, N = 40 .. 640 (40 .. 320 in the test run); every
!<    scheme's observed L1 order printed. Asserted on the finest pair above 1e-11: TENO5 and TENO6 >= 4.5 and 5.5
!<    (their order: the solution has first-order critical points, which do not degrade TENO); every scheme converges
!<    (order > 2).
!< 2. Shocked, t = 1.5/pi: L1 over all cells (first order at best: the shock), and over the cells farther than 0.1
!<    from the shock (smooth: high order expected, Shu 1998); the overshoot max(u) - 1.5 and undershoot -0.5 - min(u)
!<    below 1e-2 (no maximum principle is claimed for WENO, but large oscillations would show here). Asserted: the
!<    smooth-region L1 at N = 320 is at least 50 times smaller than the whole-domain one for the quick schemes; the
!<    overshoots bounded.
!< Report mode: the error tables, the profiles at both times, and 61 frames from t = 0 to 1.5/pi (WENO5-JS, TENO6).
use seth
use seth_test_parameters, only : ORDER_FLOOR, TOL_BISECT, TOL_KINK, T_SLACK
use seth_test_schemes
use seth_vv_scalar
use seth_vv_cases, only : vv_report, vv_open, vv_schemes, vv_label, l1_norm
use seth_vv_euler, only : PATH_SPLIT
implicit none
real(R_P), parameter      :: PI = acos(-1._R_P), CFL = 0.4_R_P
integer(I4P), allocatable :: ids(:)
integer(I4P)              :: k, j, nl, u1, u2, ns
real(R_P), allocatable    :: e1(:,:), einf(:,:), es(:,:), et(:,:), osh(:,:)
real(R_P)                 :: o
logical                   :: passed

passed = .true.
call vv_schemes(PATH_SPLIT, .not.vv_report(), ids)
nl = merge(5, 4, vv_report())
allocate(e1(size(ids),nl), einf(size(ids),nl), es(size(ids),nl), et(size(ids),nl), osh(size(ids),nl))
if (vv_report()) then
   call vv_open('burgers-smooth.dat', u1, 'scheme N L1 Linf (t = 0.5/pi)')
   call vv_open('burgers-shock.dat', u2, &
                'scheme N L1-all L1-smooth overshoot (t = 1.5/pi, smooth: > 0.1 from the shock)')
endif
print '(A16,A6,3A12,A8)', 'scheme', 'N', 'L1 smooth', 'L1 shocked', 'L1 away', 'order'
do k=1, size(ids)
   do j=1, nl
      ns = 40*2**(j-1)
      call run(ids(k), ns, 0.5_R_P/PI, e1(k,j), einf(k,j), o)
      call run(ids(k), ns, 1.5_R_P/PI, et(k,j), o, osh(k,j), es(k,j))
      o = -1._R_P ; if (j > 1) o = log(e1(k,j-1)/e1(k,j))/log(2._R_P)
      print '(A16,I6,3ES12.3,F8.2)', vv_label(ids(k)), ns, e1(k,j), et(k,j), es(k,j), o
      if (vv_report()) then
         write(u1, '(A,1X,I6,2ES24.15)') vv_label(ids(k)), ns, e1(k,j), einf(k,j)
         write(u2, '(A,1X,I6,3ES24.15)') vv_label(ids(k)), ns, et(k,j), es(k,j), osh(k,j)
      endif
   enddo
   o = finest(e1(k,:))
   select case(vv_label(ids(k)))
   case('teno5') ; if (o < 4.5_R_P) call fail('TENO5 below order 4.5 before the shock')
   case('teno6') ; if (o < 5.5_R_P) call fail('TENO6 below order 5.5 before the shock')
   end select
   if (o < 2._R_P) call fail(vv_label(ids(k))//' does not converge before the shock')
   ! TENO5-AZ (Tang 2024) is the one exception: its cutoff min(1e-5, tau) with tau = |beta_0 + beta_2 - 2 beta_1|
   ! vanishes where the shock crosses the stencil about symmetrically, every candidate is kept and the state
   ! overshoots, 3.5e-2 at 320 cells (measured, the same with the descaled cutoff); bounded at 5e-2 and printed
   if (vv_label(ids(k)) == 'teno5_az') then
      print '(A,ES10.2)', '   TENO5-AZ overshoot after the shock (documented): ', maxval(osh(k,:))
      if (maxval(osh(k,:)) > 5.e-2_R_P) call fail('teno5_az overshoot above 5e-2 after the shock')
   elseif (maxval(osh(k,:)) > 1.e-2_R_P) then
      call fail(vv_label(ids(k))//' overshoot above 1e-2 after the shock')
   endif
   if (.not.vv_report() .and. .not.(es(k,4) < et(k,4)/50._R_P)) call fail(vv_label(ids(k))//' smooth region polluted')
enddo
if (vv_report()) then
   close(u1) ; close(u2)
   call profiles
endif
if (.not.passed) error stop 'Burgers verification failed'
print '(A)', 'seth_vv_burgers: ok'

contains
   subroutine fail(what)
   character(*), intent(in) :: what
   print '(A)', '   FAILED: '//what
   passed = .false.
   endsubroutine fail

   pure function finest(e) result(o)
   !< Observed order on the finest pair whose finer error is above 1e-11 (the 9th-order errors reach round-off).
   real(R_P), intent(in) :: e(:)
   real(R_P)             :: o
   integer(I4P)          :: j
   o = -1._R_P
   do j=size(e)-1, 1, -1
      if (e(j+1) > ORDER_FLOOR) then
         o = log(e(j)/e(j+1))/log(2._R_P) ; return
      endif
   enddo
   endfunction finest

   pure function exact(x, t) result(u)
   !< The exact solution at x, t.
   real(R_P), intent(in) :: x, t
   real(R_P)             :: u, y
   y = modulo(x - 0.5_R_P*t + 1._R_P, 2._R_P) - 1._R_P ! in [-1, 1)
   if (abs(abs(y) - 1._R_P) < TOL_KINK) then
      u = 0.5_R_P ! on the shock: the mean of the two states
   elseif (y >= 0._R_P) then
      u = 0.5_R_P + v(y, t)
   else
      u = 0.5_R_P - v(-y, t)
   endif
   endfunction exact

   pure function v(y, t) result(s)
   !< v at y in [0, 1).
   real(R_P), intent(in) :: y, t
   real(R_P)             :: s, a, b, c, xm
   integer(I4P)          :: it
   xm = 1._R_P ; if (PI*t > 1._R_P) xm = acos(-1._R_P/(PI*t))/PI
   a = 0._R_P ; b = xm
   do it=1, 200
      c = 0.5_R_P*(a + b)
      if (c + t*sin(PI*c) - y > 0._R_P) then ; b = c ; else ; a = c ; endif
      if (b - a < TOL_BISECT) exit
   enddo
   s = sin(PI*0.5_R_P*(a + b))
   endfunction v

   subroutine evolve(id, n, tend, u, nframes, frames, timeacc)
   !< Advance from u0 to tend; `timeacc`: dt = CFL h (h/h0)^((o-3)/3), h0 = 2/40, so that the O(dt^3) error of RK3
   !< falls with the O(h^o) spatial one (o: the scheme's order; Henrick et al. 2005 use dt ~ h^(5/3) for 5th order).
   integer(I4P), intent(in)            :: id, n
   real(R_P),    intent(in)            :: tend
   real(R_P),    intent(out)           :: u(0:n-1)
   integer(I4P), intent(in),  optional :: nframes
   real(R_P),    intent(out), optional :: frames(0:n-1,0:*)
   logical,      intent(in),  optional :: timeacc
   real(R_P)                           :: h, t, dt, tf, sc
   integer(I4P)                        :: i, kf
   h = 2._R_P/n
   do i=0, n-1
      u(i) = 0.5_R_P + sin(PI*(-1._R_P + (i + 0.5_R_P)*h))
   enddo
   t = 0._R_P ; kf = 1
   if (present(frames)) frames(:,0) = u
   sc = 1._R_P
   if (present(timeacc)) then
      if (timeacc) sc = (h/0.05_R_P)**((SCHEMES(id)%order - 3._R_P)/3._R_P)
   endif
   do while (t < tend*(1._R_P - T_SLACK))
      dt = sc*CFL*h/maxval(abs(u)) ; dt = min(dt, tend - t)
      if (present(frames)) then
         tf = tend*kf/nframes ; dt = min(dt, tf - t)
      endif
      call scalar_step(id, FLUX_BURGERS, u, n, h, dt)
      t = t + dt
      if (present(frames)) then
         if (t >= (tend*kf/nframes)*(1._R_P - T_SLACK) .and. kf <= nframes) then
            frames(:,kf) = u ; kf = kf + 1
         endif
      endif
   enddo
   endsubroutine evolve

   subroutine run(id, n, tend, l1, linf, osh, l1s)
   !< Errors at tend: L1, Linf, overshoot (max(u) - 1.5, -0.5 - min(u)); l1s: L1 farther than 0.1 from the shock.
   integer(I4P), intent(in)            :: id, n
   real(R_P),    intent(in)            :: tend
   real(R_P),    intent(out)           :: l1, linf, osh
   real(R_P),    intent(out), optional :: l1s
   real(R_P)                           :: u(0:n-1), e(0:n-1), es_(0:n-1), h, x, ys
   integer(I4P)                        :: i
   h = 2._R_P/n
   call evolve(id, n, tend, u, timeacc=.true.)
   ys = modulo(1._R_P + 0.5_R_P*tend + 1._R_P, 2._R_P) - 1._R_P ! the shock position (after it forms)
   do i=0, n-1
      x = -1._R_P + (i + 0.5_R_P)*h
      e(i) = u(i) - exact(x, tend)
      es_(i) = 0._R_P
      if (min(abs(x - ys), 2._R_P - abs(x - ys)) > 0.1_R_P) es_(i) = e(i)
   enddo
   l1 = l1_norm(e, h) ; linf = maxval(abs(e))
   osh = max(maxval(u) - 1.5_R_P, -0.5_R_P - minval(u), 0._R_P)
   if (present(l1s)) l1s = l1_norm(es_, h)
   endsubroutine run

   subroutine profiles
   !< Report: profiles at both times for every scheme (N = 80) and frames for WENO5-JS and TENO6 (N = 80).
   integer(I4P), parameter :: N = 80, NF = 60
   real(R_P)               :: u(0:N-1), fr(0:N-1,0:NF), h
   integer(I4P)            :: k, i, uu, f
   character(4)            :: tag
   h = 2._R_P/N
   do k=1, size(ids)
      do f=1, 2
         tag = merge('pre ', 'post', f == 1)
         call evolve(ids(k), N, merge(0.5_R_P, 1.5_R_P, f == 1)/PI, u)
         call vv_open('burgers-'//trim(tag)//'-'//vv_label(ids(k))//'.dat', uu, 'x u u_exact')
         do i=0, N-1
            write(uu, '(3ES24.15)') -1._R_P + (i + 0.5_R_P)*h, u(i), &
                                    exact(-1._R_P + (i + 0.5_R_P)*h, merge(0.5_R_P, 1.5_R_P, f == 1)/PI)
         enddo
         close(uu)
      enddo
      if (vv_label(ids(k)) == 'weno5_js' .or. vv_label(ids(k)) == 'teno6') then
         call evolve(ids(k), N, 1.5_R_P/PI, u, NF, fr)
         call vv_open('burgers-frames-'//vv_label(ids(k))//'.dat', uu, &
                      'frames: line 1 x, then one line per time t_f = 1.5/pi f/60, f = 0..60: t u(1..N) u_exact(1..N)')
         write(uu, '(*(ES16.8))') ((-1._R_P + (i + 0.5_R_P)*h), i=0, N-1)
         do f=0, NF
            write(uu, '(*(ES16.8))') 1.5_R_P/PI*f/NF, fr(:,f), &
                                     (exact(-1._R_P + (i + 0.5_R_P)*h, 1.5_R_P/PI*f/NF), i=0, N-1)
         enddo
         close(uu)
      endif
   enddo
   endsubroutine profiles
endprogram seth_vv_burgers
