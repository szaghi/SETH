!< Experiment: why the centred WENO pollutes Burgers' smooth region after the shock (SETH issue #4, 6.3).
program weno_c_packet
!< Experiment: why the centred WENO pollutes Burgers' smooth region after the shock (SETH issue #4, 6.3).
!<
!< The Burgers case of seth_vv_burgers, reproduced exactly (u0 = 0.5 + sin(pi x) on [-1, 1], periodic, global
!< Lax-Friedrichs splitting with alpha = 1.1 max|u|, SSP-RK3, dt = CFL h (h/0.05)^((o-3)/3)/max|u| with CFL 0.4 and
!< o = 2S, the exact solution by bisection), with the centred WENO's blend coefficient phi in variants:
!<   adam   the library kernel seth_weno{4,6,8}_c (ADAM's scheme, eps = 1e-6, rmu = 1: the registry's parameters)
!<   copy   the same composed here: seth_weno{2S-1}_js + the central row + ADAM's detector on cells -1..2
!<   up     phi = 1: the upwind WENO-JS alone
!<   cen    phi = 0: the central row alone
!<   eps40  ADAM's detector with eps = 1e-40 in the blend (the upwind weights keep eps = 1e-6)
!<   wide   the detector over every pair of the stencil's cells 1-S..S (eps = 1e-6)
!<   wide40 both
!< Outputs (stdout): (1) copy = adam bitwise; (2) L1 before the shock (t = 0.5/pi) and, after it (t = 1.5/pi), over
!< the cells farther than 0.1 from the shock, N = 40 .. 640, per variant and S; (3) phi traced on the final state of
!< adam, N = 160 and 640: phi and the indicators against eps where the error is; (4) the smooth-region error in time
!< after the shock forms (N = 160): its L1 and where its largest value is; (5) the central row's group velocity.
!< Files: packet-trace-S<S>-N<N>.dat (per face: x, error, phi+, phi-, phi+ with eps 1e-40, I1+/eps) and
!< packet-time-S<S>.dat (t, smooth L1, distance of the largest error upstream of the shock).
use seth
use seth_z_teno_tables, only : CEN4_C, CEN6_C, CEN8_C
use seth_test_parameters, only : T_SLACK, TOL_BISECT, TOL_KINK
use seth_vv_cases, only : l1_norm
implicit none
real(R_P),    parameter :: PI = acos(-1._R_P), CFL = 0.4_R_P, EPS = 1.e-6_R_P, EPS40 = 1.e-40_R_P
integer(I4P), parameter :: V_ADAM = 1, V_COPY = 2, V_UP = 3, V_CEN = 4, V_EPS40 = 5, V_WIDE = 6, V_WIDE40 = 7, NV = 7
character(6), parameter :: VNAME(NV) = ['adam  ', 'copy  ', 'up    ', 'cen   ', 'eps40 ', 'wide  ', 'wide40']
integer(I4P), parameter :: NL = 5
integer(I4P)            :: s, v, j, n, ua
real(R_P)               :: e_pre(NL), e_post(NL), e_all(NL), osh(NL)
character(256)          :: arg
logical                 :: quick

quick = .false.
if (command_argument_count() > 0) then
   call get_command_argument(1, arg) ; quick = trim(arg) == 'quick'
endif

print '(A)', '== (1) the copy equals the library kernel'
call check_copy

print '(A)', '== (2) L1 before the shock (t = 0.5/pi), and after it farther than 0.1 from the shock (t = 1.5/pi)'
print '(A8,A7,A6,4A12,A8)', 'scheme', 'phi', 'N', 'L1 pre', 'L1 all', 'L1 away', 'overshoot', 'order'
open(newunit=ua, file='packet-adam.dat', status='replace')
write(ua, '(A)') '# scheme N L1-pre L1-all L1-away overshoot (adam, the library kernel; full precision for run.sh)'
do s=2, 4
   do v=1, NV
      if (v == V_COPY) cycle
      if (quick .and. v /= V_ADAM) cycle
      do j=1, NL
         n = 40*2**(j-1)
         if (quick .and. n > 160) exit
         call run(s, v, n, 0.5_R_P/PI, e_pre(j))
         call run(s, v, n, 1.5_R_P/PI, e_all(j), e_post(j), osh(j))
         print '(A4,I1,A3,A7,I6,4ES12.3,F8.2)', 'weno', 2*s, '_c ', VNAME(v), n, e_pre(j), e_all(j), e_post(j), osh(j), &
               merge(log(e_post(max(j-1,1))/e_post(j))/log(2._R_P), -1._R_P, j > 1)
         if (v == V_ADAM) write(ua, '(A,I0,A,1X,I6,4ES24.15)') 'weno', 2*s, '_c', n, e_pre(j), e_all(j), e_post(j), osh(j)
      enddo
   enddo
enddo
close(ua)
if (quick) stop

print '(A)', '== (3) phi on the final state of adam (t = 1.5/pi), farther than 0.1 from the shock'
do s=3, 4
   call trace(s, 160) ; call trace(s, 640)
enddo

print '(A)', '== (4) the error farther than 0.1 from the shock (or where it forms) in time, t = 0.5/pi .. 1.5/pi (N = 160)'
do s=3, 4
   call history(s, V_ADAM, 160) ; call history(s, V_EPS40, 160) ; call history(s, V_UP, 160)
enddo

print '(A)', '== (5) the central row: modified wavenumber and group velocity (FD flux form)'
call group_velocity

print '(A)', '== (6) a grid-scale packet of amplitude A on a smooth background, linear advection (N = 160, t = 0.25)'
call packet_damping

contains
   subroutine packet_damping
   !< u_t + u_x = 0 (f+ = u, f- = 0), u0 = sin(pi x) + A (-1)^i exp(-((x + 0.5)/0.1)^2); the packet's remaining
   !< amplitude at t = 0.25 over A: the grid-scale part of u(A) - u(0) (the scheme is nonlinear: the background's run
   !< is subtracted, and the smooth part of the difference left out), for A = 1e-1 .. 1e-7; and without the
   !< background (u0 = the packet alone). dt as in the Burgers runs.
   real(R_P)               :: a, r_bg, r_alone
   integer(I4P)            :: s, v, k
   integer(I4P), parameter :: VS(3) = [V_ADAM, V_UP, V_CEN]
   print '(A)', '   remaining fraction of the packet: with the background | alone'
   print '(A10,3(A9,A10,A10))', 'A', ('', VNAME(VS(v)), '', v=1,3)
   do s=3, 4
      print '(A,I1,A)', '   WENO', 2*s, '-C'
      do k=1, 7
         a = 10._R_P**(-k)
         write(*, '(ES10.1)', advance='no') a
         do v=1, 3
            r_bg = remaining(s, VS(v), a, 1._R_P) ; r_alone = remaining(s, VS(v), a, 0._R_P)
            write(*, '(A9,2ES10.2)', advance='no') '', r_bg, r_alone
         enddo
         write(*, *)
      enddo
   enddo
   endsubroutine packet_damping

   function remaining(s, v, a, bg) result(r)
   integer(I4P), intent(in) :: s, v
   real(R_P),    intent(in) :: a, bg
   real(R_P)                :: r
   integer(I4P), parameter  :: N = 160
   real(R_P)                :: ua(0:N-1), ub(0:N-1), h, x, t, dt, sc
   integer(I4P)             :: i, k
   h = 2._R_P/N
   do i=0, N-1
      x = -1._R_P + (i + 0.5_R_P)*h
      ub(i) = bg*sin(PI*x)
      ua(i) = ub(i) + a*(-1._R_P)**i*exp(-((x + 0.5_R_P)/0.1_R_P)**2)
   enddo
   sc = (h/0.05_R_P)**((2._R_P*s - 3._R_P)/3._R_P)
   dt = sc*CFL*h
   do k=1, 2
      t = 0._R_P
      do while (t < 0.25_R_P*(1._R_P - T_SLACK))
         if (k == 1) then
            call lstep(s, v, ua, N, h, min(dt, 0.25_R_P - t))
         else
            call lstep(s, v, ub, N, h, min(dt, 0.25_R_P - t))
         endif
         t = t + min(dt, 0.25_R_P - t)
      enddo
   enddo
   ub = ua - ub ! the difference made by the packet; its grid-scale part, |d_i - (d_{i-1} + d_{i+1})/2|/2, is the
   do i=0, N-1  ! amplitude of a (-1)^i component, and excludes the smooth response of the nonlinear weights
      ua(i) = 0.5_R_P*abs(ub(i) - 0.5_R_P*(ub(modulo(i-1, N)) + ub(modulo(i+1, N))))
   enddo
   r = maxval(ua)/a
   endfunction remaining

   subroutine lstep(s, v, u, n, h, dt)
   !< One SSP-RK3 step of u_t + u_x = 0 (upwind splitting: f+ = u).
   integer(I4P), intent(in)    :: s, v, n
   real(R_P),    intent(inout) :: u(0:n-1)
   real(R_P),    intent(in)    :: h, dt
   real(R_P)                   :: u1(0:n-1), u2(0:n-1), r(0:n-1)
   call lres(s, v, u, n, h, r)  ; u1 = u + dt*r
   call lres(s, v, u1, n, h, r) ; u2 = 0.75_R_P*u + 0.25_R_P*(u1 + dt*r)
   call lres(s, v, u2, n, h, r) ; u = u/3._R_P + 2._R_P/3._R_P*(u2 + dt*r)
   endsubroutine lstep

   subroutine lres(s, v, w, n, h, res)
   integer(I4P), intent(in)  :: s, v, n
   real(R_P),    intent(in)  :: w(0:n-1), h
   real(R_P),    intent(out) :: res(0:n-1)
   real(R_P)                 :: fh(0:n-1), vp(-4:5), vm(-4:5), rp, rm
   integer(I4P)              :: i, m
   vm = 0._R_P
   do i=0, n-1
      do m=-4, 5
         vp(m) = w(modulo(i+m, n))
      enddo
      call recon(s, v, vp, vm, rp, rm)
      fh(i) = rp
   enddo
   do i=0, n-1
      res(i) = -(fh(i) - fh(modulo(i-1, n)))/h
   enddo
   endsubroutine lres
   subroutine recon(s, v, vp, vm, rp, rm)
   !< The variant v of the centred WENO of order 2S, from the stencils vp, vm(-4:5).
   integer(I4P), intent(in)  :: s, v
   real(R_P),    intent(in)  :: vp(-4:5), vm(-4:5)
   real(R_P),    intent(out) :: rp, rm
   real(R_P)                 :: up, um, cp, cm, pp, pm
   if (v == V_ADAM) then
      select case(s)
      case(2) ; call seth_weno4_c(vp(-1:2), vm(-1:2), EPS, 1._R_P, rp, rm)
      case(3) ; call seth_weno6_c(vp(-2:3), vm(-2:3), EPS, 1._R_P, rp, rm)
      case(4) ; call seth_weno8_c(vp(-3:4), vm(-3:4), EPS, 1._R_P, rp, rm)
      endselect
      return
   endif
   call upwind(s, vp, vm, up, um)
   call central(s, vp, vm, cp, cm)
   select case(v)
   case(V_COPY) ; pp = phi_pairs(vp(-1:2), EPS)        ; pm = phi_pairs(vm(2:-1:-1), EPS)
   case(V_UP)   ; pp = 1._R_P                          ; pm = 1._R_P
   case(V_CEN)  ; pp = 0._R_P                          ; pm = 0._R_P
   case(V_EPS40); pp = phi_pairs(vp(-1:2), EPS40)      ; pm = phi_pairs(vm(2:-1:-1), EPS40)
   case(V_WIDE) ; pp = phi_pairs(vp(1-s:s), EPS)       ; pm = phi_pairs(vm(s:1-s:-1), EPS)
   case(V_WIDE40); pp = phi_pairs(vp(1-s:s), EPS40)    ; pm = phi_pairs(vm(s:1-s:-1), EPS40)
   endselect
   rp = cp + pp*(up - cp)
   rm = cm + pm*(um - cm)
   endsubroutine recon

   subroutine upwind(s, vp, vm, up, um)
   integer(I4P), intent(in)  :: s
   real(R_P),    intent(in)  :: vp(-4:5), vm(-4:5)
   real(R_P),    intent(out) :: up, um
   select case(s)
   case(2) ; call seth_weno3_js(vp(-1:2), vm(-1:2), EPS, 1._R_P, up, um)
   case(3) ; call seth_weno5_js(vp(-2:3), vm(-2:3), EPS, 1._R_P, up, um)
   case(4) ; call seth_weno7_js(vp(-3:4), vm(-3:4), EPS, 1._R_P, up, um)
   endselect
   endsubroutine upwind

   subroutine central(s, vp, vm, cp, cm)
   !< The central row of order 2S, in the library's summation order.
   integer(I4P), intent(in)  :: s
   real(R_P),    intent(in)  :: vp(-4:5), vm(-4:5)
   real(R_P),    intent(out) :: cp, cm
   real(R_P)                 :: c(-3:4)
   integer(I4P)              :: m
   call row(s, c)
   cp = 0._R_P ; cm = 0._R_P
   do m=1-s, s
      cp = cp + c(m)*vp(m)
      cm = cm + c(m)*vm(1-m)
   enddo
   endsubroutine central

   subroutine row(s, c)
   integer(I4P), intent(in)  :: s
   real(R_P),    intent(out) :: c(-3:4)
   c = 0._R_P
   select case(s)
   case(2) ; c(-1:2) = CEN4_C
   case(3) ; c(-2:3) = CEN6_C
   case(4) ; c(-3:4) = CEN8_C
   endselect
   endsubroutine row

   pure function phi_pairs(a, e) result(phi)
   !< ADAM's blend coefficient over consecutive cells a(1..k) (k = 4: ADAM's detector, the library's arithmetic):
   !< the largest |I_m - I_{m-1}|/(I_m + I_{m-1}) with I_m = (a(m+1) - a(m))^2 + e.
   real(R_P), intent(in) :: a(:), e
   real(R_P)             :: phi, i0, i1
   integer(I4P)          :: m
   phi = 0._R_P
   i0 = (a(2) - a(1))**2 + e
   do m=2, size(a)-1
      i1 = (a(m+1) - a(m))**2 + e
      phi = max(phi, abs(i1 - i0)/(i1 + i0))
      i0 = i1
   enddo
   endfunction phi_pairs

   subroutine check_copy
   !< Random stencils, every S: the copy's states equal the library kernel's bitwise.
   real(R_P)    :: vp(-4:5), vm(-4:5), r1p, r1m, r2p, r2m, worst
   integer(I4P) :: it, s
   worst = 0._R_P
   do s=2, 4
      do it=1, 100000
         call random_number(vp) ; call random_number(vm)
         if (mod(it, 3) == 0) vp(1:) = vp(1:) + 1._R_P
         if (mod(it, 5) == 0) vm = vm*1.e-3_R_P
         call recon(s, V_ADAM, vp, vm, r1p, r1m)
         call recon(s, V_COPY, vp, vm, r2p, r2m)
         worst = max(worst, abs(r1p - r2p), abs(r1m - r2m))
      enddo
   enddo
   print '(A,ES10.2)', '   largest |copy - library| over 3 x 100000 stencils: ', worst
   if (worst > 0._R_P) error stop 'the copy is not the library kernel'
   endsubroutine check_copy

   subroutine residual(s, v, u, n, h, res)
   !< seth_vv_scalar's residual (Burgers, global Lax-Friedrichs splitting), with the variant v.
   integer(I4P), intent(in)  :: s, v, n
   real(R_P),    intent(in)  :: u(0:n-1), h
   real(R_P),    intent(out) :: res(0:n-1)
   real(R_P)                 :: f(0:n-1), fh(0:n-1), alpha, vp(-4:5), vm(-4:5), rp, rm
   integer(I4P)              :: i, m, jj
   f = 0.5_R_P*u*u ; alpha = 1.1_R_P*maxval(abs(u))
   do i=0, n-1
      do m=-4, 5
         jj = modulo(i+m, n)
         vp(m) = 0.5_R_P*(f(jj) + alpha*u(jj))
         vm(m) = 0.5_R_P*(f(jj) - alpha*u(jj))
      enddo
      call recon(s, v, vp, vm, rp, rm)
      fh(i) = rp + rm
   enddo
   do i=0, n-1
      res(i) = -(fh(i) - fh(modulo(i-1, n)))/h
   enddo
   endsubroutine residual

   subroutine step(s, v, u, n, h, dt)
   integer(I4P), intent(in)    :: s, v, n
   real(R_P),    intent(inout) :: u(0:n-1)
   real(R_P),    intent(in)    :: h, dt
   real(R_P)                   :: u1(0:n-1), u2(0:n-1), r(0:n-1)
   call residual(s, v, u, n, h, r)  ; u1 = u + dt*r
   call residual(s, v, u1, n, h, r) ; u2 = 0.75_R_P*u + 0.25_R_P*(u1 + dt*r)
   call residual(s, v, u2, n, h, r) ; u = u/3._R_P + 2._R_P/3._R_P*(u2 + dt*r)
   endsubroutine step

   subroutine init(n, u)
   integer(I4P), intent(in)  :: n
   real(R_P),    intent(out) :: u(0:n-1)
   integer(I4P)              :: i
   do i=0, n-1
      u(i) = 0.5_R_P + sin(PI*(-1._R_P + (i + 0.5_R_P)*2._R_P/n))
   enddo
   endsubroutine init

   subroutine advance(s, v, n, u, t, tend)
   !< seth_vv_burgers' evolve (timeacc), from t to tend.
   integer(I4P), intent(in)    :: s, v, n
   real(R_P),    intent(inout) :: u(0:n-1), t
   real(R_P),    intent(in)    :: tend
   real(R_P)                   :: h, dt, sc
   h = 2._R_P/n
   sc = (h/0.05_R_P)**((2._R_P*s - 3._R_P)/3._R_P)
   do while (t < tend*(1._R_P - T_SLACK))
      dt = sc*CFL*h/maxval(abs(u)) ; dt = min(dt, tend - t)
      call step(s, v, u, n, h, dt)
      t = t + dt
   enddo
   endsubroutine advance

   subroutine errors(n, u, tend, e, away)
   !< Pointwise error, and the mask of the cells farther than 0.1 from the shock.
   integer(I4P), intent(in)  :: n
   real(R_P),    intent(in)  :: u(0:n-1), tend
   real(R_P),    intent(out) :: e(0:n-1)
   logical,      intent(out) :: away(0:n-1)
   real(R_P)                 :: h, x, ys
   integer(I4P)              :: i
   h = 2._R_P/n
   ys = modulo(1._R_P + 0.5_R_P*tend + 1._R_P, 2._R_P) - 1._R_P
   do i=0, n-1
      x = -1._R_P + (i + 0.5_R_P)*h
      e(i) = u(i) - exact(x, tend)
      away(i) = min(abs(x - ys), 2._R_P - abs(x - ys)) > 0.1_R_P
   enddo
   endsubroutine errors

   subroutine run(s, v, n, tend, l1, l1s, osh)
   integer(I4P), intent(in)            :: s, v, n
   real(R_P),    intent(in)            :: tend
   real(R_P),    intent(out)           :: l1
   real(R_P),    intent(out), optional :: l1s, osh
   real(R_P)                           :: u(0:n-1), e(0:n-1), t
   logical                             :: away(0:n-1)
   call init(n, u) ; t = 0._R_P
   call advance(s, v, n, u, t, tend)
   call errors(n, u, tend, e, away)
   l1 = l1_norm(e, 2._R_P/n)
   if (present(l1s)) l1s = l1_norm(merge(e, 0._R_P, away), 2._R_P/n)
   if (present(osh)) osh = max(maxval(u) - 1.5_R_P, -0.5_R_P - minval(u), 0._R_P)
   endsubroutine run

   subroutine trace(s, n)
   !< phi on the final state of adam, at every face farther than 0.1 from the shock.
   integer(I4P), intent(in) :: s, n
   real(R_P)                :: u(0:n-1), e(0:n-1), t, f(0:n-1), alpha, vp(-4:5), vm(-4:5), h, x
   real(R_P)                :: pp(0:n-1), pm(0:n-1), p40(0:n-1), ri(0:n-1), emax, ue(0:n-1), fe(0:n-1), pex(0:n-1)
   logical                  :: away(0:n-1), hot(0:n-1)
   integer(I4P)             :: i, m, jj, uu, nsc
   character(64)            :: fname
   h = 2._R_P/n
   call init(n, u) ; t = 0._R_P
   call advance(s, V_ADAM, n, u, t, 1.5_R_P/PI)
   call errors(n, u, 1.5_R_P/PI, e, away)
   f = 0.5_R_P*u*u ; alpha = 1.1_R_P*maxval(abs(u))
   do i=0, n-1
      ue(i) = exact(-1._R_P + (i + 0.5_R_P)*h, 1.5_R_P/PI)
   enddo
   fe = 0.5_R_P*ue*ue
   do i=0, n-1
      do m=-4, 5
         jj = modulo(i+m, n)
         vp(m) = 0.5_R_P*(fe(jj) + alpha*ue(jj))
      enddo
      pex(i) = phi_pairs(vp(-1:2), EPS) ! the same faces, the exact solution (the background without the packet)
      do m=-4, 5
         jj = modulo(i+m, n)
         vp(m) = 0.5_R_P*(f(jj) + alpha*u(jj))
         vm(m) = 0.5_R_P*(f(jj) - alpha*u(jj))
      enddo
      pp(i) = phi_pairs(vp(-1:2), EPS) ; pm(i) = phi_pairs(vm(2:-1:-1), EPS)
      p40(i) = phi_pairs(vp(-1:2), EPS40)
      ri(i) = (vp(1) - vp(0))**2/EPS
   enddo
   ! the face i+1/2 is in the smooth region when both its cells are
   do i=0, n-1
      away(i) = away(i) .and. away(modulo(i+1, n))
   enddo
   emax = maxval(abs(e), mask=away)
   hot = away .and. abs(e) > 0.1_R_P*emax ! where the packet is
   nsc = 0
   do i=0, n-2
      if (hot(i) .and. hot(i+1) .and. e(i)*e(i+1) < 0._R_P) nsc = nsc + 1
   enddo
   print '(A,I1,A,I4,A,ES9.2,A,I4,A,I4,A)', '   S = ', s, ', N = ', n, ': largest smooth error ', emax, ' (', &
         count(hot), ' cells above 0.1 of it, ', nsc, ' sign changes between neighbours there)'
   print '(A,3ES10.2)', '     phi+ where the packet is: max, mean, min        ', maxval(pp, mask=hot), &
         sum(pp, mask=hot)/max(count(hot), 1), minval(pp, mask=hot)
   print '(A,3ES10.2)', '     phi+ over the smooth region: max, mean, min     ', maxval(pp, mask=away), &
         sum(pp, mask=away)/max(count(away), 1), minval(pp, mask=away)
   print '(A,3ES10.2)', '     phi+ with eps 1e-40 where the packet is          ', maxval(p40, mask=hot), &
         sum(p40, mask=hot)/max(count(hot), 1), minval(p40, mask=hot)
   print '(A,3ES10.2)', '     phi+ of the exact solution where the packet is   ', maxval(pex, mask=hot), &
         sum(pex, mask=hot)/max(count(hot), 1), minval(pex, mask=hot)
   print '(A,ES10.2)',  '     mean |phi+ - phi+ of the exact solution| there  ', &
         sum(abs(pp - pex), mask=hot)/max(count(hot), 1)
   print '(A,3ES10.2)', '     (Delta f+)^2/eps where the packet is             ', maxval(ri, mask=hot), &
         sum(ri, mask=hot)/max(count(hot), 1), minval(ri, mask=hot)
   write(fname, '(A,I1,A,I0,A)') 'packet-trace-S', s, '-N', n, '.dat'
   open(newunit=uu, file=trim(fname), status='replace')
   write(uu, '(A)') '# x_face error(cell i) phi+ phi- phi+(eps 1e-40) phi+(exact) (Delta f+)^2/eps away hot (t = 1.5/pi, adam)'
   do i=0, n-1
      x = -1._R_P + (i + 1._R_P)*h
      write(uu, '(7ES16.7,2L3)') x, e(i), pp(i), pm(i), p40(i), pex(i), ri(i), away(i), hot(i)
   enddo
   close(uu)
   endsubroutine trace

   subroutine history(s, v, n)
   !< The error farther than 0.1 from the shock (before t = 1/pi: from the point where it forms), t = 0.5/pi .. 1.5/pi.
   integer(I4P), intent(in) :: s, v, n
   real(R_P)                :: u(0:n-1), e(0:n-1), t, tk, ys, x, d, h
   logical                  :: away(0:n-1)
   integer(I4P)             :: k, i, uu
   character(64)            :: fname
   h = 2._R_P/n
   call init(n, u) ; t = 0._R_P
   write(fname, '(A,I1,A,A,A)') 'packet-time-S', s, '-', trim(VNAME(v)), '.dat'
   open(newunit=uu, file=trim(fname), status='replace')
   write(uu, '(A)') '# t smooth-L1 max|e| distance-from-max|e|-to-the-shock-along-+x '// &
                    '(N = 160; 2 - d: along -x; before t = 1/pi the shock-to-be)'
   print '(A,I1,A,A,A)', '   S = ', s, ', ', trim(VNAME(v)), ': t*pi, smooth L1, max|e|, distance from it to the shock along +x'
   do k=0, 40
      tk = (0.5_R_P + k/40._R_P)/PI
      call advance(s, v, n, u, t, tk)
      call errors(n, u, tk, e, away)
      ys = modulo(1._R_P + 0.5_R_P*tk + 1._R_P, 2._R_P) - 1._R_P
      i = maxloc(abs(e), dim=1, mask=away) - 1
      x = -1._R_P + (i + 0.5_R_P)*h
      d = modulo(ys - x, 2._R_P) ! the flow runs to +x: the distance from x to the shock ahead of it
      write(uu, '(4ES16.7)') tk, l1_norm(merge(e, 0._R_P, away), h), maxval(abs(e), mask=away), d
      if (mod(k, 4) == 0) print '(F10.3,2ES12.3,F8.3)', tk*PI, l1_norm(merge(e, 0._R_P, away), h), &
                                maxval(abs(e), mask=away), d
   enddo
   close(uu)
   endsubroutine history

   subroutine group_velocity
   !< The FD flux-form derivative with the central row c: theta* = -i sum_m c(m) (e^{i m theta} - e^{i (m-1) theta});
   !< c is symmetric about 1/2, so theta* is real (no dissipation); its slope d theta*/d theta is the group velocity
   !< relative to the exact one.
   real(R_P)    :: c(-3:4), th, ts, ts0, g
   complex(R_P) :: z
   integer(I4P) :: s, k, m
   do s=2, 4
      call row(s, c)
      print '(A,I1,A)', '   CEN', 2*s, ': theta/pi, theta*/pi, group velocity'
      ts0 = 0._R_P
      do k=1, 20
         th = PI*k/20._R_P
         z = (0._R_P, 0._R_P)
         do m=1-s, s
            z = z + c(m)*(exp(cmplx(0._R_P, m*th, R_P)) - exp(cmplx(0._R_P, (m-1)*th, R_P)))
         enddo
         ts = aimag(z)
         g = (ts - ts0)/(PI/20._R_P) ; ts0 = ts
         if (mod(k, 2) == 0) print '(3F10.3)', th/PI, ts/PI, g
      enddo
   enddo
   endsubroutine group_velocity

   pure function exact(x, t) result(u)
   !< seth_vv_burgers' exact solution.
   real(R_P), intent(in) :: x, t
   real(R_P)             :: u, y
   y = modulo(x - 0.5_R_P*t + 1._R_P, 2._R_P) - 1._R_P
   if (abs(abs(y) - 1._R_P) < TOL_KINK) then
      u = 0.5_R_P
   elseif (y >= 0._R_P) then
      u = 0.5_R_P + vexact(y, t)
   else
      u = 0.5_R_P - vexact(-y, t)
   endif
   endfunction exact

   pure function vexact(y, t) result(sv)
   real(R_P), intent(in) :: y, t
   real(R_P)             :: sv, a, b, c, xm
   integer(I4P)          :: it
   xm = 1._R_P ; if (PI*t > 1._R_P) xm = acos(-1._R_P/(PI*t))/PI
   a = 0._R_P ; b = xm
   do it=1, 200
      c = 0.5_R_P*(a + b)
      if (c + t*sin(PI*c) - y > 0._R_P) then ; b = c ; else ; a = c ; endif
      if (b - a < TOL_BISECT) exit
   enddo
   sv = sin(PI*0.5_R_P*(a + b))
   endfunction vexact
endprogram weno_c_packet
