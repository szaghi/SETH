!< SETH verification: TENO5-DR and TENO5-THINC in the 1-D Euler harness (issue #6, steps 2.2 and 2.4).
program seth_vv_nonpoly
!< SETH verification: TENO5-DR and TENO5-THINC in the 1-D Euler harness (issue #6, steps 2.2 and 2.4).
!<
!< The splitting path of seth_vv_euler (characteristic Lax-Friedrichs splitting, SSP-RK3, CFL 0.5), transmissive ends
!< (the blast: reflective, with the positivity limiter, as seth_vv_blast).
!<
!< 1. Step 2.2, spec 7 (frozen in issue #6): Sod and Lax at N = 100: L1 of the density of TENO5-DR not above TENO5's,
!<    and strictly fewer cells across the contact, a cell counting when its density lies between 5 % and 95 % of the
!<    exact contact jump (the exact star densities, seth_vv_riemann) within 0.1 of the exact contact. TENO5-THINC is
!<    reported beside them.
!< 2. Experiment E1 (step 2.4): the self-contained kernels recompute the neighbours' masks in the characteristic basis
!<    of the face being reconstructed; the stored-mask variant (seth_vv_euler's `stored_masks`) computes each cell's
!<    masks once, in the basis of that cell's own face. On scalar data the two are equal bitwise (seth_test_nonpoly,
!<    N1); on Euler they differ where the bases differ. Sod, Lax (N = 200, against the exact solution), Shu-Osher
!<    (N = 200 on [-5, 5]) and the blast waves (N = 400): L1 against the exact solution where there is one, the largest
!<    density difference between the two variants, positivity. Asserted: both variants positive (the blast with the
!<    limiter); the rest is reported (an experiment: neither variant is better on every case, measured 2026-11-17).
!<    Report mode: nonpoly-e1.dat and nonpoly-contact.dat.
!< 3. N4 and N6 (E4) on linear advection, u_t + u_x = 0, periodic, a square wave on a sine, N = 200 cells, SSP-RK3 at
!<    CFL 0.4, one period, both composites (the overrides active: counted): N4, the sum of u conserved to NP_TOL of its
!<    magnitude (1e-13 in R8P, about 450 ulps in either kind; measured 83 ulps); N6, the domain split in two blocks, each with G ghost cells filled from the other block at every stage:
!<    with G = 6 the result equals the single block bitwise, with G = 5 (the sixth cell left unfilled, zero: a harness
!<    mutant) it differs. The halo width of the composite kernels is 6, proved. Report:
!<    nonpoly-halo.dat.
!< 5. N5, the ablation (step 3.4): TENO5, ST, RBF, THINC, DR, ST+DR, RBF+DR (the RBF kernels with s = 3 dx, Liang &
!<    Fu's alpha = 3) on Sod and Lax (N = 200, against the exact solution: L1 of the density),
!<    Shu-Osher (N = 200 on [-5, 5], t = 1.8, against WENO5-JS on 3200 cells: L1) and the blast waves (N = 400, the
!<    limiter on, against WENO5-JS on 6400 cells: L1, faces limited). Report mode only (the references cost minutes);
!<    asserted there: positivity (the blast with the limiter); reported: nonpoly-ablation.dat.
!< 6. E5 (step 3.3): Shu-Osher on [-5, 5] and the same flow rescaled by 1/8 (x and t: a power of 2, so that every
!<    operation scales exactly; with 1/10 the rounding of h/10 flips TENO cutoffs and moves the profile by 1e-4), N = 200:
!<    TENO5's density profiles equal bitwise (asserted: the scheme has no length scale); TENO5-RBF's differ by the change
!<    of s = 3 dx from 0.15 to 0.019 (reported, nonpoly-units.dat).
!< 7. E6 (step 3.3): sin(2 pi x) advected 100 periods by TENO5-RBF with lambda frozen at 3 and 6 (s = lambda dx; the
!<    sensor keeps s = s1 on this wave), N = 20 and 40, SSP-RK3 at CFL 0.4 and 0.05: the growth rate of the amplitude,
!<    log(A(t)/A(0))/t, against the prediction of Eq. 29, dx^5 k^4 (lambda^2/2 - k^2/60), and TENO5's (the same run with
!<    s = 0). Asserted (D7), inside the validated range (s <= 0.12): the RBF growth rate minus TENO5's at most twice the
!<    predicted lambda term, dx^5 k^4 lambda^2/2, and the amplitude below twice the initial one. Outside it reported
!<    only (at s = 0.3 and CFL 0.05 the amplitude nearly triples in 100 periods). Reported: nonpoly-e6.dat.
!< 4. E2 (step 2.5): a square wave of width 0.4 advected 10 periods (N = 100, CFL 0.4): its width in cells (between
!<    5 % and 95 % of the unit jump, both jumps) after each period, TENO5, TENO5-DR, TENO5-THINC. E3: a smooth
!<    periodic plateau 0.5 (tanh((x - 0.25)/(w h)) - tanh((x - 0.75)/(w h))) of w = 1, 2, 3, 5, 10 cells (its two
!<    gradients w cells wide) advected one period (N = 100): L1 against the exact
!<    profile, each composite against TENO5 (the width where the override starts to steepen a smooth gradient into a
!<    step). Reported (nonpoly-e2.dat, nonpoly-e3.dat); asserted only that every run stays in [-0.01, 1.01].
use seth
use seth_vv_euler
use seth_vv_cases
use seth_vv_riemann, only : riemann_sample
use seth_test_parameters, only : NP_TOL
use seth_test_schemes, only : SCHEMES, NALL, reconstruct, rbf_dx
implicit none
integer(I4P), parameter :: NS = 3_I4P
character(16), parameter :: NAMES(NS) = [character(16) :: 'teno5', 'teno5_dr', 'teno5_thinc']
integer(I4P)            :: ids(NS), k, j, ue, uc
real(R_P)               :: l1(NS), wid(NS), l1s(2), dq
logical                 :: passed, pos(NS), posd, poss
real(R_P), parameter    :: SOD_L(3) = [1._R_P, 0._R_P, 1._R_P], SOD_R(3) = [0.125_R_P, 0._R_P, 0.1_R_P]
real(R_P), parameter    :: LAX_L(3) = [0.445_R_P, 0.698_R_P, 3.528_R_P], LAX_R(3) = [0.5_R_P, 0._R_P, 0.571_R_P]

passed = .true.
do k=1, NS
   ids(k) = 0
   ids(k) = scheme_id(NAMES(k))
   if (ids(k) == 0) error stop 'seth_vv_nonpoly: scheme not in the registry'
enddo
if (vv_report()) call vv_open('nonpoly-contact.dat', uc, 'case scheme L1(rho) contact-cells positive (N = 100)')
if (vv_report()) call vv_open('nonpoly-e1.dat', ue, &
   'case scheme L1-default L1-stored max|rho_default-rho_stored| positive-default positive-stored')

! 1. spec 7 of step 2.2
call contact('sod', SOD_L, SOD_R, 0.2_R_P)
call contact('lax', LAX_L, LAX_R, 0.13_R_P)

! 3. N4, N6
call halo()

! 4. E2, E3
call e2_e3()

! 5. N5, the ablation; 6. E5, the length unit
if (vv_report()) call ablation() ! the references on 16 times the cells: report mode only
call units()

! 7. E6, the long-wave growth of the RBF scheme in time
call e6()

! 2. E1: default against stored masks
print '(A)', 'E1: self-contained (masks in the face''s basis) against stored masks (each in its own face''s basis)'
print '(A12,A14,3A13,2A5)', 'case', 'scheme', 'L1 default', 'L1 stored', 'max |drho|', 'pos', 'pos'
do k=2, NS
   call e1_tube('sod', ids(k), SOD_L, SOD_R, 0.2_R_P)
   call e1_tube('lax', ids(k), LAX_L, LAX_R, 0.13_R_P)
   call e1_ref('shu-osher', ids(k))
   call e1_ref('blast', ids(k))
enddo
if (vv_report()) then ; close(ue) ; close(uc) ; endif
if (.not.passed) error stop 'seth_vv_nonpoly failed'
print '(A)', 'seth_vv_nonpoly: ok'

contains
   subroutine contact(name, wl, wr, tend)
   !< L1 and the cells across the contact at N = 100, TENO5, TENO5-DR, TENO5-THINC.
   character(*), intent(in) :: name
   real(R_P),    intent(in) :: wl(3), wr(3), tend
   integer(I4P), parameter  :: N = 100_I4P
   real(R_P)                :: w0(3,N), q(3,1-NG:N+NG), h, x, xc, rl, rr, lo, hi, w(3), us
   integer(I4P)             :: k, i
   h = 1._R_P/N
   ! the exact contact: speed u*, densities either side of it
   us = star_velocity(wl, wr)
   w = riemann_sample(wl, wr, GAMMA, us - 1.e-6_R_P) ; rl = w(1)
   w = riemann_sample(wl, wr, GAMMA, us + 1.e-6_R_P) ; rr = w(1)
   xc = 0.5_R_P + us*tend
   lo = min(rl, rr) + 0.05_R_P*abs(rl - rr) ; hi = max(rl, rr) - 0.05_R_P*abs(rl - rr)
   do i=1, N
      x = (i - 0.5_R_P)*h ; w0(:,i) = merge(wl, wr, x < 0.5_R_P)
   enddo
   do k=1, NS
      call euler_run(ids(k), PATH_SPLIT, w0, N, h, BC_TRANSMISSIVE, tend, 0.5_R_P, q, pos(k))
      l1(k) = 0._R_P ; wid(k) = 0._R_P
      do i=1, N
         x = (i - 0.5_R_P)*h
         w = riemann_sample(wl, wr, GAMMA, (x - 0.5_R_P)/tend)
         l1(k) = l1(k) + abs(q(1,i) - w(1))*h
         if (abs(x - xc) < 0.1_R_P .and. q(1,i) > lo .and. q(1,i) < hi) wid(k) = wid(k) + 1._R_P
      enddo
   enddo
   print '(A,A,A)', '1. ', name, ', N = 100: L1(rho) and cells across the contact (5 %..95 % of the jump)'
   do k=1, NS
      print '(A16,ES12.4,F6.0,L3)', NAMES(k), l1(k), wid(k), pos(k)
      if (vv_report()) write(uc, '(A,1X,A,ES24.15,1X,I0,1X,L1)') name, trim(NAMES(k)), l1(k), nint(wid(k)), pos(k)
   enddo
   if (.not.all(pos)) then ; print '(A)', 'FAIL: not positive' ; passed = .false. ; endif
   if (l1(2) > l1(1)) then ; print '(A)', 'FAIL: TENO5-DR L1 above TENO5''s' ; passed = .false. ; endif
   if (wid(2) >= wid(1)) then ; print '(A)', 'FAIL: TENO5-DR contact not narrower than TENO5''s' ; passed = .false. ; endif
   endsubroutine contact

   function star_velocity(wl, wr) result(us)
   !< u*, the contact's speed: in the star region u = u* while xi varies, so u(xi) - xi vanishes at xi = u* (and
   !< nowhere else for these tubes: the outer states and the fans have u - xi of one sign). The closest of 20001
   !< samples on [-5, 5], refined by the star state's own velocity.
   real(R_P), intent(in) :: wl(3), wr(3)
   real(R_P)             :: us, w(3), xi, best
   integer(I4P)          :: i
   best = huge(1._R_P) ; us = 0._R_P
   do i=0, 20000
      xi = -5._R_P + 10._R_P*i/20000._R_P
      w = riemann_sample(wl, wr, GAMMA, xi)
      if (abs(w(2) - xi) < best) then ; best = abs(w(2) - xi) ; us = w(2) ; endif
   enddo
   endfunction star_velocity

   subroutine e1_tube(name, id, wl, wr, tend)
   !< E1 on a shock tube, N = 200: the two variants against the exact solution and against each other.
   character(*), intent(in) :: name
   integer(I4P), intent(in) :: id
   real(R_P),    intent(in) :: wl(3), wr(3), tend
   integer(I4P), parameter  :: N = 200_I4P
   real(R_P)                :: w0(3,N), qd(3,1-NG:N+NG), qs(3,1-NG:N+NG), h, x, w(3)
   integer(I4P)             :: i
   h = 1._R_P/N
   do i=1, N
      x = (i - 0.5_R_P)*h ; w0(:,i) = merge(wl, wr, x < 0.5_R_P)
   enddo
   stored_masks = .false.
   call euler_run(id, PATH_SPLIT, w0, N, h, BC_TRANSMISSIVE, tend, 0.5_R_P, qd, posd)
   stored_masks = .true.
   call euler_run(id, PATH_SPLIT, w0, N, h, BC_TRANSMISSIVE, tend, 0.5_R_P, qs, poss)
   stored_masks = .false.
   l1s = 0._R_P
   do i=1, N
      x = (i - 0.5_R_P)*h
      w = riemann_sample(wl, wr, GAMMA, (x - 0.5_R_P)/tend)
      l1s(1) = l1s(1) + abs(qd(1,i) - w(1))*h ; l1s(2) = l1s(2) + abs(qs(1,i) - w(1))*h
   enddo
   dq = maxval(abs(qd(1,1:N) - qs(1,1:N)))
   call row(name, id, l1s(1), l1s(2), dq, posd, poss)
   if (.not.(posd .and. poss)) then ; print '(A)', 'FAIL: not positive' ; passed = .false. ; endif
   endsubroutine e1_tube

   subroutine e1_ref(name, id)
   !< E1 on Shu-Osher (N = 200, [-5, 5], t = 1.8) or the blast (N = 400, reflective, limiter on, t = 0.038): the two
   !< variants against each other.
   character(*), intent(in) :: name
   integer(I4P), intent(in) :: id
   integer(I4P)             :: n, i, bc
   real(R_P), allocatable   :: w0(:,:), qd(:,:), qs(:,:)
   real(R_P)                :: xl, xr, h, tend
   if (name == 'blast') then
      n = 400 ; xl = 0._R_P ; xr = 1._R_P ; tend = 0.038_R_P ; bc = BC_REFLECTIVE ; limiter = .true.
   else
      n = 200 ; xl = -5._R_P ; xr = 5._R_P ; tend = 1.8_R_P ; bc = BC_TRANSMISSIVE
   endif
   allocate(w0(3,n), qd(3,1-NG:n+NG), qs(3,1-NG:n+NG))
   h = (xr - xl)/n
   do i=1, n
      if (name == 'blast') then
         w0(:,i) = init_blast(xl + (i - 0.5_R_P)*h)
      else
         w0(:,i) = init_shu_osher(xl + (i - 0.5_R_P)*h)
      endif
   enddo
   stored_masks = .false.
   call euler_run(id, PATH_SPLIT, w0, n, h, bc, tend, 0.5_R_P, qd, posd)
   stored_masks = .true.
   call euler_run(id, PATH_SPLIT, w0, n, h, bc, tend, 0.5_R_P, qs, poss)
   stored_masks = .false. ; limiter = .false.
   dq = maxval(abs(qd(1,1:n) - qs(1,1:n)))
   call row(name, id, -1._R_P, -1._R_P, dq, posd, poss)
   if (.not.(posd .and. poss)) then ; print '(A)', 'FAIL: not positive' ; passed = .false. ; endif
   endsubroutine e1_ref

   subroutine halo()
   !< N4 and N6 on periodic linear advection.
   integer(I4P), parameter :: N = 200_I4P
   real(R_P)               :: u0(N), u1(N), u2(N), u3(N), h, x, tot0, tot, dt
   integer(I4P)            :: k, i, nact, uh
   logical                 :: same6, same5
   h = 1._R_P/N ; dt = 0.4_R_P*h
   do i=1, N
      x = (i - 0.5_R_P)*h
      u0(i) = 0.3_R_P*sin(2._R_P*acos(-1._R_P)*x) + merge(1._R_P, 0._R_P, x > 0.3_R_P .and. x < 0.6_R_P)
   enddo
   if (vv_report()) call vv_open('nonpoly-halo.dat', uh, 'scheme overridden-states |sum-change|/sum G=6-bitwise G=5-differs')
   print '(A)', '3. periodic advection, one period: N4 conservation, N6 two blocks (6 ghost cells bitwise, 5 not)'
   do k=2, NS
      call advect(ids(k), u0, 1, 6, u1, nact)
      call advect(ids(k), u0, 2, 6, u2, i)
      call advect(ids(k), u0, 2, 5, u3, i)
      tot0 = sum(u0) ; tot = sum(u1)
      same6 = all(u1 == u2) ; same5 = all(u1 == u3)
      print '(A16,A,I7,A,ES10.3,A,L2,A,L2)', NAMES(k), ': states overridden ', nact, ', |sum change|/sum ', &
            abs(tot - tot0)/abs(tot0), ', 2 blocks G=6 bitwise', same6, ', G=5 bitwise', same5
      if (vv_report()) write(uh, '(A,1X,I0,ES14.5,2(1X,L1))') trim(NAMES(k)), nact, abs(tot - tot0)/abs(tot0), same6, &
                                                              .not.same5
      if (nact == 0) then ; print '(A)', 'FAIL: no override acted' ; passed = .false. ; endif
      if (abs(tot - tot0) > NP_TOL*abs(tot0)) then ; print '(A)', 'FAIL N4: not conservative' ; passed = .false. ; endif
      if (.not.same6) then ; print '(A)', 'FAIL N6: two blocks with 6 ghost cells differ' ; passed = .false. ; endif
      if (same5) then ; print '(A)', 'FAIL N6: 5 ghost cells do not change the result' ; passed = .false. ; endif
   enddo
   if (vv_report()) close(uh)
   endsubroutine halo

   subroutine advect(id, u0, nb, g, u, nact)
   !< One period of u_t + u_x = 0 (a = 1, the left-biased flux, vp = vm = u: the right state unused), SSP-RK3, in nb
   !< blocks (1 or 2) of equal size with g ghost cells each, refilled from the neighbouring block (periodic) at every
   !< stage; cells of a block's halo beyond g are left unfilled (zero).
   integer(I4P), intent(in)  :: id, nb, g
   real(R_P),    intent(in)  :: u0(:)
   real(R_P),    intent(out) :: u(size(u0))
   integer(I4P), intent(out) :: nact
   real(R_P)                 :: u1(size(u0)), u2(size(u0)), r(size(u0)), h, dt
   integer(I4P)              :: n, step, nsteps
   n = size(u0) ; h = 1._R_P/n ; nsteps = nint(1._R_P/(0.4_R_P*h)) ; dt = 1._R_P/nsteps
   u = u0 ; nact = 0
   do step=1, nsteps
      call adv_rhs(id, nb, g, n, h, u, r, nact) ; u1 = u + dt*r
      call adv_rhs(id, nb, g, n, h, u1, r, nact) ; u2 = 0.75_R_P*u + 0.25_R_P*(u1 + dt*r)
      call adv_rhs(id, nb, g, n, h, u2, r, nact) ; u = u/3._R_P + 2._R_P/3._R_P*(u2 + dt*r)
   enddo
   endsubroutine advect

   function scheme_id(name) result(id)
   !< The registry index of a scheme, the swept ones and the extras (NSCHEMES+1..NALL).
   character(*), intent(in) :: name
   integer(I4P)             :: id, j
   id = 0
   do j=1, NALL
      if (trim(SCHEMES(j)%name) == trim(name)) id = j
   enddo
   if (id == 0) error stop 'seth_vv_nonpoly: scheme not in the registry'
   endfunction scheme_id

   subroutine ablation()
   !< N5: seven schemes on four cases.
   integer(I4P), parameter  :: NA = 7_I4P
   character(16), parameter :: AN(NA) = [character(16) :: 'teno5', 'teno5_st', 'teno5_rbf', 'teno5_thinc', 'teno5_dr', &
                                         'teno5_st_dr', 'teno5_rbf_dr']
   character(9), parameter  :: CN(4) = [character(9) :: 'sod', 'lax', 'shu-osher', 'blast']
   integer(I4P)             :: aid(NA), k, j, ua, nlim(NA), nfb(NA)
   real(R_P)                :: l1(NA)
   logical                  :: pos(NA)
   real(R_P), allocatable   :: qref(:,:)
   do k=1, NA
      aid(k) = scheme_id(AN(k))
   enddo
   if (vv_report()) call vv_open('nonpoly-ablation.dat', ua, 'case scheme L1(rho) limited-faces positive')
   print '(A)', '5. N5, the ablation: L1(rho) (Sod, Lax: exact; Shu-Osher, blast: reference), limited faces'
   do k=1, 4
      nlim = 0
      select case(k)
      case(1)
         call riemann_problem('sod', SOD_L, SOD_R, 0.5_R_P, 0.2_R_P, 0.5_R_P, 200_I4P, aid, PATH_SPLIT, l1, pos, nfb, &
                              files=.false.)
      case(2)
         call riemann_problem('lax', LAX_L, LAX_R, 0.5_R_P, 0.13_R_P, 0.5_R_P, 200_I4P, aid, PATH_SPLIT, l1, pos, nfb, &
                              files=.false.)
      case(3)
         call reference_problem('shu-osher', init_shu_osher, -5._R_P, 5._R_P, BC_TRANSMISSIVE, 1.8_R_P, 0.5_R_P, &
                                200_I4P, 3200_I4P, 3_I4P, aid, PATH_SPLIT, l1, pos, files=.false., qref_in=qref)
         deallocate(qref)
      case default
         limiter = .true.
         call reference_problem('blast', init_blast, 0._R_P, 1._R_P, BC_REFLECTIVE, 0.038_R_P, 0.5_R_P, 400_I4P, &
                                6400_I4P, 3_I4P, aid, PATH_SPLIT, l1, pos, nlim=nlim, files=.false., qref_in=qref)
         limiter = .false.
      endselect
      do j=1, NA
         print '(A12,A16,ES12.4,I6,L3)', trim(CN(k)), AN(j), l1(j), nlim(j), pos(j)
         if (vv_report()) write(ua, '(A,1X,A,ES24.15,1X,I0,1X,L1)') trim(CN(k)), trim(AN(j)), l1(j), nlim(j), pos(j)
         if (.not.pos(j)) then ; print '(A)', 'FAIL N5: not positive' ; passed = .false. ; endif
      enddo
   enddo
   if (vv_report()) close(ua)
   endsubroutine ablation

   subroutine e6()
   !< E6: 100 periods of a sine.
   real(R_P), parameter    :: LAMS(2) = [3._R_P, 6._R_P], CFLS(2) = [0.4_R_P, 0.05_R_P]
   integer(I4P), parameter :: NS6(2) = [20_I4P, 40_I4P]
   real(R_P), allocatable  :: u0(:), u(:)
   real(R_P)               :: h, k2, rate, rate0, pred, s, a0, amax
   integer(I4P)            :: il, ic, in_, n, i, ue, idr, id5
   idr = scheme_id('teno5_rbf') ; id5 = scheme_id('teno5')
   k2 = (2._R_P*acos(-1._R_P))**2
   if (vv_report()) call vv_open('nonpoly-e6.dat', ue, 'lambda N cfl s rate-RBF rate-TENO5 predicted-lambda-term')
   print '(A)', '7. E6, sin(2 pi x), 100 periods: growth rate of the amplitude (RBF, TENO5), predicted lambda term'
   do in_=1, 2
      n = NS6(in_) ; h = 1._R_P/n
      allocate(u0(n), u(n))
      do i=1, n
         u0(i) = (cos(2._R_P*acos(-1._R_P)*(i-1)*h) - cos(2._R_P*acos(-1._R_P)*i*h))/(2._R_P*acos(-1._R_P)*h) ! averages
      enddo
      a0 = maxval(abs(u0))
      do il=1, 2
         s = LAMS(il)*h
         pred = h**5*k2*k2*LAMS(il)**2/2._R_P
         do ic=1, 2
            call long_run(idr, u0, CFLS(ic), LAMS(il)/3._R_P*h, u) ; rate = log(maxval(abs(u))/a0)/100._R_P
            amax = maxval(abs(u))
            call long_run(id5, u0, CFLS(ic), 0._R_P, u) ; rate0 = log(maxval(abs(u))/a0)/100._R_P
            print '(A,F4.1,A,I3,A,F5.2,A,F6.3,3ES12.4)', '  lambda ', LAMS(il), ' N ', n, ' CFL ', CFLS(ic), ' s ', s, &
                  rate, rate0, pred
            if (vv_report()) write(ue, '(F5.1,1X,I0,F6.2,F8.4,3ES16.7)') LAMS(il), n, CFLS(ic), s, rate, rate0, pred
            if (s <= 0.12_R_P) then ! D7: inside the validated range
               if (.not.(amax < 2._R_P*a0)) then ; print '(A)', 'FAIL E6: unbounded' ; passed = .false. ; endif
               if (rate - rate0 > 2._R_P*pred) then
                  print '(A)', 'FAIL E6: growth above twice the prediction inside the validated range' ; passed = .false.
               endif
            endif
         enddo
      enddo
      deallocate(u0, u)
   enddo
   if (vv_report()) close(ue)
   endsubroutine e6

   subroutine long_run(id, u0, cfl, dxr, u)
   !< 100 periods of u_t + u_x = 0, periodic, SSP-RK3 at the given CFL; the RBF kernels' rbf_dx = dxr.
   integer(I4P), intent(in)  :: id
   real(R_P),    intent(in)  :: u0(:), cfl, dxr
   real(R_P),    intent(out) :: u(size(u0))
   real(R_P)                 :: u1(size(u0)), u2(size(u0)), r(size(u0)), h, dt
   integer(I4P)              :: n, step, nsteps, nact
   n = size(u0) ; h = 1._R_P/n ; nsteps = nint(100._R_P/(cfl*h)) ; dt = 100._R_P/nsteps
   u = u0 ; nact = 0
   do step=1, nsteps
      rbf_dx = dxr
      call adv_rhs(id, 1_I4P, 6_I4P, n, h, u, r, nact) ; u1 = u + dt*r
      call adv_rhs(id, 1_I4P, 6_I4P, n, h, u1, r, nact) ; u2 = 0.75_R_P*u + 0.25_R_P*(u1 + dt*r)
      call adv_rhs(id, 1_I4P, 6_I4P, n, h, u2, r, nact) ; u = u/3._R_P + 2._R_P/3._R_P*(u2 + dt*r)
   enddo
   endsubroutine long_run

   subroutine units()
   !< E5: Shu-Osher on [-5, 5] and on [0, 1] (x, t / 10), N = 200.
   integer(I4P), parameter :: N = 200_I4P
   real(R_P)               :: w0(3,N), qa(3,1-NG:N+NG), qb(3,1-NG:N+NG), h, d(2)
   integer(I4P)            :: i, k, uu
   logical                 :: pa, pb
   character(16), parameter :: UN(2) = [character(16) :: 'teno5', 'teno5_rbf']
   if (vv_report()) call vv_open('nonpoly-units.dat', uu, 'scheme max|rho_[-5,5] - rho_scaled| s-on-[-5,5] s-scaled')
   print '(A)', '6. E5, Shu-Osher on [-5, 5] against the same flow scaled by 1/8: max density difference'
   do k=1, 2
      h = 10._R_P/N
      do i=1, N
         w0(:,i) = init_shu_osher(-5._R_P + (i - 0.5_R_P)*h)
      enddo
      call euler_run(scheme_id(UN(k)), PATH_SPLIT, w0, N, h, BC_TRANSMISSIVE, 1.8_R_P, 0.5_R_P, qa, pa)
      call euler_run(scheme_id(UN(k)), PATH_SPLIT, w0, N, h/8._R_P, BC_TRANSMISSIVE, 1.8_R_P/8._R_P, 0.5_R_P, qb, pb)
      d(k) = maxval(abs(qa(1,1:N) - qb(1,1:N)))
      print '(A16,ES12.4,A,F6.3,A,F6.4)', UN(k), d(k), '   s = 3 dx: ', 3._R_P*h, ' against ', 3._R_P*h/8._R_P
      if (vv_report()) write(uu, '(A,3ES16.7)') trim(UN(k)), d(k), 3._R_P*h, 3._R_P*h/8._R_P
      if (.not.(pa .and. pb)) then ; print '(A)', 'FAIL E5: not positive' ; passed = .false. ; endif
   enddo
   if (d(1) /= 0._R_P) then ; print '(A)', 'FAIL E5: TENO5 depends on the length unit' ; passed = .false. ; endif
   if (vv_report()) close(uu)
   endsubroutine units

   subroutine e2_e3()
   !< E2 and E3 on periodic linear advection, N = 100.
   integer(I4P), parameter :: N = 100_I4P, NP = 10_I4P
   real(R_P), parameter    :: WS(5) = [1._R_P, 2._R_P, 3._R_P, 5._R_P, 10._R_P]
   real(R_P)               :: u0(N), u(N), v(N), h, x, wid(NP,NS), l1(5,NS)
   integer(I4P)            :: k, i, p, j, nact, u2, u3
   h = 1._R_P/N
   do i=1, N
      x = (i - 0.5_R_P)*h ; u0(i) = merge(1._R_P, 0._R_P, x > 0.3_R_P .and. x < 0.7_R_P)
   enddo
   do k=1, NS
      u = u0
      do p=1, NP
         call advect(ids(k), u, 1, 6, v, nact) ; u = v
         wid(p,k) = count(u > 0.05_R_P .and. u < 0.95_R_P)
         if (minval(u) < -0.01_R_P .or. maxval(u) > 1.01_R_P) then
            print '(A,A)', 'FAIL E2: out of [-0.01, 1.01]: ', trim(NAMES(k)) ; passed = .false.
         endif
      enddo
   enddo
   print '(A)', '4. E2, square wave, cells across the two jumps after periods 1..10 (TENO5, DR, THINC):'
   do k=1, NS
      print '(A16,10F5.0)', NAMES(k), wid(:,k)
   enddo
   do j=1, 5
      do i=1, N
         x = (i - 0.5_R_P)*h ; u0(i) = cell_avg_tanh(x, h, WS(j))
      enddo
      do k=1, NS
         call advect(ids(k), u0, 1, 6, v, nact)
         l1(j,k) = sum(abs(v - u0))*h
         if (minval(v) < -0.01_R_P .or. maxval(v) > 1.01_R_P) then
            print '(A,A)', 'FAIL E3: out of [-0.01, 1.01]: ', trim(NAMES(k)) ; passed = .false.
         endif
      enddo
   enddo
   print '(A)', '   E3, a smooth plateau, gradients w cells wide, one period: L1 (TENO5, DR, THINC)'
   do j=1, 5
      print '(A,F5.1,3ES12.4)', '   w = ', WS(j), l1(j,:)
   enddo
   if (vv_report()) then
      call vv_open('nonpoly-e2.dat', u2, 'period cells-across-the-jumps: teno5 teno5_dr teno5_thinc')
      do p=1, NP
         write(u2, '(I0,3(1X,I0))') p, nint(wid(p,:))
      enddo
      close(u2)
      call vv_open('nonpoly-e3.dat', u3, 'w L1: teno5 teno5_dr teno5_thinc')
      do j=1, 5
         write(u3, '(F6.1,3ES16.7)') WS(j), l1(j,:)
      enddo
      close(u3)
   endif
   endsubroutine e2_e3

   pure function cell_avg_tanh(x, h, w) result(a)
   !< The average over the cell [x - h/2, x + h/2] of 0.5 (tanh((x - 0.25)/(w h)) - tanh((x - 0.75)/(w h))): a smooth
   !< periodic plateau (its values and slopes at 0 and 1 agree by symmetry), rising at 0.25 and falling at 0.75.
   real(R_P), intent(in) :: x, h, w
   real(R_P)             :: a, s
   s = w*h
   a = 0.5_R_P*s*(lcosh((x + 0.5_R_P*h - 0.25_R_P)/s) - lcosh((x - 0.5_R_P*h - 0.25_R_P)/s) &
                - lcosh((x + 0.5_R_P*h - 0.75_R_P)/s) + lcosh((x - 0.5_R_P*h - 0.75_R_P)/s))/h
   endfunction cell_avg_tanh

   pure function lcosh(y) result(r)
   !< log(cosh(y)), without overflow for large |y|.
   real(R_P), intent(in) :: y
   real(R_P)             :: r
   r = abs(y) + log(0.5_R_P*(1._R_P + exp(-2._R_P*abs(y))))
   endfunction lcosh

   subroutine adv_rhs(id, nb, g, n, h, v, res, nact)
   !< -(F_{i+1/2} - F_{i-1/2})/h block by block, each with g ghost cells from its neighbours (periodic); nact counts
   !< the faces of the first block where the composite's state differs from TENO5's.
   integer(I4P), intent(in)    :: id, nb, g, n
   real(R_P),    intent(in)    :: h, v(n)
   real(R_P),    intent(out)   :: res(n)
   integer(I4P), intent(inout) :: nact
   real(R_P)                   :: blk(1-6:n/nb+6), vp(-5:6), f(0:n/nb), rp, rm, sp, sm, a(-2:3)
   integer(I4P)                :: b, nl, i0, i, m, j
   nl = n/nb
   do b=1, nb
      i0 = (b-1)*nl
      do i=1, nl
         blk(i) = v(i0+i)
      enddo
      do j=1, g ! the halo: g cells from the neighbouring blocks (periodic)
         blk(1-j) = v(modulo(i0 - j, n) + 1)
         blk(nl+j) = v(modulo(i0 + nl + j - 1, n) + 1)
      enddo
      do j=g+1, 6 ! beyond the halo: not filled (zero, as an unfilled halo would leave it)
         blk(1-j) = 0._R_P ; blk(nl+j) = 0._R_P
      enddo
      do i=0, nl ! face i+1/2 of the block
         do m=-5, 6
            vp(m) = blk(i+m)
         enddo
         call reconstruct(id, vp, vp, rp, rm)
         do m=-2, 3
            a(m) = vp(m)
         enddo
         call seth_teno5(a, a, SETH_TENO_EPS, SETH_TENO5_CT, sp, sm)
         if (rp /= sp .and. b == 1 .and. i > 0) nact = nact + 1
         f(i) = rp
      enddo
      do i=1, nl
         res(i0+i) = -(f(i) - f(i-1))/h
      enddo
   enddo
   endsubroutine adv_rhs

   subroutine row(name, id, a, b, d, p1, p2)
   character(*), intent(in) :: name
   integer(I4P), intent(in) :: id
   real(R_P),    intent(in) :: a, b, d
   logical,      intent(in) :: p1, p2
   print '(A12,A14,3ES13.4,2L5)', name, trim(SCHEMES(id)%name), a, b, d, p1, p2
   if (vv_report()) write(ue, '(A,1X,A,3ES24.15,2(1X,L1))') name, trim(SCHEMES(id)%name), a, b, d, p1, p2
   endsubroutine row
endprogram seth_vv_nonpoly
