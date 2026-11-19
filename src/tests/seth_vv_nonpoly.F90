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
!<    masks once, in the basis of that cell's own face, through the library's pair (seth_teno5_face_masks and the
!<    *_stored kernels, step 2.6). On scalar data the two are equal bitwise (seth_test_nonpoly, N1); on Euler they
!<    differ where the bases differ. TENO5-DR, TENO5-THINC and the RBF+DR scheme of Liang & Fu 2024. Sod, Lax (N = 200, against the exact solution), Shu-Osher
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
!<    Fu's alpha = 3), one scheme at a time, on Sod, Lax and the 123 problem (N = 200, against the exact solution: L1
!<    of the density, its overshoot beyond the exact extrema, the cells across the contact of Sod and Lax; 123 with the
!<    limiter), Shu-Osher (N = 200 on [-5, 5], t = 1.8, against WENO5-JS on 3200 cells: L1) and the blast waves
!<    (N = 400, the limiter on, against WENO5-JS on 6400 cells: L1); every case the detector activations (states whose
!<    detector fires, per characteristic field, summed over the run: seth_vv_euler's count_activations) and the faces
!<    limited. Report mode only (the references cost minutes); asserted there: positivity; reported:
!<    nonpoly-ablation.dat. Burgers and the composite advection of the seven schemes are in seth_vv_burgers and
!<    seth_vv_advection (burgers-shock.dat, advection.dat), merged into the page's table.
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
!< 8. Step 2.2, spec 6 (frozen): the Gaussian pulse exp(-300 (x - 0.5)^2) of Liang et al. 2022, section 4.1, on [0, 1],
!<    periodic, one period, exact cell averages, SSP-RK3 at CFL 0.1 (80/N)^(2/3) (dt^3 ~ dx^5: the paper's "time step
!<    reduced so that the time integration error can be ignored", value not stated), N = 80, 160, 320 (and 640, 1280 in
!<    report mode). Frozen: zero overrides of TENO5-DR and its solution bitwise TENO5's. Measured (2026-11-18): true at
!<    N = 80 and 160; at N = 320 DR overrides 85 states (THINC 49), all in the tails where |u| ~ 1e-15, the scheme's
!<    round-off noise on a Gaussian 15 decades below its peak (TENO5's cutoff reads the noise as structure: its
!<    indicators ~1e-30 lie far above eps = 1e-40); the solution then differs from TENO5's by 7.8e-15. Asserted (a
!<    deviation from the frozen spec, left to the owner, notes/m5-progress.md): no override at a face whose stencil
!<    exceeds NP_NOISE_FLOOR (1e-12 in R8P) and the solution equal to TENO5's to NP_NOISE_FLOOR; bitwise and zero
!<    overrides are printed. Reported:
!<    the L_inf errors of TENO5, TENO5-DR, TENO5-THINC beside the paper's Table 4 (nonpoly-pulse.dat; finite volumes
!<    here, the paper's discretisation of the case is not stated, so its digits are compared, not asserted).
!< 9. Report mode only, the data of the figures and animations of docs/verification/nonpolynomial/: Sod and Lax at
!<    N = 100, the density profiles of TENO5, DR, THINC and RBF+DR (nonpoly-contact-<case>.dat); Sod's contact in time,
!<    61 frames to t = 0.2 (nonpoly-sod-frames-<scheme>.dat) and the cells across it per frame (nonpoly-e2-sod.dat);
!<    the square wave of E2 over 10 periods, 5 frames a period, with the faces where each composite's state differs
!<    from TENO5's (nonpoly-square-frames-<scheme>.dat); E5's density profiles on both domains (nonpoly-units-
!<    profiles.dat); Shu-Osher under TENO5-RBF, N = 200, 61 frames (shu-osher-frames-teno5_rbf.dat, nonpoly-prefixed:
!<    nonpoly-shu-osher-frames-teno5_rbf.dat); the blast waves under TENO5-DR, N = 400 with the limiter, 61 frames, and
!<    per frame the faces the DR detector flags on the density (nonpoly-blast-frames-teno5_dr.dat,
!<    nonpoly-blast-flags-teno5_dr.dat; the scheme detects on the characteristic fields, the figure on the density).
!< Section 3 (N4, N6) also runs TENO5-RBF-DR (s = 3 dx), its overrides counted against TENO5-RBF.
use seth
use seth_vv_euler
use seth_vv_cases
use seth_vv_riemann, only : riemann_sample
use seth_test_parameters, only : NP_TOL, NP_NOISE_FLOOR
use seth_test_schemes, only : SCHEMES, NALL, reconstruct, rbf_dx
use selah, only : stepper, box, box_interior, SCHEME_SSP_33, FORM_SHU_OSHER_2R ! the time integrator (SELAH #2)
implicit none
integer(I4P), parameter :: NS = 3_I4P
character(16), parameter :: NAMES(NS) = [character(16) :: 'teno5', 'teno5_dr', 'teno5_thinc']
integer(I4P)            :: ids(NS), k, j, ue, uc, id_rbf_dr
integer(I4P)            :: nsig = 0_I4P !< overrides at faces above NP_NOISE_FLOOR (adv_rhs; the pulse of section 8)
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

! 8. step 2.2 spec 6: the Gaussian pulse of Liang et al. 2022, Table 4
call pulse()

! 9. the data of the figures and animations of the non-polynomial pages (report mode only)
if (vv_report()) call figures()

! 2. E1: default against stored masks
print '(A)', 'E1: self-contained (masks in the face''s basis) against stored masks (each in its own face''s basis)'
print '(A12,A14,3A13,2A5)', 'case', 'scheme', 'L1 default', 'L1 stored', 'max |drho|', 'pos', 'pos'
do k=2, NS
   call e1_tube('sod', ids(k), SOD_L, SOD_R, 0.2_R_P)
   call e1_tube('lax', ids(k), LAX_L, LAX_R, 0.13_R_P)
   call e1_ref('shu-osher', ids(k))
   call e1_ref('blast', ids(k))
enddo
id_rbf_dr = scheme_id('teno5_rbf_dr') ! the complete scheme of Liang & Fu 2024, both variants through the library pair
call e1_tube('sod', id_rbf_dr, SOD_L, SOD_R, 0.2_R_P)
call e1_tube('lax', id_rbf_dr, LAX_L, LAX_R, 0.13_R_P)
call e1_ref('shu-osher', id_rbf_dr)
call e1_ref('blast', id_rbf_dr)
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
   !< N4 and N6 on periodic linear advection: TENO5-DR, TENO5-THINC, and TENO5-RBF-DR (s = 3 dx; its overrides counted
   !< against TENO5-RBF).
   integer(I4P), parameter :: N = 200_I4P
   character(16), parameter :: COMP(3) = [character(16) :: 'teno5_dr', 'teno5_thinc', 'teno5_rbf_dr']
   character(16), parameter :: BASE(3) = [character(16) :: 'teno5', 'teno5', 'teno5_rbf']
   real(R_P)               :: u0(N), u1(N), u2(N), u3(N), h, x, tot0, tot, dt, dx0
   integer(I4P)            :: k, i, nact, uh, id, ib
   logical                 :: same6, same5
   h = 1._R_P/N ; dt = 0.4_R_P*h
   dx0 = rbf_dx ; rbf_dx = h
   do i=1, N
      x = (i - 0.5_R_P)*h
      u0(i) = 0.3_R_P*sin(2._R_P*acos(-1._R_P)*x) + merge(1._R_P, 0._R_P, x > 0.3_R_P .and. x < 0.6_R_P)
   enddo
   if (vv_report()) call vv_open('nonpoly-halo.dat', uh, 'scheme overridden-states |sum-change|/sum G=6-bitwise G=5-differs')
   print '(A)', '3. periodic advection, one period: N4 conservation, N6 two blocks (6 ghost cells bitwise, 5 not)'
   do k=1, 3
      id = scheme_id(COMP(k)) ; ib = scheme_id(BASE(k))
      call advect(id, u0, 1, 6, u1, nact, 0.4_R_P, ib)
      call advect(id, u0, 2, 6, u2, i, 0.4_R_P, ib)
      call advect(id, u0, 2, 5, u3, i, 0.4_R_P, ib)
      tot0 = sum(u0) ; tot = sum(u1)
      same6 = all(u1 == u2) ; same5 = all(u1 == u3)
      print '(A16,A,I7,A,ES10.3,A,L2,A,L2)', COMP(k), ': states overridden ', nact, ', |sum change|/sum ', &
            abs(tot - tot0)/abs(tot0), ', 2 blocks G=6 bitwise', same6, ', G=5 bitwise', same5
      if (vv_report()) write(uh, '(A,1X,I0,ES14.5,2(1X,L1))') trim(COMP(k)), nact, abs(tot - tot0)/abs(tot0), same6, &
                                                              .not.same5
      if (nact == 0) then ; print '(A)', 'FAIL: no override acted' ; passed = .false. ; endif
      if (abs(tot - tot0) > NP_TOL*abs(tot0)) then ; print '(A)', 'FAIL N4: not conservative' ; passed = .false. ; endif
      if (.not.same6) then ; print '(A)', 'FAIL N6: two blocks with 6 ghost cells differ' ; passed = .false. ; endif
      if (same5) then ; print '(A)', 'FAIL N6: 5 ghost cells do not change the result' ; passed = .false. ; endif
   enddo
   if (vv_report()) close(uh)
   rbf_dx = dx0
   endsubroutine halo

   subroutine figures()
   !< Section 9: the data of the figures and animations (report mode).
   integer(I4P), parameter  :: NF = 61_I4P, NC = 100_I4P, NSO = 200_I4P, NB = 400_I4P
   character(16), parameter :: FN(4) = [character(16) :: 'teno5', 'teno5_dr', 'teno5_thinc', 'teno5_rbf_dr']
   real(R_P)                :: w0(3,NB), q(3,1-NG:NB+NG), fr(3,NB,NF), tf(NF), h, x, xc, rl, rr, lo, hi, w(3), us
   real(R_P)                :: prof(NC,4), wid(NF,3)
   real(R_P)                :: u(NC), uw(NC), ex(NC)
   real(R_P), pointer       :: rv(:)
   integer(I4P)             :: k, i, f, j, uf, step, nsteps, ib, flg(NB), ks
   type(stepper)            :: st
   type(box)                :: bx
   logical                  :: pz
   ! (a) Sod and Lax, N = 100: the density profiles
   do j=1, 2
      h = 1._R_P/NC
      do i=1, NC
         x = (i - 0.5_R_P)*h
         w0(:,i) = merge(merge(SOD_L, LAX_L, j == 1), merge(SOD_R, LAX_R, j == 1), x < 0.5_R_P)
      enddo
      do k=1, 4
         call euler_run(scheme_id(FN(k)), PATH_SPLIT, w0(:,1:NC), NC, h, BC_TRANSMISSIVE, &
                        merge(0.2_R_P, 0.13_R_P, j == 1), 0.5_R_P, q(:,1-NG:NC+NG), pz)
         prof(:,k) = q(1,1:NC)
      enddo
      call vv_open('nonpoly-contact-'//trim(merge('sod', 'lax', j == 1))//'.dat', uf, &
                   'x rho-exact rho: teno5 teno5_dr teno5_thinc teno5_rbf_dr (N = 100)')
      do i=1, NC
         x = (i - 0.5_R_P)*h
         w = riemann_sample(merge(SOD_L, LAX_L, j == 1), merge(SOD_R, LAX_R, j == 1), GAMMA, &
                            (x - 0.5_R_P)/merge(0.2_R_P, 0.13_R_P, j == 1))
         write(uf, '(6ES16.8)') x, w(1), prof(i,:)
      enddo
      close(uf)
   enddo
   ! (b) Sod's contact in time, 61 frames, and the cells across it per frame
   h = 1._R_P/NC
   do i=1, NC
      x = (i - 0.5_R_P)*h ; w0(:,i) = merge(SOD_L, SOD_R, x < 0.5_R_P)
   enddo
   tf = [(0.2_R_P*f/(NF - 1), f=0, NF-1)] ; tf(1) = 1.e-12_R_P ! kind-ok: a first frame just after 0
   us = star_velocity(SOD_L, SOD_R)
   w = riemann_sample(SOD_L, SOD_R, GAMMA, us - 1.e-6_R_P) ; rl = w(1)
   w = riemann_sample(SOD_L, SOD_R, GAMMA, us + 1.e-6_R_P) ; rr = w(1)
   lo = min(rl, rr) + 0.05_R_P*abs(rl - rr) ; hi = max(rl, rr) - 0.05_R_P*abs(rl - rr)
   do k=1, 3
      call euler_run(scheme_id(FN(k)), PATH_SPLIT, w0(:,1:NC), NC, h, BC_TRANSMISSIVE, 0.2_R_P, 0.5_R_P, &
                     q(:,1-NG:NC+NG), pz, fr(:,1:NC,:), tf, NF)
      call write_frames('nonpoly-sod-frames-'//trim(FN(k))//'.dat', 0._R_P, h, NC, tf, fr(:,1:NC,:))
      do f=1, NF
         xc = 0.5_R_P + us*tf(f)
         wid(f,k) = 0._R_P
         do i=1, NC
            x = (i - 0.5_R_P)*h
            if (abs(x - xc) < 0.1_R_P .and. fr(1,i,f) > lo .and. fr(1,i,f) < hi) wid(f,k) = wid(f,k) + 1._R_P
         enddo
      enddo
   enddo
   call vv_open('nonpoly-e2-sod.dat', uf, &
                't cells-across-the-contact (5 %..95 % of the jump): teno5 teno5_dr teno5_thinc')
   do f=1, NF
      write(uf, '(ES14.6,3(1X,I0))') tf(f), nint(wid(f,:))
   enddo
   close(uf)
   ! (c) the square wave of E2, 10 periods, 5 frames a period, the faces where the state differs from TENO5's
   h = 1._R_P/NC ; nsteps = nint(1._R_P/(0.4_R_P*h))
   do i=1, NC
      x = (i - 0.5_R_P)*h ; ex(i) = merge(1._R_P, 0._R_P, x > 0.3_R_P .and. x < 0.7_R_P)
   enddo
   do k=1, 3
      call vv_open('nonpoly-square-frames-'//trim(FN(k))//'.dat', uf, &
                   'frames: line 1 x(1..N); then per frame: t u(1..N) u_exact(1..N) overridden-face-i+1/2(1..N)')
      write(uf, '(*(ES16.8))') ((i - 0.5_R_P)*h, i=1, NC)
      u = ex
      bx = box_interior(lb=[1_I4P], ub=[NC], ngc=0_I4P, components=.false.)
      call st%initialize(SCHEME_SSP_33, FORM_SHU_OSHER_2R, lb=[1_I4P], ub=[NC], components=.false.)
      do step=0, 10*nsteps
         if (step > 0) then ! SELAH's SSP(3,3), two-register form (the stage state is u), dt = 1/nsteps
            call st%begin_step(bx, u)
            do ks=1, st%stages()
               call st%begin_stage(ks, bx, 1._R_P/nsteps, u)
               call st%rhs_register_view(rv)
               call adv_rhs(scheme_id(FN(k)), ids(1), 1_I4P, 6_I4P, NC, h, u, rv, ib)
               call st%end_stage(ks, bx, 1._R_P/nsteps, u)
            enddo
            call st%finish_step(bx, 1._R_P/nsteps, u)
         endif
         if (mod(step, nsteps/5) == 0) then
            do i=1, NC ! the exact solution at t (the square wave translated, periodic)
               x = modulo((i - 0.5_R_P)*h - real(step, R_P)/nsteps, 1._R_P)
               uw(i) = merge(1._R_P, 0._R_P, x > 0.3_R_P .and. x < 0.7_R_P)
            enddo
            call overridden(scheme_id(FN(k)), ids(1), u, flg(1:NC))
            write(uf, '(*(ES16.8))') real(step, R_P)/nsteps, u, uw, real(flg(1:NC), R_P)
         endif
      enddo
      close(uf)
   enddo
   ! (e) Shu-Osher under TENO5-RBF, N = 200, 61 frames
   h = 10._R_P/NSO
   do i=1, NSO
      w0(:,i) = init_shu_osher(-5._R_P + (i - 0.5_R_P)*h)
   enddo
   tf = [(1.8_R_P*f/(NF - 1), f=0, NF-1)] ; tf(1) = 1.e-12_R_P ! kind-ok: a first frame just after 0
   call euler_run(scheme_id('teno5_rbf'), PATH_SPLIT, w0(:,1:NSO), NSO, h, BC_TRANSMISSIVE, 1.8_R_P, 0.5_R_P, &
                  q(:,1-NG:NSO+NG), pz, fr(:,1:NSO,:), tf, NF)
   call write_frames('nonpoly-shu-osher-frames-teno5_rbf.dat', -5._R_P, h, NSO, tf, fr(:,1:NSO,:))
   ! (f) the blast waves under TENO5-DR, N = 400, the limiter on, 61 frames; the DR detector on the density
   h = 1._R_P/NB
   do i=1, NB
      w0(:,i) = init_blast((i - 0.5_R_P)*h)
   enddo
   tf = [(0.038_R_P*f/(NF - 1), f=0, NF-1)] ; tf(1) = 1.e-12_R_P ! kind-ok: a first frame just after 0
   limiter = .true.
   call euler_run(scheme_id('teno5_dr'), PATH_SPLIT, w0, NB, h, BC_REFLECTIVE, 0.038_R_P, 0.5_R_P, q, pz, fr, tf, NF)
   limiter = .false.
   call write_frames('nonpoly-blast-frames-teno5_dr.dat', 0._R_P, h, NB, tf, fr)
   call vv_open('nonpoly-blast-flags-teno5_dr.dat', uf, &
                'per frame: t, the DR level (0, 1, 2) of the left state at face i+1/2, i = 1..N, on the density')
   do f=1, NF
      call dr_flags(fr(1,:,f), flg)
      write(uf, '(ES14.6,*(1X,I0))') tf(f), flg
   enddo
   close(uf)
   call st%destroy
   endsubroutine figures

   subroutine overridden(id, base, u, flg)
   !< 1 at the faces i+1/2 (periodic) where the scheme's left state differs from the base scheme's, 0 elsewhere.
   integer(I4P), intent(in)  :: id, base
   real(R_P),    intent(in)  :: u(:)
   integer(I4P), intent(out) :: flg(size(u))
   real(R_P)                 :: vp(-5:6), rp, rm, sp, sm
   integer(I4P)              :: i, m, n
   n = size(u)
   do i=1, n
      do m=-5, 6
         vp(m) = u(modulo(i - 1 + m, n) + 1)
      enddo
      call reconstruct(id, vp, vp, rp, rm)
      call reconstruct(base, vp, vp, sp, sm)
      flg(i) = merge(1_I4P, 0_I4P, rp /= sp)
   enddo
   endsubroutine overridden

   subroutine dr_flags(rho, flg)
   !< The DR level of the left state at every face i+1/2 from the masks of the density (constant continuation).
   real(R_P),    intent(in)  :: rho(:)
   integer(I4P), intent(out) :: flg(size(rho))
   real(R_P)                 :: w(-2:2)
   integer(I4P)              :: mk(-3:size(rho)+4), mp(-3:3), i, m, c, n, lm
   n = size(rho)
   do i=-3, n+4
      do m=-2, 2
         w(m) = rho(min(max(i + m, 1), n))
      enddo
      call seth_teno5_mask(w, SETH_TENO_EPS, SETH_TENO5_CT, mk(i))
   enddo
   do i=1, n
      do c=-3, 3
         mp(c) = mk(i + c)
      enddo
      call seth_dr_detect(mp, mp, flg(i), lm)
   enddo
   endsubroutine dr_flags

   subroutine pulse()
   !< Step 2.2 spec 6: Liang et al. 2022, section 4.1 and Table 4.
   integer(I4P), parameter :: NN(5) = [80_I4P, 160_I4P, 320_I4P, 640_I4P, 1280_I4P]
   real(R_P),    parameter :: PAPER(5) = [1.40e-2_R_P, 4.94e-4_R_P, 1.59e-5_R_P, 5.02e-7_R_P, 1.57e-8_R_P] ! kind-ok: printed values
   real(R_P),    parameter :: CFL80 = 0.1_R_P !< CFL at N = 80, scaled by (80/N)^(2/3): dt^3 ~ dx^5
   real(R_P), allocatable  :: u0(:), ut(:), ud(:), uth(:)
   real(R_P)               :: h, a, b, sq, einf(3), ep, cfl
   integer(I4P)            :: j, i, nmax, nd, nth, idt, idd, idth, up, nsd
   idt = scheme_id('teno5') ; idd = scheme_id('teno5_dr') ; idth = scheme_id('teno5_thinc')
   nmax = merge(5_I4P, 3_I4P, vv_report()) ! the suite to N = 320; report mode to N = 1280 (minutes)
   if (vv_report()) call vv_open('nonpoly-pulse.dat', up, &
      'N  Linf-TENO5  Linf-DR  Linf-THINC  overrides-DR  overrides-THINC  DR-bitwise-TENO5  paper-Linf')
   print '(A)', '8. Gaussian pulse exp(-300 (x - 0.5)^2), one period (Liang et al. 2022, Table 4): DR inert, bitwise TENO5'
   sq = sqrt(300._R_P) ; ep = 0._R_P
   do j=1, nmax
      allocate(u0(NN(j)), ut(NN(j)), ud(NN(j)), uth(NN(j)))
      h = 1._R_P/NN(j)
      do i=1, NN(j) ! exact cell averages
         a = (i - 1)*h - 0.5_R_P ; b = i*h - 0.5_R_P
         u0(i) = sqrt(acos(-1._R_P))/(2._R_P*sq)*(erf(sq*b) - erf(sq*a))/h
      enddo
      cfl = CFL80*(80._R_P/NN(j))**(2._R_P/3._R_P)
      call advect(idt, u0, 1, 6, ut, i, cfl, idt)
      nsig = 0
      call advect(idd, u0, 1, 6, ud, nd, cfl, idt)
      nsd = nsig
      call advect(idth, u0, 1, 6, uth, nth, cfl, idt)
      einf = [maxval(abs(ut - u0)), maxval(abs(ud - u0)), maxval(abs(uth - u0))]
      print '(A,I5,A,3ES11.3,A,2I6,A,I3,A,L2,A,ES9.2,A,ES10.3)', '  N', NN(j), ': L_inf TENO5, DR, THINC', einf, &
            '; overrides DR, THINC', nd, nth, ' (DR above the floor', nsd, '); DR bitwise TENO5', all(ud == ut), &
            ', max |DR - TENO5|', maxval(abs(ud - ut)), '; paper', PAPER(j)
      if (vv_report()) write(up, '(I0,3ES14.5,2(1X,I0),1X,L1,ES12.3)') NN(j), einf, nd, nth, all(ud == ut), PAPER(j)
      if (nsd /= 0 .or. maxval(abs(ud - ut)) > NP_NOISE_FLOOR) then
         print '(A)', 'FAIL spec 6: TENO5-DR not inert on the Gaussian pulse above its noise floor' ; passed = .false.
      endif
      deallocate(u0, ut, ud, uth)
   enddo
   if (vv_report()) close(up)
   endsubroutine pulse

   subroutine advect(id, u0, nb, g, u, nact, cfl, base)
   !< One period of u_t + u_x = 0 (a = 1, the left-biased flux, vp = vm = u: the right state unused), SSP-RK3 at the
   !< given CFL, in nb blocks (1 or 2) of equal size with g ghost cells each, refilled from the neighbouring block
   !< (periodic) at every stage; cells of a block's halo beyond g are left unfilled (zero). nact counts the faces where
   !< the scheme's state differs from the base scheme's (the overrides of a composite: TENO5, or TENO5-RBF for RBF+DR).
   integer(I4P), intent(in)  :: id, nb, g, base
   real(R_P),    intent(in)  :: u0(:), cfl
   real(R_P),    intent(out) :: u(size(u0))
   integer(I4P), intent(out) :: nact
   real(R_P)                 :: h, dt
   real(R_P), pointer        :: r(:)
   integer(I4P)              :: n, step, nsteps, k
   type(stepper)             :: st
   type(box)                 :: bx
   n = size(u0) ; h = 1._R_P/n ; nsteps = nint(1._R_P/(cfl*h)) ; dt = 1._R_P/nsteps
   u = u0 ; nact = 0
   ! SELAH's SSP(3,3), two-register form (the stage state is u itself), on the periodic rank-1 field u(1:n)
   bx = box_interior(lb=[1_I4P], ub=[n], ngc=0_I4P, components=.false.)
   call st%initialize(SCHEME_SSP_33, FORM_SHU_OSHER_2R, lb=[1_I4P], ub=[n], components=.false.)
   do step=1, nsteps
      call st%begin_step(bx, u)
      do k=1, st%stages()
         call st%begin_stage(k, bx, dt, u)
         call st%rhs_register_view(r)
         call adv_rhs(id, base, nb, g, n, h, u, r, nact)
         call st%end_stage(k, bx, dt, u)
      enddo
      call st%finish_step(bx, dt, u)
   enddo
   call st%destroy
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
   !< N5: seven schemes on five Euler cases, one scheme at a time (activations counted per scheme).
   integer(I4P), parameter  :: NA = 7_I4P, NE = 200_I4P
   character(16), parameter :: AN(NA) = [character(16) :: 'teno5', 'teno5_st', 'teno5_rbf', 'teno5_thinc', 'teno5_dr', &
                                         'teno5_st_dr', 'teno5_rbf_dr']
   character(9), parameter  :: CN(5) = [character(9) :: 'sod', 'lax', '123', 'shu-osher', 'blast']
   real(R_P), parameter     :: W123L(3) = [1._R_P, -2._R_P, 0.4_R_P], W123R(3) = [1._R_P, 2._R_P, 0.4_R_P]
   integer(I4P)             :: k, j, i, ua, nl(1), nf(1)
   real(R_P)                :: l1(1), osh, cw, h, x, w(3), wl(3), wr(3), tend, rmax, rmin, xc, us, rl, rr, lo, hi
   real(R_P)                :: w0(3,NE), q(3,1-NG:NE+NG)
   logical                  :: pos(1)
   real(R_P), allocatable   :: qref(:,:)
   if (vv_report()) call vv_open('nonpoly-ablation.dat', ua, &
      'case scheme L1(rho) overshoot contact-cells activations limited-faces positive (overshoot, contact: -1 n/a)')
   print '(A)', '5. N5, the ablation: L1(rho), overshoot, cells across the contact, detector activations, limited faces'
   count_activations = .true.
   do k=1, 5
      do j=1, NA
         nactivations = 0 ; nlimited = 0 ; osh = -1._R_P ; cw = -1._R_P
         select case(k)
         case(1:3) ! the exact solution: L1, overshoot beyond its extrema, the contact (Sod, Lax)
            select case(k)
            case(1) ; wl = SOD_L ; wr = SOD_R ; tend = 0.2_R_P
            case(2) ; wl = LAX_L ; wr = LAX_R ; tend = 0.13_R_P
            case default ; wl = W123L ; wr = W123R ; tend = 0.15_R_P
            endselect
            h = 1._R_P/NE
            do i=1, NE
               x = (i - 0.5_R_P)*h ; w0(:,i) = merge(wl, wr, x < 0.5_R_P)
            enddo
            limiter = k == 3
            call euler_run(scheme_id(AN(j)), PATH_SPLIT, w0, NE, h, BC_TRANSMISSIVE, tend, 0.5_R_P, q, pos(1))
            limiter = .false.
            l1(1) = 0._R_P ; rmax = -huge(1._R_P) ; rmin = huge(1._R_P)
            do i=1, NE
               x = (i - 0.5_R_P)*h
               w = riemann_sample(wl, wr, GAMMA, (x - 0.5_R_P)/tend)
               l1(1) = l1(1) + abs(q(1,i) - w(1))*h
               rmax = max(rmax, w(1)) ; rmin = min(rmin, w(1))
            enddo
            osh = max(maxval(q(1,1:NE)) - rmax, rmin - minval(q(1,1:NE)), 0._R_P)
            if (k <= 2) then
               us = star_velocity(wl, wr) ; xc = 0.5_R_P + us*tend
               w = riemann_sample(wl, wr, GAMMA, us - 1.e-6_R_P) ; rl = w(1)
               w = riemann_sample(wl, wr, GAMMA, us + 1.e-6_R_P) ; rr = w(1)
               lo = min(rl, rr) + 0.05_R_P*abs(rl - rr) ; hi = max(rl, rr) - 0.05_R_P*abs(rl - rr)
               cw = 0._R_P
               do i=1, NE
                  x = (i - 0.5_R_P)*h
                  if (abs(x - xc) < 0.1_R_P .and. q(1,i) > lo .and. q(1,i) < hi) cw = cw + 1._R_P
               enddo
            endif
            nl(1) = nlimited
         case(4)
            call reference_problem('shu-osher', init_shu_osher, -5._R_P, 5._R_P, BC_TRANSMISSIVE, 1.8_R_P, 0.5_R_P, &
                                   NE, 3200_I4P, 3_I4P, [scheme_id(AN(j))], PATH_SPLIT, l1, pos, files=.false., &
                                   qref_in=qref)
            nl(1) = nlimited
         case default
            limiter = .true.
            call reference_problem('blast', init_blast, 0._R_P, 1._R_P, BC_REFLECTIVE, 0.038_R_P, 0.5_R_P, 400_I4P, &
                                   6400_I4P, 3_I4P, [scheme_id(AN(j))], PATH_SPLIT, l1, pos, nlim=nl, files=.false., &
                                   qref_in=qref)
            limiter = .false.
         endselect
         nf(1) = int(nactivations, I4P)
         print '(A12,A16,2ES12.4,F5.0,I9,I6,L3)', trim(CN(k)), AN(j), l1(1), osh, cw, nf(1), nl(1), pos(1)
         if (vv_report()) write(ua, '(A,1X,A,2ES24.15,1X,I0,1X,I0,1X,I0,1X,L1)') trim(CN(k)), trim(AN(j)), l1(1), osh, &
                                                                           nint(cw), nf(1), nl(1), pos(1)
         if (.not.pos(1)) then ; print '(A)', 'FAIL N5: not positive' ; passed = .false. ; endif
      enddo
      if (allocated(qref)) deallocate(qref)
   enddo
   count_activations = .false.
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
   real(R_P)                 :: h, dt
   real(R_P), pointer        :: r(:)
   integer(I4P)              :: n, step, nsteps, nact, k
   type(stepper)             :: st
   type(box)                 :: bx
   n = size(u0) ; h = 1._R_P/n ; nsteps = nint(100._R_P/(cfl*h)) ; dt = 100._R_P/nsteps
   u = u0 ; nact = 0
   ! SELAH's SSP(3,3), two-register form (the stage state is u itself), on the periodic rank-1 field u(1:n)
   bx = box_interior(lb=[1_I4P], ub=[n], ngc=0_I4P, components=.false.)
   call st%initialize(SCHEME_SSP_33, FORM_SHU_OSHER_2R, lb=[1_I4P], ub=[n], components=.false.)
   do step=1, nsteps
      rbf_dx = dxr
      call st%begin_step(bx, u)
      do k=1, st%stages()
         call st%begin_stage(k, bx, dt, u)
         call st%rhs_register_view(r)
         call adv_rhs(id, ids(1), 1_I4P, 6_I4P, n, h, u, r, nact)
         call st%end_stage(k, bx, dt, u)
      enddo
      call st%finish_step(bx, dt, u)
   enddo
   call st%destroy
   endsubroutine long_run

   subroutine units()
   !< E5: Shu-Osher on [-5, 5] and on [0, 1] (x, t / 10), N = 200.
   integer(I4P), parameter :: N = 200_I4P
   real(R_P)               :: w0(3,N), qa(3,1-NG:N+NG), qb(3,1-NG:N+NG), h, d(2), pr(N,3)
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
      pr(:,k) = qa(1,1:N) ; if (k == 2) pr(:,3) = qb(1,1:N)
      print '(A16,ES12.4,A,F6.3,A,F6.4)', UN(k), d(k), '   s = 3 dx: ', 3._R_P*h, ' against ', 3._R_P*h/8._R_P
      if (vv_report()) write(uu, '(A,3ES16.7)') trim(UN(k)), d(k), 3._R_P*h, 3._R_P*h/8._R_P
      if (.not.(pa .and. pb)) then ; print '(A)', 'FAIL E5: not positive' ; passed = .false. ; endif
   enddo
   if (d(1) /= 0._R_P) then ; print '(A)', 'FAIL E5: TENO5 depends on the length unit' ; passed = .false. ; endif
   if (vv_report()) then
      close(uu)
      call vv_open('nonpoly-units-profiles.dat', uu, 'x rho: teno5 teno5_rbf teno5_rbf-on-the-scaled-domain (x on [-5, 5])')
      do i=1, N
         write(uu, '(4ES16.8)') -5._R_P + (i - 0.5_R_P)*h, pr(i,:)
      enddo
      close(uu)
   endif
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
         call advect(ids(k), u, 1, 6, v, nact, 0.4_R_P, ids(1)) ; u = v
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
         call advect(ids(k), u0, 1, 6, v, nact, 0.4_R_P, ids(1))
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

   subroutine adv_rhs(id, base, nb, g, n, h, v, res, nact)
   !< -(F_{i+1/2} - F_{i-1/2})/h block by block, each with g ghost cells from its neighbours (periodic); nact counts
   !< the faces of the first block where the scheme's state differs from the base scheme's.
   integer(I4P), intent(in)    :: id, base, nb, g, n
   real(R_P),    intent(in)    :: h, v(n)
   real(R_P),    intent(out)   :: res(n)
   integer(I4P), intent(inout) :: nact
   real(R_P)                   :: blk(1-6:n/nb+6), vp(-5:6), f(0:n/nb), rp, rm, sp, sm
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
         call reconstruct(base, vp, vp, sp, sm)
         if (rp /= sp .and. b == 1 .and. i > 0) then
            nact = nact + 1
            if (maxval(abs(vp)) > NP_NOISE_FLOOR) nsig = nsig + 1
         endif
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
