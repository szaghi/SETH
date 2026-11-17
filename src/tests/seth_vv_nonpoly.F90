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
use seth
use seth_vv_euler
use seth_vv_cases
use seth_vv_riemann, only : riemann_sample
use seth_test_schemes, only : SCHEMES, NSCHEMES
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
   do j=1, NSCHEMES
      if (trim(SCHEMES(j)%name) == trim(NAMES(k))) ids(k) = j
   enddo
   if (ids(k) == 0) error stop 'seth_vv_nonpoly: scheme not in the registry'
enddo
if (vv_report()) call vv_open('nonpoly-contact.dat', uc, 'case scheme L1(rho) contact-cells positive (N = 100)')
if (vv_report()) call vv_open('nonpoly-e1.dat', ue, &
   'case scheme L1-default L1-stored max|rho_default-rho_stored| positive-default positive-stored')

! 1. spec 7 of step 2.2
call contact('sod', SOD_L, SOD_R, 0.2_R_P)
call contact('lax', LAX_L, LAX_R, 0.13_R_P)

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

   subroutine row(name, id, a, b, d, p1, p2)
   character(*), intent(in) :: name
   integer(I4P), intent(in) :: id
   real(R_P),    intent(in) :: a, b, d
   logical,      intent(in) :: p1, p2
   print '(A12,A14,3ES13.4,2L5)', name, trim(SCHEMES(id)%name), a, b, d, p1, p2
   if (vv_report()) write(ue, '(A,1X,A,3ES24.15,2(1X,L1))') name, trim(SCHEMES(id)%name), a, b, d, p1, p2
   endsubroutine row
endprogram seth_vv_nonpoly
