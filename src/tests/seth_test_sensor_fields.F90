!< SETH test: the stencil sensors against the Ducros sensors on analytic two-dimensional fields.
program seth_test_sensor_fields
!< SETH test: the stencil sensors against the Ducros sensors on analytic two-dimensional fields.
!<
!< Roadmap issue #4, Phase 5, the SETH half of item 5.2 (step 4.4: sensors on analytic fields here, turbulent flows in
!< ADAM). N x N cells of width h on [-1, 1]^2; the fields are sampled at the cell centres (no flow is computed).
!<
!< - Stencil sensors (seth_sensor_teno, _ren, _lcr, _tau) on the density, along x and along y, as the 1-D harness:
!<   cells i-2..i+3 around the face, normalised by max(|rho|, 1), the published thresholds. A cell is flagged when a
!<   face of it is (s < 1).
!< - Ducros sensors on the velocity gradient by second-order central differences, ADAM's policy (FLUME's defaults):
!<   eps = 1e-8 (c/h)^2 with c = 1, flagged when phi > 0.2 (seth_ducros) or psi > 0.2 (seth_ducros_compression).
!<   Without ADAM's density-jump override: that is the caller's policy, and the contact case shows why it exists.
!<
!< Fields (u, v, rho; the pressure does not enter any sensor):
!< 1. Taylor-Green vortex, u = sin(pi x) cos(pi y), v = -cos(pi x) sin(pi y), rho = 1 + (cos 2 pi x + cos 2 pi y)/4,
!<    N = 8 .. 128: smooth, theta = 0. Flagged fraction per N. Asserted: Ducros flags nothing at any N (the discrete
!<    divergence vanishes to round-off); Li, Chen & Ren's flags nothing from N = 32 (8 cells per wavelength of rho;
!<    measured: nothing from N = 8). Ren's flags the extrema of rho at every N.
!< 2. A normal shock at x = 0, Mach 2, gamma = 1.4 (Rankine-Hugoniot states, u in the shock frame), N = 64. Asserted:
!<    every sensor flags every cell next to the shock, and nothing more than 3 cells from it.
!< 3. The shock of 2 with a Gaussian vortex at (0, 0) added to the velocity, u_theta = G r exp(-r^2/r0^2), r0 = 0.25:
!<    its vorticity at the centre, 2 G, against the shock's discrete divergence, |theta| = |u_2 - u_1|/(2 h). For
!<    omega_0/|theta| = 0 .. 8: the fraction of the shock's cells within r0 of the centre that each sensor flags.
!<    Asserted: the stencil sensors flag all of them at every strength (the density does not see the vortex); Ducros
!<    misses some once omega > 2 |theta| (phi > 0.2 iff omega^2 < 4 theta^2), i.e. at omega_0/|theta| >= 4.
!< 4. A steep expansion, u = 0.5 tanh(x/h), rho = 1 + 0.2 tanh(x/h), N = 64. Asserted: seth_ducros flags it,
!<    seth_ducros_compression does not.
!< 5. A contact, rho = 1 for x < 0 and 2 beyond, u = v = 0 (uniform velocity), N = 64. Asserted: every stencil sensor
!<    flags it; neither Ducros sensor flags anything.
!< Report mode: sensor-fields.dat (case parameter sensor flagged-fraction) for the documentation's table.
use seth
use seth_test_parameters, only : EPS_FIELD_SENSOR
use seth_vv_report, only : vv_report, vv_open
implicit none
real(R_P),    parameter :: PI = acos(-1._R_P), GAM = 1.4_R_P, THR = 0.2_R_P
integer(I4P), parameter :: NS = 6
character(8), parameter :: NAMES(NS) = ['teno    ', 'ren     ', 'lcr     ', 'tau     ', 'ducros  ', 'ducros-c']
real(R_P), allocatable  :: u(:,:), v(:,:), rho(:,:)
logical, allocatable    :: fl(:,:,:)
real(R_P)               :: frac(NS), u1, u2, r1, r2, ratio, theta, om0, gv
integer(I4P)            :: n, k, i, j, m, nsh, ndet(NS)
integer                 :: ur
logical                 :: passed
passed = .true.
if (vv_report()) call vv_open('sensor-fields.dat', ur, 'case parameter sensor flagged-fraction')

! 1. Taylor-Green
print '(A)', '1. Taylor-Green vortex: fraction of cells flagged (teno ren lcr tau ducros ducros-c)'
do k=0, 4
   n = 8*2**k
   call grid(n)
   do j=1, n
      do i=1, n
         u(i,j) = sin(PI*xc(i))*cos(PI*xc(j)) ; v(i,j) = -cos(PI*xc(i))*sin(PI*xc(j))
         rho(i,j) = 1._R_P + 0.25_R_P*(cos(2._R_P*PI*xc(i)) + cos(2._R_P*PI*xc(j)))
      enddo
   enddo
   call flags(frac)
   print '(A,I5,6F8.4)', '   N ', n, frac
   call report('taylor-green', real(n, R_P), frac)
   if (any(frac(5:6) > 0._R_P)) call fail('Taylor-Green: a Ducros sensor flagged a divergence-free field')
   if (n >= 32 .and. frac(3) > 0._R_P) call fail('Taylor-Green: Li, Chen & Ren flagged a resolved vortex')
enddo

! 2. normal shock, Mach 2 (shock frame)
n = 64 ; call grid(n)
call shock_states(2._R_P, r1, u1, r2, u2)
call set_shock(0._R_P)
call flags(frac)
print '(A)', '2. normal shock, Mach 2, N = 64: fraction of the cells next to the shock flagged; cells flagged farther'
call shock_detection(0._R_P, ndet, nsh)
print '(A,6F8.4)', '   next to it ', real(ndet)/nsh
call report('shock', 2._R_P, real(ndet, R_P)/nsh)
do m=1, NS
   if (ndet(m) < nsh) call fail('shock: '//trim(NAMES(m))//' missed a shock cell')
   if (far(m) > 0) call fail('shock: '//trim(NAMES(m))//' flagged a cell more than 3 h from the shock')
enddo

! 3. the shock through a vortex
theta = abs(u2 - u1)/(2._R_P*2._R_P/n)
print '(A)', '3. shock through a vortex, N = 64: fraction of the shock cells within r0 flagged, against omega_0/|theta|'
do k=0, 4
   ratio = merge(0._R_P, 0.5_R_P*2._R_P**k, k == 0)
   om0 = ratio*theta ; gv = 0.5_R_P*om0
   call set_shock(gv)
   call flags(frac)
   call shock_detection(0.25_R_P, ndet, nsh)
   print '(A,F5.1,6F8.4)', '   omega_0/|theta| ', ratio, real(ndet)/nsh
   call report('shock-vortex', ratio, real(ndet, R_P)/nsh)
   if (any(ndet(1:4) < nsh)) call fail('shock-vortex: a stencil sensor missed a shock cell')
   if (ratio <= 1._R_P .and. ndet(5) < nsh) call fail('shock-vortex: Ducros missed a shock cell in a weak vortex')
   if (ratio >= 4._R_P .and. ndet(5) == nsh) call fail('shock-vortex: Ducros saw through a strong vortex (prediction)')
enddo

! 4. steep expansion
call grid(64)
do j=1, n
   do i=1, n
      u(i,j) = 0.5_R_P*tanh(xc(i)/(2._R_P/n)) ; v(i,j) = 0._R_P ; rho(i,j) = 1._R_P + 0.2_R_P*tanh(xc(i)/(2._R_P/n))
   enddo
enddo
call flags(frac)
print '(A,6F8.4)', '4. steep expansion, N = 64: fraction flagged ', frac
call report('expansion', 0._R_P, frac)
if (.not.(frac(5) > 0._R_P)) call fail('expansion: seth_ducros did not flag it')
if (frac(6) > 0._R_P) call fail('expansion: seth_ducros_compression flagged it')

! 5. contact
do j=1, n
   do i=1, n
      u(i,j) = 0._R_P ; v(i,j) = 0._R_P ; rho(i,j) = merge(1._R_P, 2._R_P, xc(i) < 0._R_P)
   enddo
enddo
call flags(frac)
print '(A,6F8.4)', '5. contact, N = 64: fraction flagged ', frac
call report('contact', 0._R_P, frac)
if (any(.not.(frac(1:4) > 0._R_P))) call fail('contact: a stencil sensor did not flag it')
if (any(frac(5:6) > 0._R_P)) call fail('contact: a Ducros sensor flagged a contact')

if (vv_report()) close(ur)
if (.not.passed) error stop 'sensor field test failed'
print '(A)', 'seth_test_sensor_fields: ok'

contains
   pure function xc(i) result(x)
   !< Centre of cell i.
   integer(I4P), intent(in) :: i
   real(R_P)                :: x
   x = -1._R_P + (i - 0.5_R_P)*2._R_P/n
   endfunction xc

   subroutine grid(nn)
   integer(I4P), intent(in) :: nn
   n = nn
   if (allocated(u)) deallocate(u, v, rho, fl)
   allocate(u(n,n), v(n,n), rho(n,n), fl(n,n,NS))
   endsubroutine grid

   subroutine shock_states(mach, ra, ua, rb, ub)
   !< Rankine-Hugoniot: upstream rho = 1, p = 1/gamma (c = 1), u = mach; downstream from the normal-shock relations.
   real(R_P), intent(in)  :: mach
   real(R_P), intent(out) :: ra, ua, rb, ub
   ra = 1._R_P ; ua = mach
   rb = ra*(GAM + 1._R_P)*mach**2/((GAM - 1._R_P)*mach**2 + 2._R_P)
   ub = ua*ra/rb
   endsubroutine shock_states

   subroutine set_shock(g)
   !< The shock at x = 0 (cell faces), plus the Gaussian vortex of strength g (u_theta = g r exp(-r^2/r0^2)).
   real(R_P), intent(in) :: g
   real(R_P), parameter  :: R0 = 0.25_R_P
   real(R_P)             :: e
   integer(I4P)          :: ii, jj
   do jj=1, n
      do ii=1, n
         e = g*exp(-(xc(ii)**2 + xc(jj)**2)/R0**2)
         u(ii,jj) = merge(u1, u2, xc(ii) < 0._R_P) - e*xc(jj)
         v(ii,jj) = e*xc(ii)
         rho(ii,jj) = merge(r1, r2, xc(ii) < 0._R_P)
      enddo
   enddo
   endsubroutine set_shock

   subroutine flags(fr)
   !< Flag every interior cell (3 cells from the edges) by each sensor; fr: the flagged fraction of them.
   real(R_P), intent(out) :: fr(NS)
   real(R_P)              :: st(-2:3), s, g(3,3), hh, eps, phi
   integer(I4P)           :: ii, jj, is, d, c
   hh = 2._R_P/n ; eps = EPS_FIELD_SENSOR/hh**2
   fl = .false.
   do jj=4, n-3
      do ii=4, n-3
         do is=1, 4
            do d=1, 2 ! faces of the cell along x (d = 1) and y (d = 2): ii-1/2 and ii+1/2
               do c=-1, 0
                  if (d == 1) then
                     st = rho(ii+c-2:ii+c+3,jj)
                  else
                     st = rho(ii,jj+c-2:jj+c+3)
                  endif
                  call stencil_sensor(is, st, s)
                  if (s < 1._R_P) fl(ii,jj,is) = .true.
               enddo
            enddo
         enddo
         g = 0._R_P
         g(1,1) = (u(ii+1,jj) - u(ii-1,jj))/(2._R_P*hh) ; g(1,2) = (u(ii,jj+1) - u(ii,jj-1))/(2._R_P*hh)
         g(2,1) = (v(ii+1,jj) - v(ii-1,jj))/(2._R_P*hh) ; g(2,2) = (v(ii,jj+1) - v(ii,jj-1))/(2._R_P*hh)
         call seth_ducros(eps, g, phi) ; fl(ii,jj,5) = phi > THR
         call seth_ducros_compression(eps, g, phi) ; fl(ii,jj,6) = phi > THR
      enddo
   enddo
   do is=1, NS
      fr(is) = count(fl(4:n-3,4:n-3,is))/real((n-6)**2, R_P)
   enddo
   endsubroutine flags

   subroutine shock_detection(rad, nd, nt)
   !< The cells next to the shock (the two columns around x = 0), within rad of the centre (rad = 0: all of them):
   !< nt of them, nd(is) flagged by sensor is.
   real(R_P),    intent(in)  :: rad
   integer(I4P), intent(out) :: nd(NS), nt
   integer(I4P)              :: ii, jj
   nd = 0 ; nt = 0
   do jj=4, n-3
      do ii=n/2, n/2+1
         if (rad > 0._R_P .and. xc(ii)**2 + xc(jj)**2 > rad**2) cycle
         nt = nt + 1
         where (fl(ii,jj,:)) nd = nd + 1
      enddo
   enddo
   endsubroutine shock_detection

   integer(I4P) function far(is)
   !< Cells flagged by sensor is farther than 3 h from the shock.
   integer(I4P), intent(in) :: is
   integer(I4P)             :: ii, jj
   far = 0
   do jj=4, n-3
      do ii=4, n-3
         if (fl(ii,jj,is) .and. abs(xc(ii)) > 3._R_P*2._R_P/n) far = far + 1
      enddo
   enddo
   endfunction far

   subroutine stencil_sensor(is, st, s)
   integer(I4P), intent(in)  :: is
   real(R_P),    intent(in)  :: st(-2:3)
   real(R_P),    intent(out) :: s
   real(R_P)                 :: rmu
   rmu = 1._R_P/max(maxval(abs(st)), 1._R_P)
   select case(is)
   case(1) ; call seth_sensor_teno(st, rmu, SETH_SENSOR_TENO_CT, s)
   case(2) ; call seth_sensor_ren(st, rmu, SETH_SENSOR_REN_PSIC, s)
   case(3) ; call seth_sensor_lcr(st, rmu, SETH_SENSOR_LCR_PSIC, s)
   case(4) ; call seth_sensor_tau(st, rmu, SETH_SENSOR_TAU_ETA, s)
   endselect
   endsubroutine stencil_sensor

   subroutine report(case, par, fr)
   character(*), intent(in) :: case
   real(R_P),    intent(in) :: par, fr(NS)
   integer(I4P)             :: is
   if (.not.vv_report()) return
   do is=1, NS
      write(ur, '(A,1X,ES12.4,1X,A,F10.5)') case, par, trim(NAMES(is)), fr(is)
   enddo
   endsubroutine report

   subroutine fail(what)
   character(*), intent(in) :: what
   print '(A)', '   FAILED: '//what
   passed = .false.
   endsubroutine fail
endprogram seth_test_sensor_fields
