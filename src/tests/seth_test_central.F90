!< SETH test: the KEEP-PE central fluxes (seth_keep_pe*) and the Ducros sensors (seth_ducros*), the hybrid's other half.
program seth_test_central
!< SETH test: the KEEP-PE central fluxes (seth_keep_pe*) and the Ducros sensors (seth_ducros*), the hybrid's other half.
!<
!< On periodic lines of N cells (h = 1/N), in each direction d = 1, 2, 3 (the normal velocity the d-th component):
!< 1. Order: (F_{i+1/2} - F_{i-1/2})/h of a smooth state against the exact derivative of the Euler flux, N = 32 .. 256:
!<    the observed L_inf order of every component >= 2L - 0.3 (2, 4, 6, 8; the optimised eight-point a_l: 6).
!< 2. Consistency: a uniform state gives F = f(q), the Euler flux, to round-off.
!< 3. Kinetic energy: on a random state of uniform pressure, sum_i (u_i . R(rho u)_i - |u_i|^2/2 R(rho)_i) = 0 to 1e-12
!<    of its scale (R the flux differences): the convective terms neither create nor destroy kinetic energy.
!< 4. Pressure equilibrium: a uniform velocity and pressure across a random density (ratios up to 100): R(rho u) =
!<    u R(rho) and R(E) = |u|^2/2 R(rho) to 1e-12 (velocity and pressure stay uniform).
!< 5. Ducros: pure dilatation phi = psi = 1 (compression), expansion phi = 1, psi = 0, pure rotation 0, and
!<    psi = sqrt(phi) in compression.
!< 6. Device: seth_keep_pe6 in an OpenACC loop (or its OpenMP-offload twin) equals the host to 1e-13 (built with
!<    tests-nvf-oac or tests-amd-omp; on the host otherwise).
!< The comparison with ADAM's compute_kep_flux and compute_ducros is seth_test_central_adam (R8P only).
use seth
use seth_test_parameters, only : CENTRAL_FLOOR_N, COMPLEX_STEP, ORDERS_AT_FLOOR_OK, ORDER_FLOOR, &
                                 TOL_DEVICE_HOST, TOL_DUCROS, TOL_DUCROS_SQRT, TOL_FLUX, TOL_KE_PEP
use seth_test_rng, only : seed_rng
implicit none
real(R_P),    parameter :: GAM = 1.4_R_P
real(R_P),    parameter :: PI = acos(-1._R_P)
integer(I4P), parameter :: NK = 5_I4P
character(12), parameter :: NAMES(NK) = ['keep_pe2    ', 'keep_pe4    ', 'keep_pe6    ', 'keep_pe8    ', 'keep_pe8_opt']
integer(I4P), parameter :: HALF(NK) = [1, 2, 3, 4, 4]
real(R_P),    parameter :: DESIGN(NK) = [2._R_P, 4._R_P, 6._R_P, 8._R_P, 6._R_P]
real(R_P)               :: err(4), order, w(5,-3:4), f(5), fe(5), x, s, phi, psi, g(3,3), worst
integer(I4P)            :: k, d, n, ig, i
logical                 :: passed

passed = .true.
call seed_rng(1013_I4P)
do k=1, NK
   do d=1, 3
      ! 1. order
      do ig=1, 4
         n = 32*2**(ig-1)
         err(ig) = derivative_error(k, d, n)
      enddo
      order = -1._R_P ! the finest pair whose finer error is above round-off: ORDER_FLOOR + CENTRAL_FLOOR_N n (R8P 1e-12)
      do ig=3, 1, -1
         if (err(ig+1) > ORDER_FLOOR + CENTRAL_FLOOR_N*32*2**ig) then
            order = log(err(ig)/err(ig+1))/log(2._R_P) ; exit
         endif
      enddo
      print '(A12,A,I1,A,4ES10.2,A,F6.2)', NAMES(k), ' d=', d, ' Linf error ', err, ' order ', order
      if (order < DESIGN(k) - 0.3_R_P .and. .not.(ORDERS_AT_FLOOR_OK .and. order < 0._R_P)) passed = .false.
      ! 2. consistency
      do i=-3, 4
         w(:,i) = [1.3_R_P, 0.4_R_P, -0.2_R_P, 0.7_R_P, 2.1_R_P]
      enddo
      call flux(k, d, w, f)
      call euler_flux(d, w(:,0), fe)
      if (maxval(abs(f - fe)) > TOL_FLUX*maxval(abs(fe))) then
         print '(A,5ES12.4)', NAMES(k)//': uniform state, F - f(q) = ', f - fe ; passed = .false.
      endif
   enddo
   ! 3. kinetic energy and 4. pressure equilibrium
   call ke_pep(k, x, s)
   print '(A12,A,ES10.2,A,ES10.2)', NAMES(k), ' kinetic-energy identity ', x, ', pressure equilibrium ', s
   if (x > TOL_KE_PEP .or. s > TOL_KE_PEP) passed = .false.
enddo
! 6. the device: seth_keep_pe6 in an OpenACC loop against the host, on random states (without OpenACC both on the
!    host); within 1e-13 of the flux scale (contractions differ)
block
integer(I4P), parameter :: NS = 4096_I4P
real(R_P)               :: ws(5,-2:3,NS), fdev(5,NS), fhost(5,NS), wl(5,-2:3)
integer(I4P)            :: j
do j=1, NS
   do i=-2, 3
      call random_number(ws(:,i,j))
      ws(1,i,j) = 0.1_R_P + ws(1,i,j) ; ws(5,i,j) = 0.1_R_P + ws(5,i,j) ; ws(2:4,i,j) = 2._R_P*ws(2:4,i,j) - 1._R_P
   enddo
enddo
!$acc parallel loop gang vector copyin(ws) copyout(fdev) private(wl)
!$omp target teams distribute parallel do map(to:ws) map(from:fdev) private(wl)
do j=1, NS
   wl = ws(:,:,j)
   call seth_keep_pe6(1 + mod(j, 3), GAM, wl, fdev(:,j))
enddo
do j=1, NS
   wl = ws(:,:,j)
   call seth_keep_pe6(1 + mod(j, 3), GAM, wl, fhost(:,j))
enddo
worst = maxval(abs(fdev - fhost))/maxval(abs(fhost))
print '(A,ES10.2)', 'seth_keep_pe6, device against host, worst scaled difference ', worst
if (worst > TOL_DEVICE_HOST) passed = .false.
endblock
! 5. Ducros
g = 0._R_P ; g(1,1) = -1._R_P ; g(2,2) = -1._R_P ; g(3,3) = -1._R_P
call seth_ducros(0._R_P, g, phi) ; call seth_ducros_compression(0._R_P, g, psi)
if (abs(phi - 1._R_P) > TOL_DUCROS .or. abs(psi - 1._R_P) > TOL_DUCROS) then
   print '(A)', 'Ducros: compression not 1' ; passed = .false.
endif
g = -g
call seth_ducros(0._R_P, g, phi) ; call seth_ducros_compression(0._R_P, g, psi)
if (abs(phi - 1._R_P) > TOL_DUCROS .or. psi /= 0._R_P) then ; print '(A)', 'Ducros: expansion' ; passed = .false. ; endif
g = 0._R_P ; g(1,2) = 1._R_P ; g(2,1) = -1._R_P
call seth_ducros(0._R_P, g, phi) ; call seth_ducros_compression(0._R_P, g, psi)
if (phi /= 0._R_P .or. psi /= 0._R_P) then ; print '(A)', 'Ducros: rotation not 0' ; passed = .false. ; endif
worst = 0._R_P
do n=1, 10000
   call random_number(g) ; g = 2._R_P*g - 1._R_P
   g(1,1) = -abs(g(1,1)) ; g(2,2) = -abs(g(2,2)) ; g(3,3) = -abs(g(3,3)) ! compression
   call seth_ducros(1.e-6_R_P, g, phi) ; call seth_ducros_compression(1.e-6_R_P, g, psi)
   worst = max(worst, abs(psi - sqrt(phi)))
enddo
print '(A,ES10.2)', 'Ducros: |psi - sqrt(phi)| in compression, worst ', worst
if (worst > TOL_DUCROS_SQRT) passed = .false.
if (.not.passed) error stop 'central flux test failed'
print '(A)', 'seth_test_central: ok'

contains
   subroutine flux(k, d, w, f)
   integer(I4P), intent(in)  :: k, d
   real(R_P),    intent(in)  :: w(5,-3:4)
   real(R_P),    intent(out) :: f(5)
   select case(k)
   case(1) ; call seth_keep_pe2(d, GAM, w(:,0:1), f)
   case(2) ; call seth_keep_pe4(d, GAM, w(:,-1:2), f)
   case(3) ; call seth_keep_pe6(d, GAM, w(:,-2:3), f)
   case(4) ; call seth_keep_pe8(d, GAM, w(:,-3:4), f)
   case(5) ; call seth_keep_pe8_opt(d, GAM, w(:,-3:4), f)
   endselect
   endsubroutine flux

   pure subroutine euler_flux(d, q, f)
   !< The Euler flux in direction d of the primitive state q = (rho, u, v, w, p).
   integer(I4P), intent(in)  :: d
   real(R_P),    intent(in)  :: q(5)
   real(R_P),    intent(out) :: f(5)
   real(R_P)                 :: e
   e = q(5)/(GAM - 1._R_P) + 0.5_R_P*q(1)*(q(2)**2 + q(3)**2 + q(4)**2)
   f(1) = q(1)*q(1+d)
   f(2:4) = q(1)*q(1+d)*q(2:4)
   f(1+d) = f(1+d) + q(5)
   f(5) = (e + q(5))*q(1+d)
   endsubroutine euler_flux

   pure function state(x) result(q)
   !< A smooth periodic primitive state on [0, 1].
   real(R_P), intent(in) :: x
   real(R_P)             :: q(5)
   q(1) = 1._R_P + 0.3_R_P*sin(2._R_P*PI*x)
   q(2) = 0.5_R_P + 0.2_R_P*cos(2._R_P*PI*x)
   q(3) = -0.3_R_P + 0.1_R_P*sin(4._R_P*PI*x + 0.3_R_P)
   q(4) = 0.2_R_P*cos(2._R_P*PI*x + 1._R_P)
   q(5) = 1._R_P + 0.2_R_P*cos(4._R_P*PI*x)
   endfunction state

   function derivative_error(k, d, n) result(e)
   !< Linf over the cells of (F_{i+1/2} - F_{i-1/2})/h - d f(q)/dx, all components, the state with its d-th velocity
   !< component as the normal one; the exact derivative by the complex step (exact to round-off).
   integer(I4P), intent(in) :: k, d, n
   real(R_P)                :: e, h, ff(5,0:n), fd(5), q(5,-3:4)
   integer(I4P)             :: i, m
   h = 1._R_P/n
   do i=0, n ! face i+1/2, cells i-3 .. i+4 (periodic)
      do m=-3, 4
         q(:,m) = state((modulo(i+m, n) + 0.5_R_P)*h)
         call swap(d, q(:,m))
      enddo
      call flux(k, d, q, ff(:,i))
   enddo
   e = 0._R_P
   do i=1, n
      fd = aimag(cflux(d, cmplx((modulo(i, n) + 0.5_R_P)*h, COMPLEX_STEP, R_P)))/COMPLEX_STEP
      e = max(e, maxval(abs((ff(:,i) - ff(:,i-1))/h - fd)))
   enddo
   endfunction derivative_error

   pure function cflux(d, z) result(f)
   !< The Euler flux of the smooth state at the complex point z (for the complex-step derivative).
   integer(I4P), intent(in) :: d
   complex(R_P), intent(in) :: z
   complex(R_P)             :: f(5), q(5), e, t
   q(1) = 1._R_P + 0.3_R_P*sin(2._R_P*PI*z)
   q(2) = 0.5_R_P + 0.2_R_P*cos(2._R_P*PI*z)
   q(3) = -0.3_R_P + 0.1_R_P*sin(4._R_P*PI*z + 0.3_R_P)
   q(4) = 0.2_R_P*cos(2._R_P*PI*z + 1._R_P)
   q(5) = 1._R_P + 0.2_R_P*cos(4._R_P*PI*z)
   t = q(1+d) ; q(1+d) = q(2) ; q(2) = t
   e = q(5)/(GAM - 1._R_P) + 0.5_R_P*q(1)*(q(2)**2 + q(3)**2 + q(4)**2)
   f(1) = q(1)*q(1+d)
   f(2:4) = q(1)*q(1+d)*q(2:4)
   f(1+d) = f(1+d) + q(5)
   f(5) = (e + q(5))*q(1+d)
   endfunction cflux

   pure subroutine swap(d, q)
   !< Put the first velocity component of the smooth state into the d-th slot (the normal one).
   integer(I4P), intent(in)    :: d
   real(R_P),    intent(inout) :: q(5)
   real(R_P)                   :: t
   t = q(1+d) ; q(1+d) = q(2) ; q(2) = t
   endsubroutine swap

   pure function sw(d, q) result(r)
   integer(I4P), intent(in) :: d
   real(R_P),    intent(in) :: q(5)
   real(R_P)                :: r(5)
   r = q ; call swap(d, r)
   endfunction sw

   subroutine ke_pep(k, ke, pep)
   !< The kinetic-energy identity on a random state of uniform pressure, and the pressure equilibrium of a random
   !< density with uniform velocity and pressure; both on 64 periodic cells, direction 1 (the others by symmetry).
   integer(I4P), intent(in)  :: k
   real(R_P),    intent(out) :: ke, pep
   integer(I4P), parameter   :: N = 64_I4P
   real(R_P)                 :: q(5,0:N-1), ff(5,0:N-1), r(5,0:N-1), sc, x
   integer(I4P)              :: i
   ! kinetic energy: random rho, u, v, w; p = 1
   do i=0, N-1
      call random_number(q(:,i))
      q(1,i) = 0.5_R_P + q(1,i) ; q(2:4,i) = 2._R_P*q(2:4,i) - 1._R_P ; q(5,i) = 1._R_P
   enddo
   call faces(k, q, ff) ; call divergence(ff, r)
   ke = 0._R_P ; sc = 0._R_P
   do i=0, N-1
      ke = ke + sum(q(2:4,i)*r(2:4,i)) - 0.5_R_P*sum(q(2:4,i)**2)*r(1,i)
      sc = sc + abs(sum(q(2:4,i)*r(2:4,i))) + abs(0.5_R_P*sum(q(2:4,i)**2)*r(1,i))
   enddo
   ke = abs(ke)/sc
   ! pressure equilibrium: random rho in [0.01, 1]; u = (0.7, -0.4, 0.2), p = 1
   do i=0, N-1
      call random_number(x)
      q(:,i) = [0.01_R_P + 0.99_R_P*x, 0.7_R_P, -0.4_R_P, 0.2_R_P, 1._R_P]
   enddo
   call faces(k, q, ff) ; call divergence(ff, r)
   pep = 0._R_P ; sc = maxval(abs(r))
   do i=0, N-1
      pep = max(pep, maxval(abs(r(2:4,i) - q(2:4,i)*r(1,i))), abs(r(5,i) - 0.5_R_P*sum(q(2:4,i)**2)*r(1,i)))
   enddo
   pep = pep/sc
   endsubroutine ke_pep

   subroutine faces(k, q, ff)
   !< The face fluxes of a periodic line of NP cells, direction 1.
   integer(I4P), intent(in)  :: k
   real(R_P),    intent(in)  :: q(:,0:)
   real(R_P),    intent(out) :: ff(:,0:)
   real(R_P)                 :: wl(5,-3:4)
   integer(I4P)              :: i, m, np
   np = size(q, 2)
   do i=0, np-1
      do m=-3, 4
         wl(:,m) = q(:,modulo(i+m, np))
      enddo
      call flux(k, 1_I4P, wl, ff(:,i))
   enddo
   endsubroutine faces

   subroutine divergence(ff, r)
   !< R_i = -(F_{i+1/2} - F_{i-1/2}), periodic.
   real(R_P), intent(in)  :: ff(:,0:)
   real(R_P), intent(out) :: r(:,0:)
   integer(I4P)           :: i, np
   np = size(ff, 2)
   do i=0, np-1
      r(:,i) = -(ff(:,i) - ff(:,modulo(i-1, np)))
   enddo
   endsubroutine divergence
endprogram seth_test_central
