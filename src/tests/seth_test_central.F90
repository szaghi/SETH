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
!< 5. ADAM: seth_keep_pe6 equals ADAM's compute_kep_flux bitwise on random states (gfortran; other compilers within
!<    2e-15 of the flux scale, fused multiply-adds), and seth_ducros equals compute_ducros with its override off.
!< 6. Ducros: pure dilatation phi = psi = 1 (compression), expansion phi = 1, psi = 0, pure rotation 0, and
!<    psi = sqrt(phi) in compression.
!< 7. Device: seth_keep_pe6 in an OpenACC loop equals the host to 1e-13 (built with tests-nvf-oac; on the host
!<    otherwise).
use seth
use seth_test_rng, only : seed_rng
use seth_adam_reference, only : compute_kep_flux, compute_ducros
implicit none
real(R8P),    parameter :: GAM = 1.4_R8P
real(R8P),    parameter :: PI = acos(-1._R8P)
integer(I4P), parameter :: NK = 5_I4P
character(12), parameter :: NAMES(NK) = ['keep_pe2    ', 'keep_pe4    ', 'keep_pe6    ', 'keep_pe8    ', 'keep_pe8_opt']
integer(I4P), parameter :: HALF(NK) = [1, 2, 3, 4, 4]
real(R8P),    parameter :: DESIGN(NK) = [2._R8P, 4._R8P, 6._R8P, 8._R8P, 6._R8P]
real(R8P)               :: err(4), order, w(5,-3:4), f(5), fa(5), fe(5), x, s, phi, psi, g(3,3), worst
integer(I4P)            :: k, d, n, ig, i, nbit, ntot
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
      order = log(err(3)/err(4))/log(2._R8P)
      if (err(4) < 1.e-11_R8P) order = log(err(2)/err(3))/log(2._R8P) ! round-off floor near 1e-12
      print '(A12,A,I1,A,4ES10.2,A,F6.2)', NAMES(k), ' d=', d, ' Linf error ', err, ' order ', order
      if (order < DESIGN(k) - 0.3_R8P) passed = .false.
      ! 2. consistency
      do i=-3, 4
         w(:,i) = [1.3_R8P, 0.4_R8P, -0.2_R8P, 0.7_R8P, 2.1_R8P]
      enddo
      call flux(k, d, w, f)
      call euler_flux(d, w(:,0), fe)
      if (maxval(abs(f - fe)) > 1.e-14_R8P*maxval(abs(fe))) then
         print '(A,5ES12.4)', NAMES(k)//': uniform state, F - f(q) = ', f - fe ; passed = .false.
      endif
   enddo
   ! 3. kinetic energy and 4. pressure equilibrium
   call ke_pep(k, x, s)
   print '(A12,A,ES10.2,A,ES10.2)', NAMES(k), ' kinetic-energy identity ', x, ', pressure equilibrium ', s
   if (x > 1.e-12_R8P .or. s > 1.e-12_R8P) passed = .false.
enddo
! 5. ADAM
nbit = 0 ; ntot = 0 ; worst = 0._R8P
do n=1, 20000
   do i=-2, 3
      call random_number(w(:,i))
      w(1,i) = 0.1_R8P + w(1,i) ; w(5,i) = 0.1_R8P + w(5,i) ; w(2:4,i) = 2._R8P*w(2:4,i) - 1._R8P
   enddo
   d = 1 + mod(n, 3)
   call seth_keep_pe6(d, GAM, w(:,-2:3), f)
   call compute_kep_flux(d, GAM, w(:,-2:3), fa)
   ntot = ntot + 5 ; nbit = nbit + count(f == fa)
   worst = max(worst, maxval(abs(f - fa))/maxval(abs(fa)))
   call random_number(g) ; g = 2._R8P*g - 1._R8P
   call seth_ducros(1.e-3_R8P, g, phi)
   call compute_ducros(1.e-3_R8P, 0._R8P, 1._R8P, g, psi)
   if (phi /= psi) then ; print '(A)', 'seth_ducros differs from ADAM''s compute_ducros' ; passed = .false. ; exit ; endif
enddo
print '(A,I0,A,I0,A,ES10.2)', 'seth_keep_pe6 against ADAM''s compute_kep_flux: ', nbit, ' of ', ntot, &
      ' bitwise, worst scaled difference ', worst
if (worst > 2.e-15_R8P) passed = .false.
! 7. the device: seth_keep_pe6 in an OpenACC loop against the host, on the states of item 5 (without OpenACC both on
!    the host); within 1e-13 of the flux scale (contractions differ)
block
integer(I4P), parameter :: NS = 4096_I4P
real(R8P)               :: ws(5,-2:3,NS), fdev(5,NS), fhost(5,NS), wl(5,-2:3)
integer(I4P)            :: j
do j=1, NS
   do i=-2, 3
      call random_number(ws(:,i,j))
      ws(1,i,j) = 0.1_R8P + ws(1,i,j) ; ws(5,i,j) = 0.1_R8P + ws(5,i,j) ; ws(2:4,i,j) = 2._R8P*ws(2:4,i,j) - 1._R8P
   enddo
enddo
!$acc parallel loop gang vector copyin(ws) copyout(fdev) private(wl)
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
if (worst > 1.e-13_R8P) passed = .false.
endblock
! 6. Ducros
g = 0._R8P ; g(1,1) = -1._R8P ; g(2,2) = -1._R8P ; g(3,3) = -1._R8P
call seth_ducros(0._R8P, g, phi) ; call seth_ducros_compression(0._R8P, g, psi)
if (abs(phi - 1._R8P) > 1.e-15_R8P .or. abs(psi - 1._R8P) > 1.e-15_R8P) then
   print '(A)', 'Ducros: compression not 1' ; passed = .false.
endif
g = -g
call seth_ducros(0._R8P, g, phi) ; call seth_ducros_compression(0._R8P, g, psi)
if (abs(phi - 1._R8P) > 1.e-15_R8P .or. psi /= 0._R8P) then ; print '(A)', 'Ducros: expansion' ; passed = .false. ; endif
g = 0._R8P ; g(1,2) = 1._R8P ; g(2,1) = -1._R8P
call seth_ducros(0._R8P, g, phi) ; call seth_ducros_compression(0._R8P, g, psi)
if (phi /= 0._R8P .or. psi /= 0._R8P) then ; print '(A)', 'Ducros: rotation not 0' ; passed = .false. ; endif
worst = 0._R8P
do n=1, 10000
   call random_number(g) ; g = 2._R8P*g - 1._R8P
   g(1,1) = -abs(g(1,1)) ; g(2,2) = -abs(g(2,2)) ; g(3,3) = -abs(g(3,3)) ! compression
   call seth_ducros(1.e-6_R8P, g, phi) ; call seth_ducros_compression(1.e-6_R8P, g, psi)
   worst = max(worst, abs(psi - sqrt(phi)))
enddo
print '(A,ES10.2)', 'Ducros: |psi - sqrt(phi)| in compression, worst ', worst
if (worst > 1.e-14_R8P) passed = .false.
if (.not.passed) error stop 'central flux test failed'
print '(A)', 'seth_test_central: ok'

contains
   subroutine flux(k, d, w, f)
   integer(I4P), intent(in)  :: k, d
   real(R8P),    intent(in)  :: w(5,-3:4)
   real(R8P),    intent(out) :: f(5)
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
   real(R8P),    intent(in)  :: q(5)
   real(R8P),    intent(out) :: f(5)
   real(R8P)                 :: e
   e = q(5)/(GAM - 1._R8P) + 0.5_R8P*q(1)*(q(2)**2 + q(3)**2 + q(4)**2)
   f(1) = q(1)*q(1+d)
   f(2:4) = q(1)*q(1+d)*q(2:4)
   f(1+d) = f(1+d) + q(5)
   f(5) = (e + q(5))*q(1+d)
   endsubroutine euler_flux

   pure function state(x) result(q)
   !< A smooth periodic primitive state on [0, 1].
   real(R8P), intent(in) :: x
   real(R8P)             :: q(5)
   q(1) = 1._R8P + 0.3_R8P*sin(2._R8P*PI*x)
   q(2) = 0.5_R8P + 0.2_R8P*cos(2._R8P*PI*x)
   q(3) = -0.3_R8P + 0.1_R8P*sin(4._R8P*PI*x + 0.3_R8P)
   q(4) = 0.2_R8P*cos(2._R8P*PI*x + 1._R8P)
   q(5) = 1._R8P + 0.2_R8P*cos(4._R8P*PI*x)
   endfunction state

   function derivative_error(k, d, n) result(e)
   !< Linf over the cells of (F_{i+1/2} - F_{i-1/2})/h - d f(q)/dx, all components, the state with its d-th velocity
   !< component as the normal one; the exact derivative by the complex step (exact to round-off).
   integer(I4P), intent(in) :: k, d, n
   real(R8P)                :: e, h, ff(5,0:n), fd(5), q(5,-3:4)
   integer(I4P)             :: i, m
   h = 1._R8P/n
   do i=0, n ! face i+1/2, cells i-3 .. i+4 (periodic)
      do m=-3, 4
         q(:,m) = state((modulo(i+m, n) + 0.5_R8P)*h)
         call swap(d, q(:,m))
      enddo
      call flux(k, d, q, ff(:,i))
   enddo
   e = 0._R8P
   do i=1, n
      fd = aimag(cflux(d, cmplx((modulo(i, n) + 0.5_R8P)*h, 1.e-30_R8P, R8P)))/1.e-30_R8P
      e = max(e, maxval(abs((ff(:,i) - ff(:,i-1))/h - fd)))
   enddo
   endfunction derivative_error

   pure function cflux(d, z) result(f)
   !< The Euler flux of the smooth state at the complex point z (for the complex-step derivative).
   integer(I4P), intent(in) :: d
   complex(R8P), intent(in) :: z
   complex(R8P)             :: f(5), q(5), e, t
   q(1) = 1._R8P + 0.3_R8P*sin(2._R8P*PI*z)
   q(2) = 0.5_R8P + 0.2_R8P*cos(2._R8P*PI*z)
   q(3) = -0.3_R8P + 0.1_R8P*sin(4._R8P*PI*z + 0.3_R8P)
   q(4) = 0.2_R8P*cos(2._R8P*PI*z + 1._R8P)
   q(5) = 1._R8P + 0.2_R8P*cos(4._R8P*PI*z)
   t = q(1+d) ; q(1+d) = q(2) ; q(2) = t
   e = q(5)/(GAM - 1._R8P) + 0.5_R8P*q(1)*(q(2)**2 + q(3)**2 + q(4)**2)
   f(1) = q(1)*q(1+d)
   f(2:4) = q(1)*q(1+d)*q(2:4)
   f(1+d) = f(1+d) + q(5)
   f(5) = (e + q(5))*q(1+d)
   endfunction cflux

   pure subroutine swap(d, q)
   !< Put the first velocity component of the smooth state into the d-th slot (the normal one).
   integer(I4P), intent(in)    :: d
   real(R8P),    intent(inout) :: q(5)
   real(R8P)                   :: t
   t = q(1+d) ; q(1+d) = q(2) ; q(2) = t
   endsubroutine swap

   pure function sw(d, q) result(r)
   integer(I4P), intent(in) :: d
   real(R8P),    intent(in) :: q(5)
   real(R8P)                :: r(5)
   r = q ; call swap(d, r)
   endfunction sw

   subroutine ke_pep(k, ke, pep)
   !< The kinetic-energy identity on a random state of uniform pressure, and the pressure equilibrium of a random
   !< density with uniform velocity and pressure; both on 64 periodic cells, direction 1 (the others by symmetry).
   integer(I4P), intent(in)  :: k
   real(R8P),    intent(out) :: ke, pep
   integer(I4P), parameter   :: N = 64_I4P
   real(R8P)                 :: q(5,0:N-1), ff(5,0:N-1), r(5,0:N-1), sc, x
   integer(I4P)              :: i
   ! kinetic energy: random rho, u, v, w; p = 1
   do i=0, N-1
      call random_number(q(:,i))
      q(1,i) = 0.5_R8P + q(1,i) ; q(2:4,i) = 2._R8P*q(2:4,i) - 1._R8P ; q(5,i) = 1._R8P
   enddo
   call faces(k, q, ff) ; call divergence(ff, r)
   ke = 0._R8P ; sc = 0._R8P
   do i=0, N-1
      ke = ke + sum(q(2:4,i)*r(2:4,i)) - 0.5_R8P*sum(q(2:4,i)**2)*r(1,i)
      sc = sc + abs(sum(q(2:4,i)*r(2:4,i))) + abs(0.5_R8P*sum(q(2:4,i)**2)*r(1,i))
   enddo
   ke = abs(ke)/sc
   ! pressure equilibrium: random rho in [0.01, 1]; u = (0.7, -0.4, 0.2), p = 1
   do i=0, N-1
      call random_number(x)
      q(:,i) = [0.01_R8P + 0.99_R8P*x, 0.7_R8P, -0.4_R8P, 0.2_R8P, 1._R8P]
   enddo
   call faces(k, q, ff) ; call divergence(ff, r)
   pep = 0._R8P ; sc = maxval(abs(r))
   do i=0, N-1
      pep = max(pep, maxval(abs(r(2:4,i) - q(2:4,i)*r(1,i))), abs(r(5,i) - 0.5_R8P*sum(q(2:4,i)**2)*r(1,i)))
   enddo
   pep = pep/sc
   endsubroutine ke_pep

   subroutine faces(k, q, ff)
   !< The face fluxes of a periodic line of NP cells, direction 1.
   integer(I4P), intent(in)  :: k
   real(R8P),    intent(in)  :: q(:,0:)
   real(R8P),    intent(out) :: ff(:,0:)
   real(R8P)                 :: wl(5,-3:4)
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
   real(R8P), intent(in)  :: ff(:,0:)
   real(R8P), intent(out) :: r(:,0:)
   integer(I4P)           :: i, np
   np = size(ff, 2)
   do i=0, np-1
      r(:,i) = -(ff(:,i) - ff(:,modulo(i-1, np)))
   enddo
   endsubroutine divergence
endprogram seth_test_central
