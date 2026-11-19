!< SETH test: whole-domain convergence, the spatial operator and SSP-RK3 advection, critical points included.
program seth_test_convergence
!< SETH test: whole-domain convergence, the spatial operator and SSP-RK3 advection, critical points included.
!<
!< u = sin(pi x - sin(pi x)/pi) on [-1, 1], periodic (Henrick et al. 2005's field: first-order critical points that
!< are not symmetry points, where sin(pi x) would let the degraded errors cancel), point values on N nodes
!< (finite-difference WENO, a = 1). Unlike the pointwise order tests every node counts, and the literature predicts
!< which schemes keep their order in Linf across the critical points (eps = 1e-40, as Borges et al.):
!<
!< 1. the operator: Linf of D_i - u'(x_i), D_i = (F_{i+1/2} - F_{i-1/2})/h, N = 20 .. 640, no time error;
!< 2. SSP-RK3 advection over one period (t = 2), dt = 0.4 h (h/h0)^((o-3)/3) so that the time error O(dt^3) falls
!<    with the spatial one (Henrick et al. 2005 use dt ~ h^(5/3) for fifth order), N = 20 .. 160, Linf of u - u0.
!<
!< Asserted on the finest pair of N whose errors stay above round-off (1e3 eps/h): WENO-Z with p = 2 and TENO at
!< their order (first-order critical points do not
!< degrade them: Borges 2008, Tang 2024, Zhang 2024b); WENO-JS at least Henrick et al.'s critical-point rate
!< max(2r - 3, 1) (an Linf over all nodes is bounded by its worst node); WENO-Z with p = 1 at least 4 for WENO-Z5;
!< the centred WENO the rate of its upwind WENO-JS (where a difference of the blend's cells vanishes, phi = O(1) and
!< the state is the upwind one: measured 1.55, 3.56, 4.76 for orders 4, 6, 8).
!< Every scheme and order is printed.
!< Per kind (seth_test_parameters): the round-off floor of the operator, and in R4P an order loss and the acceptance of
!< a scheme that reaches round-off before any pair of levels is above it (it converged at least as fast as the floor
!< allows: a degraded scheme would stay above it). The R8P values keep the test as it was.
use seth
use seth_test_parameters, only : CONV_FLOOR_OP, CONV_ORDER_LOSS, EPS_ZT, ORDERS_AT_FLOOR_OK
use seth_test_schemes
use seth_vv_report, only : vv_report, vv_open
implicit none
real(R_P),    parameter :: PI = acos(-1._R_P)
integer(I4P), parameter :: NLEV = 6_I4P
real(R_P)               :: eop(NLEV), eadv(4), oop, oadv, need
integer(I4P)            :: id, p, j, n
logical                 :: passed, advect_it
integer                 :: ue, uo

passed = .true.
if (vv_report()) then
   call vv_open('convergence.dat', ue, 'scheme p N Linf(D - du/dx) (Henrick field)')
   call vv_open('convergence-orders.dat', uo, 'scheme p operator-order advection-order (-1: not run)')
endif
print '(A16,A3,A12,A12)', 'scheme', 'p', 'operator', 'RK3 advect'
do id=1, NSCHEMES
   if (SCHEMES(id)%interp) cycle ! the finite-difference flux setting reconstructs: interpolations elsewhere
   do p=1, merge(2, 1, SCHEMES(id)%family == FAM_Z)
      do j=1, NLEV
         n = 20*2**(j-1)
         eop(j) = operator_error(id, p, n)
         if (vv_report()) write(ue, '(A,1X,I1,1X,I5,ES24.15)') trim(SCHEMES(id)%name), p, n, eop(j)
      enddo
      oop = finest_order(eop, NLEV, 1)
      advect_it = (SCHEMES(id)%family == FAM_TENO .or. (SCHEMES(id)%family == FAM_Z .and. p == 2)) .and. &
                  SCHEMES(id)%s <= 3
      oadv = -1._R_P
      if (advect_it) then
         do j=1, 4
            n = 20*2**(j-1)
            eadv(j) = advection_error(id, p, n)
         enddo
         oadv = finest_order(eadv, 4, 0)
      endif
      print '(A16,I3,2F12.2)', SCHEMES(id)%name, p, oop, oadv
      if (vv_report()) write(uo, '(A,1X,I1,2F8.3)') trim(SCHEMES(id)%name), p, oop, oadv
      ! the expected order
      if (SCHEMES(id)%family == FAM_JS .or. SCHEMES(id)%family == FAM_C) then
         need = real(max(2*SCHEMES(id)%s - 3, 1), R_P)
      elseif (SCHEMES(id)%family == FAM_Z .and. p == 1) then
         need = merge(4._R_P, -1._R_P, trim(SCHEMES(id)%name) == 'wenoz5')
      elseif (SCHEMES(id)%family == FAM_Z .and. SCHEMES(id)%s == 2) then
         need = -1._R_P ! WENO-Z3 drops at critical points
      else
         need = SCHEMES(id)%order
      endif
      if (need > 0._R_P .and. oop < need - 0.5_R_P - CONV_ORDER_LOSS .and. .not.(ORDERS_AT_FLOOR_OK .and. oop < 0._R_P)) then
         print '(A,F5.1)', '   ^ operator below ', need ; passed = .false.
      endif
      if (advect_it .and. oadv < need - 0.5_R_P - CONV_ORDER_LOSS .and. .not.(ORDERS_AT_FLOOR_OK .and. oadv < 0._R_P)) then
         print '(A,F5.1)', '   ^ advection below ', need ; passed = .false.
      endif
   enddo
enddo
if (vv_report()) then ; close(ue) ; close(uo) ; endif
if (.not.passed) error stop 'convergence test failed'
print '(A)', 'seth_test_convergence: ok'

contains
   subroutine fluxes(id, p, n, u, f)
   !< F_{i+1/2} (left state, a = 1) on every face of the periodic grid.
   integer(I4P), intent(in)  :: id, p, n
   real(R_P),    intent(in)  :: u(0:n-1)
   real(R_P),    intent(out) :: f(0:n-1)
   real(R_P)                 :: v(VLO:VHI), rp, rm, par
   integer(I4P)              :: i, m
   par = SCHEMES(id)%par
   if (SCHEMES(id)%family == FAM_Z) par = real(p, R_P)
   do i=0, n-1
      do m=VLO, VHI
         v(m) = u(modulo(i+m, n))
      enddo
      call reconstruct(id, v, v, rp, rm, eps=EPS_ZT, par=par)
      f(i) = rp
   enddo
   endsubroutine fluxes

   function operator_error(id, p, n) result(e)
   !< Linf of D - u' over all nodes, u = sin(pi x).
   integer(I4P), intent(in) :: id, p, n
   real(R_P)                :: e, h, u(0:n-1), f(0:n-1), x
   integer(I4P)             :: i
   h = 2._R_P/n
   do i=0, n-1
      u(i) = henrick(-1._R_P + i*h)
   enddo
   call fluxes(id, p, n, u, f)
   e = 0._R_P
   do i=0, n-1
      x = -1._R_P + i*h
      e = max(e, abs((f(i) - f(modulo(i-1, n)))/h - cos(PI*x - sin(PI*x)/PI)*(PI - cos(PI*x))))
   enddo
   endfunction operator_error

   function advection_error(id, p, n) result(e)
   !< Linf of u(t = 2) - u0 after SSP-RK3 with dt = 0.4 h (h/h0)^((o-3)/3).
   integer(I4P), intent(in) :: id, p, n
   real(R_P)                :: e, h, dt, u(0:n-1), u0(0:n-1), u1(0:n-1), u2(0:n-1), l(0:n-1), f(0:n-1)
   integer(I4P)             :: i, s, nsteps
   h = 2._R_P/n
   do i=0, n-1
      u0(i) = henrick(-1._R_P + i*h)
   enddo
   dt = 0.4_R_P*h*(h/0.1_R_P)**((SCHEMES(id)%order - 3._R_P)/3._R_P)
   nsteps = ceiling(2._R_P/dt) ; dt = 2._R_P/nsteps
   u = u0
   do s=1, nsteps
      call fluxes(id, p, n, u, f)  ; call dif(n, h, f, l) ; u1 = u + dt*l
      call fluxes(id, p, n, u1, f) ; call dif(n, h, f, l) ; u2 = 0.75_R_P*u + 0.25_R_P*(u1 + dt*l)
      call fluxes(id, p, n, u2, f) ; call dif(n, h, f, l) ; u = u/3._R_P + 2._R_P/3._R_P*(u2 + dt*l)
   enddo
   e = maxval(abs(u - u0))
   endfunction advection_error

   pure function henrick(x) result(y)
   !< sin(pi x - sin(pi x)/pi).
   real(R_P), intent(in) :: x
   real(R_P)             :: y
   y = sin(PI*x - sin(PI*x)/PI)
   endfunction henrick

   function finest_order(e, nl, kind) result(o)
   !< Observed order on the finest pair of levels whose finer error is above round-off: CONV_FLOOR_OP eps/h for the
   !< derivative (kind 1, h = 2/(20 2^(j-1)); R8P 1e3), 1e3 eps for the solution (kind 0); -1 if none.
   integer(I4P), intent(in) :: nl, kind
   real(R_P),    intent(in) :: e(nl)
   real(R_P)                :: o, flo
   integer(I4P)             :: j
   o = -1._R_P
   do j=nl-1, 1, -1
      flo = 1.e3_R_P*epsilon(1._R_P)
      if (kind == 1) flo = CONV_FLOOR_OP*epsilon(1._R_P)*(20._R_P*2._R_P**j)/2._R_P
      if (e(j+1) > flo) then
         o = log(e(j)/e(j+1))/log(2._R_P) ; return
      endif
   enddo
   endfunction finest_order

   subroutine dif(n, h, f, l)
   !< -(F_{i+1/2} - F_{i-1/2})/h, periodic.
   integer(I4P), intent(in)  :: n
   real(R_P),    intent(in)  :: h, f(0:n-1)
   real(R_P),    intent(out) :: l(0:n-1)
   integer(I4P)              :: i
   do i=0, n-1
      l(i) = -(f(i) - f(modulo(i-1, n)))/h
   enddo
   endsubroutine dif
endprogram seth_test_convergence
