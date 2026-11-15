!< SETH test: the interpolations from point values, over a whole periodic domain and at critical points.
program seth_test_interpolation
!< SETH test: the interpolations from point values, over a whole periodic domain and at critical points.
!<
!< Every interpolation of the registry (WENO-JS, WENO-Z with p = 1 and 2, TENO): point values at the cell centres,
!< the two states at the faces. The finite-difference flux tests (critical points, convergence, advection) apply to
!< the reconstructions only; this is their counterpart for the interpolations, where the face value itself is the
!< result (no division by h).
!<
!< 1. Whole domain: u = sin(pi x - sin(pi x)/pi) on [-1, 1], periodic (Henrick et al. 2005's field: first-order
!<    critical points that are not symmetry points), N = 20 .. 640 cells, Linf of the face values against the exact
!<    ones over every face, both states; the order on the finest pair above round-off (100 eps).
!< 2. Pointwise at a critical point placed at the face x = 0: f = g - g'(0) x (first order: f' = 0, f'' /= 0) and
!<    f = g - g'(0) x - g''(0) x^2/2 (second order), g the generic function of seth_test_functions, h halved from h0.
!<
!< Asserted (eps = 1e-40 for WENO-Z and TENO, the registry's 1e-6 for WENO-JS): WENO-Z with p = 2 and TENO keep the
!< order of their linear scheme minus 0.5 over the whole domain and at first-order critical points; WENO-JS reaches at
!< least r + 1 over the whole domain (2r - 1 smooth, its weights deviate by O(h) at critical points, so the face value
!< errs by O(h) times the O(h^r) candidate errors), WENO-Z with p = 1 at least that of WENO-JS. The pointwise orders at
!< second-order critical points are printed (measured, not asserted: no published interpolation analysis), and so is
!< WENO-Z3, which drops at critical points as its reconstruction does (seth_test_convergence): its worst face sits at
!< the critical point of the field, and the Linf ratio jumps with the position of the grid relative to it (measured:
!< 4.2, 2.0, 4.3, 4.9 from N = 160 to 2560). Tolerance on an order: 0.5.
use seth
use seth_test_schemes
use seth_test_functions, only : gfun, gder, gder2
use seth_vv_cases, only : vv_report, vv_open
implicit none
real(R8P),    parameter :: PI = acos(-1._R8P)
integer(I4P), parameter :: NLEV = 6_I4P, NH = 7_I4P
real(R8P)               :: e(2,NLEV), ep(2,NH), o(2), ocp1(2), ocp2(2), need, par, epsv
integer(I4P)            :: id, p, j
logical                 :: passed
integer                 :: ue, uo

passed = .true.
if (vv_report()) then
   call vv_open('interpolation.dat', ue, 'scheme p N Linf (both states, Henrick field)')
   call vv_open('interpolation-orders.dat', uo, 'scheme p domain(+) domain(-) cp1(+) cp1(-) cp2(+) cp2(-)')
endif
print '(A16,A3,6A10)', 'scheme', 'p', 'Linf +', 'Linf -', 'cp1 +', 'cp1 -', 'cp2 +', 'cp2 -'
do id=1, NSCHEMES
   if (.not.SCHEMES(id)%interp) cycle
   do p=1, merge(2, 1, SCHEMES(id)%family == FAM_Z)
      par = SCHEMES(id)%par ; if (SCHEMES(id)%family == FAM_Z) par = real(p, R8P)
      epsv = merge(1.e-40_R8P, SCHEMES(id)%eps, SCHEMES(id)%family == FAM_Z .or. SCHEMES(id)%family == FAM_TENO)
      do j=1, NLEV
         call domain_error(id, 20*2**(j-1), par, epsv, e(:,j))
         if (vv_report()) write(ue, '(A,1X,I1,1X,I5,ES24.15)') trim(SCHEMES(id)%name), p, 20*2**(j-1), maxval(e(:,j))
      enddo
      o = [finest(e(1,:), NLEV), finest(e(2,:), NLEV)]
      call critical(id, 1, par, epsv, ep) ; ocp1 = [finest(ep(1,:), NH), finest(ep(2,:), NH)]
      call critical(id, 2, par, epsv, ep) ; ocp2 = [finest(ep(1,:), NH), finest(ep(2,:), NH)]
      print '(A16,I3,6F10.2)', SCHEMES(id)%name, p, o, ocp1, ocp2
      if (vv_report()) write(uo, '(A,1X,I1,6F8.3)') trim(SCHEMES(id)%name), p, o, ocp1, ocp2
      if (SCHEMES(id)%family == FAM_JS_INTERP) then
         need = real(merge(SCHEMES(id)%s + 1, 1, SCHEMES(id)%s >= 2), R8P)
      elseif (SCHEMES(id)%family == FAM_Z .and. SCHEMES(id)%s == 2) then
         need = -1._R8P ! WENO-Z3 drops at critical points (as in seth_test_convergence): printed, not asserted
      elseif (SCHEMES(id)%family == FAM_Z .and. p == 1) then
         need = real(SCHEMES(id)%s + 1, R8P)
      else
         need = SCHEMES(id)%order
      endif
      if (need > 0._R8P .and. any(o < need - 0.5_R8P)) then
         print '(A,F5.1)', '   ^ whole-domain order below ', need ; passed = .false.
      endif
      if (SCHEMES(id)%family == FAM_TENO .or. (SCHEMES(id)%family == FAM_Z .and. p == 2 .and. SCHEMES(id)%s >= 3)) then
         if (any(ocp1 < need - 0.5_R8P)) then
            print '(A,F5.1)', '   ^ order at first-order critical points below ', need ; passed = .false.
         endif
      endif
   enddo
enddo
if (vv_report()) then ; close(ue) ; close(uo) ; endif
if (.not.passed) error stop 'interpolation test failed'
print '(A)', 'seth_test_interpolation: ok'

contains
   subroutine domain_error(id, n, par, epsv, e)
   !< Linf of the two face states over the periodic grid of n cells.
   integer(I4P), intent(in)  :: id, n
   real(R8P),    intent(in)  :: par, epsv
   real(R8P),    intent(out) :: e(2)
   real(R8P)                 :: h, u(0:n-1), v(VLO:VHI), rp, rm, xf, ex
   integer(I4P)              :: i, m
   h = 2._R8P/n
   do i=0, n-1
      u(i) = henrick(-1._R8P + (i + 0.5_R8P)*h) ! cell i = [-1 + i h, -1 + (i+1) h], its centre
   enddo
   e = 0._R8P
   do i=0, n-1 ! face i+1/2 = right face of cell i: stencil cell m is cell i+m
      do m=VLO, VHI
         v(m) = u(modulo(i+m, n))
      enddo
      call reconstruct(id, v, v, rp, rm, eps=epsv, par=par)
      xf = -1._R8P + (i + 1)*h ; ex = henrick(xf)
      e = max(e, abs([rp, rm] - ex))
   enddo
   endsubroutine domain_error

   subroutine critical(id, kind, par, epsv, e)
   !< Face errors at a critical point of order `kind` placed at the face x = 0, h halved from h0.
   integer(I4P), intent(in)  :: id, kind
   real(R8P),    intent(in)  :: par, epsv
   real(R8P),    intent(out) :: e(2,NH)
   real(R8P)                 :: h, v(VLO:VHI), rp, rm, d1, d2, x
   integer(I4P)              :: j, m
   d1 = gder(0._R8P)
   d2 = gder2(0._R8P)
   do j=1, NH
      h = 0.4_R8P*2._R8P**max(0, SCHEMES(id)%s - 3)/2._R8P**(j-1)
      do m=VLO, VHI
         x = (m - 0.5_R8P)*h
         v(m) = gfun(x) - d1*x
         if (kind == 2) v(m) = v(m) - 0.5_R8P*d2*x*x
      enddo
      call reconstruct(id, v, v, rp, rm, eps=epsv, par=par)
      e(:,j) = abs([rp, rm] - gfun(0._R8P))
   enddo
   endsubroutine critical

   pure function henrick(x) result(y)
   real(R8P), intent(in) :: x
   real(R8P)             :: y
   y = sin(PI*x - sin(PI*x)/PI)
   endfunction henrick

   pure function finest(e, nl) result(o)
   !< Observed order on the finest pair of levels whose finer error is above round-off (100 eps); -1 if none.
   integer(I4P), intent(in) :: nl
   real(R8P),    intent(in) :: e(nl)
   real(R8P)                :: o
   integer(I4P)             :: j
   o = -1._R8P
   do j=nl-1, 1, -1
      if (e(j+1) > 1.e2_R8P*epsilon(1._R8P)) then
         o = log(e(j)/e(j+1))/log(2._R8P) ; return
      endif
   enddo
   endfunction finest
endprogram seth_test_interpolation
