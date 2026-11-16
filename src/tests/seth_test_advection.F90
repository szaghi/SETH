!< SETH test: the kernels inside a finite-difference scheme: derivative convergence and SSP-RK3 linear advection.
program seth_test_advection
!< SETH test: the kernels inside a finite-difference scheme: derivative convergence and SSP-RK3 linear advection.
!<
!< Finite-difference WENO for u_t + a u_x = 0 (point values on nodes, the flux reconstructed at the faces by the
!< kernels, du_i/dt = -(F_{i+1/2} - F_{i-1/2})/h; a > 0 uses the left state, a < 0 the right one).
!<
!< 1. Derivative convergence (no time error, no critical point): D = (F_{i+1/2} - F_{i-1/2})/h of u = g(x), the generic
!<    function of seth_test_functions, at x = 0.3 against g'(0.3), h halved from h0; every reconstruction scheme, both signs of a; the observed order on the finest
!<    pair above round-off must reach the order of the scheme's linear background minus 0.5.
!< 2. TENO keeps its linear scheme on smooth data (Fu, Hu & Adams 2016, via llm-wiki `targeted-eno.md`): one period
!<    of a Gaussian pulse exp(-50 x^2) on [-1, 1], 200 nodes, SSP-RK3, CFL 0.4: the TENO solution equals the solution
!<    of the linear scheme sum_k d_k q_k to round-off (1e-13).
!< 3. The Jiang-Shu composite profile (Jiang & Shu 1996; llm-wiki `jiang-shu-composite-advection.md`): 200 nodes,
!<    t = 8 (four periods), SSP-RK3, CFL 0.4. Asserted: finite, over- and undershoot below 0.1 (an oscillating or
!<    unstable scheme exceeds it; the values are printed), and L1(WENO5-Z) < L1(WENO5-JS), the ranking of Borges et
!<    al. 2008 (Table 4).
use seth
use seth_test_parameters, only : TOL_TENO_LINEAR
use seth_test_functions, only : gfun, gder
use seth_test_schemes
use seth_z_teno_tables
implicit none
integer(I4P), parameter :: NH = 6_I4P
real(R_P),    parameter :: PI = acos(-1._R_P)
real(R_P)               :: o(2), l1_js5, l1_z5, l1, umax, umin
integer(I4P)            :: id, k
logical                 :: passed
character(16)           :: name

passed = .true.
! 1. derivative convergence
print '(A16,2A10)', 'scheme', 'a > 0', 'a < 0'
do id=1, NSCHEMES
   if (SCHEMES(id)%interp) cycle ! the finite-difference flux setting reconstructs: interpolations elsewhere
   call derivative_orders(id, o)
   print '(A16,2F10.2)', SCHEMES(id)%name, o
   if (any(o < SCHEMES(id)%order - 0.5_R_P)) then
      print '(A,F5.1)', '   ^ below the order of the linear scheme ', SCHEMES(id)%order ; passed = .false.
   endif
enddo
! 2. TENO = its linear scheme on a smooth pulse
do id=1, NSCHEMES
   if (SCHEMES(id)%family /= FAM_TENO .or. SCHEMES(id)%interp) cycle
   l1 = teno_vs_linear(id)
   print '(A16,A,ES9.2)', SCHEMES(id)%name, ': max |TENO - linear| after one period of a Gaussian: ', l1
   if (l1 > TOL_TENO_LINEAR) passed = .false.
enddo
! 3. composite profile
l1_js5 = -1._R_P ; l1_z5 = -1._R_P
print '(A16,3A12)', 'scheme', 'L1', 'max u - 1', 'min u'
do id=1, NSCHEMES
   name = SCHEMES(id)%name
   if (.not.any(name == [character(16) :: 'weno5_js', 'weno7_js', 'wenoz5', 'wenoz7', 'teno5', 'teno6', &
                                          'teno6_opt'])) cycle
   call composite(id, l1, umax, umin)
   print '(A16,3ES12.3)', name, l1, umax - 1._R_P, umin
   if (.not.(l1 < 1._R_P .and. umax - 1._R_P < 0.1_R_P .and. umin > -0.1_R_P)) then
      print '(A)', '   ^ unbounded, oscillating or not finite' ; passed = .false.
   endif
   if (name == 'weno5_js') l1_js5 = l1
   if (name == 'wenoz5')   l1_z5 = l1
enddo
if (.not.(l1_z5 < l1_js5)) then
   print '(A)', 'WENO5-Z is not less dissipative than WENO5-JS on the composite profile (Borges et al. 2008)'
   passed = .false.
endif
if (.not.passed) error stop 'advection test failed'
print '(A)', 'seth_test_advection: ok'

contains
   subroutine derivative_orders(id, o)
   !< Orders of the finite-difference derivative of g (seth_test_functions) at 0.3, both signs of a.
   integer(I4P), intent(in)  :: id
   real(R_P),    intent(out) :: o(2)
   real(R_P)                 :: up(VLO:VHI), um(VLO:VHI), h, h0, err(2,NH), flo(NH), rp1, rm1, rp0, rm0
   integer(I4P)              :: j, m, side
   h0 = 0.2_R_P*2._R_P**max(0, SCHEMES(id)%s - 3)
   do j=1, NH
      h = h0/2._R_P**(j-1)
      do m=VLO, VHI ! face i+1/2: cell m is node i+m; face i-1/2: node i+m-1 (x_i = 0.3)
         up(m) = gfun(0.3_R_P + m*h)
         um(m) = gfun(0.3_R_P + (m-1)*h)
      enddo
      call reconstruct(id, up, up, rp1, rm1)
      call reconstruct(id, um, um, rp0, rm0)
      err(1,j) = abs((rp1 - rp0)/h - gder(0.3_R_P))
      err(2,j) = abs((rm1 - rm0)/h - gder(0.3_R_P))
      flo(j) = 2.e3_R_P*epsilon(1._R_P)/h ! |g| < 2 near 0.3
   enddo
   do side=1, 2
      o(side) = 99._R_P
      do j=NH-1, 1, -1
         if (err(side,j+1) > flo(j+1)) then
            o(side) = log(err(side,j)/err(side,j+1))/log(2._R_P) ; exit
         endif
      enddo
   enddo
   endsubroutine derivative_orders

   subroutine rhs(id, linear, n, h, u, dudt)
   !< du/dt = -(F_{i+1/2} - F_{i-1/2})/h, a = 1, periodic; `linear`: the TENO scheme's linear background.
   integer(I4P), intent(in)  :: id, n
   logical,      intent(in)  :: linear
   real(R_P),    intent(in)  :: h, u(0:n-1)
   real(R_P),    intent(out) :: dudt(0:n-1)
   real(R_P)                 :: f(0:n-1), v(VLO:VHI), rp, rm
   integer(I4P)              :: i, m
   do i=0, n-1 ! face i+1/2
      do m=VLO, VHI
         v(m) = u(modulo(i+m, n))
      enddo
      if (linear) then
         f(i) = linear_teno(id, v)
      else
         call reconstruct(id, v, v, rp, rm)
         f(i) = rp
      endif
   enddo
   do i=0, n-1
      dudt(i) = -(f(i) - f(modulo(i-1, n)))/h
   enddo
   endsubroutine rhs

   subroutine advect(id, linear, n, h, tend, u)
   !< SSP-RK3 (Shu & Osher) to tend, CFL 0.4 (the last step shortened).
   integer(I4P), intent(in)    :: id, n
   logical,      intent(in)    :: linear
   real(R_P),    intent(in)    :: h, tend
   real(R_P),    intent(inout) :: u(0:n-1)
   real(R_P)                   :: t, dt, u1(0:n-1), u2(0:n-1), l(0:n-1)
   integer(I4P)                :: nsteps, s
   nsteps = ceiling(tend/(0.4_R_P*h))
   dt = tend/nsteps
   t = 0._R_P
   do s=1, nsteps
      call rhs(id, linear, n, h, u, l)  ; u1 = u + dt*l
      call rhs(id, linear, n, h, u1, l) ; u2 = 0.75_R_P*u + 0.25_R_P*(u1 + dt*l)
      call rhs(id, linear, n, h, u2, l) ; u = u/3._R_P + 2._R_P/3._R_P*(u2 + dt*l)
   enddo
   endsubroutine advect

   function teno_vs_linear(id) result(d)
   !< max |u_TENO - u_linear| after one period of exp(-50 x^2).
   integer(I4P), intent(in) :: id
   real(R_P)                :: d
   integer(I4P), parameter  :: n = 200_I4P
   real(R_P)                :: u(0:n-1), w(0:n-1), h, x
   integer(I4P)             :: i
   h = 2._R_P/n
   do i=0, n-1
      x = -1._R_P + i*h
      u(i) = exp(-50._R_P*x*x)
   enddo
   w = u
   call advect(id, .false., n, h, 2._R_P, u)
   call advect(id, .true.,  n, h, 2._R_P, w)
   d = maxval(abs(u - w))
   endfunction teno_vs_linear

   function linear_teno(id, v) result(r)
   !< The linear background of the TENO scheme id, sum_k d_k q_k / sum d, left-biased.
   integer(I4P), intent(in) :: id
   real(R_P),    intent(in) :: v(VLO:VHI)
   real(R_P)                :: r
   select case(trim(SCHEMES(id)%name))
   case('teno5')       ; r = lin_of(TENO5_NC, TENO5_MAXW, TENO5_CLO, TENO5_CW, TENO5_P, TENO5_D, v)
   case('teno6')       ; r = lin_of(TENO6_NC, TENO6_MAXW, TENO6_CLO, TENO6_CW, TENO6_P, TENO6_D, v)
   case('teno6_opt')   ; r = lin_of(TENO6_OPT_NC, TENO6_OPT_MAXW, TENO6_OPT_CLO, TENO6_OPT_CW, TENO6_OPT_P, &
                                    TENO6_OPT_D, v)
   case('teno6_neno6') ; r = lin_of(TENO6_NENO6_NC, TENO6_NENO6_MAXW, TENO6_NENO6_CLO, TENO6_NENO6_CW, &
                                    TENO6_NENO6_P, TENO6_NENO6_D, v)
   case('teno5_opt')    ; r = lin_of(TENO5_OPT_NC, TENO5_OPT_MAXW, TENO5_OPT_CLO, TENO5_OPT_CW, TENO5_OPT_P, &
                                    TENO5_OPT_D, v)
   case('teno5_az')     ; r = lin_of(TENO5_AZ_NC, TENO5_AZ_MAXW, TENO5_AZ_CLO, TENO5_AZ_CW, TENO5_AZ_P, &
                                    TENO5_AZ_D, v)
   case('teno8')        ; r = lin_of(TENO8_NC, TENO8_MAXW, TENO8_CLO, TENO8_CW, TENO8_P, &
                                    TENO8_D, v)
   case('teno5_a')     ; r = lin_of(TENO5_NC, TENO5_MAXW, TENO5_CLO, TENO5_CW, TENO5_P, TENO5_D, v)
   case('teno8_a')     ; r = lin_of(TENO8_NC, TENO8_MAXW, TENO8_CLO, TENO8_CW, TENO8_P, TENO8_D, v)
   case default        ; error stop 'linear_teno: not a TENO scheme'
   endselect
   endfunction linear_teno

   function lin_of(nc, maxw, clo, cw, cp, d, v) result(r)
   integer(I4P), intent(in) :: nc, maxw, clo(0:nc-1), cw(0:nc-1)
   real(R_P),    intent(in) :: cp(0:maxw-1,0:nc-1), d(0:nc-1), v(VLO:VHI)
   real(R_P)                :: r, q, ds
   integer(I4P)             :: k, a
   ds = 0._R_P
   do k=0, nc-1
      ds = ds + d(k)
   enddo
   r = 0._R_P
   do k=0, nc-1
      q = 0._R_P
      do a=0, cw(k)-1
         q = q + cp(a,k)*v(clo(k)+a)
      enddo
      r = r + (d(k)/ds)*q
   enddo
   endfunction lin_of

   subroutine composite(id, l1, umax, umin)
   !< The Jiang-Shu composite profile after four periods: L1 error, max, min.
   integer(I4P), intent(in)  :: id
   real(R_P),    intent(out) :: l1, umax, umin
   integer(I4P), parameter   :: n = 200_I4P
   real(R_P)                 :: u(0:n-1), u0(0:n-1), h
   integer(I4P)              :: i
   h = 2._R_P/n
   do i=0, n-1
      u0(i) = jiang_shu(-1._R_P + i*h)
   enddo
   u = u0
   call advect(id, .false., n, h, 8._R_P, u)
   l1 = h*sum(abs(u - u0))
   umax = maxval(u) ; umin = minval(u)
   endsubroutine composite

   pure function jiang_shu(x) result(u)
   !< The composite profile, constants as llm-wiki `jiang-shu-composite-advection.md` prints them.
   real(R_P), intent(in) :: x
   real(R_P)             :: u
   real(R_P), parameter  :: A = 0.5_R_P, Z = -0.7_R_P, DG = 0.005_R_P, AF = 10._R_P
   real(R_P)             :: bg
   bg = log(2._R_P)/(36._R_P*DG*DG)
   if (x >= -0.8_R_P .and. x <= -0.6_R_P) then
      u = (g(x, bg, Z - DG) + g(x, bg, Z + DG) + 4._R_P*g(x, bg, Z))/6._R_P
   elseif (x >= -0.4_R_P .and. x <= -0.2_R_P) then
      u = 1._R_P
   elseif (x >= 0._R_P .and. x <= 0.2_R_P) then
      u = 1._R_P - abs(10._R_P*(x - 0.1_R_P))
   elseif (x >= 0.4_R_P .and. x <= 0.6_R_P) then
      u = (ef(x, AF, A - DG) + ef(x, AF, A + DG) + 4._R_P*ef(x, AF, A))/6._R_P
   else
      u = 0._R_P
   endif
   endfunction jiang_shu

   pure function g(x, b, c) result(y)
   real(R_P), intent(in) :: x, b, c
   real(R_P)             :: y
   y = exp(-b*(x - c)**2)
   endfunction g

   pure function ef(x, al, c) result(y)
   real(R_P), intent(in) :: x, al, c
   real(R_P)             :: y
   y = sqrt(max(1._R_P - al*al*(x - c)**2, 0._R_P))
   endfunction ef
endprogram seth_test_advection
