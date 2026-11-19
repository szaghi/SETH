!< Tutorial 8: positivity. A double rarefaction that opens a vacuum, with and without the positivity limiter.
program t08_positivity
use seth
use selah, only : stepper, box, box_interior, SCHEME_SSP_33, FORM_SHU_OSHER_2R   ! the time integrator: SELAH's
implicit none
real(R_P), parameter :: GAMMA = 1.4_R_P
integer,   parameter :: N = 400, NG = 3
real(R_P)            :: q(3,N), h, x, t, rmin, pmin, wp(3)
integer              :: i, k, nlim
logical              :: limit, ok

h = 1._R_P/N
print '(A)', 'Woodward and Colella blast waves: p = 1000, 0.01, 100 on [0, 0.1], [0.1, 0.9], [0.9, 1], reflecting walls'
print '(A)', 'TENO5, characteristic splitting, 400 points, to t = 0.038'
do k=1, 2
   limit = k == 2
   do i=1, N
      x = (i - 0.5_R_P)*h
      q(:,i) = prim2cons([1._R_P, 0._R_P, merge(1000._R_P, merge(0.01_R_P, 100._R_P, x < 0.9_R_P), x < 0.1_R_P)])
   enddo
   call evolve(q, 0.038_R_P, limit, t, ok, nlim)
   if (.not.ok) then
      print '(A,F8.5)', '  without the limiter: a negative density or pressure at t = ', t
   else
      rmin = huge(1._R_P) ; pmin = huge(1._R_P)
      do i=1, N
         wp = cons2prim(q(:,i)) ; rmin = min(rmin, wp(1)) ; pmin = min(pmin, wp(3))
      enddo
      print '(A,A,ES9.2,A,ES9.2)', merge('  with the limiter:   ', '  without the limiter:', limit), &
                                   ' reached t = 0.038, min rho ', rmin, ', min p ', pmin
      print '(A,I0)', '  face fluxes limited (theta < 1), over all stages: ', nlim
   endif
enddo

contains
   pure function prim2cons(w) result(q)
   real(R_P), intent(in) :: w(3)                ! rho, u, p
   real(R_P)             :: q(3)                ! rho, rho u, E
   q = [w(1), w(1)*w(2), w(3)/(GAMMA - 1) + 0.5_R_P*w(1)*w(2)**2]
   endfunction prim2cons

   pure function cons2prim(q) result(w)
   real(R_P), intent(in) :: q(3)
   real(R_P)             :: w(3)
   w = [q(1), q(2)/q(1), (GAMMA - 1)*(q(3) - 0.5_R_P*q(2)**2/q(1))]
   endfunction cons2prim

   pure function flux(q) result(f)
   real(R_P), intent(in) :: q(3)
   real(R_P)             :: f(3), w(3)
   w = cons2prim(q)
   f = [q(2), q(2)*w(2) + w(3), (q(3) + w(3))*w(2)]
   endfunction flux

   pure subroutine eigenvectors(ql, qr, l, r)
   !< The left (rows of l) and right (columns of r) eigenvectors of the flux Jacobian at the Roe average of ql, qr.
   real(R_P), intent(in)  :: ql(3), qr(3)
   real(R_P), intent(out) :: l(3,3), r(3,3)
   real(R_P)              :: wl(3), wr(3), sl, sr, u, hh, c, b1, b2
   wl = cons2prim(ql) ; wr = cons2prim(qr)
   sl = sqrt(wl(1)) ; sr = sqrt(wr(1))
   u  = (sl*wl(2) + sr*wr(2))/(sl + sr)
   hh = (sl*(ql(3) + wl(3))/wl(1) + sr*(qr(3) + wr(3))/wr(1))/(sl + sr)   ! enthalpy
   c  = sqrt((GAMMA - 1)*(hh - 0.5_R_P*u*u))
   r(:,1) = [1._R_P, u - c, hh - u*c]
   r(:,2) = [1._R_P, u, 0.5_R_P*u*u]
   r(:,3) = [1._R_P, u + c, hh + u*c]
   b1 = (GAMMA - 1)/(c*c) ; b2 = 0.5_R_P*b1*u*u
   l(1,:) = 0.5_R_P*[b2 + u/c, -b1*u - 1/c, b1]
   l(2,:) = [1 - b2, b1*u, -b1]
   l(3,:) = 0.5_R_P*[b2 - u/c, -b1*u + 1/c, b1]
   endsubroutine eigenvectors

   subroutine residual(q, dt, limit, res, nlim)
   !< res = -(F_{i+1/2} - F_{i-1/2})/h; with limit, every face flux limited by seth_positivity_euler.
   real(R_P), intent(in)    :: q(3,N), dt
   logical,   intent(in)    :: limit
   real(R_P), intent(out)   :: res(3,N)
   integer,   intent(inout) :: nlim
   real(R_P)                :: qg(3,1-NG:N+NG), fg(3,1-NG:N+NG), wc(3,-2:3), gc(3,-2:3), l(3,3), r(3,3)
   real(R_P)                :: alpha(3), fc(3), fhat(3,0:N), flow(3), wp(3), rp, rm, c, theta, a
   integer                  :: i, k, m
   qg(:,1:N) = q
   do i=1, NG                                   ! reflecting walls: mirrored cells, the momentum reversed
      qg(:,1-i) = q(:,i) ; qg(2,1-i) = -q(2,i)
      qg(:,N+i) = q(:,N+1-i) ; qg(2,N+i) = -q(2,N+1-i)
   enddo
   do i=1-NG, N+NG
      fg(:,i) = flux(qg(:,i))
   enddo
   alpha = 0
   do i=1, N
      wp = cons2prim(q(:,i)) ; c = sqrt(GAMMA*wp(3)/wp(1))
      alpha = max(alpha, abs([wp(2) - c, wp(2), wp(2) + c]))
   enddo
   a = maxval(alpha)                            ! the first-order flux's alpha: max(|u| + c)
   alpha = 1.1_R_P*alpha
   do i=0, N
      call eigenvectors(qg(:,i), qg(:,i+1), l, r)
      do m=-2, 3
         wc(:,m) = matmul(l, qg(:,i+m)) ; gc(:,m) = matmul(l, fg(:,i+m))
      enddo
      do k=1, 3
         call seth_teno5(0.5_R_P*(gc(k,:) + alpha(k)*wc(k,:)), 0.5_R_P*(gc(k,:) - alpha(k)*wc(k,:)), &
                         SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)
         fc(k) = rp + rm
      enddo
      fhat(:,i) = matmul(r, fc)
      if (limit) then
         ! the first-order Lax-Friedrichs flux, positive for 2 dt/h a <= 1; theta blends the two
         flow = 0.5_R_P*(fg(:,i) + fg(:,i+1)) - 0.5_R_P*a*(qg(:,i+1) - qg(:,i))
         call seth_positivity_euler(3, GAMMA, 2*dt/h, SETH_POSITIVITY_EPS, SETH_POSITIVITY_EPS, qg(:,i), qg(:,i+1), &
                                    fhat(:,i), flow, fc, theta)
         fhat(:,i) = fc
         if (theta < 1) nlim = nlim + 1
      endif
   enddo
   res = -(fhat(:,1:N) - fhat(:,0:N-1))/h
   endsubroutine residual

   subroutine evolve(q, tend, limit, t, ok, nlim)
   !< SSP-RK3 to tend, CFL 0.4; ok false (and t the time) at the first stage with rho <= 0 or p <= 0.
   real(R_P), intent(inout) :: q(3,N)
   real(R_P), intent(in)    :: tend
   logical,   intent(in)    :: limit
   real(R_P), intent(out)   :: t
   logical,   intent(out)   :: ok
   integer,   intent(out)   :: nlim
   real(R_P)                :: q0(3,N), dt, smax, wp(3)
   integer                  :: i
   type(stepper)            :: st
   type(box)                :: bx
   real(R_P), pointer       :: r(:,:)
   integer                  :: k
   bx = box_interior(lb=[1, 1], ub=[3, N], ngc=0, components=.true.)
   call st%initialize(SCHEME_SSP_33, FORM_SHU_OSHER_2R, lb=[1, 1], ub=[3, N], components=.true.)
   t = 0 ; nlim = 0 ; ok = .true.
   do while (t < tend)
      smax = 0
      do i=1, N
         wp = cons2prim(q(:,i)) ; smax = max(smax, abs(wp(2)) + sqrt(GAMMA*wp(3)/wp(1)))
      enddo
      dt = min(0.4_R_P*h/smax, tend - t)
      q0 = q                                     ! the state to keep if a stage is not admissible
      call st%begin_step(bx, q)
      do k=1, st%stages()                        ! SELAH's SSP(3,3): the stage state is q itself
         call st%begin_stage(k, bx, dt, q)
         call st%rhs_register_view(r)
         call residual(q, dt, limit, r, nlim)
         call st%end_stage(k, bx, dt, q)
         ok = admissible(q)
         if (.not.ok) then
            q = q0 ; call st%destroy ; return
         endif
      enddo
      call st%finish_step(bx, dt, q)
      t = t + dt
   enddo
   call st%destroy
   endsubroutine evolve

   pure logical function admissible(q)
   !< Every cell with rho > 0 and p > 0.
   real(R_P), intent(in) :: q(3,N)
   integer               :: i
   admissible = .true.
   do i=1, N
      admissible = admissible .and. seth_euler_admissible(3, GAMMA, q(:,i), tiny(1._R_P), tiny(1._R_P))
   enddo
   endfunction admissible
endprogram t08_positivity
