!< SETH verification: a minimal 1-D Euler solver around the SETH kernels (test harness, not library API).
module seth_vv_euler
!< SETH verification: a minimal 1-D Euler solver around the SETH kernels (test harness, not library API).
!<
!< Ideal gas, gamma = 1.4, conservative variables q = (rho, rho u, E) on N cells of width h with NG = 5 ghost cells
!< (WENO-9), point values at the cell centres (conservative finite differences: the fluxes at the faces, du/dt =
!< -(F_{i+1/2} - F_{i-1/2})/h). Two space discretizations, the two ways SETH's kernels are used:
!<
!< - PATH_SPLIT, the reconstruction kernels: flux splitting (Shu & Osher 1989; Jiang & Shu 1996) on characteristic
!<   variables. At each face the Roe-averaged left and right eigenvectors L, R project the stencil's states and fluxes,
!<   w = L q, g = L f; global Lax-Friedrichs splitting per field, g+- = (g +- alpha_k w)/2, alpha_k = max |lambda_k|
!<   over the domain times 1.1; g+ is reconstructed with the left bias (vp), g- with the right bias (vm), and
!<   F = R (rp + rm). The 10% margin keeps the split fluxes free of degenerate critical points: with alpha = max|lambda|
!<   exactly, d g-/dx = (lambda - alpha) dw/dx vanishes to second order where |lambda| peaks, and TENO drops to second
!<   order at second-order critical points (seth_test_critical_points);
!< - PATH_INTERP, the interpolation kernels: the alternative finite-difference formulation (Shu & Osher 1988; Jiang,
!<   Shu & Zhang 2013). The characteristic variables w = L q are interpolated to the face (left and right states), a
!<   Riemann solver (HLLC, Toro et al. 1994 with Davis-Einfeldt speeds) gives F0, and the flux is corrected to high
!<   order with the even derivatives of the flux, F = F0 - h^2/24 f_xx + 7 h^4/5760 f_xxxx, f_xx h^2 and f_xxxx h^4
!<   from the six node fluxes around the face (exact coefficients: O(h^4) and O(h^2), so the correction errs by
!<   O(h^6)). A face state with rho <= 0 or p <= 0 falls back to the first-order states (counted in `nfallback`).
!<
!< Positivity: with `limiter` on, the flux of every face is blended with the first-order Lax-Friedrichs flux,
!< F = theta F_high + (1 - theta) F_LF, theta in [0, 1] as large as keeps rho and p of the two one-sided states
!< q_i - 2 lambda F and q_{i+1} + 2 lambda F above 1e-13 (lambda = dt/h): the positivity-preserving flux limiter of
!< Hu, Adams & Shu (2013) for finite-difference WENO, theta by bisection. Each stage is then a convex combination of
!< positive states when dt max(|u| + c) <= h/2 (CFL <= 0.5), and SSP-RK3 keeps that. Every face with theta < 1 is
!< counted in `nlimited`: the cases assert that the limiter never acts on smooth flows, and report how often each
!< scheme needs it on the stress tests.
!<
!< Boundaries: periodic, transmissive (zeroth-order extrapolation) or reflecting walls. Time: SSP-RK3 (Shu & Osher
!< 1988).
use seth
use seth_test_schemes, only : reconstruct, VLO, VHI
implicit none
private
public :: GAMMA, NG, PATH_SPLIT, PATH_INTERP, BC_PERIODIC, BC_TRANSMISSIVE, BC_REFLECTIVE
public :: prim2cons, cons2prim, euler_dt, euler_step, nfallback, corrections, limiter, nlimited, positive_state

real(R8P),    parameter :: GAMMA = 1.4_R8P
integer(I4P), parameter :: NG = 5_I4P
integer(I4P), parameter :: PATH_SPLIT = 1, PATH_INTERP = 2
integer(I4P), parameter :: BC_PERIODIC = 1, BC_TRANSMISSIVE = 2, BC_REFLECTIVE = 3
! f'' h^2 and f'''' h^4 at the face from the node fluxes of cells -2..3 (exact rationals)
real(R8P), parameter :: C2(-2:3) = [-5._R8P/48._R8P, 13._R8P/16._R8P, -17._R8P/24._R8P, -17._R8P/24._R8P, &
                                    13._R8P/16._R8P, -5._R8P/48._R8P]
real(R8P), parameter :: C4(-2:3) = [0.5_R8P, -1.5_R8P, 1._R8P, 1._R8P, -1.5_R8P, 0.5_R8P]
integer(I4P) :: nfallback = 0 !< Interpolated face states replaced by first-order ones (non-positive rho or p).
logical      :: corrections = .true. !< PATH_INTERP: add the flux corrections (off: the second-order flux F0).
logical      :: limiter = .false.    !< Positivity-preserving flux limiter (Hu, Adams & Shu 2013).
integer(I4P) :: nlimited = 0         !< Faces where the limiter blended in the first-order flux (theta < 1).

contains
   pure function prim2cons(w) result(q)
   !< (rho, u, p) -> (rho, rho u, E).
   real(R8P), intent(in) :: w(3)
   real(R8P)             :: q(3)
   q = [w(1), w(1)*w(2), w(3)/(GAMMA - 1._R8P) + 0.5_R8P*w(1)*w(2)**2]
   endfunction prim2cons

   pure function cons2prim(q) result(w)
   !< (rho, rho u, E) -> (rho, u, p).
   real(R8P), intent(in) :: q(3)
   real(R8P)             :: w(3)
   w(1) = q(1) ; w(2) = q(2)/q(1) ; w(3) = (GAMMA - 1._R8P)*(q(3) - 0.5_R8P*q(2)*w(2))
   endfunction cons2prim

   pure function flux(q) result(f)
   real(R8P), intent(in) :: q(3)
   real(R8P)             :: f(3), w(3)
   w = cons2prim(q)
   f = [q(2), q(2)*w(2) + w(3), (q(3) + w(3))*w(2)]
   endfunction flux

   pure subroutine eigen(ql, qr, l, r, lam)
   !< Left and right eigenvectors of the flux Jacobian at the Roe average of ql and qr, and its eigenvalues.
   real(R8P), intent(in)  :: ql(3), qr(3)
   real(R8P), intent(out) :: l(3,3), r(3,3), lam(3)
   real(R8P)              :: wl(3), wr(3), sl, sr, u, hh, c, b1, b2
   wl = cons2prim(ql) ; wr = cons2prim(qr)
   sl = sqrt(wl(1)) ; sr = sqrt(wr(1))
   u = (sl*wl(2) + sr*wr(2))/(sl + sr)
   hh = (sl*(ql(3) + wl(3))/wl(1) + sr*(qr(3) + wr(3))/wr(1))/(sl + sr)
   c = sqrt(max((GAMMA - 1._R8P)*(hh - 0.5_R8P*u*u), tiny(1._R8P)))
   lam = [u - c, u, u + c]
   r(:,1) = [1._R8P, u - c, hh - u*c]
   r(:,2) = [1._R8P, u,     0.5_R8P*u*u]
   r(:,3) = [1._R8P, u + c, hh + u*c]
   b2 = (GAMMA - 1._R8P)/(c*c) ; b1 = 0.5_R8P*b2*u*u
   l(1,:) = 0.5_R8P*[b1 + u/c, -(b2*u + 1._R8P/c), b2]
   l(2,:) = [1._R8P - b1, b2*u, -b2]
   l(3,:) = 0.5_R8P*[b1 - u/c, -(b2*u - 1._R8P/c), b2]
   endsubroutine eigen

   pure function hllc(ql, qr) result(f)
   !< HLLC flux (Toro, Spruce & Speares 1994), Davis-Einfeldt wave speed estimates.
   real(R8P), intent(in) :: ql(3), qr(3)
   real(R8P)             :: f(3), wl(3), wr(3), cl, cr, sl, sr, sm, qs(3)
   wl = cons2prim(ql) ; wr = cons2prim(qr)
   cl = sqrt(GAMMA*wl(3)/wl(1)) ; cr = sqrt(GAMMA*wr(3)/wr(1))
   sl = min(wl(2) - cl, wr(2) - cr) ; sr = max(wl(2) + cl, wr(2) + cr)
   sm = (wr(3) - wl(3) + ql(2)*(sl - wl(2)) - qr(2)*(sr - wr(2)))/(wl(1)*(sl - wl(2)) - wr(1)*(sr - wr(2)))
   if (sl >= 0._R8P) then
      f = flux(ql)
   elseif (sr <= 0._R8P) then
      f = flux(qr)
   elseif (sm >= 0._R8P) then
      qs = wl(1)*(sl - wl(2))/(sl - sm)*[1._R8P, sm, ql(3)/wl(1) + (sm - wl(2))*(sm + wl(3)/(wl(1)*(sl - wl(2))))]
      f = flux(ql) + sl*(qs - ql)
   else
      qs = wr(1)*(sr - wr(2))/(sr - sm)*[1._R8P, sm, qr(3)/wr(1) + (sm - wr(2))*(sm + wr(3)/(wr(1)*(sr - wr(2))))]
      f = flux(qr) + sr*(qs - qr)
   endif
   endfunction hllc

   subroutine ghosts(q, n, bc)
   !< Fill the NG ghost cells of each end.
   integer(I4P), intent(in)    :: n, bc
   real(R8P),    intent(inout) :: q(3,1-NG:n+NG)
   integer(I4P)                :: g
   do g=1, NG
      select case(bc)
      case(BC_PERIODIC)
         q(:,1-g) = q(:,n+1-g) ; q(:,n+g) = q(:,g)
      case(BC_TRANSMISSIVE)
         q(:,1-g) = q(:,1) ; q(:,n+g) = q(:,n)
      case(BC_REFLECTIVE)
         q(:,1-g) = q(:,g) ; q(2,1-g) = -q(2,g)
         q(:,n+g) = q(:,n+1-g) ; q(2,n+g) = -q(2,n+1-g)
      endselect
   enddo
   endsubroutine ghosts

   function euler_dt(q, n, h, cfl) result(dt)
   !< dt = cfl h / max(|u| + c).
   integer(I4P), intent(in) :: n
   real(R8P),    intent(in) :: q(3,1-NG:n+NG), h, cfl
   real(R8P)                :: dt, w(3), smax
   integer(I4P)             :: i
   smax = 0._R8P
   do i=1, n
      w = cons2prim(q(:,i))
      smax = max(smax, abs(w(2)) + sqrt(GAMMA*w(3)/w(1)))
   enddo
   dt = cfl*h/smax
   endfunction euler_dt

   subroutine residual(id, path, q, n, h, dt, bc, res, par)
   !< res = -(F_{i+1/2} - F_{i-1/2})/h for i = 1..n (dt: for the positivity limiter only).
   integer(I4P), intent(in)           :: id, path, n, bc
   real(R8P),    intent(inout)        :: q(3,1-NG:n+NG)
   real(R8P),    intent(in)           :: h, dt
   real(R8P),    intent(out)          :: res(3,n)
   real(R8P),    intent(in), optional :: par
   real(R8P)                          :: fhat(3,0:n), fn(3,1-NG:n+NG), l(3,3), r(3,3), lam(3), alpha(3)
   real(R8P)                          :: w(3,VLO:VHI), g(3,VLO:VHI), vp(VLO:VHI), vm(VLO:VHI), rp, rm, fc(3)
   real(R8P)                          :: wl(3), wr(3), ql(3), qr(3), pl(3), pr(3)
   integer(I4P)                       :: i, m, k
   call ghosts(q, n, bc)
   do i=1-NG, n+NG
      fn(:,i) = flux(q(:,i))
   enddo
   alpha = 0._R8P
   if (path == PATH_SPLIT) then
      do i=0, n
         call eigen(q(:,i), q(:,i+1), l, r, lam)
         alpha = max(alpha, abs(lam))
      enddo
      do i=1, n ! the cell eigenvalues too: the splitting must dominate every |lambda| on the stencils
         pl = cons2prim(q(:,i))
         alpha = max(alpha, abs(pl(2)) + sqrt(GAMMA*pl(3)/pl(1)))
      enddo
      alpha = 1.1_R8P*alpha ! margin: alpha = max|lambda| exactly makes g- = (lambda - alpha) w/2 a second-order
                            ! critical point where |lambda| peaks, where TENO drops to second order
   endif
   do i=0, n ! face i+1/2: stencil cell m is cell i+m
      call eigen(q(:,i), q(:,i+1), l, r, lam)
      do m=VLO, VHI
         w(:,m) = matmul(l, q(:,i+m))
         g(:,m) = matmul(l, fn(:,i+m))
      enddo
      if (path == PATH_SPLIT) then
         do k=1, 3
            vp = 0.5_R8P*(g(k,:) + alpha(k)*w(k,:))
            vm = 0.5_R8P*(g(k,:) - alpha(k)*w(k,:))
            call reconstruct(id, vp, vm, rp, rm, par=par)
            fc(k) = rp + rm
         enddo
         fhat(:,i) = matmul(r, fc)
      else
         do k=1, 3
            vp = w(k,:)
            call reconstruct(id, vp, vp, wl(k), wr(k), par=par)
         enddo
         ql = matmul(r, wl) ; qr = matmul(r, wr)
         pl = cons2prim(ql) ; pr = cons2prim(qr)
         if (.not.(pl(1) > 0._R8P .and. pl(3) > 0._R8P .and. pr(1) > 0._R8P .and. pr(3) > 0._R8P)) then
            ql = q(:,i) ; qr = q(:,i+1) ; nfallback = nfallback + 1
         endif
         fhat(:,i) = hllc(ql, qr)
         if (corrections) then
         do m=-2, 3
            fhat(:,i) = fhat(:,i) - (C2(m)/24._R8P - 7._R8P*C4(m)/5760._R8P)*fn(:,i+m)
         enddo
         endif
      endif
   enddo
   if (limiter) call limit(q, n, h, dt, fhat)
   do i=1, n
      res(:,i) = -(fhat(:,i) - fhat(:,i-1))/h
   enddo
   endsubroutine residual

   subroutine limit(q, n, h, dt, fhat)
   !< Hu, Adams & Shu (2013): blend each face flux with the Lax-Friedrichs flux so that q_i - 2 (dt/h) F and
   !< q_{i+1} + 2 (dt/h) F keep rho, p >= 1e-13.
   integer(I4P), intent(in)    :: n
   real(R8P),    intent(in)    :: q(3,1-NG:n+NG), h, dt
   real(R8P),    intent(inout) :: fhat(3,0:n)
   real(R8P),    parameter     :: FLOOR = 1.e-13_R8P
   real(R8P)                   :: lam, a, w(3), flf(3), lo, hi, th
   integer(I4P)                :: i, it
   lam = 2._R8P*dt/h
   a = 0._R8P
   do i=0, n+1
      w = cons2prim(q(:,i))
      a = max(a, abs(w(2)) + sqrt(max(GAMMA*w(3)/w(1), 0._R8P)))
   enddo
   do i=0, n
      if (ok(fhat(:,i))) cycle
      flf = 0.5_R8P*(flux(q(:,i)) + flux(q(:,i+1))) - 0.5_R8P*a*(q(:,i+1) - q(:,i))
      lo = 0._R8P ; hi = 1._R8P
      do it=1, 50
         th = 0.5_R8P*(lo + hi)
         if (ok(th*fhat(:,i) + (1._R8P - th)*flf)) then ; lo = th ; else ; hi = th ; endif
      enddo
      fhat(:,i) = lo*fhat(:,i) + (1._R8P - lo)*flf
      nlimited = nlimited + 1
   enddo
   contains
      pure logical function ok(f)
      real(R8P), intent(in) :: f(3)
      real(R8P)             :: ql(3), qr(3)
      ql = q(:,i) - lam*f ; qr = q(:,i+1) + lam*f
      ok = ql(1) >= FLOOR .and. qr(1) >= FLOOR .and. &
           (GAMMA - 1._R8P)*(ql(3) - 0.5_R8P*ql(2)**2/ql(1)) >= FLOOR .and. &
           (GAMMA - 1._R8P)*(qr(3) - 0.5_R8P*qr(2)**2/qr(1)) >= FLOOR
      endfunction ok
   endsubroutine limit

   subroutine euler_step(id, path, q, n, h, dt, bc, ok, par)
   !< One SSP-RK3 step; `ok` false (and q left at the failed stage) when a stage has a cell with rho <= 0 or p <= 0:
   !< the next residual would take the square root of a negative pressure (a floating-point trap in the debug builds).
   integer(I4P), intent(in)           :: id, path, n, bc
   real(R8P),    intent(inout)        :: q(3,1-NG:n+NG)
   real(R8P),    intent(in)           :: h, dt
   logical,      intent(out)          :: ok
   real(R8P),    intent(in), optional :: par
   real(R8P)                          :: q1(3,1-NG:n+NG), q2(3,1-NG:n+NG), res(3,n)
   call residual(id, path, q, n, h, dt, bc, res, par)
   q1 = q ; q1(:,1:n) = q(:,1:n) + dt*res
   ok = positive_state(q1, n) ; if (.not.ok) then ; q = q1 ; return ; endif
   call residual(id, path, q1, n, h, dt, bc, res, par)
   q2 = q ; q2(:,1:n) = 0.75_R8P*q(:,1:n) + 0.25_R8P*(q1(:,1:n) + dt*res)
   ok = positive_state(q2, n) ; if (.not.ok) then ; q = q2 ; return ; endif
   call residual(id, path, q2, n, h, dt, bc, res, par)
   q(:,1:n) = q(:,1:n)/3._R8P + 2._R8P/3._R8P*(q2(:,1:n) + dt*res)
   ok = positive_state(q, n)
   endsubroutine euler_step

   pure function positive_state(q, n) result(ok)
   !< rho > 0 and p > 0 in every interior cell (no division: rho is tested first, p as rho E - (rho u)^2/2 > 0).
   integer(I4P), intent(in) :: n
   real(R8P),    intent(in) :: q(3,1-NG:n+NG)
   logical                  :: ok
   integer(I4P)             :: i
   ok = .true.
   do i=1, n
      if (.not.(q(1,i) > 0._R8P)) then ; ok = .false. ; return ; endif
      if (.not.(q(1,i)*q(3,i) - 0.5_R8P*q(2,i)**2 > 0._R8P)) then ; ok = .false. ; return ; endif
   enddo
   endfunction positive_state
endmodule seth_vv_euler
