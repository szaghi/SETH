!< SETH verification: a minimal 2-D Euler harness, the 1-D characteristic splitting swept in x and y (issue #6, 4.3).
module seth_vv_euler2d
!< SETH verification: a minimal 2-D Euler harness, the 1-D characteristic splitting swept in x and y (issue #6, 4.3).
!<
!< Test code, fenced by owner decision D8: a uniform Cartesian grid, no AMR, no MPI, no curvilinear grid, no viscous
!< terms. Conservative variables q = (rho, rho u, rho v, E), gamma = 1.4, NG = 6 ghost cells. The residual is the sum
!< of two line residuals: along x, the scheme of seth_vv_euler's PATH_SPLIT (Roe-averaged eigenvectors of the x flux,
!< characteristic variables, global Lax-Friedrichs splitting per field, alpha = 1.1 max |lambda| over the domain,
!< the registry's reconstruct on -5:6); along y, the same code on the data with the two momenta swapped. The registry's
!< rbf_dx is set to h (the RBF kernels' s = alpha dx). SSP-RK3. Periodic or transmissive ends. With the x and y passes
!< the same code, data symmetric about the diagonal stay symmetric bit for bit (the two contributions are summed in
!< either order, and a + b = b + a in floating point).
use penf, only : I4P, R_P
use seth_vv_euler, only : GAMMA, NG
use seth_test_schemes, only : reconstruct, VLO, VHI, rbf_dx
use selah, only : stepper, box, box_interior, SCHEME_SSP_33, FORM_SHU_OSHER_2R ! the time integrator (SELAH #2)
implicit none
private
public :: BC2_PERIODIC, BC2_TRANSMISSIVE, prim2cons2, cons2prim2, euler2_step, euler2_dt, eigen2, positive2

integer(I4P), parameter :: BC2_PERIODIC = 1, BC2_TRANSMISSIVE = 2

contains
   pure function prim2cons2(w) result(q)
   !< (rho, u, v, p) -> (rho, rho u, rho v, E).
   real(R_P), intent(in) :: w(4)
   real(R_P)             :: q(4)
   q = [w(1), w(1)*w(2), w(1)*w(3), w(4)/(GAMMA - 1._R_P) + 0.5_R_P*w(1)*(w(2)**2 + w(3)**2)]
   endfunction prim2cons2

   pure function cons2prim2(q) result(w)
   !< (rho, rho u, rho v, E) -> (rho, u, v, p).
   real(R_P), intent(in) :: q(4)
   real(R_P)             :: w(4)
   w(1) = q(1) ; w(2) = q(2)/q(1) ; w(3) = q(3)/q(1)
   w(4) = (GAMMA - 1._R_P)*(q(4) - 0.5_R_P*(q(2)*w(2) + q(3)*w(3)))
   endfunction cons2prim2

   pure function flux_x(q) result(f)
   !< The x flux.
   real(R_P), intent(in) :: q(4)
   real(R_P)             :: f(4), w(4)
   w = cons2prim2(q)
   f = [q(2), q(2)*w(2) + w(4), q(3)*w(2), (q(4) + w(4))*w(2)]
   endfunction flux_x

   pure subroutine eigen2(ql, qr, l, r, lam)
   !< Roe-averaged left and right eigenvectors of the x flux (rows of l, columns of r) and its eigenvalues
   !< u - c, u, u, u + c.
   real(R_P), intent(in)  :: ql(4), qr(4)
   real(R_P), intent(out) :: l(4,4), r(4,4), lam(4)
   real(R_P)              :: wl(4), wr(4), sl, sr, u, v, hh, c, b1, b2, q2
   wl = cons2prim2(ql) ; wr = cons2prim2(qr)
   sl = sqrt(wl(1)) ; sr = sqrt(wr(1))
   u = (sl*wl(2) + sr*wr(2))/(sl + sr)
   v = (sl*wl(3) + sr*wr(3))/(sl + sr)
   hh = (sl*(ql(4) + wl(4))/wl(1) + sr*(qr(4) + wr(4))/wr(1))/(sl + sr)
   q2 = 0.5_R_P*(u*u + v*v)
   c = sqrt((GAMMA - 1._R_P)*(hh - q2))
   lam = [u - c, u, u, u + c]
   r(:,1) = [1._R_P, u - c, v, hh - u*c]
   r(:,2) = [1._R_P, u, v, q2]
   r(:,3) = [0._R_P, 0._R_P, 1._R_P, v]
   r(:,4) = [1._R_P, u + c, v, hh + u*c]
   b1 = (GAMMA - 1._R_P)/(c*c) ; b2 = b1*q2
   l(1,:) = 0.5_R_P*[b2 + u/c, -b1*u - 1._R_P/c, -b1*v, b1]
   l(2,:) = [1._R_P - b2, b1*u, b1*v, -b1]
   l(3,:) = [-v, 0._R_P, 1._R_P, 0._R_P]
   l(4,:) = 0.5_R_P*[b2 - u/c, -b1*u + 1._R_P/c, -b1*v, b1]
   endsubroutine eigen2

   subroutine line_residual(id, ql, n, h, alpha, res)
   !< -(F_{i+1/2} - F_{i-1/2})/h along one line of the x direction (cells 1..n, ghosts already filled).
   integer(I4P), intent(in)  :: id, n
   real(R_P),    intent(in)  :: ql(4,1-NG:n+NG), h, alpha(4)
   real(R_P),    intent(out) :: res(4,n)
   real(R_P)                 :: fn(4,1-NG:n+NG), fhat(4,0:n), l(4,4), r(4,4), lam(4), w(4,VLO:VHI), g(4,VLO:VHI)
   real(R_P)                 :: vp(VLO:VHI), vm(VLO:VHI), rp, rm, fc(4)
   integer(I4P)              :: i, m, k
   do i=1-NG, n+NG
      fn(:,i) = flux_x(ql(:,i))
   enddo
   do i=0, n
      call eigen2(ql(:,i), ql(:,i+1), l, r, lam)
      do m=VLO, VHI
         w(:,m) = matmul(l, ql(:,i+m))
         g(:,m) = matmul(l, fn(:,i+m))
      enddo
      do k=1, 4
         vp = 0.5_R_P*(g(k,:) + alpha(k)*w(k,:))
         vm = 0.5_R_P*(g(k,:) - alpha(k)*w(k,:))
         call reconstruct(id, vp, vm, rp, rm)
         fc(k) = rp + rm
      enddo
      fhat(:,i) = matmul(r, fc)
   enddo
   do i=1, n
      res(:,i) = -(fhat(:,i) - fhat(:,i-1))/h
   enddo
   endsubroutine line_residual

   subroutine ghosts2(q, n, bc)
   !< Fill the NG ghost cells of each side (corners unused).
   integer(I4P), intent(in)    :: n, bc
   real(R_P),    intent(inout) :: q(4,1-NG:n+NG,1-NG:n+NG)
   integer(I4P)                :: g
   do g=1, NG
      if (bc == BC2_PERIODIC) then
         q(:,1-g,1:n) = q(:,n+1-g,1:n) ; q(:,n+g,1:n) = q(:,g,1:n)
         q(:,1:n,1-g) = q(:,1:n,n+1-g) ; q(:,1:n,n+g) = q(:,1:n,g)
      else
         q(:,1-g,1:n) = q(:,1,1:n) ; q(:,n+g,1:n) = q(:,n,1:n)
         q(:,1:n,1-g) = q(:,1:n,1) ; q(:,1:n,n+g) = q(:,1:n,n)
      endif
   enddo
   endsubroutine ghosts2

   subroutine residual2(id, q, n, h, bc, res)
   !< The 2-D residual: the x lines plus the y lines (the y lines as x lines of the data with the momenta swapped).
   integer(I4P), intent(in)    :: id, n, bc
   real(R_P),    intent(inout) :: q(4,1-NG:n+NG,1-NG:n+NG)
   real(R_P),    intent(in)    :: h
   real(R_P),    intent(out)   :: res(4,n,n)
   real(R_P)                   :: ql(4,1-NG:n+NG), rl(4,n), ax(4), ay(4), w(4), c
   integer(I4P)                :: i, j
   rbf_dx = h
   call ghosts2(q, n, bc)
   ax = 0._R_P ; ay = 0._R_P
   do j=1, n
      do i=1, n
         w = cons2prim2(q(:,i,j)) ; c = sqrt(GAMMA*w(4)/w(1))
         ax = max(ax, abs([w(2) - c, w(2), w(2), w(2) + c]))
         ay = max(ay, abs([w(3) - c, w(3), w(3), w(3) + c]))
      enddo
   enddo
   ax = 1.1_R_P*ax ; ay = 1.1_R_P*ay
   do j=1, n
      ql = q(:,:,j)
      call line_residual(id, ql, n, h, ax, rl)
      res(:,:,j) = rl
   enddo
   do i=1, n
      ql(1,:) = q(1,i,:) ; ql(2,:) = q(3,i,:) ; ql(3,:) = q(2,i,:) ; ql(4,:) = q(4,i,:)
      call line_residual(id, ql, n, h, ay, rl)
      res(1,i,:) = res(1,i,:) + rl(1,:) ; res(2,i,:) = res(2,i,:) + rl(3,:)
      res(3,i,:) = res(3,i,:) + rl(2,:) ; res(4,i,:) = res(4,i,:) + rl(4,:)
   enddo
   endsubroutine residual2

   function euler2_dt(q, n, h, cfl) result(dt)
   !< dt = cfl h / max(|u| + c, |v| + c).
   integer(I4P), intent(in) :: n
   real(R_P),    intent(in) :: q(4,1-NG:n+NG,1-NG:n+NG), h, cfl
   real(R_P)                :: dt, w(4), smax, c
   integer(I4P)             :: i, j
   smax = 0._R_P
   do j=1, n
      do i=1, n
         w = cons2prim2(q(:,i,j)) ; c = sqrt(GAMMA*w(4)/w(1))
         smax = max(smax, abs(w(2)) + c, abs(w(3)) + c)
      enddo
   enddo
   dt = cfl*h/smax
   endfunction euler2_dt

   subroutine euler2_step(id, q, n, h, dt, bc, ok)
   !< One SSP-RK3 step, SELAH's SSP(3,3) in its two-register form on the rank-3 field q(v, i, j) (the stage state is q
   !< itself: the residual fills its ghost cells at each stage); ok false, q left at the failed stage, when a stage has
   !< a cell with rho <= 0 or p <= 0.
   integer(I4P), intent(in)    :: id, n, bc
   real(R_P),    intent(inout) :: q(4,1-NG:n+NG,1-NG:n+NG)
   real(R_P),    intent(in)    :: h, dt
   logical,      intent(out)   :: ok
   type(stepper)               :: st
   type(box)                   :: bx
   real(R_P), pointer          :: r(:,:,:)
   integer(I4P)                :: k
   bx = box_interior(lb=[1_I4P, 1-NG, 1-NG], ub=[4_I4P, n+NG, n+NG], ngc=NG, components=.true.)
   call st%initialize(SCHEME_SSP_33, FORM_SHU_OSHER_2R, lb=[1_I4P, 1-NG, 1-NG], ub=[4_I4P, n+NG, n+NG], &
                      components=.true.)
   call st%begin_step(bx, q)
   do k=1, st%stages()
      call st%begin_stage(k, bx, dt, q)
      call st%rhs_register_view(r)
      call residual2(id, q, n, h, bc, r(:,1:n,1:n))
      call st%end_stage(k, bx, dt, q)
      ok = positive2(q, n)
      if (.not.ok) then
         call st%destroy
         return
      endif
   enddo
   call st%finish_step(bx, dt, q)
   call st%destroy
   endsubroutine euler2_step

   pure function positive2(q, n) result(ok)
   !< rho > 0 and p > 0 in every interior cell.
   integer(I4P), intent(in) :: n
   real(R_P),    intent(in) :: q(4,1-NG:n+NG,1-NG:n+NG)
   logical                  :: ok
   integer(I4P)             :: i, j
   ok = .true.
   do j=1, n
      do i=1, n
         if (.not.(q(1,i,j) > 0._R_P)) then ; ok = .false. ; return ; endif
         if (.not.(q(1,i,j)*q(4,i,j) - 0.5_R_P*(q(2,i,j)**2 + q(3,i,j)**2) > 0._R_P)) then ; ok = .false. ; return ; endif
      enddo
   enddo
   endfunction positive2
endmodule seth_vv_euler2d
