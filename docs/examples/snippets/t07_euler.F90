!< Tutorial 7: a system. The Euler equations on Sod's shock tube, by characteristic and by component-wise splitting.
program t07_euler
use seth
use selah, only : stepper, box, box_interior, SCHEME_SSP_33, FORM_SHU_OSHER_2R   ! the time integrator: SELAH's
implicit none
real(R_P), parameter :: GAMMA = 1.4_R_P
integer,   parameter :: N = 200, NG = 3
real(R_P)            :: q(3,N), h, x, w(3), err, umax, pstar, ustar, tv, u_old
integer              :: i, k
logical              :: characteristic

h = 1._R_P/N
call star_state([1._R_P, 0._R_P, 1._R_P], [0.125_R_P, 0._R_P, 0.1_R_P], pstar, ustar)
print '(A)', "Sod's shock tube, 200 points, t = 0.2, TENO5"
print '(A,F8.5)', '  exact velocity behind the shock, u* = ', ustar
print '(A)', '  splitting         L1(rho)   max u - u*   TV(u) - 2 u*'
do k=1, 2
   characteristic = k == 1
   do i=1, N                                    ! point values at x_i = (i - 1/2) h: finite differences
      x = (i - 0.5_R_P)*h
      q(:,i) = prim2cons(merge([1._R_P, 0._R_P, 1._R_P], [0.125_R_P, 0._R_P, 0.1_R_P], x < 0.5_R_P))
   enddo
   call evolve(q, 0.2_R_P, characteristic)
   err = 0 ; umax = 0 ; tv = 0
   do i=1, N
      w = cons2prim(q(:,i))
      err = err + abs(w(1) - exact_rho((i - 0.5_R_P)*h, 0.2_R_P, pstar, ustar))*h
      umax = max(umax, w(2))
      if (i > 1) tv = tv + abs(w(2) - u_old)
      u_old = w(2)
   enddo
   print '(2X,A16,ES10.3,F13.5,F15.5)', merge('characteristic', 'component-wise', characteristic), err, &
                                        max(umax - ustar, 0._R_P), tv - 2*ustar
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

   subroutine residual(q, res, characteristic)
   !< res = -(F_{i+1/2} - F_{i-1/2})/h: Lax-Friedrichs splitting per field, TENO5, transmissive ends.
   real(R_P), intent(in)  :: q(3,N)
   real(R_P), intent(out) :: res(3,N)
   logical,   intent(in)  :: characteristic
   real(R_P)              :: qg(3,1-NG:N+NG), fg(3,1-NG:N+NG), wc(3,-2:3), gc(3,-2:3), l(3,3), r(3,3)
   real(R_P)              :: alpha(3), fc(3), fhat(3,0:N), wp(3), rp, rm, c
   integer                :: i, k, m
   qg(:,1:N) = q
   do i=1, NG
      qg(:,1-i) = q(:,1) ; qg(:,N+i) = q(:,N)   ! transmissive: copies of the end cells
   enddo
   do i=1-NG, N+NG
      fg(:,i) = flux(qg(:,i))
   enddo
   alpha = 0                                    ! alpha_k = 1.1 max |lambda_k| over the domain
   do i=1, N
      wp = cons2prim(q(:,i)) ; c = sqrt(GAMMA*wp(3)/wp(1))
      alpha = max(alpha, abs([wp(2) - c, wp(2), wp(2) + c]))
   enddo
   alpha = 1.1_R_P*alpha
   if (.not.characteristic) alpha = maxval(alpha)
   do i=0, N
      if (characteristic) then                  ! project the stencil on the face's characteristic fields
         call eigenvectors(qg(:,i), qg(:,i+1), l, r)
         do m=-2, 3
            wc(:,m) = matmul(l, qg(:,i+m)) ; gc(:,m) = matmul(l, fg(:,i+m))
         enddo
      else                                      ! the conservative variables themselves
         wc = qg(:,i-2:i+3) ; gc = fg(:,i-2:i+3)
      endif
      do k=1, 3                                 ! each field: split, reconstruct, sum
         call seth_teno5(0.5_R_P*(gc(k,:) + alpha(k)*wc(k,:)), 0.5_R_P*(gc(k,:) - alpha(k)*wc(k,:)), &
                         SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)
         fc(k) = rp + rm
      enddo
      if (characteristic) then
         fhat(:,i) = matmul(r, fc)              ! back to the conservative variables
      else
         fhat(:,i) = fc
      endif
   enddo
   res = -(fhat(:,1:N) - fhat(:,0:N-1))/h
   endsubroutine residual

   subroutine evolve(q, tend, characteristic)
   !< SSP-RK3 to tend, CFL 0.4.
   real(R_P), intent(inout) :: q(3,N)
   real(R_P), intent(in)    :: tend
   logical,   intent(in)    :: characteristic
   real(R_P)                :: t, dt, smax, wp(3)
   integer                  :: i
   type(stepper)            :: st
   type(box)                :: bx
   real(R_P), pointer       :: r(:,:)
   integer                  :: k
   bx = box_interior(lb=[1, 1], ub=[3, N], ngc=0, components=.true.)
   call st%initialize(SCHEME_SSP_33, FORM_SHU_OSHER_2R, lb=[1, 1], ub=[3, N], components=.true.)
   t = 0
   do while (t < tend)
      smax = 0
      do i=1, N
         wp = cons2prim(q(:,i)) ; smax = max(smax, abs(wp(2)) + sqrt(GAMMA*wp(3)/wp(1)))
      enddo
      dt = min(0.4_R_P*h/smax, tend - t)
      call st%begin_step(bx, q)
      do k=1, st%stages()                          ! SELAH's SSP(3,3): the stage state is q itself
         call st%begin_stage(k, bx, dt, q)
         call st%rhs_register_view(r)
         call residual(q, r, characteristic)
         call st%end_stage(k, bx, dt, q)
      enddo
      call st%finish_step(bx, dt, q)
      t = t + dt
   enddo
   call st%destroy
   endsubroutine evolve

   subroutine star_state(wl, wr, pstar, ustar)
   !< The pressure and velocity between the waves of the Riemann problem (Toro, chapter 4: Newton on p).
   real(R_P), intent(in)  :: wl(3), wr(3)
   real(R_P), intent(out) :: pstar, ustar
   real(R_P)              :: fl, fr, dl, dr
   integer                :: it
   pstar = 0.5_R_P*(wl(3) + wr(3))
   do it=1, 50
      call pfun(pstar, wl, fl, dl) ; call pfun(pstar, wr, fr, dr)
      pstar = pstar - (fl + fr + wr(2) - wl(2))/(dl + dr)
   enddo
   call pfun(pstar, wl, fl, dl) ; call pfun(pstar, wr, fr, dr)
   ustar = 0.5_R_P*(wl(2) + wr(2) + fr - fl)
   endsubroutine star_state

   pure subroutine pfun(p, w, f, d)
   real(R_P), intent(in)  :: p, w(3)
   real(R_P), intent(out) :: f, d
   real(R_P)              :: a, b, c
   c = sqrt(GAMMA*w(3)/w(1))
   if (p > w(3)) then                           ! shock
      a = 2/((GAMMA + 1)*w(1)) ; b = (GAMMA - 1)/(GAMMA + 1)*w(3)
      f = (p - w(3))*sqrt(a/(p + b)) ; d = sqrt(a/(p + b))*(1 - 0.5_R_P*(p - w(3))/(p + b))
   else                                         ! rarefaction
      f = 2*c/(GAMMA - 1)*((p/w(3))**((GAMMA - 1)/(2*GAMMA)) - 1)
      d = (p/w(3))**(-(GAMMA + 1)/(2*GAMMA))/(w(1)*c)
   endif
   endsubroutine pfun

   pure function exact_rho(x, t, pstar, ustar) result(rho)
   !< The exact density of Sod's problem: rarefaction left, contact, shock right.
   real(R_P), intent(in) :: x, t, pstar, ustar
   real(R_P)             :: rho, s, cl, cstar, rhol, rhor, rho_sl, rho_sr, sh
   rhol = 1 ; rhor = 0.125_R_P ; cl = sqrt(GAMMA)
   s = (x - 0.5_R_P)/t
   rho_sl = rhol*(pstar/1)**(1/GAMMA)
   cstar = cl*(pstar/1)**((GAMMA - 1)/(2*GAMMA))
   rho_sr = rhor*((pstar/0.1_R_P + (GAMMA - 1)/(GAMMA + 1))/((GAMMA - 1)/(GAMMA + 1)*pstar/0.1_R_P + 1))
   sh = sqrt(GAMMA*0.1_R_P/rhor)*sqrt((GAMMA + 1)/(2*GAMMA)*pstar/0.1_R_P + (GAMMA - 1)/(2*GAMMA))
   if (s < -cl) then
      rho = rhol
   elseif (s < ustar - cstar) then              ! inside the rarefaction fan
      rho = rhol*(2/(GAMMA + 1) + (GAMMA - 1)/((GAMMA + 1)*cl)*(-s))**(2/(GAMMA - 1))
   elseif (s < ustar) then
      rho = rho_sl
   elseif (s < sh) then
      rho = rho_sr
   else
      rho = rhor
   endif
   endfunction exact_rho
endprogram t07_euler
