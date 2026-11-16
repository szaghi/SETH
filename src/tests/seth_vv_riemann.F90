!< SETH verification: the exact solution of the Riemann problem of the 1-D Euler equations (ideal gas).
module seth_vv_riemann
!< SETH verification: the exact solution of the Riemann problem of the 1-D Euler equations (ideal gas).
!<
!< Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics, 3rd ed. (2009), chapter 4: the star pressure from
!< the pressure functions f_L + f_R + (u_R - u_L) = 0 by Newton's method from the two-rarefaction guess (section 4.3.2,
!< iterated to a relative change of 1e-14), the star velocity (4.9), and the sampling of the self-similar solution at
!< xi = (x - x0)/t (section 4.5). Vacuum generation is not handled (none of the cases here creates it).
use penf, only : I4P, R_P
use seth_test_parameters, only : TOL_RIEMANN
implicit none
private
public :: riemann_sample

contains
   pure subroutine pressure_function(p, rk, pk, ck, g, f, fd)
   real(R_P), intent(in)  :: p, rk, pk, ck, g
   real(R_P), intent(out) :: f, fd
   real(R_P)              :: ak, bk, qrt
   if (p > pk) then ! shock
      ak = 2._R_P/((g + 1._R_P)*rk) ; bk = (g - 1._R_P)/(g + 1._R_P)*pk
      qrt = sqrt(ak/(bk + p))
      f = (p - pk)*qrt
      fd = (1._R_P - 0.5_R_P*(p - pk)/(bk + p))*qrt
   else ! rarefaction
      f = 2._R_P*ck/(g - 1._R_P)*((p/pk)**((g - 1._R_P)/(2._R_P*g)) - 1._R_P)
      fd = (p/pk)**(-(g + 1._R_P)/(2._R_P*g))/(rk*ck)
   endif
   endsubroutine pressure_function

   pure subroutine star(wl, wr, g, ps, us)
   !< Star pressure and velocity.
   real(R_P), intent(in)  :: wl(3), wr(3), g
   real(R_P), intent(out) :: ps, us
   real(R_P)              :: cl, cr, fl, fr, fdl, fdr, p0, z
   integer(I4P)           :: it
   cl = sqrt(g*wl(3)/wl(1)) ; cr = sqrt(g*wr(3)/wr(1))
   z = (g - 1._R_P)/(2._R_P*g)
   ps = ((cl + cr - 0.5_R_P*(g - 1._R_P)*(wr(2) - wl(2)))/(cl/wl(3)**z + cr/wr(3)**z))**(1._R_P/z) ! two-rarefaction
   do it=1, 100
      call pressure_function(ps, wl(1), wl(3), cl, g, fl, fdl)
      call pressure_function(ps, wr(1), wr(3), cr, g, fr, fdr)
      p0 = ps
      ps = max(ps - (fl + fr + wr(2) - wl(2))/(fdl + fdr), TOL_RIEMANN*p0)
      if (abs(ps - p0) <= TOL_RIEMANN*0.5_R_P*(ps + p0)) exit
   enddo
   call pressure_function(ps, wl(1), wl(3), cl, g, fl, fdl)
   call pressure_function(ps, wr(1), wr(3), cr, g, fr, fdr)
   us = 0.5_R_P*(wl(2) + wr(2) + fr - fl)
   endsubroutine star

   pure function riemann_sample(wl, wr, g, xi) result(w)
   !< The primitive state (rho, u, p) at xi = (x - x0)/t.
   real(R_P), intent(in) :: wl(3), wr(3), g, xi
   real(R_P)             :: w(3), ps, us, cl, cr, sh, st, sl, sr, c, gm, gp
   call star(wl, wr, g, ps, us)
   cl = sqrt(g*wl(3)/wl(1)) ; cr = sqrt(g*wr(3)/wr(1))
   gm = (g - 1._R_P)/(g + 1._R_P) ; gp = (g + 1._R_P)/(2._R_P*g)
   if (xi <= us) then ! left of the contact
      if (ps > wl(3)) then ! left shock
         sl = wl(2) - cl*sqrt(gp*ps/wl(3) + (g - 1._R_P)/(2._R_P*g))
         if (xi <= sl) then
            w = wl
         else
            w = [wl(1)*(ps/wl(3) + gm)/(gm*ps/wl(3) + 1._R_P), us, ps]
         endif
      else ! left rarefaction
         sh = wl(2) - cl ; st = us - cl*(ps/wl(3))**((g - 1._R_P)/(2._R_P*g))
         if (xi <= sh) then
            w = wl
         elseif (xi > st) then
            w = [wl(1)*(ps/wl(3))**(1._R_P/g), us, ps]
         else
            c = 2._R_P/(g + 1._R_P)*(cl + 0.5_R_P*(g - 1._R_P)*(wl(2) - xi))
            w(2) = 2._R_P/(g + 1._R_P)*(cl + 0.5_R_P*(g - 1._R_P)*wl(2) + xi)
            w(1) = wl(1)*(c/cl)**(2._R_P/(g - 1._R_P))
            w(3) = wl(3)*(c/cl)**(2._R_P*g/(g - 1._R_P))
         endif
      endif
   else ! right of the contact
      if (ps > wr(3)) then ! right shock
         sr = wr(2) + cr*sqrt(gp*ps/wr(3) + (g - 1._R_P)/(2._R_P*g))
         if (xi >= sr) then
            w = wr
         else
            w = [wr(1)*(ps/wr(3) + gm)/(gm*ps/wr(3) + 1._R_P), us, ps]
         endif
      else ! right rarefaction
         sh = wr(2) + cr ; st = us + cr*(ps/wr(3))**((g - 1._R_P)/(2._R_P*g))
         if (xi >= sh) then
            w = wr
         elseif (xi < st) then
            w = [wr(1)*(ps/wr(3))**(1._R_P/g), us, ps]
         else
            c = 2._R_P/(g + 1._R_P)*(cr - 0.5_R_P*(g - 1._R_P)*(wr(2) - xi))
            w(2) = 2._R_P/(g + 1._R_P)*(-cr + 0.5_R_P*(g - 1._R_P)*wr(2) + xi)
            w(1) = wr(1)*(c/cr)**(2._R_P/(g - 1._R_P))
            w(3) = wr(3)*(c/cr)**(2._R_P*g/(g - 1._R_P))
         endif
      endif
   endif
   endfunction riemann_sample
endmodule seth_vv_riemann
