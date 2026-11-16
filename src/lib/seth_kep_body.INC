   ! Body of a KEEP-PE flux kernel of SETH, included after the routine declares NL (half the stencil) and A(NL) (the
   ! central first-derivative coefficients). The arithmetic, and its order, are those of ADAM's compute_kep_flux, so
   ! that the sixth-order kernel reproduces it bitwise.
   integer(I4P), intent(in)  :: d             !< Normal direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: gamma         !< Ratio of specific heats.
   real(R8P),    intent(in)  :: w(5,1-NL:NL)    !< rho, u, v, w, p of the cells around the face, cell 0 on its left.
   real(R8P),    intent(out) :: f(5)          !< Flux of rho, rho u, rho v, rho w, E.
   real(R8P)                 :: uu            !< u_c . u_r.
   real(R8P)                 :: mf, pf        !< Two-point mass flux and pressure.
   real(R8P)                 :: rg            !< 1/(gamma - 1).
   integer(I4P)              :: l, m, n, c, r !< Counters.
   !$acc routine seq
   !$omp declare target

   rg = 1._R8P/(gamma - 1._R8P)
   do n=1, 5
      f(n) = 0._R8P
   enddo
   do l=1, NL
      do m=0, l-1
         c = -m ; r = -m + l
         mf = 0.25_R8P*(w(1,c) + w(1,r))*(w(1+d,c) + w(1+d,r))
         pf = 0.5_R8P*(w(5,c) + w(5,r))
         f(1) = f(1) + 2._R8P*A(l)*mf
         do n=1, 3
            f(1+n) = f(1+n) + 2._R8P*A(l)*mf*0.5_R8P*(w(1+n,c) + w(1+n,r))
         enddo
         f(1+d) = f(1+d) + 2._R8P*A(l)*pf
         uu = w(2,c)*w(2,r) + w(3,c)*w(3,r) + w(4,c)*w(4,r)
         f(5) = f(5) + 2._R8P*A(l)*(0.5_R8P*mf*uu + 0.5_R8P*rg*pf*(w(1+d,c) + w(1+d,r)) + &
                                    0.5_R8P*(w(1+d,c)*w(5,r) + w(1+d,r)*w(5,c)))
      enddo
   enddo
