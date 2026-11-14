   ! Body of a WENO-Z or TENO kernel of SETH, included after the routine declares its scheme's named constants
   ! (seth_z_teno_tables): NC, LO, UMAX, CLO(0:NC-1), CW(0:NC-1), CP(0:MAXW-1,0:NC-1) (face coefficients),
   ! CL(0:MAXW-1,0:MAXW-1,0:NC-1), CDQ(0:MAXW-1,0:NC-1) (indicators as sums of squares), D(0:NC-1), TC(0:NC-1), TG,
   ! GLO, GW, CLG(0:GW-1,0:GW-1), CDG(0:GW-1) (the global indicator, the same form).
   !
   ! The indicators are evaluated as beta = sum_i CDQ(i) (sum_{a>=i} CL(a,i) u_a)**2 with CDQ >= 0 (an exact LDL^T of
   ! their quadratic forms): non-negative in floating point by construction. The expanded quadratic form cancels and
   ! can round to -1e-14 on constant data, which with eps = 1e-40 turned tau/(beta + eps) into -inf and the state into
   ! NaN (found by seth_test_properties).
   ! SETH_TENO selects the TENO weights (cutoff), otherwise the WENO-Z ones.
   !
   ! Cells by offset from the face i+1/2 (cell 0 left of it). The left state is the left-biased reconstruction of vp;
   ! the right state is the same reconstruction of the mirrored right-biased values, st(m) = vm(1-m).
   real(R8P),    intent(in)  :: vp(LO:1-LO) !< Values reconstructed with a left bias.
   real(R8P),    intent(in)  :: vm(LO:1-LO) !< Values reconstructed with a right bias.
   real(R8P),    intent(in)  :: eps         !< Regulariser of the smoothness indicators.
#ifdef SETH_TENO
   real(R8P),    intent(in)  :: ct          !< Cutoff C_T of the normalised smoothness measures.
#else
   integer(I4P), intent(in)  :: p           !< Exponent of tau/beta in the weights (1 or 2).
#endif
   real(R8P),    intent(out) :: rp          !< Left state at the face, from vp.
   real(R8P),    intent(out) :: rm          !< Right state at the face, from vm.
   real(R8P)                 :: st(LO:UMAX) !< Stencil of the current side, left-biased.
   real(R8P)                 :: q(0:NC-1)   !< Candidate face values.
   real(R8P)                 :: beta(0:NC-1)!< Candidate smoothness indicators.
   real(R8P)                 :: wgt(0:NC-1) !< Weights.
   real(R8P)                 :: r(1:2)      !< States (1 left, 2 right).
   real(R8P)                 :: tau         !< Global smoothness measure.
   real(R8P)                 :: betag       !< Global indicator (TENO6).
   real(R8P)                 :: x, wsum     !< Work.
   real(R8P)                 :: lsum        !< One term of an indicator, before squaring.
#ifdef SETH_TENO
   real(R8P)                 :: g(0:NC-1)   !< C + tau/(beta + eps), then its normalised q-th power.
   real(R8P)                 :: gmax        !< Largest g.
#endif
   integer(I4P)              :: side, k, a, i, m
   !$acc routine seq
   !$omp declare target

   do side=1, 2
      if (side == 1) then
         do m=LO, UMAX
            st(m) = vp(m)
         enddo
      else
         do m=LO, UMAX
            st(m) = vm(1-m)
         enddo
      endif
      ! candidates
      do k=0, NC-1
         q(k) = 0._R8P
         do a=0, CW(k)-1
            q(k) = q(k) + CP(a,k)*st(CLO(k)+a)
         enddo
         beta(k) = 0._R8P
         do i=0, CW(k)-1
            lsum = 0._R8P
            do a=i, CW(k)-1
               lsum = lsum + CL(a,i,k)*st(CLO(k)+a)
            enddo
            beta(k) = beta(k) + CDQ(i,k)*(lsum*lsum)
         enddo
      enddo
      ! global smoothness measure
      tau = 0._R8P
      do k=0, NC-1
         tau = tau + TC(k)*beta(k)
      enddo
      if (TG /= 0._R8P) then ! a constant: folds
         betag = 0._R8P
         do i=0, GW-1
            lsum = 0._R8P
            do a=i, GW-1
               lsum = lsum + CLG(a,i)*st(GLO+a)
            enddo
            betag = betag + CDG(i)*(lsum*lsum)
         enddo
         tau = tau + TG*betag
      endif
      tau = abs(tau)
#ifdef SETH_TENO
      ! chi_k = gamma_k / sum gamma, gamma_k = g_k**6, evaluated as (g_k/g_max)**6 / sum (g/g_max)**6: the same chi,
      ! without the overflow of g**6 when tau/beta is large (a big jump on piecewise-constant data)
      gmax = 0._R8P
      do k=0, NC-1
         g(k) = 1._R8P + tau/(beta(k) + eps)
         gmax = max(gmax, g(k))
      enddo
      x = 0._R8P
      do k=0, NC-1
         g(k) = g(k)/gmax
         g(k) = (g(k)*g(k))*(g(k)*g(k))*(g(k)*g(k))
         x = x + g(k)
      enddo
      wsum = 0._R8P
      do k=0, NC-1
         wgt(k) = 0._R8P
         if (g(k)/x >= ct) wgt(k) = D(k)
         wsum = wsum + wgt(k)
      enddo
#else
      wsum = 0._R8P
      do k=0, NC-1
         x = tau/(beta(k) + eps)
         select case(p)
         case(1)
         case(2)
            x = x*x
         case default
            x = x**p
         endselect
         wgt(k) = D(k)*(1._R8P + x)
         wsum = wsum + wgt(k)
      enddo
#endif
      r(side) = 0._R8P
      do k=0, NC-1
         r(side) = r(side) + (wgt(k)/wsum)*q(k)
      enddo
   enddo
   rp = r(1)
   rm = r(2)
