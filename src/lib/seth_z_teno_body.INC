   ! Body of a WENO-Z or TENO kernel of SETH, included after the routine declares its scheme's named constants
   ! (seth_z_teno_tables): NC, LO, UMAX, CLO(0:NC-1), CW(0:NC-1), CP(0:MAXW-1,0:NC-1) (face coefficients),
   ! CL(0:MAXW-1,0:MAXW-1,0:NC-1), CDQ(0:MAXW-1,0:NC-1) (indicators as sums of squares), D(0:NC-1), TC(0:NC-1), TG,
   ! GLO, GW, CLG(0:GW-1,0:GW-1), CDG(0:GW-1) (the global indicator, the same form).
   !
   ! The indicators are evaluated as beta = sum_i CDQ(i) (sum_{a>=i} CL(a,i) u_a)**2 with CDQ >= 0 (an exact LDL^T of
   ! their quadratic forms): non-negative in floating point by construction. The expanded quadratic form cancels and
   ! can round to -1e-14 on constant data, which with eps = 1e-40 turned tau/(beta + eps) into -inf and the state into
   ! NaN (found by seth_test_properties).
   ! SETH_TENO selects the TENO weights (cutoff), otherwise the WENO-Z ones. SETH_TENO_AZ (with SETH_TENO) makes the
   ! cutoff adaptive as in Tang 2024, Eq. (20): min(C_T, tau), zero where the data is smooth (tau -> 0). tau has the
   ! dimension of the data squared, so the kernel takes the inverse descaler rmu and compares tau rmu^2 (rmu = 1: Tang's
   ! form; rmu = 1/|data|: scale invariant), evaluated as (tau rmu) rmu so that zero data (rmu huge) gives 0, not NaN.
   ! SETH_TENO_A (with SETH_TENO) computes the cutoff from the data instead of taking it (Fu, Hu & Adams 2018, as
   ! written by Fu 2021, Eqs. (17)-(19), and Liang & Fu 2024b, Eqs. (47)-(49)): eta_j = (|2 Da Db| + e)/(Da^2 + Db^2 + e)
   ! at the nodes j = -1..AHI, e = 0.9 cr/(1 - 0.9 cr) (xi/rmu)^2, m = 1 - min(1, min_j eta_j / cr),
   ! g = (1 - m)^4 (1 + 4 m), C_T = 10^-floor(a1 - a2 (1 - g)). The kernel declares AHI (1 or 2: the nodes of the min).
   !
   ! Cells by offset from the face i+1/2 (cell 0 left of it). The left state is the left-biased reconstruction of vp;
   ! the right state is the same reconstruction of the mirrored right-biased values, st(m) = vm(1-m).
   real(R_P),    intent(in)  :: vp(LO:1-LO) !< Values reconstructed with a left bias.
   real(R_P),    intent(in)  :: vm(LO:1-LO) !< Values reconstructed with a right bias.
   real(R_P),    intent(in)  :: eps         !< Regulariser of the smoothness indicators.
#ifdef SETH_TENO
#ifdef SETH_TENO_A
   real(R_P),    intent(in)  :: cr          !< Threshold C_r of the discontinuity sensor.
   real(R_P),    intent(in)  :: xi          !< Regulariser scale xi of the sensor (e = 0.9 cr/(1-0.9 cr) xi^2).
   real(R_P),    intent(in)  :: a1, a2      !< Cutoff exponents: C_T from 10^-a1 (smooth) to 10^-(a1-a2) (jump).
   real(R_P),    intent(in)  :: rmu         !< Inverse descaler of the sensor's regulariser (1: the papers' form).
   real(R_P)                 :: ct          !< Cutoff C_T of this state, from the sensor.
   real(R_P)                 :: da, db, eta, etamin, epss, mm, gg !< Sensor work.
   real(R_P),    parameter   :: POW10(0:20) = [1.e0_R_P, 1.e-1_R_P, 1.e-2_R_P, 1.e-3_R_P, 1.e-4_R_P, 1.e-5_R_P, &
                                               1.e-6_R_P, 1.e-7_R_P, 1.e-8_R_P, 1.e-9_R_P, 1.e-10_R_P, 1.e-11_R_P, &
                                               1.e-12_R_P, 1.e-13_R_P, 1.e-14_R_P, 1.e-15_R_P, 1.e-16_R_P, &
                                               1.e-17_R_P, 1.e-18_R_P, 1.e-19_R_P, 1.e-20_R_P] !< 10^-k, k = 0..20.
#else
   real(R_P),    intent(in)  :: ct          !< Cutoff C_T of the normalised smoothness measures.
#endif
#ifdef SETH_TENO_AZ
   real(R_P),    intent(in)  :: rmu         !< Inverse descaler of tau in the adaptive cutoff (1: Tang 2024).
#endif
#else
   integer(I4P), intent(in)  :: p           !< Exponent of tau/beta in the weights (1 or 2).
#endif
   real(R_P),    intent(out) :: rp          !< Left state at the face, from vp.
   real(R_P),    intent(out) :: rm          !< Right state at the face, from vm.
   real(R_P)                 :: st(LO:UMAX) !< Stencil of the current side, left-biased.
   real(R_P)                 :: q(0:NC-1)   !< Candidate face values.
   real(R_P)                 :: beta(0:NC-1)!< Candidate smoothness indicators.
   real(R_P)                 :: wgt(0:NC-1) !< Weights.
   real(R_P)                 :: r(1:2)      !< States (1 left, 2 right).
   real(R_P)                 :: tau         !< Global smoothness measure.
   real(R_P)                 :: betag       !< Global indicator (TENO6).
   real(R_P)                 :: x, wsum     !< Work.
   real(R_P)                 :: lsum        !< One term of an indicator, before squaring.
#ifdef SETH_TENO
   real(R_P)                 :: g(0:NC-1)   !< C + tau/(beta + eps), then its normalised q-th power.
   real(R_P)                 :: gmax        !< Largest g.
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
#ifdef SETH_TENO_A
      ! the cutoff of this state, from the first differences at the nodes -1..AHI
      epss = (0.9_R_P*cr/(1._R_P - 0.9_R_P*cr))*(xi/rmu)**2
      etamin = 1._R_P
      do m=-1, AHI
         da = st(m+1) - st(m) ; db = st(m) - st(m-1)
         eta = (abs(2._R_P*da*db) + epss)/(da*da + db*db + epss)
         etamin = min(etamin, eta)
      enddo
      mm = 1._R_P - min(1._R_P, etamin/cr)
      gg = (1._R_P - mm)**4*(1._R_P + 4._R_P*mm)
      ct = POW10(min(max(int(a1 - a2*(1._R_P - gg)), 0), 20))
#endif
      ! candidates
      do k=0, NC-1
         q(k) = 0._R_P
         do a=0, CW(k)-1
            q(k) = q(k) + CP(a,k)*st(CLO(k)+a)
         enddo
         beta(k) = 0._R_P
         do i=0, CW(k)-1
            lsum = 0._R_P
            do a=i, CW(k)-1
               lsum = lsum + CL(a,i,k)*st(CLO(k)+a)
            enddo
            beta(k) = beta(k) + CDQ(i,k)*(lsum*lsum)
         enddo
      enddo
      ! global smoothness measure
      tau = 0._R_P
      do k=0, NC-1
         tau = tau + TC(k)*beta(k)
      enddo
      if (TG /= 0._R_P) then ! a constant: folds
         betag = 0._R_P
         do i=0, GW-1
            lsum = 0._R_P
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
      gmax = 0._R_P
      do k=0, NC-1
         g(k) = 1._R_P + tau/(beta(k) + eps)
         gmax = max(gmax, g(k))
      enddo
      x = 0._R_P
      do k=0, NC-1
         g(k) = g(k)/gmax
         g(k) = (g(k)*g(k))*(g(k)*g(k))*(g(k)*g(k))
         x = x + g(k)
      enddo
      wsum = 0._R_P
      do k=0, NC-1
         wgt(k) = 0._R_P
#ifdef SETH_TENO_AZ
         if (g(k)/x >= min(ct, (tau*rmu)*rmu)) wgt(k) = D(k)
#else
         if (g(k)/x >= ct) wgt(k) = D(k)
#endif
         wsum = wsum + wgt(k)
      enddo
#else
      wsum = 0._R_P
      do k=0, NC-1
         x = tau/(beta(k) + eps)
         select case(p)
         case(1)
         case(2)
            x = x*x
         case default
            x = x**p
         endselect
         wgt(k) = D(k)*(1._R_P + x)
         wsum = wsum + wgt(k)
      enddo
#endif
      r(side) = 0._R_P
      do k=0, NC-1
         r(side) = r(side) + (wgt(k)/wsum)*q(k)
      enddo
   enddo
   rp = r(1)
   rm = r(2)
