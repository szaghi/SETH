   ! Bench body of the positivity limiter with explicit loops and local arrays of dimension NVD, loops to NVL (macros of
   ! the including module: NVD a constant; NVL the dummy nv or a constant). The same algorithm as
   ! src/lib/seth_positivity.F90, whose array expressions over (nv) and automatic arrays it replaces.
   pure subroutine positivity_euler(nv, gamma, lam, eps_rho, eps_p, ql, qr, f_high, f_low, f, theta)
   integer(I4P), intent(in)  :: nv
   real(R8P),    intent(in)  :: gamma, lam, eps_rho, eps_p
   real(R8P),    intent(in)  :: ql(NVD), qr(NVD), f_high(NVD), f_low(NVD)
   real(R8P),    intent(out) :: f(NVD), theta
   real(R8P)                 :: a0(NVD), a1(NVD), b0(NVD), b1(NVD), t
   integer(I4P)              :: it, k
   !$acc routine seq
   !$omp declare target

   do k=1, NVL
      a0(k) = ql(k) - lam*f_low(k) ; a1(k) = ql(k) - lam*f_high(k)
      b0(k) = qr(k) + lam*f_low(k) ; b1(k) = qr(k) + lam*f_high(k)
   enddo
   if (admissible(nv, gamma, a1, eps_rho, eps_p) .and. admissible(nv, gamma, b1, eps_rho, eps_p)) then
      theta = 1._R8P
      do k=1, NVL
         f(k) = f_high(k)
      enddo
      return
   endif
   if (.not.(admissible(nv, gamma, a0, eps_rho, eps_p) .and. admissible(nv, gamma, b0, eps_rho, eps_p))) then
      theta = 0._R8P
      do k=1, NVL
         f(k) = f_low(k)
      enddo
      return
   endif
   t = min(segment_bound(nv, gamma, eps_rho, eps_p, a0, a1), segment_bound(nv, gamma, eps_rho, eps_p, b0, b1))
   theta = t
   do it=0, 52
      do k=1, NVL
         f(k) = theta*f_high(k) + (1._R8P - theta)*f_low(k)
         a1(k) = ql(k) - lam*f(k) ; b1(k) = qr(k) + lam*f(k)
      enddo
      if (admissible(nv, gamma, a1, eps_rho, eps_p) .and. admissible(nv, gamma, b1, eps_rho, eps_p)) return
      theta = max(t*(1._R8P - 1.e-15_R8P*2._R8P**it), 0._R8P)
   enddo
   theta = 0._R8P
   do k=1, NVL
      f(k) = f_low(k)
   enddo
   endsubroutine positivity_euler

   pure function segment_bound(nv, gamma, eps_rho, eps_p, s0, s1) result(t)
   integer(I4P), intent(in) :: nv
   real(R8P),    intent(in) :: gamma, eps_rho, eps_p, s0(NVD), s1(NVD)
   real(R8P)                :: t, d(NVD), kappa, qa, qb, qc, disc, r1, r2, dd, sd, ss
   integer(I4P)             :: k
   !$acc routine seq
   !$omp declare target

   do k=1, NVL
      d(k) = s1(k) - s0(k)
   enddo
   t = 1._R8P
   if (s1(1) < eps_rho) t = (s0(1) - eps_rho)/(s0(1) - s1(1))
   kappa = eps_p/(gamma - 1._R8P)
   dd = 0._R8P ; sd = 0._R8P ; ss = 0._R8P
   do k=2, NVL-1
      dd = dd + d(k)*d(k) ; sd = sd + s0(k)*d(k) ; ss = ss + s0(k)*s0(k)
   enddo
   qa = d(1)*d(NVL) - 0.5_R8P*dd
   qb = s0(1)*d(NVL) + s0(NVL)*d(1) - sd - kappa*d(1)
   qc = s0(1)*s0(NVL) - 0.5_R8P*ss - kappa*s0(1)
   if (qa*t*t + qb*t + qc >= 0._R8P) return
   if (abs(qa) <= epsilon(1._R8P)*(abs(qb) + abs(qc))) then
      t = min(t, -qc/qb)
   else
      disc = max(qb*qb - 4._R8P*qa*qc, 0._R8P)
      r1 = -0.5_R8P*(qb + sign(sqrt(disc), qb))
      r2 = qc/r1
      r1 = r1/qa
      disc = t
      if (r1 >= 0._R8P .and. r1 <= disc) disc = r1
      if (r2 >= 0._R8P .and. r2 <= disc) disc = r2
      t = min(t, disc)
   endif
   t = max(t, 0._R8P)
   endfunction segment_bound

   pure logical function admissible(nv, gamma, q, eps_rho, eps_p) result(ok)
   integer(I4P), intent(in) :: nv
   real(R8P),    intent(in) :: gamma, q(NVD), eps_rho, eps_p
   real(R8P)                :: ke
   integer(I4P)             :: k
   !$acc routine seq
   !$omp declare target
   ok = .false.
   if (.not.(q(1) >= eps_rho)) return
   ke = 0._R8P
   do k=2, NVL-1
      ke = ke + q(k)*q(k)
   enddo
   ok = (gamma - 1._R8P)*(q(NVL) - 0.5_R8P*ke/q(1)) >= eps_p
   endfunction admissible
