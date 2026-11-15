!< SETH test: the positivity-preserving flux limiter (Hu, Adams & Shu 2013) on random states, 1-D to 3-D layouts.
program seth_test_positivity
!< SETH test: the positivity-preserving flux limiter (Hu, Adams & Shu 2013) on random states, 1-D to 3-D layouts.
!<
!< For nv = 3, 4, 5 (q = rho, nv-2 momentum components, E; ideal gas, gamma = 1.4), on random cell pairs with densities
!< and pressures from 1e-8 to 1 and velocities up to 10 (near vacuum: kinetic energy far above the internal one), the
!< low-order flux is the Lax-Friedrichs flux with alpha = max(|u_n| + c) of the two cells and lam = 1/alpha (CFL 1/2,
!< the limit of HAS), and the high-order flux is that flux plus a random perturbation, often large enough to break
!< positivity. Asserted:
!<
!< 1. theta in [0, 1] and f = theta f_high + (1 - theta) f_low (to round-off);
!< 2. the limited one-sided states q_l - lam f, q_r + lam f are admissible (rho >= eps, p >= eps, eps = 1e-13);
!< 3. theta is maximal: for theta < 1 it lies within 1e-9 of the largest admissible factor found by bisection (60
!<    steps) on the exact condition;
!< 4. theta = 1 and f = f_high bitwise when f_high keeps both one-sided states admissible (no perturbation);
!< 5. theta = 0 and f = f_low when even f_low gives an inadmissible state (lam ten times too large);
!< 6. the low-order flux under its CFL limit is itself admissible (the premise of the limiter, HAS eq. 13).
!<
!< Random numbers from a fixed, printed seed (SETH_TEST_SEED overrides it).
use seth
use seth_test_rng, only : seed_rng
implicit none
real(R8P),    parameter :: G = 1.4_R8P, EPS = SETH_POSITIVITY_EPS
integer(I4P), parameter :: NR = 20000_I4P
real(R8P)               :: ql(5), qr(5), fh(5), fl(5), f(5), theta, lam, alpha, x, worst_comb, worst_opt, ref
integer(I4P)            :: nv, n_, nlim, nfull, failures
logical                 :: passed

passed = .true.
call seed_rng(1006_I4P)
do nv=3, 5
   failures = 0 ; nlim = 0 ; nfull = 0 ; worst_comb = 0._R8P ; worst_opt = 0._R8P
   do n_=1, NR
      call random_state(nv, ql(1:nv)) ; call random_state(nv, qr(1:nv))
      alpha = max(speed(nv, ql(1:nv)), speed(nv, qr(1:nv)))
      lam = 1._R8P/alpha
      fl(1:nv) = 0.5_R8P*(flux(nv, ql(1:nv)) + flux(nv, qr(1:nv))) - 0.5_R8P*alpha*(qr(1:nv) - ql(1:nv))
      ! 6. the premise: the low-order flux is admissible under its CFL limit
      if (.not.(seth_euler_admissible(nv, G, ql(1:nv) - lam*fl(1:nv), 0._R8P, 0._R8P) .and. &
                seth_euler_admissible(nv, G, qr(1:nv) + lam*fl(1:nv), 0._R8P, 0._R8P))) &
         call fail('the Lax-Friedrichs flux is not admissible at CFL 1/2')
      ! 4. no perturbation: theta = 1 bitwise when f_high = f_low is admissible
      if (seth_euler_admissible(nv, G, ql(1:nv) - lam*fl(1:nv), EPS, EPS) .and. &
          seth_euler_admissible(nv, G, qr(1:nv) + lam*fl(1:nv), EPS, EPS)) then
         call seth_positivity_euler(nv, G, lam, EPS, EPS, ql(1:nv), qr(1:nv), fl(1:nv), fl(1:nv), f(1:nv), theta)
         if (theta /= 1._R8P .or. any(f(1:nv) /= fl(1:nv))) call fail('theta /= 1 on an admissible flux')
      else
         cycle ! the floors exclude this pair: the limiter has nothing admissible to blend towards
      endif
      ! a perturbed high-order flux
      call random_number(x)
      fh(1:nv) = fl(1:nv) + (10._R8P**(4._R8P*x - 3._R8P))*perturbation(nv, ql(1:nv), qr(1:nv))
      call seth_positivity_euler(nv, G, lam, EPS, EPS, ql(1:nv), qr(1:nv), fh(1:nv), fl(1:nv), f(1:nv), theta)
      ! 1.
      if (.not.(theta >= 0._R8P .and. theta <= 1._R8P)) call fail('theta outside [0, 1]')
      worst_comb = max(worst_comb, maxval(abs(f(1:nv) - (theta*fh(1:nv) + (1._R8P - theta)*fl(1:nv)))) / &
                                   (maxval(abs(fh(1:nv))) + maxval(abs(fl(1:nv)))))
      ! 2.
      if (.not.(seth_euler_admissible(nv, G, ql(1:nv) - lam*f(1:nv), EPS, EPS) .and. &
                seth_euler_admissible(nv, G, qr(1:nv) + lam*f(1:nv), EPS, EPS))) &
         call fail('a limited one-sided state is not admissible')
      ! 3.
      if (theta < 1._R8P) then
         nlim = nlim + 1
         ref = bisection(nv, lam, ql(1:nv), qr(1:nv), fh(1:nv), fl(1:nv))
         worst_opt = max(worst_opt, ref - theta)
         if (ref - theta > 1.e-9_R8P) call fail('theta not maximal')
      else
         nfull = nfull + 1
      endif
      if (failures > 5) exit
   enddo
   ! 5. lam ten times too large: f_low itself fails on some pair; theta = 0, f = f_low there
   call five(nv)
   print '(A,I1,A,I6,A,I6,A,ES9.2,A,ES9.2,A,I0)', 'nv=', nv, ': limited ', nlim, ', unlimited ', nfull, &
         ', worst |f - blend| ', worst_comb, ', worst theta below the maximum ', worst_opt, ', failures ', failures
   if (failures > 0) passed = .false.
enddo
if (.not.passed) error stop 'positivity limiter test failed'
print '(A)', 'seth_test_positivity: ok'

contains
   subroutine fail(what)
   character(*), intent(in) :: what
   if (failures < 3) print '(A)', '   FAILED: '//what
   failures = failures + 1
   endsubroutine fail

   subroutine random_state(nv, q)
   !< rho, p log-uniform in [1e-8, 1], each velocity component uniform in [-10, 10].
   integer(I4P), intent(in)  :: nv
   real(R8P),    intent(out) :: q(nv)
   real(R8P)                 :: r(nv), rho, p
   call random_number(r)
   rho = 10._R8P**(-8._R8P*r(1)) ; p = 10._R8P**(-8._R8P*r(nv))
   q(1) = rho
   q(2:nv-1) = rho*(20._R8P*r(2:nv-1) - 10._R8P)
   q(nv) = p/(G - 1._R8P) + 0.5_R8P*dot_product(q(2:nv-1), q(2:nv-1))/rho
   endsubroutine random_state

   pure function flux(nv, q) result(fx)
   !< The Euler flux normal to direction 1.
   integer(I4P), intent(in) :: nv
   real(R8P),    intent(in) :: q(nv)
   real(R8P)                :: fx(nv), u, p
   u = q(2)/q(1)
   p = (G - 1._R8P)*(q(nv) - 0.5_R8P*dot_product(q(2:nv-1), q(2:nv-1))/q(1))
   fx = u*q
   fx(2) = fx(2) + p
   fx(nv) = fx(nv) + p*u
   endfunction flux

   pure function speed(nv, q) result(s)
   integer(I4P), intent(in) :: nv
   real(R8P),    intent(in) :: q(nv)
   real(R8P)                :: s, p
   p = (G - 1._R8P)*(q(nv) - 0.5_R8P*dot_product(q(2:nv-1), q(2:nv-1))/q(1))
   s = abs(q(2)/q(1)) + sqrt(G*p/q(1))
   endfunction speed

   function perturbation(nv, a, b) result(d)
   !< A random vector of the magnitude of the two states' fluxes.
   integer(I4P), intent(in) :: nv
   real(R8P),    intent(in) :: a(nv), b(nv)
   real(R8P)                :: d(nv), r(nv)
   call random_number(r)
   d = (2._R8P*r - 1._R8P)*(abs(flux(nv, a)) + abs(flux(nv, b)) + abs(a) + abs(b))
   endfunction perturbation

   function bisection(nv, lam, a, b, fh, fl) result(t)
   !< The largest admissible factor by 60 bisection steps on the exact condition (theta = 1 checked first).
   integer(I4P), intent(in) :: nv
   real(R8P),    intent(in) :: lam, a(nv), b(nv), fh(nv), fl(nv)
   real(R8P)                :: t, lo, hi, m, ff(nv)
   integer(I4P)             :: k
   ff = fh
   if (seth_euler_admissible(nv, G, a - lam*ff, EPS, EPS) .and. seth_euler_admissible(nv, G, b + lam*ff, EPS, EPS)) then
      t = 1._R8P ; return
   endif
   lo = 0._R8P ; hi = 1._R8P
   do k=1, 60
      m = 0.5_R8P*(lo + hi) ; ff = m*fh + (1._R8P - m)*fl
      if (seth_euler_admissible(nv, G, a - lam*ff, EPS, EPS) .and. &
          seth_euler_admissible(nv, G, b + lam*ff, EPS, EPS)) then
         lo = m
      else
         hi = m
      endif
   enddo
   t = lo
   endfunction bisection

   subroutine five(nv)
   !< lam ten times too large: wherever f_low itself fails, theta = 0 and f = f_low.
   integer(I4P), intent(in) :: nv
   integer(I4P)             :: k, nfound
   nfound = 0
   do k=1, 2000
      call random_state(nv, ql(1:nv)) ; call random_state(nv, qr(1:nv))
      alpha = max(speed(nv, ql(1:nv)), speed(nv, qr(1:nv)))
      lam = 10._R8P/alpha
      fl(1:nv) = 0.5_R8P*(flux(nv, ql(1:nv)) + flux(nv, qr(1:nv))) - 0.5_R8P*alpha*(qr(1:nv) - ql(1:nv))
      if (seth_euler_admissible(nv, G, ql(1:nv) - lam*fl(1:nv), EPS, EPS) .and. &
          seth_euler_admissible(nv, G, qr(1:nv) + lam*fl(1:nv), EPS, EPS)) cycle
      nfound = nfound + 1
      fh(1:nv) = 2._R8P*fl(1:nv)
      call seth_positivity_euler(nv, G, lam, EPS, EPS, ql(1:nv), qr(1:nv), fh(1:nv), fl(1:nv), f(1:nv), theta)
      if (theta /= 0._R8P .or. any(f(1:nv) /= fl(1:nv))) call fail('inadmissible f_low: theta /= 0')
   enddo
   if (nfound == 0) call fail('no inadmissible low-order flux found at lam = 10/alpha')
   endsubroutine five
endprogram seth_test_positivity
