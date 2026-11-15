!< SETH, positivity: the positivity-preserving flux limiter of Hu, Adams & Shu (2013) for the Euler equations.
module seth_positivity
!< SETH, positivity: the positivity-preserving flux limiter of Hu, Adams & Shu (2013) for the Euler equations.
!<
!< A high-order finite-difference (or finite-volume) scheme can drive density or pressure negative near vacuum and at
!< strong shocks. The cut-off flux limiter of Hu, Adams & Shu (2013, J. Comput. Phys. 242) keeps the scheme and limits
!< only the face fluxes that would: at each face
!<
!<    F = theta F_high + (1 - theta) F_low,   0 <= theta <= 1,
!<
!< with F_low a first-order positivity-preserving flux, and theta the largest value for which the two one-sided
!< states the face flux enters, q_l - lam F and q_r + lam F, keep rho >= eps_rho and p >= eps_p. Writing the
!< forward-Euler update of a cell as the mean of its two one-sided states, positive one-sided states make the update
!< positive; SSP Runge-Kutta stages are convex combinations of forward-Euler steps and keep it. The limited flux is a
!< convex combination of two consistent fluxes, hence consistent and conservative (HAS, section 2.3).
!<
!<    call seth_positivity_euler(nv, gamma, lam, eps_rho, eps_p, ql, qr, f_high, f_low, f, theta)
!<
!< - q = (rho, rho u_1, .., rho u_d, E), nv = d + 2 (1-D: 3, 2-D: 4, 3-D: 5), ideal gas p = (gamma - 1)(E - |m|^2/2rho);
!< - ql, qr: the states of the cells left and right of the face; f_high, f_low: the high-order and the low-order face
!<   fluxes, in the same variables;
!< - lam: 2 dt/h in 1-D. In d dimensions the update splits into the directions with partition weights w_k summing to 1
!<   (HAS section 2.5: equal halves as Zhang & Shu, or speed-weighted): lam = 2 dt/(w_k h_k) for the faces of direction
!<   k. The low-order flux must be positivity preserving for that lam: the Lax-Friedrichs flux
!<   (f(ql) + f(qr))/2 - alpha (qr - ql)/2 with alpha >= max(|u_n| + c) is, when lam alpha <= 1 (CFL <= 1/2 in 1-D);
!< - eps_rho, eps_p: the floors, HAS take min(1e-13, the initial minimum of rho, of p) (SETH_POSITIVITY_EPS = 1e-13);
!< - f: the limited flux; theta: its factor (1: f_high unchanged; 0: f_low). If even f_low gives a one-sided state
!<   below the floors (lam too large, or an inadmissible cell), theta = 0 and f = f_low: the caller can check its
!<   cells with seth_euler_admissible.
!<
!< theta in closed form: the density of a one-sided state is linear in theta, the bound is a ratio; for the pressure,
!< rho p/(gamma - 1) = rho E - |m|^2/2 is quadratic in theta, and since the admissible set is convex the admissible
!< thetas form an interval [0, theta*]: theta* is the root of rho E - |m|^2/2 - eps_p rho/(gamma - 1) in
!< (0, theta_rho] when that quadratic is negative at theta_rho (HAS solve the same quadratic). A rounding guard
!< re-checks the one-sided states of the flux actually returned and, if one is a few ulp below its floor, lowers theta
!< by relative steps from 1e-15 (doubling). Every routine is pure, `acc routine seq` and `omp declare target`.
use penf, only : I4P, R8P
implicit none
private
public :: SETH_POSITIVITY_EPS
public :: seth_positivity_euler
public :: seth_euler_admissible

real(R8P), parameter :: SETH_POSITIVITY_EPS = 1.e-13_R8P !< The floor of Hu, Adams & Shu (2013).

contains
   pure subroutine seth_positivity_euler(nv, gamma, lam, eps_rho, eps_p, ql, qr, f_high, f_low, f, theta)
   !< The positivity-preserving limited flux at one face (Hu, Adams & Shu 2013).
   integer(I4P), intent(in)  :: nv              !< Variables: rho, nv-2 momentum components, E.
   real(R8P),    intent(in)  :: gamma           !< Ratio of specific heats.
   real(R8P),    intent(in)  :: lam             !< 2 dt/h (2 dt/(w h) per direction in multi-D).
   real(R8P),    intent(in)  :: eps_rho, eps_p  !< Floors of density and pressure.
   real(R8P),    intent(in)  :: ql(nv), qr(nv)  !< States of the cells left and right of the face.
   real(R8P),    intent(in)  :: f_high(nv)      !< High-order face flux.
   real(R8P),    intent(in)  :: f_low(nv)       !< Low-order, positivity-preserving face flux.
   real(R8P),    intent(out) :: f(nv)           !< Limited face flux.
   real(R8P),    intent(out) :: theta           !< Limiting factor: f = theta f_high + (1 - theta) f_low.
   real(R8P)                 :: a0(nv), a1(nv), b0(nv), b1(nv), t
   integer(I4P)              :: it
   !$acc routine seq
   !$omp declare target

   ! the one-sided states with the low-order (0) and the high-order (1) flux: left q_l - lam F, right q_r + lam F
   a0 = ql - lam*f_low  ; a1 = ql - lam*f_high
   b0 = qr + lam*f_low  ; b1 = qr + lam*f_high
   if (seth_euler_admissible(nv, gamma, a1, eps_rho, eps_p) .and. &
       seth_euler_admissible(nv, gamma, b1, eps_rho, eps_p)) then
      theta = 1._R8P ; f = f_high ; return
   endif
   if (.not.(seth_euler_admissible(nv, gamma, a0, eps_rho, eps_p) .and. &
             seth_euler_admissible(nv, gamma, b0, eps_rho, eps_p))) then
      theta = 0._R8P ; f = f_low ; return ! the low-order flux itself fails: nothing to blend
   endif
   t = min(segment_bound(nv, gamma, eps_rho, eps_p, a0, a1), segment_bound(nv, gamma, eps_rho, eps_p, b0, b1))
   ! rounding guard on the flux actually returned: the closed form puts a state on the floor, where rounding can leave
   ! it a few ulp below; back off by relative steps 1e-15, 2e-15, 4e-15, ... (theta stays within rounding of its
   ! maximum), down to theta = 0 (f_low, admissible) at worst
   theta = t
   do it=0, 52
      f = theta*f_high + (1._R8P - theta)*f_low
      if (seth_euler_admissible(nv, gamma, ql - lam*f, eps_rho, eps_p) .and. &
          seth_euler_admissible(nv, gamma, qr + lam*f, eps_rho, eps_p)) return
      theta = max(t*(1._R8P - 1.e-15_R8P*2._R8P**it), 0._R8P)
   enddo
   theta = 0._R8P ; f = f_low
   endsubroutine seth_positivity_euler

   pure function segment_bound(nv, gamma, eps_rho, eps_p, s0, s1) result(t)
   !< The largest t in [0, 1] with s0 + tau (s1 - s0) admissible for every tau in [0, t] (s0 admissible).
   integer(I4P), intent(in) :: nv
   real(R8P),    intent(in) :: gamma, eps_rho, eps_p, s0(nv), s1(nv)
   real(R8P)                :: t, d(nv), kappa, qa, qb, qc, disc, r1, r2
   !$acc routine seq
   !$omp declare target

   d = s1 - s0
   ! density: linear
   t = 1._R8P
   if (s1(1) < eps_rho) t = (s0(1) - eps_rho)/(s0(1) - s1(1))
   ! pressure: G(tau) = rho E - |m|^2/2 - kappa rho, kappa = eps_p/(gamma - 1), G = qa tau^2 + qb tau + qc, qc >= 0
   kappa = eps_p/(gamma - 1._R8P)
   qa = d(1)*d(nv) - 0.5_R8P*dot_product(d(2:nv-1), d(2:nv-1))
   qb = s0(1)*d(nv) + s0(nv)*d(1) - dot_product(s0(2:nv-1), d(2:nv-1)) - kappa*d(1)
   qc = s0(1)*s0(nv) - 0.5_R8P*dot_product(s0(2:nv-1), s0(2:nv-1)) - kappa*s0(1)
   if (qa*t*t + qb*t + qc >= 0._R8P) return ! admissible up to the density bound (the set is convex)
   ! a sign change in (0, t): the smallest positive root
   if (abs(qa) <= epsilon(1._R8P)*(abs(qb) + abs(qc))) then
      t = min(t, -qc/qb)
   else
      disc = max(qb*qb - 4._R8P*qa*qc, 0._R8P)
      r1 = -0.5_R8P*(qb + sign(sqrt(disc), qb)) ! stable: no cancellation
      r2 = qc/r1                                 ! the two roots: r1/qa and qc/r1
      r1 = r1/qa
      t = min(t, smallest_in(r1, r2, t))
   endif
   t = max(t, 0._R8P)
   endfunction segment_bound

   pure function smallest_in(r1, r2, tmax) result(t)
   !< The smallest of r1, r2 in [0, tmax] (tmax if neither is).
   real(R8P), intent(in) :: r1, r2, tmax
   real(R8P)             :: t
   !$acc routine seq
   !$omp declare target
   t = tmax
   if (r1 >= 0._R8P .and. r1 <= t) t = r1
   if (r2 >= 0._R8P .and. r2 <= t) t = r2
   endfunction smallest_in

   pure logical function seth_euler_admissible(nv, gamma, q, eps_rho, eps_p) result(ok)
   !< rho >= eps_rho and p >= eps_p for the conservative state q = (rho, m(1:nv-2), E) of an ideal gas.
   integer(I4P), intent(in) :: nv
   real(R8P),    intent(in) :: gamma, q(nv), eps_rho, eps_p
   !$acc routine seq
   !$omp declare target
   ok = .false.
   if (.not.(q(1) >= eps_rho)) return
   ok = (gamma - 1._R8P)*(q(nv) - 0.5_R8P*dot_product(q(2:nv-1), q(2:nv-1))/q(1)) >= eps_p
   endfunction seth_euler_admissible
endmodule seth_positivity
