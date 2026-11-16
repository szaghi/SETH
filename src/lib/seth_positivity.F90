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
!<
!< Round-off: the floors hold for the one-sided states as this routine evaluates them (q -+ lam f, in a loop). A
!< caller that evaluates them differently (an array expression, a fused multiply-add) gets states that differ by
!< epsilon (|q| + |lam f|), and pressures that differ by about that times (1 + |u|^2): choose eps_rho and eps_p above
!< it. With HAS's 1e-13, a state of density 1e-7 at |u| = 400 next to a cell 1e9 times denser can come out with a
!< pressure of -1e-12 in the caller's arithmetic (seth_test_positivity, adm_ro).
!<
!< Device rule (notes/gpu-kernel-audit.md): no local array has a bound that depends on a dummy argument, and no array
!< expression makes a temporary of such a size. Automatic arrays in a device routine are allocated from the device heap
!< at every call: on an RTX 4070 that made this kernel 300 times slower, and above the default heap size the launch
!< failed. The locals are dimensioned SETH_POSITIVITY_NV_MAX = 5 and the loops run to nv; nv > 5 is not supported (the
!< routine then returns theta = 0 and f = f_low).
use penf, only : I4P, R_P
use seth_kind_parameters, only : SETH_POSITIVITY_BACKOFF
implicit none
private
public :: SETH_POSITIVITY_EPS
public :: SETH_POSITIVITY_NV_MAX
public :: seth_positivity_euler
public :: seth_euler_admissible

real(R_P),    parameter :: SETH_POSITIVITY_EPS = 1.e-13_R_P !< The floor of Hu, Adams & Shu (2013).
integer(I4P), parameter :: SETH_POSITIVITY_NV_MAX = 5_I4P   !< Largest nv (3-D Euler): the dimension of the locals.
integer(I4P), parameter :: NVM = SETH_POSITIVITY_NV_MAX

contains
   pure subroutine seth_positivity_euler(nv, gamma, lam, eps_rho, eps_p, ql, qr, f_high, f_low, f, theta)
   !< The positivity-preserving limited flux at one face (Hu, Adams & Shu 2013).
   integer(I4P), intent(in)  :: nv              !< Variables: rho, nv-2 momentum components, E (nv <= 5).
   real(R_P),    intent(in)  :: gamma           !< Ratio of specific heats.
   real(R_P),    intent(in)  :: lam             !< 2 dt/h (2 dt/(w h) per direction in multi-D).
   real(R_P),    intent(in)  :: eps_rho, eps_p  !< Floors of density and pressure.
   real(R_P),    intent(in)  :: ql(nv), qr(nv)  !< States of the cells left and right of the face.
   real(R_P),    intent(in)  :: f_high(nv)      !< High-order face flux.
   real(R_P),    intent(in)  :: f_low(nv)       !< Low-order, positivity-preserving face flux.
   real(R_P),    intent(out) :: f(nv)           !< Limited face flux.
   real(R_P),    intent(out) :: theta           !< Limiting factor: f = theta f_high + (1 - theta) f_low.
   real(R_P)                 :: a0(NVM), a1(NVM), b0(NVM), b1(NVM), t
   integer(I4P)              :: it, k
   !$acc routine seq
   !$omp declare target

   if (nv > NVM) then ! not supported: the low-order flux
      theta = 0._R_P
      do k=1, nv
         f(k) = f_low(k)
      enddo
      return
   endif
   ! the one-sided states with the low-order (0) and the high-order (1) flux: left q_l - lam F, right q_r + lam F
   do k=1, nv
      a0(k) = ql(k) - lam*f_low(k) ; a1(k) = ql(k) - lam*f_high(k)
      b0(k) = qr(k) + lam*f_low(k) ; b1(k) = qr(k) + lam*f_high(k)
   enddo
   if (seth_euler_admissible(nv, gamma, a1, eps_rho, eps_p) .and. &
       seth_euler_admissible(nv, gamma, b1, eps_rho, eps_p)) then
      theta = 1._R_P
      do k=1, nv
         f(k) = f_high(k)
      enddo
      return
   endif
   if (.not.(seth_euler_admissible(nv, gamma, a0, eps_rho, eps_p) .and. &
             seth_euler_admissible(nv, gamma, b0, eps_rho, eps_p))) then
      theta = 0._R_P ! the low-order flux itself fails: nothing to blend
      do k=1, nv
         f(k) = f_low(k)
      enddo
      return
   endif
   t = min(segment_bound(nv, gamma, eps_rho, eps_p, a0, a1), segment_bound(nv, gamma, eps_rho, eps_p, b0, b1))
   ! rounding guard on the flux actually returned: the closed form puts a state on the floor, where rounding can leave
   ! it a few ulp below; back off by relative steps 1e-15, 2e-15, 4e-15, ... (theta stays within rounding of its
   ! maximum), down to theta = 0 (f_low, admissible) at worst
   theta = t
   do it=0, 52
      do k=1, nv
         f(k) = theta*f_high(k) + (1._R_P - theta)*f_low(k)
      enddo
      ! a loop of its own, on the stored f: fused with the one above, nvfortran -fast contracts the two statements and
      ! checks states of an unrounded flux, not of the one returned (seth_test_positivity then fails)
      do k=1, nv
         a1(k) = ql(k) - lam*f(k) ; b1(k) = qr(k) + lam*f(k)
      enddo
      if (seth_euler_admissible(nv, gamma, a1, eps_rho, eps_p) .and. &
          seth_euler_admissible(nv, gamma, b1, eps_rho, eps_p)) return
      theta = max(t*(1._R_P - SETH_POSITIVITY_BACKOFF*2._R_P**it), 0._R_P)
   enddo
   theta = 0._R_P
   do k=1, nv
      f(k) = f_low(k)
   enddo
   endsubroutine seth_positivity_euler

   pure function segment_bound(nv, gamma, eps_rho, eps_p, s0, s1) result(t)
   !< The largest t in [0, 1] with s0 + tau (s1 - s0) admissible for every tau in [0, t] (s0 admissible).
   integer(I4P), intent(in) :: nv
   real(R_P),    intent(in) :: gamma, eps_rho, eps_p, s0(NVM), s1(NVM)
   real(R_P)                :: t, d(NVM), kappa, qa, qb, qc, disc, r1, r2, dd, sd, ss
   integer(I4P)             :: k
   !$acc routine seq
   !$omp declare target

   do k=1, nv
      d(k) = s1(k) - s0(k)
   enddo
   ! density: linear
   t = 1._R_P
   if (s1(1) < eps_rho) t = (s0(1) - eps_rho)/(s0(1) - s1(1))
   ! pressure: G(tau) = rho E - |m|^2/2 - kappa rho, kappa = eps_p/(gamma - 1), G = qa tau^2 + qb tau + qc, qc >= 0
   kappa = eps_p/(gamma - 1._R_P)
   dd = 0._R_P ; sd = 0._R_P ; ss = 0._R_P
   do k=2, nv-1
      dd = dd + d(k)*d(k) ; sd = sd + s0(k)*d(k) ; ss = ss + s0(k)*s0(k)
   enddo
   qa = d(1)*d(nv) - 0.5_R_P*dd
   qb = s0(1)*d(nv) + s0(nv)*d(1) - sd - kappa*d(1)
   qc = s0(1)*s0(nv) - 0.5_R_P*ss - kappa*s0(1)
   if (qa*t*t + qb*t + qc >= 0._R_P) return ! admissible up to the density bound (the set is convex)
   ! a sign change in (0, t): the smallest positive root
   if (abs(qa) <= epsilon(1._R_P)*(abs(qb) + abs(qc))) then
      t = min(t, -qc/qb)
   else
      disc = max(qb*qb - 4._R_P*qa*qc, 0._R_P)
      r1 = -0.5_R_P*(qb + sign(sqrt(disc), qb)) ! stable: no cancellation
      r2 = qc/r1                                 ! the two roots: r1/qa and qc/r1
      r1 = r1/qa
      t = min(t, smallest_in(r1, r2, t))
   endif
   t = max(t, 0._R_P)
   endfunction segment_bound

   pure function smallest_in(r1, r2, tmax) result(t)
   !< The smallest of r1, r2 in [0, tmax] (tmax if neither is).
   real(R_P), intent(in) :: r1, r2, tmax
   real(R_P)             :: t
   !$acc routine seq
   !$omp declare target
   t = tmax
   if (r1 >= 0._R_P .and. r1 <= t) t = r1
   if (r2 >= 0._R_P .and. r2 <= t) t = r2
   endfunction smallest_in

   pure logical function seth_euler_admissible(nv, gamma, q, eps_rho, eps_p) result(ok)
   !< rho >= eps_rho and p >= eps_p for the conservative state q = (rho, m(1:nv-2), E) of an ideal gas.
   integer(I4P), intent(in) :: nv
   real(R_P),    intent(in) :: gamma, q(nv), eps_rho, eps_p
   real(R_P)                :: ke
   integer(I4P)             :: k
   !$acc routine seq
   !$omp declare target
   ok = .false.
   if (.not.(q(1) >= eps_rho)) return
   ke = 0._R_P
   do k=2, nv-1
      ke = ke + q(k)*q(k)
   enddo
   ok = (gamma - 1._R_P)*(q(nv) - 0.5_R_P*ke/q(1)) >= eps_p
   endfunction seth_euler_admissible
endmodule seth_positivity
