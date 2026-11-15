!< SETH benchmark: the positivity limiter as first written (commit ea4959b), with automatic arrays.
module seth_bench_positivity_auto
!< SETH benchmark: the positivity limiter as first written (commit ea4959b), with automatic local arrays (a0(nv), ..)
!< and array expressions over nv. Kept verbatim (names apart) as the reference of the measurement that led to its
!< rewrite: on the device every call allocates them from the device heap (notes/gpu-kernel-audit.md).
use penf, only : I4P, R8P
implicit none
private
public :: positivity_euler_auto
public :: admissible_auto


contains
   pure subroutine positivity_euler_auto(nv, gamma, lam, eps_rho, eps_p, ql, qr, f_high, f_low, f, theta)
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
   if (admissible_auto(nv, gamma, a1, eps_rho, eps_p) .and. &
       admissible_auto(nv, gamma, b1, eps_rho, eps_p)) then
      theta = 1._R8P ; f = f_high ; return
   endif
   if (.not.(admissible_auto(nv, gamma, a0, eps_rho, eps_p) .and. &
             admissible_auto(nv, gamma, b0, eps_rho, eps_p))) then
      theta = 0._R8P ; f = f_low ; return ! the low-order flux itself fails: nothing to blend
   endif
   t = min(segment_bound(nv, gamma, eps_rho, eps_p, a0, a1), segment_bound(nv, gamma, eps_rho, eps_p, b0, b1))
   ! rounding guard on the flux actually returned: the closed form puts a state on the floor, where rounding can leave
   ! it a few ulp below; back off by relative steps 1e-15, 2e-15, 4e-15, ... (theta stays within rounding of its
   ! maximum), down to theta = 0 (f_low, admissible) at worst
   theta = t
   do it=0, 52
      f = theta*f_high + (1._R8P - theta)*f_low
      if (admissible_auto(nv, gamma, ql - lam*f, eps_rho, eps_p) .and. &
          admissible_auto(nv, gamma, qr + lam*f, eps_rho, eps_p)) return
      theta = max(t*(1._R8P - 1.e-15_R8P*2._R8P**it), 0._R8P)
   enddo
   theta = 0._R8P ; f = f_low
   endsubroutine positivity_euler_auto

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

   pure logical function admissible_auto(nv, gamma, q, eps_rho, eps_p) result(ok)
   !< rho >= eps_rho and p >= eps_p for the conservative state q = (rho, m(1:nv-2), E) of an ideal gas.
   integer(I4P), intent(in) :: nv
   real(R8P),    intent(in) :: gamma, q(nv), eps_rho, eps_p
   !$acc routine seq
   !$omp declare target
   ok = .false.
   if (.not.(q(1) >= eps_rho)) return
   ok = (gamma - 1._R8P)*(q(nv) - 0.5_R8P*dot_product(q(2:nv-1), q(2:nv-1))/q(1)) >= eps_p
   endfunction admissible_auto
endmodule seth_bench_positivity_auto
