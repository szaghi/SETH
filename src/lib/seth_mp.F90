!< SETH, the monotonicity-preserving bound of Suresh & Huynh (1997) on reconstructed face states.
module seth_mp
!< SETH, the monotonicity-preserving bound of Suresh & Huynh (1997) on reconstructed face states.
!<
!< A WENO or TENO state at a face is not bounded: near a discontinuity it can overshoot, and at high order (WENO-7/9,
!< WENO-Z9) it does. Suresh & Huynh (1997, J. Comput. Phys. 136:83) bound the left state at x_{j+1/2} by an interval
!< built from the five values v_{j-2} .. v_{j+2} (their Eqs. (2.8), (2.12), (2.16), (2.19), (2.24)-(2.29)):
!<
!<    d_k     = v_{k-1} - 2 v_k + v_{k+1}                     (k = j-1, j, j+1)
!<    dM4_+   = minmod(4 d_j - d_{j+1}, 4 d_{j+1} - d_j, d_j, d_{j+1})
!<    dM4_-   = minmod(4 d_j - d_{j-1}, 4 d_{j-1} - d_j, d_j, d_{j-1})
!<    v_UL    = v_j + alpha (v_j - v_{j-1})                   upper limit
!<    v_MD    = (v_j + v_{j+1})/2 - dM4_+/2                   median
!<    v_LC    = v_j + (v_j - v_{j-1})/2 + 4/3 dM4_-          large curvature
!<    v_min   = max(min(v_j, v_{j+1}, v_MD), min(v_j, v_UL, v_LC))
!<    v_max   = min(max(v_j, v_{j+1}, v_MD), max(v_j, v_UL, v_LC))
!<    v^L    <- median(v^L, v_min, v_max)
!<
!< On monotone data the interval is the monotone one (the state stays between v_j and v_{j+1}, after one stage between
!< v_{j-1} and v_j under the CFL bound 1/(1 + alpha)); near a smooth extremum the curvature terms widen it so that an
!< accurate state is left alone. Suresh & Huynh skip the bound where (v^L - v_j)(v^L - v_MP) <= 1e-10 (Eq. (2.30)),
!< v_MP = v_j + minmod(v_{j+1} - v_j, alpha (v_j - v_{j-1})): an absolute tolerance and a branch, for speed only (there
!< the median returns v^L unchanged). This kernel always applies the median: branch-free, and scale invariant.
!<
!<    call seth_mp_bound(vp, vm, alpha, rp, rm)
!<
!< vp, vm(-2:3): the data the states were reconstructed from (SETH's face convention, cells by offset from the face,
!< the right state mirrored internally); rp, rm: the states, bounded in place. alpha = SETH_MP_ALPHA = 4, Suresh &
!< Huynh's value for Runge-Kutta stepping (CFL <= 1/(1 + alpha) = 0.2 for the monotonicity proof). Pure,
!< `acc routine seq`, `omp declare target`.
use penf, only : I4P, R_P
implicit none
private
public :: SETH_MP_ALPHA
public :: seth_mp_bound

real(R_P), parameter :: SETH_MP_ALPHA = 4._R_P !< alpha of Suresh & Huynh 1997 for Runge-Kutta stepping (section 2).

contains
   pure subroutine seth_mp_bound(vp, vm, alpha, rp, rm)
   !< Bound the two states of a face by the monotonicity-preserving interval of Suresh & Huynh.
   real(R_P), intent(in)    :: vp(-2:3) !< Values reconstructed with a left bias.
   real(R_P), intent(in)    :: vm(-2:3) !< Values reconstructed with a right bias.
   real(R_P), intent(in)    :: alpha    !< Upper-limit factor (SETH_MP_ALPHA).
   real(R_P), intent(inout) :: rp       !< Left state at the face, from vp: bounded.
   real(R_P), intent(inout) :: rm       !< Right state at the face, from vm: bounded.
   !$acc routine seq
   !$omp declare target

   call seth_mp_side(vp(-2), vp(-1), vp(0), vp(1), vp(2), alpha, rp)
   call seth_mp_side(vm(3), vm(2), vm(1), vm(0), vm(-1), alpha, rm) ! mirrored: vm(1-m) plays vp(m)
   endsubroutine seth_mp_bound

   pure subroutine seth_mp_side(a, b, c, d, e, alpha, v)
   !< The interval from v_{j-2} .. v_{j+2} = a .. e, and the median of v with its ends.
   real(R_P), intent(in)    :: a, b, c, d, e !< v_{j-2}, v_{j-1}, v_j, v_{j+1}, v_{j+2}.
   real(R_P), intent(in)    :: alpha         !< Upper-limit factor.
   real(R_P), intent(inout) :: v             !< The state, bounded.
   real(R_P)                :: dm, d0, dp    !< Second differences at j-1, j, j+1.
   real(R_P)                :: m4p, m4m      !< dM4 at j+1/2 and j-1/2.
   real(R_P)                :: ul, md, lc    !< v_UL, v_MD, v_LC.
   real(R_P)                :: vmin, vmax    !< The interval.
   !$acc routine seq
   !$omp declare target

   dm = a - 2._R_P*b + c
   d0 = b - 2._R_P*c + d
   dp = c - 2._R_P*d + e
   m4p = minmod4(4._R_P*d0 - dp, 4._R_P*dp - d0, d0, dp)
   m4m = minmod4(4._R_P*d0 - dm, 4._R_P*dm - d0, d0, dm)
   ul = c + alpha*(c - b)
   md = 0.5_R_P*(c + d) - 0.5_R_P*m4p
   lc = c + 0.5_R_P*(c - b) + (4._R_P/3._R_P)*m4m
   vmin = max(min(c, d, md), min(c, ul, lc))
   vmax = min(max(c, d, md), max(c, ul, lc))
   v = max(vmin, min(v, vmax)) ! the median of v, vmin, vmax (vmin <= c <= vmax: both intervals contain c)
   endsubroutine seth_mp_side

   pure function minmod4(w, x, y, z) result(r)
   !< minmod of four numbers: the one of least magnitude if all have one sign, else 0 (Suresh & Huynh, Eq. (2.6)).
   real(R_P), intent(in) :: w, x, y, z
   real(R_P)             :: r
   !$acc routine seq
   !$omp declare target

   r = 0.125_R_P*(sign(1._R_P, w) + sign(1._R_P, x))*abs((sign(1._R_P, w) + sign(1._R_P, y)) &
       *(sign(1._R_P, w) + sign(1._R_P, z)))*min(abs(w), abs(x), abs(y), abs(z))
   endfunction minmod4
endmodule seth_mp
