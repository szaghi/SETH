!< SETH, WENO-6 point interpolation: the value at an off-node point of a line of 6 equispaced point values.
module seth_point_weno6
!< SETH, WENO-6 point interpolation: the value at an off-node point of a line of 6 equispaced point values.
!<
!< Unlike SETH's face kernels, the target is any quarter of the node spacing between the nodes (and a quarter beyond
!< the end nodes), on one array of point values:
!<
!<    v = seth_weno6_point(u, t4)
!<
!< `u(1:6)` the values at 6 equispaced nodes, node k at x = k; `t4` the target at x = t4/4 node spacings from node 0,
!< t4 = 3..25 (odd: the quarter positions, x = 3/4..25/4; 2 mod 4: the midpoints between two nodes; 0 mod 4: a node,
!< whose value is returned unchanged). Other values of t4 are outside the contract (a bounds-checked build stops on the
!< table index). The value is a convex combination of the cubic interpolants of the 4-point sub-stencils s = nodes
!< s..s+3 that contain the target's interval (the nodes m, m+1, m = clamp(floor(t4/4), 1, 5)) and no jump of
!< seth_weno6_jumps: WENO-Z weights (Borges, Carmona, Costa & Don 2008, power 2) on the Jiang-Shu indicators (Jiang &
!< Shu 1996) of each cubic over its central interval, tau the highest undivided difference of the candidates' union
!< squared. On smooth data the weights tend to the linear ones, which make the combination the Lagrange interpolant of
!< the candidates' union: order 6 with three candidates (x = 13/4, 14/4, 15/4), 5 with two, 4 with one (toward the
!< ends; at t4 = 3 and 25 the one candidate is the one-sided cubic of the end interval). With no candidate left (a jump
!< inside the target's interval) the value is the troubled-cell fallback: at a midpoint the mean of the interval's two
!< values, elsewhere the linear interpolation on the interval clamped to its two values.
!<
!< ADAM's 2:1 AMR seam transfer is this kernel (ADAM issue 68, C3 and D2): its restriction at the midpoints, its
!< coarse-to-fine fill at the quarters; the R8P results equal ADAM's routines bit for bit
!< (seth_test_weno6_point_adam). The weights are in seth_point_weno6_tables (scripts/gen_weno6_tables.py), the guard
!< SETH_WENO6_EPS in seth_kind_parameters. Every procedure is pure, `acc routine seq` and `omp declare target`.
use penf, only : I4P, R_P
use seth_kind_parameters, only : SETH_WENO6_EPS
use seth_point_weno6_tables
implicit none
private
public :: SETH_WENO6_JUMP_RATIO
public :: seth_weno6_jumps
public :: seth_weno6_point

real(R_P), parameter :: SETH_WENO6_JUMP_RATIO = 2._R_P !< An interval whose undivided difference exceeds this ratio
                                                       !< times its neighbours' is a jump (ADAM's value). Resolved
                                                       !< smooth data never trips it: the neighbours' differences
                                                       !< match the interval's up to O(h^3).

contains
   pure subroutine seth_weno6_jumps(u, jump)
   !< Flag the jumps of 6 values: interval i (nodes i, i+1) is a jump when its undivided difference exceeds
   !< SETH_WENO6_JUMP_RATIO times the larger of its neighbours' (at the ends, the inner neighbour and its linear
   !< extrapolation).
   real(R_P),    intent(in)  :: u(1:6)    !< Values at the nodes.
   logical,      intent(out) :: jump(1:5) !< Jump flag per interval.
   real(R_P)                 :: d(1:5)    !< Undivided first differences.
   real(R_P)                 :: ref       !< Neighbours' difference.
   integer(I4P)              :: i         !< Counter.
   !$acc routine seq
   !$omp declare target

   do i=1, 5
      d(i) = u(i+1) - u(i)
   enddo
   do i=1, 5
      if (i == 1) then
         ref = max(abs(d(2)), abs(2._R_P * d(2) - d(3)))
      elseif (i == 5) then
         ref = max(abs(d(4)), abs(2._R_P * d(4) - d(3)))
      else
         ref = max(abs(d(i-1)), abs(d(i+1)))
      endif
      jump(i) = abs(d(i)) > SETH_WENO6_JUMP_RATIO * ref
   enddo
   endsubroutine seth_weno6_jumps

   pure function seth_weno6_point(u, t4) result(v)
   !< The WENO-6 interpolation of the 6 values u at x = t4/4 (node k at x = k, t4 = 3..25); a node returns its value.
   !< The operations are those of ADAM's seam_weno6_core, in its order: jumps, candidates, the one-candidate shortcut,
   !< tau, beta = c^2 + 13/12 a^2 + 781/720 t^2 (c, a, t the cubic's central first, mean second and third undivided
   !< differences), alpha = d (1 + (tau/(beta + eps))^2), the normalised sum. R4P: the weights are scaled by the
   !< smallest beta + eps when tau exceeds it (seth_kind_parameters, SETH_WENO6_EPS).
   real(R_P),    intent(in) :: u(1:6) !< Values at the nodes.
   integer(I4P), intent(in) :: t4     !< Target at x = t4/4 node spacings from node 0 (3..25).
   real(R_P)                :: v      !< Interpolated value.
   real(R_P)                :: w(1:6) !< Work array of the highest difference.
   real(R_P)                :: tau    !< WENO-Z global indicator.
   real(R_P)                :: c, a, t !< Central first, mean second and third differences of a sub-stencil.
   real(R_P)                :: alpha  !< Unnormalised weight.
   real(R_P)                :: asum   !< Sum of the unnormalised weights.
#ifdef _R_P_IS_R4P
   real(R_P)                :: b(1:3) !< beta + eps of each candidate.
   real(R_P)                :: bmin   !< Its smallest value.
#endif
   logical                  :: jump(1:5) !< Jump flag per interval.
   integer(I4P)             :: i, m, s   !< Counters.
   integer(I4P)             :: lo, hi    !< First and last candidate.
   !$acc routine seq
   !$omp declare target

   if (mod(t4, 4_I4P) == 0_I4P) then ! a node
      v = u(t4 / 4_I4P)
      return
   endif
   call seth_weno6_jumps(u=u, jump=jump)
   lo = 0_I4P
   hi = 0_I4P
   do s=1, 3 ! the candidates are contiguous: a sub-stencil between two jump-free ones is jump-free
      if (WENO6_POINT_D(s,t4) > 0._R_P .and. .not.(jump(s) .or. jump(s+1) .or. jump(s+2))) then
         if (lo == 0_I4P) lo = s
         hi = s
      endif
   enddo
   if (lo == 0_I4P) then ! the troubled-cell fallback on the target's interval
      m = min(max(t4 / 4_I4P, 1_I4P), 5_I4P)
      if (mod(t4, 4_I4P) == 2_I4P) then ! a midpoint: the mean (ADAM's restriction)
         v = 0.5_R_P * (u(m) + u(m+1))
      else                              ! a quarter: the clamped linear interpolation (ADAM's fill)
         v = u(m) + 0.25_R_P * real(t4 - 4_I4P * m, R_P) * (u(m+1) - u(m))
         v = min(max(v, min(u(m), u(m+1))), max(u(m), u(m+1)))
      endif
      return
   endif
   if (lo == hi) then
      v = WENO6_POINT_W(1,lo,t4) * u(lo) + WENO6_POINT_W(2,lo,t4) * u(lo+1) + WENO6_POINT_W(3,lo,t4) * u(lo+2) + &
          WENO6_POINT_W(4,lo,t4) * u(lo+3)
      return
   endif
   do i=1, hi - lo + 4 ! nodes lo..hi+3
      w(i) = u(lo+i-1)
   enddo
   do m=1, hi - lo + 3
      do i=1, hi - lo + 4 - m
         w(i) = w(i+1) - w(i)
      enddo
   enddo
   tau = w(1) * w(1)
   v = 0._R_P
   asum = 0._R_P
#ifdef _R_P_IS_R4P
   ! R4P: (tau/(beta + eps))**2 overflows once tau/(beta + eps) passes 1.8e19, i.e. next to a flat sub-stencil
   ! (beta = 0, eps = 1e-30) once tau passes 2e-11: where tau exceeds the smallest beta + eps, bmin, every alpha is
   ! multiplied by (bmin/tau)**2, so each factor lies in [0, 1] and the normalised weights are unchanged
   bmin = huge(1._R_P)
   do s=lo, hi
      c = u(s+2) - u(s+1)
      a = 0.5_R_P * (u(s+3) - u(s+2) - u(s+1) + u(s))
      t = u(s+3) - 3._R_P * u(s+2) + 3._R_P * u(s+1) - u(s)
      b(s) = c * c + 13._R_P / 12._R_P * a * a + 781._R_P / 720._R_P * t * t + SETH_WENO6_EPS
      bmin = min(bmin, b(s))
   enddo
   do s=lo, hi
      if (tau > bmin) then
         alpha = WENO6_POINT_D(s,t4) * ((bmin / tau)**2 + (bmin / b(s))**2)
      else
         alpha = WENO6_POINT_D(s,t4) * (1._R_P + (tau / b(s))**2)
      endif
      v = v + alpha * (WENO6_POINT_W(1,s,t4) * u(s) + WENO6_POINT_W(2,s,t4) * u(s+1) + &
                       WENO6_POINT_W(3,s,t4) * u(s+2) + WENO6_POINT_W(4,s,t4) * u(s+3))
      asum = asum + alpha
   enddo
#else
   do s=lo, hi
      c = u(s+2) - u(s+1)
      a = 0.5_R_P * (u(s+3) - u(s+2) - u(s+1) + u(s))
      t = u(s+3) - 3._R_P * u(s+2) + 3._R_P * u(s+1) - u(s)
      alpha = WENO6_POINT_D(s,t4) * (1._R_P + (tau / (c * c + 13._R_P / 12._R_P * a * a + &
                                                   781._R_P / 720._R_P * t * t + SETH_WENO6_EPS))**2)
      v = v + alpha * (WENO6_POINT_W(1,s,t4) * u(s) + WENO6_POINT_W(2,s,t4) * u(s+1) + &
                       WENO6_POINT_W(3,s,t4) * u(s+2) + WENO6_POINT_W(4,s,t4) * u(s+3))
      asum = asum + alpha
   enddo
#endif
   v = v / asum
   endfunction seth_weno6_point
endmodule seth_point_weno6
