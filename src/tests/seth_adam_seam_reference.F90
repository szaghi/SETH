!< SETH tests: ADAM's seam WENO-6 restriction and fill, verbatim, as the reference of seth_weno6_point (frozen history).
module seth_adam_seam_reference
!< SETH tests: ADAM's seam WENO-6 restriction and fill, verbatim, as the reference of seth_weno6_point (frozen history).
!<
!< From ADAM commit 2ea71524, src/lib/common/adam_seam_interpolation_library.F90 (ADAM issue 68, C3 and D2), copied by
!< a script without any edit: the two constants the routines read (lines 74 and 158), the WENO tables with their
!< comments and the module's `contains` (lines 173-235: SEAM_W_CUBIC_MID, SEAM_D_WENO6, SEAM_WENO6_JUMP_RATIO,
!< SEAM_WENO6_EPS, SEAM_W_CUBIC_Q, SEAM_D_FILL_WENO6) and the routines seam_weno6_jumps, seam_weno6_core,
!< seam_restriction_weno6_1d and seam_fill_weno6_1d (lines 512-663); only the module around them is new. The reference
!< of seth_test_weno6_point_adam (SETH issue 7, E2). ADAM's step S7 (ADAM issue 139) replaces these routines with
!< SETH's seth_weno6_point and seth_weno6_jumps; the copy stays as it was copied, the fixed reference of the
!< equivalence test. R8P only.
use penf, only : I4P, R8P
implicit none
private
public :: SEAM_SUB_COINCIDENT, SEAM_RESTRICT_IDENTITY
public :: SEAM_W_CUBIC_MID, SEAM_D_WENO6, SEAM_WENO6_JUMP_RATIO, SEAM_WENO6_EPS, SEAM_W_CUBIC_Q, SEAM_D_FILL_WENO6
public :: seam_weno6_jumps, seam_weno6_core, seam_restriction_weno6_1d, seam_fill_weno6_1d

integer(I4P), parameter :: SEAM_SUB_COINCIDENT  = 3_I4P !< Sub-position of an unrefined direction (eta = 0, issue #46).
integer(I4P), parameter :: SEAM_RESTRICT_IDENTITY  = 5_I4P !< Per-direction code of an unrefined direction (issue #46).

! The WENO restriction (issue #68 C3): per direction a nonlinear convex combination of the cubic interpolants of the
! 4-point sub-stencils of the 6 nodes (sub-stencil s = nodes s..s+3, s = 1..3) that contain the coarse centre, i.e. the
! interval between nodes code+1 and code+2. Sub-stencils that would extrapolate are never used: across a jump they
! import the far side's state. The linear weights reproduce the Lagrange interpolant of the candidates' union, so the
! design order is 6 for code 2 (three candidates, the 6-point Lagrange of `SEAM_W_LAGRANGE6`), 5 for codes 1 and 3
! (two candidates, the 5-point Lagrange) and 4 for codes 0 and 4 (one candidate, at a block edge).
real(R8P), parameter :: SEAM_W_CUBIC_MID(1:4,1:3) = reshape([ &
    5._R8P,  15._R8P,  -5._R8P,   1._R8P,  & ! the midpoint of the sub-stencil's interval 1 (nodes 1, 2)
   -1._R8P,   9._R8P,   9._R8P,  -1._R8P,  & ! interval 2 (nodes 2, 3), centred
    1._R8P,  -5._R8P,  15._R8P,   5._R8P], & ! interval 3 (nodes 3, 4)
   [4,3]) / 16._R8P
   !< Cubic Lagrange weights at the midpoint of each interval of a 4-point sub-stencil (node 1:4, interval 1:3).

real(R8P), parameter :: SEAM_D_WENO6(1:3,0:4) = reshape([ &
   16._R8P,   0._R8P,   0._R8P,  & ! code 0: sub-stencil 1 only
   10._R8P,   6._R8P,   0._R8P,  & ! code 1: sub-stencils 1, 2 -> 5-point Lagrange on nodes 1..5
    3._R8P,  10._R8P,   3._R8P,  & ! code 2: sub-stencils 1, 2, 3 -> 6-point Lagrange (centred)
    0._R8P,   6._R8P,  10._R8P,  & ! code 3: sub-stencils 2, 3 -> 5-point Lagrange on nodes 2..6
    0._R8P,   0._R8P,  16._R8P], & ! code 4: sub-stencil 3 only
   [3,5]) / 16._R8P
   !< Linear weights of the sub-stencils (sub-stencil 1:3, code 0:4); zero marks a sub-stencil that is not a candidate.

real(R8P), parameter :: SEAM_WENO6_JUMP_RATIO = 2._R8P
   !< Jump detector of the WENO restriction: an interval whose undivided difference exceeds this ratio times its
   !< neighbours' is a jump. Smooth data never trips it once resolved (the neighbours average to the interval's
   !< difference up to O(h^3)); a candidate containing a jump is dropped, and with no candidate left the restriction
   !< falls back to the mean of the fine pair (the troubled-cell fallback).
real(R8P), parameter :: SEAM_WENO6_EPS = 1.e-40_R8P !< Guard of the WENO-Z weights against 0/0 on constant data.

! The WENO coarse->fine fill (issue #68 D2): the target, anchor node p + eta (eta = -+1/4), sits at t4/4 in node units,
! t4 = 4 p -+ 1 (3..25, target index j = (t4 - 1)/2 = 1..12); its candidates are the sub-stencils containing the
! interval m = clamp(floor(t4/4), 1, 5) (at the block edge, j = 1 or 12, the one candidate extrapolates a quarter cell,
! the tricubic's own one-sided cubic).
real(R8P), parameter :: SEAM_W_CUBIC_Q(1:4,1:8) = reshape([ &
   195._R8P, -117._R8P,   65._R8P,  -15._R8P,  & ! r = 1: a quarter node spacing before node 1 (block edge)
    77._R8P,   77._R8P,  -33._R8P,    7._R8P,  & ! r = 2: node 1 + 1/4
    15._R8P,  135._R8P,  -27._R8P,    5._R8P,  & ! r = 3: node 1 + 3/4
    -7._R8P,  105._R8P,   35._R8P,   -5._R8P,  & ! r = 4: node 2 + 1/4
    -5._R8P,   35._R8P,  105._R8P,   -7._R8P,  & ! r = 5: node 2 + 3/4
     5._R8P,  -27._R8P,  135._R8P,   15._R8P,  & ! r = 6: node 3 + 1/4
     7._R8P,  -33._R8P,   77._R8P,   77._R8P,  & ! r = 7: node 3 + 3/4
   -15._R8P,   65._R8P, -117._R8P,  195._R8P], & ! r = 8: node 4 + 1/4 (block edge)
   [4,8]) / 128._R8P
   !< Cubic Lagrange weights of a 4-point sub-stencil at its quarter positions (node 1:4, r = (t4 - 4 s + 3)/2 = 1:8).

real(R8P), parameter :: SEAM_D_FILL_WENO6(1:3,1:12) = reshape([ &
   320._R8P,   0._R8P,   0._R8P,  & ! j = 1, t = 3/4 (block edge): sub-stencil 1 only
   320._R8P,   0._R8P,   0._R8P,  & ! j = 2, t = 5/4
   320._R8P,   0._R8P,   0._R8P,  & ! j = 3, t = 7/4
   220._R8P, 100._R8P,   0._R8P,  & ! j = 4, t = 9/4: sub-stencils 1, 2 -> 5-point Lagrange
   180._R8P, 140._R8P,   0._R8P,  & ! j = 5, t = 11/4
    77._R8P, 198._R8P,  45._R8P,  & ! j = 6, t = 13/4: sub-stencils 1, 2, 3 -> 6-point Lagrange (centred, eta = +1/4)
    45._R8P, 198._R8P,  77._R8P,  & ! j = 7, t = 15/4: (centred, eta = -1/4)
     0._R8P, 140._R8P, 180._R8P,  & ! j = 8, t = 17/4
     0._R8P, 100._R8P, 220._R8P,  & ! j = 9, t = 19/4
     0._R8P,   0._R8P, 320._R8P,  & ! j = 10, t = 21/4
     0._R8P,   0._R8P, 320._R8P,  & ! j = 11, t = 23/4
     0._R8P,   0._R8P, 320._R8P], & ! j = 12, t = 25/4 (block edge)
   [3,12]) / 320._R8P
   !< Linear weights of the fill's sub-stencils (sub-stencil 1:3, target 1:12), all positive (computed with exact
   !< fractions); zero marks a sub-stencil that does not contain the target's interval.

contains
   pure subroutine seam_weno6_jumps(u, jump)
   !< Flag the jumps of 6 values (issue #68 C3): interval i (nodes i, i+1) is a jump when its undivided difference
   !< exceeds `SEAM_WENO6_JUMP_RATIO` times the larger of its neighbours' (at the footprint ends, the inner neighbour
   !< and its linear extrapolation). Smooth data never trips it once resolved.
   real(R8P), intent(in)  :: u(1:6)    !< Values at the footprint nodes.
   logical,   intent(out) :: jump(1:5) !< Jump flag per interval.
   real(R8P)              :: d(1:5)    !< Undivided first differences.
   real(R8P)              :: ref       !< Neighbours' difference.
   integer(I4P)           :: i         !< Counter.
   !$acc routine seq
   !$omp declare target

   do i=1, 5
      d(i) = u(i+1) - u(i)
   enddo
   do i=1, 5
      if (i == 1) then
         ref = max(abs(d(2)), abs(2._R8P * d(2) - d(3)))
      elseif (i == 5) then
         ref = max(abs(d(4)), abs(2._R8P * d(4) - d(3)))
      else
         ref = max(abs(d(i-1)), abs(d(i+1)))
      endif
      jump(i) = abs(d(i)) > SEAM_WENO6_JUMP_RATIO * ref
   enddo
   endsubroutine seam_weno6_jumps

   pure function seam_weno6_core(u, dlin, wsub, fallback) result(value_)
   !< The WENO kernel shared by the restriction and the fill (issue #68 C3, D2): a convex combination of the cubic
   !< interpolants of the candidate 4-point sub-stencils (sub-stencil s = nodes s..s+3 of the 6) at one target.
   !< 1. Jumps: interval i (nodes i, i+1) is a jump when its undivided difference exceeds `SEAM_WENO6_JUMP_RATIO` times
   !<    the larger of its neighbours' (at the footprint ends, the inner neighbour and its linear extrapolation).
   !< 2. Candidates: the sub-stencils with a positive linear weight `dlin(s)` and free of jumps; none left: `fallback`
   !<    (the troubled-cell value, chosen by the caller).
   !< 3. WENO-Z weights (Borges et al. 2008, p = 2) on the Jiang-Shu indicators of each cubic over its central
   !<    interval, beta = c^2 + 13/12 a^2 + 781/720 t^2 (c, a, t its central first, mean second and third differences),
   !<    with tau the highest undivided difference of the candidates' union squared, O(h^8..h^12) on smooth data, so the
   !<    weights reach the linear ones fast enough to keep the design order at critical points.
   real(R8P), intent(in) :: u(1:6)        !< Values at the footprint nodes.
   real(R8P), intent(in) :: dlin(1:3)     !< Linear weight per sub-stencil (0: not a candidate).
   real(R8P), intent(in) :: wsub(1:4,1:3) !< Cubic weights of each candidate sub-stencil at the target.
   real(R8P), intent(in) :: fallback      !< Value returned when no candidate is left.
   real(R8P)             :: value_        !< WENO value.
   real(R8P)             :: w(1:6)        !< Work array of the highest difference.
   real(R8P)             :: tau           !< WENO-Z global indicator.
   real(R8P)             :: c, a, t       !< Central first, mean second and third differences of a sub-stencil.
   real(R8P)             :: alpha         !< Unnormalised weight.
   real(R8P)             :: asum          !< Sum of the unnormalised weights.
   logical               :: jump(1:5)     !< Jump flag per interval.
   integer(I4P)          :: i, m, s       !< Counters.
   integer(I4P)          :: lo, hi        !< First and last candidate.
   !$acc routine seq
   !$omp declare target

   call seam_weno6_jumps(u=u, jump=jump)
   lo = 0_I4P
   hi = 0_I4P
   do s=1, 3 ! the candidates are contiguous: a sub-stencil between two jump-free ones is jump-free
      if (dlin(s) > 0._R8P .and. .not.(jump(s) .or. jump(s+1) .or. jump(s+2))) then
         if (lo == 0_I4P) lo = s
         hi = s
      endif
   enddo
   if (lo == 0_I4P) then
      value_ = fallback
      return
   endif
   if (lo == hi) then
      value_ = wsub(1,lo) * u(lo) + wsub(2,lo) * u(lo+1) + wsub(3,lo) * u(lo+2) + wsub(4,lo) * u(lo+3)
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
   value_ = 0._R8P
   asum = 0._R8P
   do s=lo, hi
      c = u(s+2) - u(s+1)
      a = 0.5_R8P * (u(s+3) - u(s+2) - u(s+1) + u(s))
      t = u(s+3) - 3._R8P * u(s+2) + 3._R8P * u(s+1) - u(s)
      alpha = dlin(s) * (1._R8P + (tau / (c * c + 13._R8P / 12._R8P * a * a + &
                                          781._R8P / 720._R8P * t * t + SEAM_WENO6_EPS))**2)
      value_ = value_ + alpha * (wsub(1,s) * u(s) + wsub(2,s) * u(s+1) + wsub(3,s) * u(s+2) + wsub(4,s) * u(s+3))
      asum = asum + alpha
   enddo
   value_ = value_ / asum
   endfunction seam_weno6_core

   pure function seam_restriction_weno6_1d(u, code) result(value_)
   !< Return the WENO restriction of 6 fine values `u` to the coarse centre between nodes code+1 and code+2 (#68 C3):
   !< `seam_weno6_core` on the sub-stencils of `SEAM_D_WENO6(:,code)` with their midpoint cubics, the troubled-cell
   !< fallback the mean of the fine pair. `SEAM_RESTRICT_IDENTITY` (an unrefined direction) returns u(1).
   real(R8P),    intent(in) :: u(1:6)        !< Fine values (footprint nodes).
   integer(I4P), intent(in) :: code          !< Restriction code (0..4, or SEAM_RESTRICT_IDENTITY).
   real(R8P)                :: value_        !< Restricted value.
   real(R8P)                :: wsub(1:4,1:3) !< Cubic weights of the candidates at the coarse centre.
   integer(I4P)             :: s             !< Sub-stencil counter.
   !$acc routine seq
   !$omp declare target

   if (code == SEAM_RESTRICT_IDENTITY) then
      value_ = u(1)
      return
   endif
   wsub = 0._R8P
   do s=1, 3
      if (SEAM_D_WENO6(s,code) > 0._R8P) wsub(:,s) = SEAM_W_CUBIC_MID(:,code-s+2)
   enddo
   value_ = seam_weno6_core(u=u, dlin=SEAM_D_WENO6(:,code), wsub=wsub, fallback=0.5_R8P * (u(code+1) + u(code+2)))
   endfunction seam_restriction_weno6_1d

   pure function seam_fill_weno6_1d(u, sub, p) result(value_)
   !< Return the WENO coarse->fine fill of 6 coarse values `u` at the fine ghost's centre, anchor node `p` + eta(sub)
   !< (issue #68 D2): `seam_weno6_core` on the sub-stencils that contain the target's interval (none extrapolates but
   !< at the block edge, where the one candidate is the tricubic's own one-sided cubic), linear weights
   !< `SEAM_D_FILL_WENO6` reproducing the Lagrange interpolant of their union (order 6 at the centred anchor
   !< positions, 5 one cell from the block edge, 4 at it); the troubled-cell fallback is the linear interpolation on
   !< the target's interval clamped to its two values. An unrefined direction (sub = 3) returns u(p).
   real(R8P),    intent(in) :: u(1:6)        !< Coarse values (footprint nodes).
   integer(I4P), intent(in) :: sub           !< Octant sub-position (1 => eta=-1/4, 2 => eta=+1/4, 3 => unrefined).
   integer(I4P), intent(in) :: p             !< Anchor position (1..6).
   real(R8P)                :: value_        !< Fine-ghost value.
   real(R8P)                :: wsub(1:4,1:3) !< Cubic weights of the candidates at the target.
   real(R8P)                :: fallback      !< Clamped linear interpolation on the target's interval.
   integer(I4P)             :: t4            !< Target position in quarter units of the node spacing (3..25).
   integer(I4P)             :: j             !< Target index (1..12).
   integer(I4P)             :: m             !< Interval holding the target (1..5, the block edge clamped).
   integer(I4P)             :: s             !< Sub-stencil counter.
   !$acc routine seq
   !$omp declare target

   if (sub == SEAM_SUB_COINCIDENT) then
      value_ = u(p)
      return
   endif
   t4 = 4_I4P * p + 2_I4P * sub - 3_I4P
   j = (t4 - 1_I4P) / 2_I4P
   m = min(max(t4 / 4_I4P, 1_I4P), 5_I4P)
   wsub = 0._R8P
   do s=1, 3
      if (SEAM_D_FILL_WENO6(s,j) > 0._R8P) wsub(:,s) = SEAM_W_CUBIC_Q(:,(t4 - 4_I4P * s + 3_I4P) / 2_I4P)
   enddo
   fallback = u(m) + 0.25_R8P * real(t4 - 4_I4P * m, R8P) * (u(m+1) - u(m))
   fallback = min(max(fallback, min(u(m), u(m+1))), max(u(m), u(m+1)))
   value_ = seam_weno6_core(u=u, dlin=SEAM_D_FILL_WENO6(:,j), wsub=wsub, fallback=fallback)
   endfunction seam_fill_weno6_1d
endmodule seth_adam_seam_reference
