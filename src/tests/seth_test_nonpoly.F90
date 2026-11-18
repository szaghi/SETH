!< SETH test: the pieces of the non-polynomial reconstructions (issue #6), against their frozen specifications.
program seth_test_nonpoly
!< SETH test: the pieces of the non-polynomial reconstructions (issue #6), against their frozen specifications.
!<
!< Step 1.1, the mask kernel seth_teno5_mask (frozen before the kernel was written, issue #6, section 4):
!< 1. The mask equals the cutoff outcome inside seth_teno5 on 10^6 random stencils: the test twin teno5_flags (the
!<    library's body with SETH_TENO_FLAGS) returns seth_teno5's states bitwise and the masks of its two states; the
!<    mask kernel on vp(-2:2), and on the mirrored vm(1-m), m = -2..2, returns them.
!< 2. Liang, Xiao, Shyy & Fu 2022, Tables 1-3: three profiles on cells 0..10 (two of them also mirrored in value),
!<    C_T = 1e-5, constant continuation beyond the ends, the flags of cells 2..8 asserted digit for digit. eps is
!<    NP_FLAGS_EPS: the papers' 1e-40 in R8P; 1e-12 in R4P, whose indicators of flat data are round-off (~4e-15 at
!<    f = 1) and turn all-ones data into the mask 5 with SETH_TENO_EPS = 1e-30 (a finding, notes/m5-progress.md).
!< 3. Mirroring the data, v(m) -> v(-m), swaps bits 1 and 2 and leaves bit 0 (10^5 random stencils, exact).
!< Step 1.2, the detectors seth_dr_detect and seth_thinc_detect (issue #6, section 4, with the equation numbers of
!< notes/m5-sources.md: Liang et al. 2022, Eqs. 17, 24, 32; Takagi et al. 2022, Eqs. 19, 20):
!< 5. Exhaustive: every combination of the masks the detectors read (8^6 for DR, 8^7 for THINC), both states, against
!<    a reference written here from the printed conditions as set membership (zeta from a table of the four cases).
!< 6. Tables 1-3 (masks of seth_teno5_mask, constant continuation): profile 1 -> cell 5 at DR level 2, profile 2 ->
!<    level 1, profile 3 -> level 2, the mirrored profiles the same, every other cell of 2..8 at level 0; the THINC
!<    detector flags cell 5 of profile 1 and nothing of profiles 2 and 3 (Liang et al. 2022, p. 7).
!< Steps 2.1 and 2.2, the overrides seth_thinc_override and seth_dr_override (frozen in issue #6, section 5; the
!< THINC guard added from Takagi's Eq. 17, notes/m5-sources.md):
!< 7. Level 0: both states bitwise unchanged, for both overrides.
!< 8. N3 on 10^6 random strictly monotone triples, both overrides, every level: the state lies in [f_{i-1}, f_{i+1}]
!<    (the papers') and in [f_i, f_{i+1}] (derived), to round-off; the right state of the mirrored data equals the
!<    left state bitwise.
!< 9. THINC: the cell average of the step function through the face value is f_i to NP_TOL (1e-13 in R8P) of the
!<    data's magnitude (the conservation constraint), by quadrature as frozen (ten-point Gauss-Legendre on each quarter
!<    of the cell, exact to round-off for tanh of slope beta = 1.6); the closed-form integral of tanh (log cosh) is printed beside it.
!< 10. DR: the state equals m + h tanh(A + atanh B) to NP_TOL; on linear data f_i +- d it is f_i + d tanh(A) to NP_TOL;
!<    the denominator 1 + tanh(A) B of Eq. 31 is positive on every monotone triple (spec 4); no NaN or Inf on 10^6
!<    random triples of any sign; the flat and the non-monotone triples give f_i (the guard).
!< 11. Scaling: f -> sigma f with rmu -> rmu/sigma, sigma = 2^k, bitwise sigma times the state; f -> f + mu moves the
!<    state by mu to NP_OFFSET_TOL (1e-12 in R8P, frozen) of |mu|, |mu| in [1, 10] (below 1 the rounding of f + mu,
!<    eps |f|, is not a fraction of |mu|: the spec's bound is relative to the offset).
!< Step 2.3, the self-contained kernels seth_teno5_dr and seth_teno5_thinc (issue #6, section 5):
!< 12. N1: on periodic scalar fields (smooth sines, a square wave, Liang et al.'s static profiles (c)-(e), random
!<    piecewise cubics with jumps, noise) the kernels equal a two-pass reference written here, bitwise: every cell's
!<    mask computed once and stored (left-biased, and on the mirrored data for the right state), the stored masks read
!<    by the detector, TENO5's states overridden. The reference shares only seth_teno5, seth_teno5_mask and the
!<    detectors and overrides with the kernels (no early exit, no neighbour recomputation).
!< 13. N2: on the smooth fields no face is detected and the states equal seth_teno5's bitwise.
!<     Specs 12-13 also for seth_teno5_rbf_dr (step 3.4; s1 = s2 = 3 dx): its reference is seth_teno5_rbf's states, the
!<     stored masks, the DR detector and override; on the smooth fields it equals seth_teno5_rbf bitwise.
!< Steps 3.1 and 3.2, seth_teno5_st and seth_teno5_rbf (frozen in issue #6, section 6):
!< 14. ST: order 5.0 +- 0.1 at a point (the generic g of seth_test_functions, cell averages, face at 0.3), both states;
!<     equal to seth_linear5 to 64 eps of the data on smooth data; beside a jump that leaves the mask (1,1,0) the order
!<     of the smooth side's error is 4.0 +- 0.2 where TENO5 gives 3.0 +- 0.2; across a unit jump at every cell boundary
!<     no state outside [0, 1] by more than 1e-12 (R8P; TOL scaled in R4P).
!< 15. RBF: s1 = s2 = 0 gives seth_teno5_st (bitwise with gfortran and ifx; nvfortran -fast reassociates the row
!<     sum, so 4 eps of the data is asserted, a deviation from the frozen "bitwise" recorded in notes/m5-progress.md);
!<     with alpha_1 = alpha_2 = 3 fixed and dx -> 0 the order is
!<     5.0 +- 0.1 (s = 3 dx .. 6 dx shrinks with the grid: the paper's degeneration to the polynomial scheme).
!< 16. RBF spec 3 (G1, G2): the generated rows: the polynomial rows of the masks 1, 2, 4 are TENO5's candidates, of 7
!<     linear5's, of 6 TENO5's renormalised weights, exactly; every s^2 row sums to zero; the leading error of each of
!<     the paper's six rows (masks 7, 3, 5, 1, 2, 4; the row of 6 is D4's, not the paper's) with lambda = 3 fixed
!<     equals Liang & Fu's Eq. 29 within 1 %, e.g. the five-point row -dx^5 (lambda^2 g'''/2 + g^(5)/60) (Richardson on
!<     dx, the generic g; R4P: rows whose error is near round-off are reported, not checked). And the kernel equals a reference written here (the sensor of Eqs. 30-31,
!<     the row of the state's mask, P + s^2 C) to 64 eps on 10^5 random stencils with s1 = s2 = 0.3 (the sensor active).
!< 17. E7 (step 3.3): how often each mask occurs, by data: smooth cubics, a jump of random place and height on a
!<     linear function, a kink, independent values; 10^5 stencils each (10^7 in report mode, nonpoly-census.dat). The
!<     mask (0,1,1), which no row of Liang & Fu covers (D4), is reported, not asserted.
!< 18. G6, E6 (step 3.2 spec 5, step 3.3): the linear RBF scheme, the five-point row P + s^2 C with s fixed: the
!<     imaginary part of its modified wavenumber, Im theta~ = -Re sum_m c_m (e^{i theta m} - e^{i theta (m-1)}), over
!<     theta in (0, pi], at s^2 = 9e-4, 3.6e-3, 2.25e-2, 9e-2 (frozen in issue #6). Predicted from Eq. 29 (issue #6,
!<     section 1.3): Im theta~ = theta^4 (s^2/2 - theta^2/60) at leading order, a band theta^2 < 30 s^2 of growth with a
!<     peak of (200/3) s^6 at theta^2 = 20 s^2. Asserted at the two smaller s (where the leading order holds): the band's
!<     edge within 10 % of sqrt(30) s, the peak within a factor 2 of the prediction (the bound of D7). Reported: the four
!<     bands and peaks, and the largest s whose peak stays at or below 2e-4, TENO6's (D7's validated range of s).
!< 19. Step 1.3, Liang et al. 2022, section 3.1: the detectors on point values at the cell centres of [0, 2], constant
!<     continuation beyond the ends, the left state of each face (the cell containing the discontinuity, positive
!<     speed). (a) sin(10 pi x), (b) sin(24 pi x), N = 120: no activation of either detector (the paper's statement).
!<     (c), (d), N = 100 and 400: every jump of the function (above 0.1) has a DR activation within 1.5 cells, and every
!<     DR activation lies within 1.5 cells of a jump of the function or of its derivative. (e), N = 200 and 400: counts
!<     only (its frozen spec contradicts the paper's own Fig. 5: left to the owner). Report: nonpoly-static.dat.
!< 4. Constant data and a jump of BIG_JUMP: the mask is defined and never 0 (one chi_k >= 1/3 > C_T). Not 7 on
!<    constant data in general (printed): the indicators of a non-dyadic constant are round-off, far above eps.
use seth
use seth_test_parameters, only : BIG_JUMP, NP_FLAGS_EPS, NP_TOL, NP_OFFSET_TOL, ORDERS_AT_FLOOR_OK, POINT_FLOOR, POINT_ORDER_LOSS
use seth_test_rng, only : seed_rng
use seth_test_teno_flags, only : teno5_flags
use seth_vv_cases, only : vv_report, vv_open
use seth_test_functions, only : gfun, gavg, GL_X, GL_W
use seth_nonpoly_tables, only : NP_RBF_P, NP_RBF_C
use seth_z_teno_tables, only : TENO5_P, LIN5_C
implicit none
integer(I4P), parameter :: NR = 1000000_I4P, NM = 100000_I4P
real(R_P)               :: vp(-2:3), vm(-2:3), w(-2:2), rp, rm, sp, sm, y
integer(I4P)            :: n_, m, mp, mm, k, mk, mk2, nfail, count_(0:7), j
logical                 :: passed
integer                 :: uf, um

passed = .true.
call seed_rng(1101_I4P)
if (vv_report()) then
   call vv_open('nonpoly-flags.dat', uf, 'profile cell f delta0 delta1 delta2 printed0 printed1 printed2 (cells 0..10)')
   call vv_open('nonpoly-masks.dat', um, 'population mask0..mask7 (counts)')
endif

! 1. the mask is the kernel's cutoff
nfail = 0 ; count_ = 0
do n_=1, NR
   call random_stencil(vp, vm)
   call seth_teno5(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)
   call teno5_flags(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, sp, sm, mp, mm)
   if (rp /= sp .or. rm /= sm) then
      nfail = nfail + 1 ; if (nfail <= 3) print '(A,2ES24.16)', 'twin states differ from seth_teno5: ', rp - sp, rm - sm
   endif
   call seth_teno5_mask(vp(-2:2), SETH_TENO_EPS, SETH_TENO5_CT, mk)
   do m=-2, 2
      w(m) = vm(1-m)
   enddo
   call seth_teno5_mask(w, SETH_TENO_EPS, SETH_TENO5_CT, mk2)
   if (mk /= mp .or. mk2 /= mm) then
      nfail = nfail + 1 ; if (nfail <= 3) print '(A,4I3)', 'mask differs from the kernel''s cutoff: ', mk, mp, mk2, mm
   endif
   count_(mp) = count_(mp) + 1
enddo
print '(A,I0,A,8I8)', 'masks equal to the cutoff of seth_teno5 on ', NR, ' stencils; masks 0..7 seen: ', count_
if (vv_report()) write(um, '(A,8(1X,I0))') 'random', count_
if (nfail > 0) then ; print '(A,I0)', 'FAIL 1: ', nfail ; passed = .false. ; endif

! 2. Liang et al. 2022, Tables 1-3
call table('profile1', 'profile 1 (Table 1)', &
           [0._R_P, 0._R_P, 0._R_P, 0._R_P, 0._R_P, 0.5_R_P, 1._R_P, 1._R_P, 1._R_P, 1._R_P, 1._R_P], &
           [1,1,1,1,0,1,0,1,1,1,1], [1,1,1,0,0,1,1,1,1,1,1], [1,1,1,1,1,1,0,0,1,1,1])
call table('profile2', 'profile 2 (Table 2)', &
           [0._R_P, 0._R_P, 0._R_P, 0._R_P, 0._R_P, 0._R_P, 1._R_P, 1._R_P, 1._R_P, 1._R_P, 1._R_P], &
           [1,1,1,1,1,0,0,1,1,1,1], [1,1,1,1,0,0,1,1,1,1,1], [1,1,1,1,1,1,0,0,1,1,1])
call table('profile2m', 'profile 2 mirrored in value', &
           [1._R_P, 1._R_P, 1._R_P, 1._R_P, 1._R_P, 1._R_P, 0._R_P, 0._R_P, 0._R_P, 0._R_P, 0._R_P], &
           [1,1,1,1,1,0,0,1,1,1,1], [1,1,1,1,0,0,1,1,1,1,1], [1,1,1,1,1,1,0,0,1,1,1])
call table('profile3', 'profile 3 (Table 3)', &
           [0._R_P, 0._R_P, 0._R_P, 0._R_P, 0._R_P, 0.1_R_P, 1._R_P, 1._R_P, 1._R_P, 1._R_P, 1._R_P], &
           [1,1,1,1,0,0,0,1,1,1,1], [1,1,1,0,0,0,1,1,1,1,1], [1,1,1,1,1,1,0,0,1,1,1])
call table('profile3m', 'profile 3 mirrored in value', &
           [1._R_P, 1._R_P, 1._R_P, 1._R_P, 1._R_P, 0.9_R_P, 0._R_P, 0._R_P, 0._R_P, 0._R_P, 0._R_P], &
           [1,1,1,1,0,0,0,1,1,1,1], [1,1,1,0,0,0,1,1,1,1,1], [1,1,1,1,1,1,0,0,1,1,1])

! 3. mirror symmetry
nfail = 0
do n_=1, NM
   call random_stencil(vp, vm)
   call seth_teno5_mask(vp(-2:2), SETH_TENO_EPS, SETH_TENO5_CT, mk)
   do m=-2, 2
      w(m) = vp(-m)
   enddo
   call seth_teno5_mask(w, SETH_TENO_EPS, SETH_TENO5_CT, mk2)
   if (mk2 /= iand(mk, 1) + 2*ibits(mk, 2, 1) + 4*ibits(mk, 1, 1)) then
      nfail = nfail + 1 ; if (nfail <= 3) print '(A,2I3)', 'mirror: ', mk, mk2
   endif
enddo
print '(A,I0,A)', 'mirroring swaps bits 1 and 2 on ', NM, ' stencils'
if (nfail > 0) then ; print '(A,I0)', 'FAIL 3: ', nfail ; passed = .false. ; endif

! 4. constant data and a big jump: defined, never 0. Constant data does not give 7 in general: the indicators of a
! non-dyadic constant are round-off, of the order of (eps_mach c)^2, far above eps = 1e-40, so tau/beta_k is O(1) and
! the cutoff drops candidates (the masks seen are printed: a finding, issue #6 step 1.1)
count_ = 0
do k=1, 1000
   call random_number(y) ; y = 20._R_P*(y - 0.5_R_P)
   w = y
   call seth_teno5_mask(w, SETH_TENO_EPS, SETH_TENO5_CT, mk)
   count_(mk) = count_(mk) + 1
   if (mk < 1 .or. mk > 7) then ; print '(A,ES24.16,I3)', 'FAIL 4: constant data, mask ', y, mk ; passed = .false. ; endif
enddo
print '(A,8I5)', 'constant data, 1000 random values, masks 0..7 seen:', count_
if (vv_report()) write(um, '(A,8(1X,I0))') 'constant', count_
do j=-1, 2 ! a jump of BIG_JUMP between cells j-1 and j, both orientations
   do k=0, 1
      do m=-2, 2
         w(m) = merge(BIG_JUMP, 0._R_P, (m >= j) .neqv. (k == 1))
      enddo
      call seth_teno5_mask(w, SETH_TENO_EPS, SETH_TENO5_CT, mk)
      if (mk < 1 .or. mk > 7) then ; print '(A,2I3,I4)', 'FAIL 4: big jump, mask ', j, k, mk ; passed = .false. ; endif
   enddo
enddo
print '(A)', 'constant data and a big jump: mask in 1..7'

! 12.-13. the self-contained kernels against the two-pass reference
call self_contained()

! 14.-15. stencil targeting and the RBF fluxes
call targeting()

! 17. E7, the census of the masks
call census()

! 18. G6 and E6, the long-wave growth of the linear RBF scheme
call rbf_adr()

! 19. the published static tests
call static_tests()

! 7.-11. the overrides
call overrides()

! 5. the detectors, exhaustively against the reference
call exhaustive()

! 6. the detector levels on Tables 1-3
call levels('profile 1', [0._R_P, 0._R_P, 0._R_P, 0._R_P, 0._R_P, 0.5_R_P, 1._R_P, 1._R_P, 1._R_P, 1._R_P, 1._R_P], 2, 1)
call levels('profile 2', [0._R_P, 0._R_P, 0._R_P, 0._R_P, 0._R_P, 0._R_P, 1._R_P, 1._R_P, 1._R_P, 1._R_P, 1._R_P], 1, 0)
call levels('profile 2 mirrored', [1._R_P, 1._R_P, 1._R_P, 1._R_P, 1._R_P, 1._R_P, 0._R_P, 0._R_P, 0._R_P, 0._R_P, &
            0._R_P], 1, 0)
call levels('profile 3', [0._R_P, 0._R_P, 0._R_P, 0._R_P, 0._R_P, 0.1_R_P, 1._R_P, 1._R_P, 1._R_P, 1._R_P, 1._R_P], 2, 0)
call levels('profile 3 mirrored', [1._R_P, 1._R_P, 1._R_P, 1._R_P, 1._R_P, 0.9_R_P, 0._R_P, 0._R_P, 0._R_P, 0._R_P, &
            0._R_P], 2, 0)

if (vv_report()) then ; close(uf) ; close(um) ; endif
if (.not.passed) error stop 'seth_test_nonpoly failed'
print '(A)', 'seth_test_nonpoly: ok'

contains
   pure function zeta_ref(mk) result(z)
   !< Eq. 17 as a table of (delta_1, delta_2): (0,0) -> 2, (0,1) -> 1, (1,0) -> -1, (1,1) -> 0.
   integer(I4P), intent(in) :: mk
   integer(I4P)             :: z
   integer(I4P), parameter  :: TAB(0:1,0:1) = reshape([2, -1, 1, 0], [2, 2]) ! TAB(delta_1, delta_2)
   z = TAB(ibits(mk, 1, 1), ibits(mk, 2, 1))
   endfunction zeta_ref

   pure function dr_ref(mk) result(l)
   !< Eqs. 24 and 32, the zeta of cells -3..2 by set membership.
   integer(I4P), intent(in) :: mk(-3:3)
   integer(I4P)             :: l, z(-3:2), c
   logical                  :: e24, e32
   do c=-3, 2
      z(c) = zeta_ref(mk(c))
   enddo
   e24 = (any(z(-3:-1) == 1 .or. z(-3:-1) == 2) .and. z(0) == 0 .and. any(z(1:2) == -1 .or. z(1:2) == 2)) .or. &
         (any(z(-2:-1) == 1 .or. z(-2:-1) == 2) .and. z(0) == 1 .and. z(1) == -1 .and. (z(2) == -1 .or. z(2) == 2))
   e32 = (all(z(-2:-1) == 1 .or. z(-2:-1) == 2) .and. z(0) == 0 .and. all(z(1:2) == -1 .or. z(1:2) == 2)) .or. &
         (all(z(-2:-1) == 1 .or. z(-2:-1) == 2) .and. z(0) == 1 .and. z(1) == -1 .and. (z(2) == -1 .or. z(2) == 2))
   l = 0
   if (e24) l = 1
   if (e24 .and. e32) l = 2
   endfunction dr_ref

   pure function thinc_ref(mk) result(l)
   !< Takagi's Eq. 20, the zeta of cells -3..3 by set membership.
   integer(I4P), intent(in) :: mk(-3:3)
   integer(I4P)             :: l, z(-3:3), c
   do c=-3, 3
      z(c) = zeta_ref(mk(c))
   enddo
   l = merge(1, 0, any(z(-3:-1) == 1 .or. z(-3:-1) == 2) .and. z(0) == 0 .and. any(z(1:3) == -1 .or. z(1:3) == 2))
   endfunction thinc_ref

   subroutine rbf_adr()
   !< Spec 18.
   real(R_P), parameter    :: S2(4) = [9.e-4_R_P, 3.6e-3_R_P, 2.25e-2_R_P, 9.e-2_R_P]
   integer(I4P), parameter :: NT = 20000_I4P
   real(R_P)               :: th, g, gmax, thmax, edge, pk, pe, s, smax
   integer(I4P)            :: k, j, ua
   logical                 :: ok
   ok = .true.
   if (vv_report()) call vv_open('nonpoly-adr.dat', ua, 's^2 band-edge predicted-edge peak predicted-peak theta-of-peak')
   print '(A)', 'RBF linear five-point row: the band of growth (Im theta~ > 0) and its peak'
   do k=1, 4
      gmax = 0._R_P ; thmax = 0._R_P ; edge = 0._R_P
      do j=1, NT
         th = acos(-1._R_P)*j/NT
         g = growth(th, S2(k))
         if (g > gmax) then ; gmax = g ; thmax = th ; endif
         if (g > 0._R_P) edge = th
      enddo
      pe = sqrt(30._R_P*S2(k)) ; pk = 200._R_P/3._R_P*S2(k)**3
      print '(A,ES9.2,A,F8.4,A,F8.4,A,ES10.3,A,ES10.3)', '  s^2 ', S2(k), ': edge ', edge, ' (pred. ', pe, &
            '), peak ', gmax, ' (pred. ', pk
      if (vv_report()) write(ua, '(6ES16.7)') S2(k), edge, pe, gmax, pk, thmax
      if (k <= 2) then
         if (abs(edge/pe - 1._R_P) > 0.1_R_P .or. gmax > 2._R_P*pk .or. gmax < 0.5_R_P*pk) ok = .false.
      endif
   enddo
   ! D7: the largest s whose peak stays at or below 2e-4
   smax = 0._R_P
   do k=1, 600
      s = 0.001_R_P*k
      gmax = 0._R_P
      do j=1, 2000
         gmax = max(gmax, growth(acos(-1._R_P)*j/2000._R_P, s*s))
      enddo
      if (gmax <= 2.e-4_R_P) smax = s
      if (gmax > 2.e-4_R_P) exit
   enddo
   print '(A,F6.3,A)', '  largest s with a peak <= 2e-4 (TENO6''s): ', smax, ' (D7: the validated range)'
   if (vv_report()) then
      write(ua, '(A,F8.4)') '# validated-s ', smax
      close(ua)
   endif
   if (.not.ok) then ; print '(A)', 'FAIL 18: band or peak off the prediction' ; passed = .false. ; endif
   endsubroutine rbf_adr

   pure function growth(th, s2) result(g)
   !< Im theta~ of the linear five-point RBF row with s^2 = s2.
   real(R_P), intent(in) :: th, s2
   real(R_P)             :: g
   integer(I4P)          :: m
   g = 0._R_P
   do m=-2, 2
      g = g - (NP_RBF_P(m,7) + s2*NP_RBF_C(m,7))*(cos(th*m) - cos(th*(m-1)))
   enddo
   endfunction growth

   subroutine static_tests()
   !< Spec 19.
   integer              :: us, up, j
   logical              :: ok
   ok = .true.
   if (vv_report()) call vv_open('nonpoly-static.dat', us, 'case N detector x-of-each-activation...')
   if (vv_report()) then ! the three functions, for the figure: 2000 points on [0, 2]
      call vv_open('nonpoly-static-profile.dat', up, 'x f_c f_d f_e (Liang et al. 2022, static tests (c), (d), (e))')
      do j=0, 1999
         write(up, '(4ES16.8)') 2._R_P*(j + 0.5_R_P)/2000._R_P, static_c(2._R_P*(j + 0.5_R_P)/2000._R_P), &
                                static_d(2._R_P*(j + 0.5_R_P)/2000._R_P), static_e(2._R_P*(j + 0.5_R_P)/2000._R_P)
      enddo
      close(up)
   endif
   call static_one('a', 120, 0, us, ok)
   call static_one('b', 120, 0, us, ok)
   call static_one('c', 100, 1, us, ok) ; call static_one('c', 400, 1, us, ok)
   call static_one('d', 100, 2, us, ok) ; call static_one('d', 400, 2, us, ok)
   call static_one('e', 200, 3, us, ok) ; call static_one('e', 400, 3, us, ok)
   if (vv_report()) close(us)
   if (.not.ok) then ; print '(A)', 'FAIL 19: the static tests' ; passed = .false. ; endif
   endsubroutine static_tests

   subroutine static_one(case_, n, kind_, us, ok)
   character(*), intent(in)    :: case_
   integer(I4P), intent(in)    :: n, kind_
   integer,      intent(in)    :: us
   real(R_P), parameter :: JC(2) = [0.6_R_P, 0.8_R_P], KC(8) = [0.2_R_P, 0.4_R_P, 0.6_R_P, 0.8_R_P, 1._R_P, 1.1_R_P, &
                                                            1.2_R_P, 1.4_R_P]
   real(R_P), parameter :: JD(9) = [0.4_R_P, 0.5_R_P, 0.6_R_P, 0.9_R_P, 1._R_P, 1.1_R_P, 1.2_R_P, 1.5_R_P, 1.8_R_P]
   real(R_P), parameter :: JE(7) = [0.25_R_P, 0.27_R_P, 0.29_R_P, 0.6_R_P, 0.63_R_P, 0.66_R_P, 1.14_R_P]
   logical,      intent(inout) :: ok
   real(R_P)                   :: f(-3:n+2), w(-2:2), x, h, xs(n)
   integer(I4P)                :: mk(-3:n+2), i, m, l, lm, ndr, nth, j
   real(R_P)                   :: xdr(n), xth(n)
   h = 2._R_P/n
   do i=0, n-1
      x = (i + 0.5_R_P)*h ; xs(i+1) = x
      select case(case_)
      case('a') ; f(i) = sin(10._R_P*acos(-1._R_P)*x)
      case('b') ; f(i) = sin(24._R_P*acos(-1._R_P)*x)
      case('c') ; f(i) = static_c(x)
      case('d') ; f(i) = static_d(x)
      case default ; f(i) = static_e(x)
      endselect
   enddo
   f(-3:-1) = f(0) ; f(n:n+2) = f(n-1)
   do i=-1, n
      do m=-2, 2
         w(m) = f(min(max(i+m, 0), n-1))
      enddo
      call seth_teno5_mask(w, SETH_TENO_EPS, SETH_TENO5_CT, mk(i))
   enddo
   mk(-3:-2) = mk(-1) ; mk(n+1:n+2) = mk(n)
   ndr = 0 ; nth = 0
   do i=0, n-1
      call seth_dr_detect(mk(i-3:i+3), mk(i-3:i+3), l, lm)
      if (l > 0) then ; ndr = ndr + 1 ; xdr(ndr) = xs(i+1) ; endif
      call seth_thinc_detect(mk(i-3:i+3), mk(i-3:i+3), l, lm)
      if (l > 0) then ; nth = nth + 1 ; xth(nth) = xs(i+1) ; endif
   enddo
   print '(A,A,A,I4,A,I3,A,I3)', '  static (', case_, '), N = ', n, ': DR activations ', ndr, ', THINC ', nth
   if (vv_report()) then
      write(us, '(A,1X,I0,1X,A,*(1X,F8.4))') case_, n, 'dr', xdr(1:ndr)
      write(us, '(A,1X,I0,1X,A,*(1X,F8.4))') case_, n, 'thinc', xth(1:nth)
   endif
   select case(kind_)
   case(0)
      if (ndr > 0 .or. nth > 0) then ; print '(A)', '  FAIL: activations on a sine' ; ok = .false. ; endif
   case(1)
      call complete_local(JC, KC, ndr, xdr, h, ok)
   case(2)
      call complete_local(JD, JD, ndr, xdr, h, ok)
   endselect
   endsubroutine static_one

   subroutine complete_local(jumps, marks, nf, xf, h, ok)
   !< Every jump has an activation within 1.5 cells; every activation lies within 1.5 cells of a mark (a jump of the
   !< function or of its derivative).
   real(R_P),    intent(in)    :: jumps(:), marks(:), xf(:), h
   integer(I4P), intent(in)    :: nf
   logical,      intent(inout) :: ok
   integer(I4P)                :: j
   do j=1, size(jumps)
      if (.not.any(abs(xf(1:nf) - jumps(j)) <= 1.5_R_P*h)) then
         print '(A,F6.3)', '  FAIL: jump not flagged at x = ', jumps(j) ; ok = .false.
      endif
   enddo
   do j=1, nf
      if (.not.any(abs(marks - xf(j)) <= 1.5_R_P*h)) then
         print '(A,F6.3)', '  FAIL: activation away from every jump at x = ', xf(j) ; ok = .false.
      endif
   enddo
   endsubroutine complete_local

   subroutine census()
   !< Spec 17 (E7).
   character(8), parameter :: DNAME(4) = ['smooth  ', 'jump    ', 'kink    ', 'noise   ']
   real(R_P)                :: v(-2:2), cf(0:3), y, hgt
   integer(I4P)             :: d, n_, m, mk, jp, cnt(0:7), uc
   integer(I8P)             :: nn
   nn = merge(10000000_I8P, 100000_I8P, vv_report())
   if (vv_report()) call vv_open('nonpoly-census.dat', uc, 'data stencils count-of-mask-0..7')
   print '(A,I0,A)', 'E7: masks by data, ', nn, ' stencils each (mask = d0 + 2 d1 + 4 d2)'
   do d=1, 4
      cnt = 0
      do n_=1, int(nn, I4P)
         call random_number(cf) ; cf = cf - 0.5_R_P
         call random_number(y) ; jp = -1 + int(4._R_P*y)
         call random_number(hgt) ; hgt = 10._R_P**(8._R_P*hgt - 6._R_P)
         do m=-2, 2
            select case(d)
            case(1) ; v(m) = 1._R_P + cf(1)*m + cf(2)*m*m*0.1_R_P + cf(3)*m**3*0.01_R_P
            case(2) ; v(m) = cf(0) + cf(1)*m*0.1_R_P + merge(hgt, 0._R_P, m >= jp)
            case(3) ; v(m) = cf(0) + cf(1)*abs(m - jp + 0.5_R_P)
            case default ; call random_number(v(m))
            endselect
         enddo
         call seth_teno5_mask(v, SETH_TENO_EPS, SETH_TENO5_CT, mk)
         cnt(mk) = cnt(mk) + 1
      enddo
      print '(2X,A8,8I10)', DNAME(d), cnt
      if (vv_report()) write(uc, '(A,1X,I0,8(1X,I0))') trim(DNAME(d)), nn, cnt
   enddo
   if (vv_report()) close(uc)
   endsubroutine census

   subroutine targeting()
   !< Specs 14-15.
   integer(I4P), parameter :: NH = 6_I4P
   real(R_P), parameter    :: XF = 0.3_R_P
   real(R_P)               :: vp(-2:3), vm(-2:3), h, rp, rm, sp, sm, ex, e(4,NH), o(4), worst, lo, hi, x
   integer(I4P)            :: j, m, k, jump, n_, us
   logical                 :: ok
   ok = .true.
   if (vv_report()) call vv_open('nonpoly-st.dat', us, 'quantity value')
   ! order on smooth data: ST (both states), RBF with s = 3 dx .. 6 dx; ST against linear5
   worst = 0._R_P
   do j=1, NH
      h = 0.2_R_P/2._R_P**(j-1)
      do m=-2, 3
         vp(m) = gavg(XF + (m-1)*h, XF + m*h) ; vm(m) = gavg(XF + (m-1)*h, XF + m*h)
      enddo
      ex = gfun(XF)
      call seth_teno5_st(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)
      call seth_linear5(vp, vm, sp, sm)
      worst = max(worst, abs(rp - sp)/maxval(abs(vp)), abs(rm - sm)/maxval(abs(vp)))
      e(1,j) = abs(rp - ex) ; e(2,j) = abs(rm - ex)
      call seth_teno5_rbf(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, SETH_RBF_CR, SETH_RBF_ALPHA1*h, SETH_RBF_ALPHA2*h, &
                          rp, rm)
      e(3,j) = abs(rp - ex) ; e(4,j) = abs(rm - ex)
   enddo
   do k=1, 4
      o(k) = pair_order(e(k,:))
   enddo
   print '(A,4F7.3,A,ES9.2)', 'ST order (left, right), RBF order (left, right): ', o, '; ST - linear5: ', worst
   if (vv_report()) write(us, '(A,4(1X,F8.4))') 'orders-smooth', o
   if (vv_report()) write(us, '(A,1X,ES12.4)') 'st-minus-linear5', worst
   if (any(abs(o - 5._R_P) > 0.1_R_P + POINT_ORDER_LOSS .and. .not.(ORDERS_AT_FLOOR_OK .and. o < 0._R_P))) then
      print '(A)', 'FAIL 14-15: order' ; ok = .false.
   endif
   if (worst > 64._R_P*epsilon(1._R_P)) then ; print '(A)', 'FAIL 14: ST differs from linear5' ; ok = .false. ; endif
   ! beside a jump between the cells -2 and -1 (mask (1,1,0)): the error of the smooth side, ST against TENO5
   do j=1, NH
      h = 0.1_R_P/2._R_P**(j-1)
      do m=-2, 3
         vp(m) = gavg(XF + (m-1)*h, XF + m*h) + merge(1._R_P, 0._R_P, m <= -2)
      enddo
      ex = gfun(XF)
      call seth_teno5_st(vp, vp, SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)
      call seth_teno5(vp, vp, SETH_TENO_EPS, SETH_TENO5_CT, sp, sm)
      e(1,j) = abs(rp - ex) ; e(2,j) = abs(sp - ex)
   enddo
   o(1) = pair_order(e(1,:)) ; o(2) = pair_order(e(2,:))
   print '(A,2F7.3)', 'beside a jump, mask (1,1,0): order of ST, of TENO5: ', o(1:2)
   if (vv_report()) write(us, '(A,2(1X,F8.4))') 'orders-jump', o(1:2)
   if ((abs(o(1) - 4._R_P) > 0.2_R_P + POINT_ORDER_LOSS .or. abs(o(2) - 3._R_P) > 0.2_R_P + POINT_ORDER_LOSS) .and. &
       .not.(ORDERS_AT_FLOOR_OK .and. min(o(1), o(2)) < 0._R_P)) then
      print '(A)', 'FAIL 14: orders beside the jump' ; ok = .false.
   endif
   ! a unit jump at every cell boundary, both orientations: no overshoot
   worst = 0._R_P
   do jump=-1, 3
      do k=0, 1
         do m=-2, 3
            vp(m) = merge(1._R_P, 0._R_P, (m >= jump) .neqv. (k == 1)) ; vm(m) = vp(m)
         enddo
         call seth_teno5_st(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)
         worst = max(worst, -rp, rp - 1._R_P, -rm, rm - 1._R_P)
      enddo
   enddo
   print '(A,ES9.2)', 'ST across a unit jump: worst overshoot ', worst
   if (vv_report()) write(us, '(A,1X,ES12.4)') 'overshoot', worst
   if (worst > NP_TOL*1.e1_R_P) then ; print '(A)', 'FAIL 14: overshoot' ; ok = .false. ; endif
   ! RBF with s1 = s2 = 0 is ST: bitwise with gfortran and ifx; nvfortran -fast reassociates the row sum (P + 0 C) f,
   ! so the frozen "bitwise" is asserted as 4 eps of the data, the bitwise count printed (deviation recorded in
   ! notes/m5-progress.md)
   k = 0 ; worst = 0._R_P
   do n_=1, 100000
      call random_stencil(vp, vm)
      call seth_teno5_st(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)
      call seth_teno5_rbf(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, SETH_RBF_CR, 0._R_P, 0._R_P, sp, sm)
      if (rp == sp .and. rm == sm) k = k + 1
      worst = max(worst, abs(rp - sp)/max(maxval(abs(vp)), tiny(1._R_P)), abs(rm - sm)/max(maxval(abs(vm)), tiny(1._R_P)))
   enddo
   print '(A,I0,A,ES9.2)', 'RBF at s1 = s2 = 0 against ST on 10^5 stencils: bitwise ', k, ', worst ', worst
   if (worst > 4._R_P*epsilon(1._R_P)) then ; print '(A)', 'FAIL 15: RBF at s = 0 differs from ST' ; ok = .false. ; endif
   call rbf_rows(ok)
   if (vv_report()) close(us)
   if (.not.ok) passed = .false.
   endsubroutine targeting

   subroutine rbf_rows(ok)
   !< Spec 16.
   logical, intent(inout) :: ok
   logical                :: rows_ok
   real(R_P), parameter   :: XF = 0.3_R_P, LAM = 3._R_P
   integer(I4P), parameter :: RMASK(6) = [7_I4P, 3_I4P, 5_I4P, 1_I4P, 2_I4P, 4_I4P] !< the paper's six rows
   integer(I4P), parameter :: RLEAD(6) = [5_I4P, 4_I4P, 4_I4P, 3_I4P, 3_I4P, 3_I4P] !< their leading dx power
   real(R_P),    parameter :: RA2(6) = [-1._R_P/2._R_P, -2._R_P/3._R_P, 1._R_P, 1._R_P, -1._R_P, -3._R_P] !< Eq. 29
   real(R_P),    parameter :: RA0(6) = [-1._R_P/60._R_P, -1._R_P/30._R_P, 1._R_P/20._R_P, 1._R_P/12._R_P, &
                                        -1._R_P/12._R_P, -1._R_P/4._R_P]
   real(R_P)              :: vp(-2:3), vm(-2:3), w(-2:3), rp, rm, sp, sm, h, e(2), c, pred, worst
   integer(I4P)           :: k, m, j, n_, r_, mk, lead
   ! the polynomial rows, exactly
   rows_ok = .true.
   do m=-2, 2
      if (NP_RBF_P(m,7) /= LIN5_C(m)) rows_ok = .false.
   enddo
   do k=0, 2
      do m=-2, 2
         if (m >= TENO5_CLO_(k) .and. m <= TENO5_CLO_(k) + 2) then
            if (NP_RBF_P(m,2**k) /= TENO5_P(m - TENO5_CLO_(k),k)) rows_ok = .false.
         elseif (NP_RBF_P(m,2**k) /= 0._R_P) then
            rows_ok = .false.
         endif
      enddo
   enddo
   do m=-2, 2
      c = 0._R_P
      if (m >= 0) c = c + 0.75_R_P*TENO5_P(m,1)
      if (m <= 0) c = c + 0.25_R_P*TENO5_P(m+2,2)
      if (abs(NP_RBF_P(m,6) - c) > 4._R_P*epsilon(1._R_P)) rows_ok = .false.
   enddo
   do k=0, 7
      if (abs(sum(NP_RBF_C(:,k))) > 16._R_P*epsilon(1._R_P)) rows_ok = .false.
   enddo
   if (.not.rows_ok) then ; print '(A)', 'FAIL 16: a generated row' ; ok = .false. ; endif
   ! Eq. 29: the leading error of each of the paper's six rows, lambda fixed, P + s^2 C applied to cell averages of g
   ! (cell m = [XF + (m-1) h, XF + m h], the face at XF): error / dx^L -> a2 lambda^2 g^(L-2) + a0 g^(L), Richardson
   ! on h, h/2 removing the dx^(L+1) term; the coefficients as printed (scripts/gen_nonpoly_tables.py, (b))
   do r_=1, 6
      mk = RMASK(r_) ; lead = RLEAD(r_)
      pred = RA2(r_)*LAM*LAM*gdn(lead-2, XF) + RA0(r_)*gdn(lead, XF)
      do j=1, 2
         h = 0.025_R_P/2._R_P**(j-1)
         rp = 0._R_P
         do m=-2, 2
            rp = rp + (NP_RBF_P(m,mk) + (LAM*h)**2*NP_RBF_C(m,mk))*gavg(XF + (m-1)*h, XF + m*h)
         enddo
         e(j) = (rp - gfun(XF))/h**lead
      enddo
      c = 2._R_P*e(2) - e(1)
      if (ORDERS_AT_FLOOR_OK .and. 8._R_P*epsilon(1._R_P)/0.0125_R_P**lead > 0.002_R_P*abs(pred)) then
         print '(A,I0,A)', 'RBF row of the mask ', mk, ': its error is near round-off at this precision (not checked)'
      else
         print '(A,I0,A,I0,A,2ES13.5,A,F7.3,A)', 'RBF row of the mask ', mk, ', lambda = 3: error / dx^', lead, ' ', &
               c, pred, ' (', 100._R_P*abs(c/pred - 1._R_P), ' % from Eq. 29)'
         if (abs(c/pred - 1._R_P) > 0.01_R_P) then ; print '(A,I0)', 'FAIL 16: Eq. 29, mask ', mk ; ok = .false. ; endif
      endif
   enddo
   ! the kernel against the reference, sensor active
   worst = 0._R_P
   do n_=1, 100000
      call random_stencil(vp, vm)
      call seth_teno5_rbf(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, SETH_RBF_CR, 0.3_R_P, 0.3_R_P, rp, rm)
      sp = rbf_ref(vp, 0.3_R_P, 0.3_R_P)
      do m=-2, 3
         w(m) = vm(1-m)
      enddo
      sm = rbf_ref(w, 0.3_R_P, 0.3_R_P)
      worst = max(worst, abs(rp - sp)/maxval(abs(vp)), abs(rm - sm)/maxval(abs(vm)))
   enddo
   print '(A,ES9.2)', 'RBF kernel against the reference on 10^5 stencils (s1 = s2 = 0.3): worst ', worst
   if (worst > 64._R_P*epsilon(1._R_P)) then ; print '(A)', 'FAIL 16: RBF kernel differs from its reference' ; ok = .false. ; endif
   endsubroutine rbf_rows

   function rbf_ref(v, s1, s2) result(r)
   !< Liang & Fu 2024, Eqs. 25 and 30-32 for the left-biased state of v(-2:3): the mask of the cell's TENO5 cutoff
   !< selects the row, s = s1 + s2 (1 - g(m)) from the minimum of eta over the nodes -1..1.
   real(R_P), intent(in) :: v(-2:3), s1, s2
   real(R_P)             :: r, eta(-1:1), mm, g, s, u(-2:2)
   integer(I4P)          :: j, mk
   do j=-2, 2
      u(j) = v(j)
   enddo
   call seth_teno5_mask(u, SETH_TENO_EPS, SETH_TENO5_CT, mk)
   do j=-1, 1
      eta(j) = (abs(2._R_P*(v(j) - v(j-1))*(v(j+1) - v(j))) + SETH_TENO_EPS) &
               /((v(j) - v(j-1))**2 + (v(j+1) - v(j))**2 + SETH_TENO_EPS)
   enddo
   mm = 1._R_P - min(1._R_P, minval(eta)/SETH_RBF_CR)
   g = (1._R_P - mm)**4*(1._R_P + 4._R_P*mm)
   s = s1 + s2*(1._R_P - g)
   r = 0._R_P
   do j=-2, 2
      r = r + (NP_RBF_P(j,mk) + s*s*NP_RBF_C(j,mk))*v(j)
   enddo
   endfunction rbf_ref

   subroutine overrides()
   !< Specs 7-11 of the two overrides.
   real(R_P)    :: f(-1:2), g(-1:2), w(-1:2), rp, rm, sp, sm, x(3), lo, hi, worst(7), fs, fd, y, avg, b, sc, mu, ex
   real(R_P)    :: rp0, rm0, a1, avgq
   integer(I4P) :: n_, lv, k, nfail, ex2, q, nsat
   logical      :: thinc
   nfail = 0 ; worst = 0._R_P ; nsat = 0
   do n_=1, NR
      call random_number(x) ; x = 6._R_P*(x - 0.5_R_P)
      call random_number(y)
      ! a strictly monotone triple, increasing or decreasing
      if (x(1) == x(2) .or. x(2) == x(3) .or. x(1) == x(3)) cycle
      call sort3(x)
      if (y < 0.5_R_P) x = x(3:1:-1)
      f(-1:1) = x ; f(2) = x(3) + 1._R_P
      if ((f(1) - f(0))*(f(0) - f(-1)) <= SETH_NONPOLY_GUARD) cycle
      do m=-1, 2
         g(m) = f(1-m) ! the mirrored data: the right state of g is the left state of f
      enddo
      call random_number(y) ; rp0 = y ; rm0 = -y
      do k=0, 1
         thinc = k == 0
         do lv=0, merge(1, 2, thinc)
            rp = rp0 ; rm = rm0 ; sp = rp0 ; sm = rm0
            if (thinc) then
               call seth_thinc_override(f, g, lv, lv, SETH_THINC_BETA, 1._R_P, rp, rm)
            else
               call seth_dr_override(f, g, lv, lv, SETH_DR_T1, SETH_DR_T2, 1._R_P, rp, rm)
            endif
            if (lv == 0) then ! 7.
               if (rp /= rp0 .or. rm /= rm0) then ; nfail = nfail + 1 ; print '(A)', 'FAIL 7: level 0 changed a state' ; endif
               cycle
            endif
            ! 8. bounded, mirror
            lo = min(f(-1), f(1)) ; hi = max(f(-1), f(1))
            sc = NP_TOL*max(abs(lo), abs(hi), 1._R_P)
            if (rp < lo - sc .or. rp > hi + sc) then
               nfail = nfail + 1 ; if (nfail <= 3) print '(A,4ES16.8)', 'FAIL 8: outside [f_i-1, f_i+1] ', rp, lo, hi, sc
            endif
            lo = min(f(0), f(1)) ; hi = max(f(0), f(1))
            worst(1) = max(worst(1), max(lo - rp, rp - hi, 0._R_P)/max(abs(hi), abs(lo), 1._R_P))
            if (rp /= rm) then
               nfail = nfail + 1 ; if (nfail <= 3) print '(A,2ES24.16)', 'FAIL 8: mirror ', rp, rm
            endif
            fs = 0.5_R_P*(f(-1) + f(1)) ; fd = 0.5_R_P*(f(1) - f(-1))
            if (thinc) then ! 9. conservation: d from the face, the average of h(X) = fs + fd tanh(a + beta (X - 1))
               y = (rp - fs)/fd ! over the cell X in [0, 1], by quadrature (frozen) and in closed form (printed)
               if (abs(y) >= 1._R_P - 4._R_P*epsilon(1._R_P)) then ! the face saturated: d not recoverable from it
                  nsat = nsat + 1 ; cycle                           ! (R4P only: the face rounds to fs +- fd)
               endif
               a1 = atanh(y)
               avgq = 0._R_P ! four sub-cells: tanh's poles lie pi/(2 beta) off the axis, which limits ten points
               do ex2=0, 3   ! over the whole cell to ~4e-13; on a quarter of it ~1e-24
                  do q=1, 5
                     avgq = avgq + GL_W(q)*(tanh(a1 + SETH_THINC_BETA*(0.125_R_P*(2*ex2 + 1 - GL_X(q)) - 1._R_P)) + &
                                            tanh(a1 + SETH_THINC_BETA*(0.125_R_P*(2*ex2 + 1 + GL_X(q)) - 1._R_P)))
                  enddo
               enddo
               avgq = fs + 0.125_R_P*fd*avgq
               avg = fs + fd*(log(cosh(a1)) - log(cosh(a1 - SETH_THINC_BETA)))/SETH_THINC_BETA
               worst(2) = max(worst(2), abs(avgq - f(0))/maxval(abs(f(-1:1))))
               worst(7) = max(worst(7), abs(avg - f(0))/maxval(abs(f(-1:1))))
            else ! 10. the closed form, and the denominator of Eq. 31
               b = (2._R_P*f(0) - (f(1) + f(-1)))/(f(1) - f(-1))
               a1 = atanh(merge(SETH_DR_T2, SETH_DR_T1, lv == 2))
               ex = fs + fd*tanh(a1 + atanh(b))
               worst(3) = max(worst(3), abs(rp - ex)/maxval(abs(f(-1:1))))
               if (.not.(1._R_P + merge(SETH_DR_T2, SETH_DR_T1, lv == 2)*b > 0._R_P)) then
                  nfail = nfail + 1 ; if (nfail <= 3) print '(A,ES12.4)', 'FAIL 10: denominator not positive, B = ', b
               endif
            endif
            ! 11. scaling by 2^k with rmu/2^k, and an offset
            call random_number(y) ; ex2 = int(40._R_P*y) - 20 ; sc = 2._R_P**ex2
            sp = rp0 ; sm = rm0
            if (thinc) then
               call seth_thinc_override(sc*f, sc*g, lv, lv, SETH_THINC_BETA, 1._R_P/sc, sp, sm)
            else
               call seth_dr_override(sc*f, sc*g, lv, lv, SETH_DR_T1, SETH_DR_T2, 1._R_P/sc, sp, sm)
            endif
            if (sp /= sc*rp) then
               nfail = nfail + 1 ; if (nfail <= 3) print '(A,2ES24.16)', 'FAIL 11: scaling ', sp/sc, rp
            endif
            call random_number(mu) ; call random_number(y) ; mu = sign(1._R_P + 9._R_P*mu, y - 0.5_R_P) ! |mu| in [1, 10]
            sp = rp0 ; sm = rm0
            if (thinc) then
               call seth_thinc_override(f + mu, g + mu, lv, lv, SETH_THINC_BETA, 1._R_P, sp, sm)
            else
               call seth_dr_override(f + mu, g + mu, lv, lv, SETH_DR_T1, SETH_DR_T2, 1._R_P, sp, sm)
            endif
            worst(4) = max(worst(4), abs(sp - mu - rp)/abs(mu))
         enddo
      enddo
   enddo
   if (worst(1) > NP_TOL .or. worst(2) > NP_TOL .or. worst(3) > NP_TOL .or. worst(4) > NP_OFFSET_TOL) nfail = nfail + 1
   ! 10. linear data f_i +- d: f_i + d tanh(A)
   do n_=1, 1000
      call random_number(x) ; x = 4._R_P*(x - 0.5_R_P)
      if (abs(x(2)) < 1.e-3_R_P) cycle
      f = [x(1) - x(2), x(1), x(1) + x(2), x(1) + 2._R_P*x(2)]
      do lv=1, 2
         rp = 0._R_P ; rm = 0._R_P
         call seth_dr_override(f, f, lv, 0, SETH_DR_T1, SETH_DR_T2, 1._R_P, rp, rm)
         ex = x(1) + x(2)*merge(SETH_DR_T2, SETH_DR_T1, lv == 2)
         worst(5) = max(worst(5), abs(rp - ex)/(abs(x(1)) + abs(x(2))))
      enddo
   enddo
   if (worst(5) > NP_TOL) nfail = nfail + 1
   ! 10. any sign: finite; the guard
   do n_=1, NR
      call random_number(x) ; x = 2._R_P*(x - 0.5_R_P)
      call random_number(y) ; x = x*10._R_P**int(12._R_P*y - 6._R_P)
      f(-1:1) = x ; f(2) = 0._R_P
      rp = 0._R_P ; rm = 0._R_P ; sp = 0._R_P ; sm = 0._R_P
      call seth_dr_override(f, f, 2, 1, SETH_DR_T1, SETH_DR_T2, 1._R_P, rp, rm)
      call seth_thinc_override(f, f, 1, 1, SETH_THINC_BETA, 1._R_P, sp, sm)
      if (.not.(abs(rp) <= huge(1._R_P) .and. abs(rm) <= huge(1._R_P) .and. abs(sp) <= huge(1._R_P) .and. &
                abs(sm) <= huge(1._R_P))) then
         nfail = nfail + 1 ; if (nfail <= 3) print '(A,3ES12.4)', 'FAIL 10: not finite on ', x
      endif
      if ((f(1) - f(0))*(f(0) - f(-1)) <= SETH_NONPOLY_GUARD .and. (rp /= f(0) .or. sp /= f(0))) then
         nfail = nfail + 1 ; if (nfail <= 3) print '(A,3ES12.4)', 'FAIL 10: guard not applied on ', x
      endif
   enddo
   f = [1._R_P, 1._R_P, 1._R_P, 1._R_P] ; rp = 0._R_P ; rm = 0._R_P
   call seth_dr_override(f, f, 2, 2, SETH_DR_T1, SETH_DR_T2, 1._R_P, rp, rm)
   call seth_thinc_override(f, f, 1, 1, SETH_THINC_BETA, 1._R_P, sp, sm)
   if (rp /= 1._R_P .or. rm /= 1._R_P) nfail = nfail + 1
   print '(A,I0,A)', 'overrides on ', NR, ' monotone triples: level 0 unchanged, mirror bitwise, scaling bitwise'
   print '(A,I0)', '  THINC faces saturated (d not recoverable from the face; conservation and scaling not checked): ', nsat
   print '(A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2,A)', '  worst: outside [f_i, f_i+1] ', worst(1), &
         ', THINC conservation ', worst(2), ' (closed form ', worst(7), '), DR closed form ', worst(3), &
         ', offset / |mu| ', worst(4), ', linear ', worst(5)
   if (nfail > 0) then ; print '(A,I0)', 'FAIL 7-11: ', nfail ; passed = .false. ; endif
   endsubroutine overrides

   subroutine self_contained()
   !< Specs 12-13: N1 and N2 on periodic scalar fields of NF cells.
   integer(I4P), parameter :: NF = 200_I4P, NFIELD = 9_I4P
   character(14), parameter :: FNAME(NFIELD) = ['sine k=2      ', 'sine k=8      ', 'square wave   ', 'static (c)    ', &
                                                'static (d)    ', 'static (e)    ', 'jumps, cubics ', 'noise         ', &
                                                'gaussian      ']
   logical,       parameter :: SMOOTH(NFIELD) = [.true., .true., .false., .false., .false., .false., .false., .false., &
                                                 .true.]
   real(R_P)    :: u(0:NF-1), vp(-5:6), vm(-5:6), rp, rm, sp, sm, tp, tm, a(-2:3), b(-2:3), w(-2:2), x, hgt
   real(R_P), parameter :: S_RBF = 3._R_P*2._R_P/NF !< s = alpha dx, alpha = 3 (SETH_RBF_ALPHA1 = ALPHA2), dx = 2/NF
   integer(I4P) :: ml(0:NF-1), mr(0:NF-1), mp(-3:3), mm(-3:3), f_, j, m, c, lp, lm, nact(3), nfail, k, jp
   nfail = 0
   do f_=1, NFIELD
      do j=0, NF-1
         x = 2._R_P*(j + 0.5_R_P)/NF ! cell centres on [0, 2]
         select case(f_)
         case(1) ; u(j) = sin(2._R_P*acos(-1._R_P)*x)
         case(2) ; u(j) = sin(8._R_P*acos(-1._R_P)*x)
         case(3) ; u(j) = merge(1._R_P, 0._R_P, x > 0.5_R_P .and. x < 1.3_R_P)
         case(4) ; u(j) = static_c(x)
         case(5) ; u(j) = static_d(x)
         case(6) ; u(j) = static_e(x)
         case(7)
            if (j == 0) then ; hgt = 0._R_P ; jp = 0 ; endif
            if (j == jp) then ! a new piece every 7..30 cells
               call random_number(x) ; jp = j + 7 + int(23._R_P*x) ; call random_number(hgt) ; hgt = 4._R_P*hgt - 2._R_P
            endif
            u(j) = hgt + 0.1_R_P*sin(0.3_R_P*j)
         case(8) ; call random_number(u(j))
         case default ; u(j) = exp(-300._R_P*(x/2._R_P - 0.5_R_P)**2)
         endselect
      enddo
      ! pass 1: every cell's masks, stored
      do j=0, NF-1
         do m=-2, 2
            w(m) = u(modulo(j+m, NF))
         enddo
         call seth_teno5_mask(w, SETH_TENO_EPS, SETH_TENO5_CT, ml(j))
         do m=-2, 2
            w(m) = u(modulo(j-m, NF)) ! the mirrored data: the right state's mask of cell j
         enddo
         call seth_teno5_mask(w, SETH_TENO_EPS, SETH_TENO5_CT, mr(j))
      enddo
      nact = 0
      do k=1, 3 ! 1 DR, 2 THINC, 3 RBF+DR (s = 3 dx, dx = 2/NF; its base the RBF states, not TENO5's)
         do j=0, NF-1 ! face j+1/2
            do m=-5, 6
               vp(m) = u(modulo(j+m, NF)) ; vm(m) = vp(m)
            enddo
            ! pass 2 (the reference): stored masks, detector, override
            do m=-2, 3
               a(m) = vp(m) ; b(m) = vm(m)
            enddo
            if (k == 3) then
               call seth_teno5_rbf(a, b, SETH_TENO_EPS, SETH_TENO5_CT, SETH_RBF_CR, S_RBF, S_RBF, sp, sm)
            else
               call seth_teno5(a, b, SETH_TENO_EPS, SETH_TENO5_CT, sp, sm)
            endif
            tp = sp ; tm = sm
            do c=-3, 3
               mp(c) = ml(modulo(j+c, NF)) ; mm(c) = mr(modulo(j+1-c, NF))
            enddo
            if (k == 1) then
               call seth_dr_detect(mp, mm, lp, lm)
               call seth_dr_override(vp(-1:2), vm(-1:2), lp, lm, SETH_DR_T1, SETH_DR_T2, 1._R_P, sp, sm)
               call seth_teno5_dr(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, SETH_DR_T1, SETH_DR_T2, 1._R_P, rp, rm)
            elseif (k == 2) then
               call seth_thinc_detect(mp, mm, lp, lm)
               call seth_thinc_override(vp(-1:2), vm(-1:2), lp, lm, SETH_THINC_BETA, 1._R_P, sp, sm)
               call seth_teno5_thinc(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, SETH_THINC_BETA, 1._R_P, rp, rm)
            else
               call seth_dr_detect(mp, mm, lp, lm)
               call seth_dr_override(vp(-1:2), vm(-1:2), lp, lm, SETH_DR_T1, SETH_DR_T2, 1._R_P, sp, sm)
               call seth_teno5_rbf_dr(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, SETH_RBF_CR, S_RBF, S_RBF, SETH_DR_T1, &
                                      SETH_DR_T2, 1._R_P, rp, rm)
            endif
            if (lp > 0) nact(k) = nact(k) + 1
            if (lm > 0) nact(k) = nact(k) + 1
            if (rp /= sp .or. rm /= sm) then
               nfail = nfail + 1
               if (nfail <= 3) print '(A,A,I2,I5,4ES24.16)', 'FAIL 12: ', FNAME(f_), k, j, rp, sp, rm, sm
            endif
            if (SMOOTH(f_) .and. (lp > 0 .or. lm > 0 .or. rp /= tp .or. rm /= tm)) then
               nfail = nfail + 1 ; if (nfail <= 3) print '(A,A,I2,I5)', 'FAIL 13: detected on smooth data ', FNAME(f_), k, j
            endif
         enddo
      enddo
      print '(A,A,A,I4,A,I4,A,I4)', '  ', FNAME(f_), ' states overridden (of 400): DR ', nact(1), ', THINC ', nact(2), &
            ', RBF+DR ', nact(3)
   enddo
   print '(A)', 'self-contained kernels equal to the two-pass stored-mask reference, bitwise (N1); smooth fields inert (N2)'
   if (nfail > 0) then ; print '(A,I0)', 'FAIL 12-13: ', nfail ; passed = .false. ; endif
   endsubroutine self_contained

   pure function gdn(n, x) result(y)
   !< The n-th derivative of the generic g of seth_test_functions.
   integer(I4P), intent(in) :: n
   real(R_P),    intent(in) :: x
   real(R_P)                :: y
   real(R_P), parameter     :: HALF_PI = 2._R_P*atan(1._R_P)
   y = 0.9_R_P**n*sin(0.9_R_P*x + 0.4_R_P + n*HALF_PI) + 0.5_R_P*0.6_R_P**n*exp(0.6_R_P*x) &
       + 0.3_R_P*1.7_R_P**n*cos(1.7_R_P*x - 1.1_R_P + n*HALF_PI)
   endfunction gdn

   pure function static_c(x) result(y)
   !< Liang et al. 2022, static test (c): Jiang and Shu's multiple waves on [0, 2] (shifted from [-1, 1]).
   real(R_P), intent(in) :: x
   real(R_P)             :: y, xx, b
   real(R_P), parameter  :: A_ = 0.5_R_P, Z = -0.7_R_P, DT = 0.005_R_P, AL = 10._R_P
   xx = x - 1._R_P
   b = log(2._R_P)/(36._R_P*DT*DT)
   if (xx >= -0.8_R_P .and. xx <= -0.6_R_P) then
      y = (g(xx, b, Z - DT) + g(xx, b, Z + DT) + 4._R_P*g(xx, b, Z))/6._R_P
   elseif (xx >= -0.4_R_P .and. xx <= -0.2_R_P) then
      y = 1._R_P
   elseif (xx >= 0._R_P .and. xx <= 0.2_R_P) then
      y = 1._R_P - abs(10._R_P*(xx - 0.1_R_P))
   elseif (xx >= 0.4_R_P .and. xx <= 0.6_R_P) then
      y = (ff(xx, AL, A_ - DT) + ff(xx, AL, A_ + DT) + 4._R_P*ff(xx, AL, A_))/6._R_P
   else
      y = 0._R_P
   endif
   endfunction static_c

   pure function g(x, b, z) result(r)
   real(R_P), intent(in) :: x, b, z
   real(R_P)             :: r
   r = exp(-b*(x - z)**2)
   endfunction g

   pure function ff(x, al, a) result(r)
   real(R_P), intent(in) :: x, al, a
   real(R_P)             :: r
   r = sqrt(max(1._R_P - al*al*(x - a)**2, 0._R_P))
   endfunction ff

   pure function static_d(x) result(y)
   !< Liang et al. 2022, static test (d) (notes/m5-sources.md).
   real(R_P), intent(in) :: x
   real(R_P)             :: y
   if (x <= 0.4_R_P) then ; y = sin(8._R_P*x)
   elseif (x <= 0.5_R_P) then ; y = -1._R_P
   elseif (x <= 0.6_R_P) then ; y = 1._R_P
   elseif (x > 0.9_R_P .and. x <= 1.0_R_P) then ; y = 0._R_P
   elseif (x > 1.0_R_P .and. x <= 1.1_R_P) then ; y = 1._R_P
   elseif (x > 1.1_R_P .and. x <= 1.2_R_P) then ; y = 0.5_R_P
   elseif (x > 1.5_R_P .and. x <= 1.8_R_P) then ; y = exp(x) - 5._R_P
   else ; y = sin(9.42_R_P*x)
   endif
   endfunction static_d

   pure function static_e(x) result(y)
   !< Liang et al. 2022, static test (e) (notes/m5-sources.md).
   real(R_P), intent(in) :: x
   real(R_P)             :: y
   if ((x > 0.27_R_P .and. x <= 0.29_R_P) .or. (x > 0.63_R_P .and. x <= 0.66_R_P) .or. &
       (x > 1.14_R_P .and. x <= 1.17_R_P) .or. (x > 1.53_R_P .and. x <= 1.56_R_P)) then
      y = 1.25_R_P
   elseif ((x >= 0.25_R_P .and. x <= 0.27_R_P) .or. (x > 0.6_R_P .and. x <= 0.63_R_P) .or. &
           (x > 1.5_R_P .and. x <= 1.53_R_P)) then
      y = -1.25_R_P
   else
      y = sin(2._R_P*acos(-1._R_P)*x*exp(x + 1._R_P))
   endif
   endfunction static_e

   pure subroutine sort3(x)
   real(R_P), intent(inout) :: x(3)
   real(R_P)                :: t
   if (x(1) > x(2)) then ; t = x(1) ; x(1) = x(2) ; x(2) = t ; endif
   if (x(2) > x(3)) then ; t = x(2) ; x(2) = x(3) ; x(3) = t ; endif
   if (x(1) > x(2)) then ; t = x(1) ; x(1) = x(2) ; x(2) = t ; endif
   endsubroutine sort3

   subroutine exhaustive()
   !< Every combination of masks: DR on cells -3..2 (cell 3 set to the combination's cell -3, unread), THINC on -3..3;
   !< the left state takes the combination, the right state its reverse (another combination, the same coverage).
   integer(I4P) :: mk(-3:3), mr(-3:3), idx, t, c, lp, lm, nfail, n1, n2, nt
   nfail = 0 ; n1 = 0 ; n2 = 0 ; nt = 0
   do c=0, 7
      if (seth_teno5_zeta(c) /= zeta_ref(c)) then
         nfail = nfail + 1 ; print '(A,3I3)', 'zeta: ', c, seth_teno5_zeta(c), zeta_ref(c)
      endif
   enddo
   do idx=0, 8**6-1
      t = idx
      do c=-3, 2
         mk(c) = mod(t, 8) ; t = t/8
      enddo
      mk(3) = mk(-3)
      mr = mk(3:-3:-1)
      call seth_dr_detect(mk, mr, lp, lm)
      if (lp /= dr_ref(mk) .or. lm /= dr_ref(mr)) then
         nfail = nfail + 1 ; if (nfail <= 3) print '(A,7I2,2I3)', 'DR detector: ', mk, lp, dr_ref(mk)
      endif
      if (lp == 1) n1 = n1 + 1
      if (lp == 2) n2 = n2 + 1
   enddo
   do idx=0, 8**7-1
      t = idx
      do c=-3, 3
         mk(c) = mod(t, 8) ; t = t/8
      enddo
      mr = mk(3:-3:-1)
      call seth_thinc_detect(mk, mr, lp, lm)
      if (lp /= thinc_ref(mk) .or. lm /= thinc_ref(mr)) then
         nfail = nfail + 1 ; if (nfail <= 3) print '(A,7I2,2I3)', 'THINC detector: ', mk, lp, thinc_ref(mk)
      endif
      if (lp == 1) nt = nt + 1
   enddo
   print '(A,I0,A,I0,A,I0,A,I0,A)', 'detectors equal to the printed conditions on every combination: DR ', 8**6, &
         ' (level 1: ', n1, ', level 2: ', n2, '), THINC ', 8**7, ''
   print '(A,I0)', '  THINC level 1: ', nt
   if (nfail > 0) then ; print '(A,I0)', 'FAIL 5: ', nfail ; passed = .false. ; endif
   endsubroutine exhaustive

   subroutine levels(label, f, l5dr, l5th)
   !< DR and THINC levels of the left state at the faces c+1/2, c = 2..8, of the profile f(0:10): cell 5 at the given
   !< levels, every other cell at 0. Masks of cells -1..11 on the constant continuation.
   character(*), intent(in) :: label
   real(R_P),    intent(in) :: f(0:10)
   integer,      intent(in) :: l5dr, l5th
   integer(I4P)             :: mk(-1:11), c, m, lp, lm, ldr(2:8), lth(2:8)
   real(R_P)                :: w(-2:2)
   do c=-1, 11
      do m=-2, 2
         w(m) = f(min(max(c+m, 0), 10))
      enddo
      call seth_teno5_mask(w, NP_FLAGS_EPS, SETH_TENO5_CT, mk(c))
   enddo
   do c=2, 8
      call seth_dr_detect(mk(c-3:c+3), mk(c-3:c+3), ldr(c), lm)
      call seth_thinc_detect(mk(c-3:c+3), mk(c-3:c+3), lth(c), lm)
   enddo
   lp = merge(1, 0, all(ldr(2:4) == 0) .and. ldr(5) == l5dr .and. all(ldr(6:8) == 0) .and. &
                    all(lth(2:4) == 0) .and. lth(5) == l5th .and. all(lth(6:8) == 0))
   print '(A,T24,A,7I1,A,7I1,A)', label, 'DR levels (cells 2..8) ', ldr, ', THINC ', lth, merge(' ok  ', ' FAIL', lp == 1)
   if (lp /= 1) passed = .false.
   endsubroutine levels

   pure function pair_order(e) result(o)
   !< The order from the finest pair of grids whose finer error is above round-off (POINT_FLOOR eps); -1: none.
   real(R_P), intent(in) :: e(:)
   real(R_P)             :: o
   integer(I4P)          :: j
   o = -1._R_P
   do j=size(e)-1, 1, -1
      if (e(j+1) > POINT_FLOOR*epsilon(1._R_P)) then
         o = log(e(j)/e(j+1))/log(2._R_P) ; exit
      endif
   enddo
   endfunction pair_order

   pure function TENO5_CLO_(k) result(c)
   !< The first cell of TENO5's candidate k (S0 central, S1 downwind, S2 upwind).
   integer(I4P), intent(in) :: k
   integer(I4P)             :: c
   integer(I4P), parameter  :: CLO(0:2) = [-1, 0, -2]
   c = CLO(k)
   endfunction TENO5_CLO_

   subroutine random_stencil(vp, vm)
   !< A random stencil pair: a smooth cubic, independent values, a jump (random place and height) on a cubic, or a kink.
   real(R_P), intent(out) :: vp(-2:3), vm(-2:3)
   real(R_P)              :: y, cf(0:3), hgt
   integer(I4P)           :: m, kind_, jp
   call random_number(y) ; kind_ = int(4._R_P*y)
   call random_number(cf) ; cf = cf - 0.5_R_P
   call random_number(y) ; jp = -2 + int(6._R_P*y)
   call random_number(hgt) ; hgt = 10._R_P**(8._R_P*hgt - 6._R_P)
   select case(kind_)
   case(0)
      do m=-2, 3
         vp(m) = 1._R_P + cf(1)*m + cf(2)*m*m*0.1_R_P + cf(3)*m**3*0.01_R_P
         vm(m) = 2._R_P - cf(1)*m + cf(0)*m*m*0.1_R_P
      enddo
   case(1)
      call random_number(vp) ; call random_number(vm)
      vp = vp - 0.5_R_P ; vm = vm - 0.3_R_P
   case(2)
      do m=-2, 3
         vp(m) = cf(0) + cf(1)*m*0.1_R_P + merge(hgt, 0._R_P, m >= jp)
         vm(m) = cf(2) + cf(3)*m*0.1_R_P - merge(hgt, 0._R_P, m >= jp)
      enddo
   case default
      do m=-2, 3
         vp(m) = cf(0) + cf(1)*abs(m - jp + 0.5_R_P)
         vm(m) = cf(2) - cf(3)*abs(m - jp + 0.5_R_P)
      enddo
   endselect
   endsubroutine random_stencil

   subroutine table(key, label, f, d0, d1, d2)
   !< The masks of cells 2..8 of the profile f(0:10), constant continuation, against the printed flags; every cell
   !< 0..10 written to the report data.
   character(*), intent(in) :: key, label
   real(R_P),    intent(in) :: f(0:10)
   integer,      intent(in) :: d0(0:10), d1(0:10), d2(0:10)
   real(R_P)                :: w(-2:2)
   integer(I4P)             :: c, m, mk
   character(7)             :: got(0:2)
   logical                  :: ok
   ok = .true.
   got = ''
   if (vv_report()) then
      do c=0, 10
         do m=-2, 2
            w(m) = f(min(max(c+m, 0), 10))
         enddo
         call seth_teno5_mask(w, NP_FLAGS_EPS, SETH_TENO5_CT, mk)
         write(uf, '(A,1X,I0,1X,ES12.4,6(1X,I0))') key, c, f(c), ibits(mk, 0, 1), ibits(mk, 1, 1), ibits(mk, 2, 1), &
                                                   d0(c), d1(c), d2(c)
      enddo
   endif
   do c=2, 8
      do m=-2, 2
         w(m) = f(min(max(c+m, 0), 10))
      enddo
      call seth_teno5_mask(w, NP_FLAGS_EPS, SETH_TENO5_CT, mk)
      got(0)(c-1:c-1) = achar(48 + ibits(mk, 0, 1))
      got(1)(c-1:c-1) = achar(48 + ibits(mk, 1, 1))
      got(2)(c-1:c-1) = achar(48 + ibits(mk, 2, 1))
      if (ibits(mk, 0, 1) /= d0(c) .or. ibits(mk, 1, 1) /= d1(c) .or. ibits(mk, 2, 1) /= d2(c)) ok = .false.
   enddo
   print '(A,T34,A,1X,A,1X,A,A)', label, got(0), got(1), got(2), merge(' ok  ', ' FAIL', ok)
   if (.not.ok) then
      print '(A,T34,7I1,1X,7I1,1X,7I1)', '  printed (cells 2..8)', d0(2:8), d1(2:8), d2(2:8)
      passed = .false.
   endif
   endsubroutine table
endprogram seth_test_nonpoly
