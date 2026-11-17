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
!<    data's magnitude (the conservation constraint, by the closed-form integral of tanh).
!< 10. DR: the state equals m + h tanh(A + atanh B) to NP_TOL; on linear data f_i +- d it is f_i + d tanh(A) to NP_TOL;
!<    no NaN or Inf on 10^6 random triples of any sign; the flat and the non-monotone triples give f_i (the guard).
!< 11. Scaling: f -> sigma f with rmu -> rmu/sigma, sigma = 2^k, bitwise sigma times the state; f -> f + mu equal to
!<    NP_TOL of |mu| + |f|.
!< Step 2.3, the self-contained kernels seth_teno5_dr and seth_teno5_thinc (issue #6, section 5):
!< 12. N1: on periodic scalar fields (smooth sines, a square wave, Liang et al.'s static profiles (c)-(e), random
!<    piecewise cubics with jumps, noise) the kernels equal a two-pass reference written here, bitwise: every cell's
!<    mask computed once and stored (left-biased, and on the mirrored data for the right state), the stored masks read
!<    by the detector, TENO5's states overridden. The reference shares only seth_teno5, seth_teno5_mask and the
!<    detectors and overrides with the kernels (no early exit, no neighbour recomputation).
!< 13. N2: on the smooth fields no face is detected and the states equal seth_teno5's bitwise.
!< 4. Constant data and a jump of BIG_JUMP: the mask is defined and never 0 (one chi_k >= 1/3 > C_T). Not 7 on
!<    constant data in general (printed): the indicators of a non-dyadic constant are round-off, far above eps.
use seth
use seth_test_parameters, only : BIG_JUMP, NP_FLAGS_EPS, NP_TOL
use seth_test_rng, only : seed_rng
use seth_test_teno_flags, only : teno5_flags
use seth_vv_cases, only : vv_report, vv_open
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

   subroutine overrides()
   !< Specs 7-11 of the two overrides.
   real(R_P)    :: f(-1:2), g(-1:2), w(-1:2), rp, rm, sp, sm, x(3), lo, hi, worst(6), fs, fd, y, avg, b, sc, mu, ex
   real(R_P)    :: rp0, rm0, a1
   integer(I4P) :: n_, lv, k, nfail, ex2
   logical      :: thinc
   nfail = 0 ; worst = 0._R_P
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
            if (thinc) then ! 9. conservation: d from the face, the average of h over the cell
               y = (rp - fs)/fd
               avg = fs + fd*(log(cosh(atanh(y))) - log(cosh(atanh(y) - SETH_THINC_BETA)))/SETH_THINC_BETA
               worst(2) = max(worst(2), abs(avg - f(0))/maxval(abs(f(-1:1))))
            else ! 10. the closed form
               b = (2._R_P*f(0) - (f(1) + f(-1)))/(f(1) - f(-1))
               a1 = atanh(merge(SETH_DR_T2, SETH_DR_T1, lv == 2))
               ex = fs + fd*tanh(a1 + atanh(b))
               worst(3) = max(worst(3), abs(rp - ex)/maxval(abs(f(-1:1))))
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
            call random_number(mu) ; mu = 20._R_P*(mu - 0.5_R_P)
            sp = rp0 ; sm = rm0
            if (thinc) then
               call seth_thinc_override(f + mu, g + mu, lv, lv, SETH_THINC_BETA, 1._R_P, sp, sm)
            else
               call seth_dr_override(f + mu, g + mu, lv, lv, SETH_DR_T1, SETH_DR_T2, 1._R_P, sp, sm)
            endif
            worst(4) = max(worst(4), abs(sp - mu - rp)/(abs(mu) + maxval(abs(f(-1:1)))))
         enddo
      enddo
   enddo
   if (worst(1) > NP_TOL .or. worst(2) > NP_TOL .or. worst(3) > NP_TOL .or. worst(4) > NP_TOL) nfail = nfail + 1
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
   print '(A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2)', '  worst: outside [f_i, f_i+1] ', worst(1), &
         ', THINC conservation ', worst(2), ', DR closed form ', worst(3), ', offset ', worst(4), ', linear ', worst(5)
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
   integer(I4P) :: ml(0:NF-1), mr(0:NF-1), mp(-3:3), mm(-3:3), f_, j, m, c, lp, lm, nact(2), nfail, k, jp
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
      do k=1, 2 ! 1 DR, 2 THINC
         do j=0, NF-1 ! face j+1/2
            do m=-5, 6
               vp(m) = u(modulo(j+m, NF)) ; vm(m) = vp(m)
            enddo
            ! pass 2 (the reference): stored masks, detector, override
            do m=-2, 3
               a(m) = vp(m) ; b(m) = vm(m)
            enddo
            call seth_teno5(a, b, SETH_TENO_EPS, SETH_TENO5_CT, sp, sm)
            tp = sp ; tm = sm
            do c=-3, 3
               mp(c) = ml(modulo(j+c, NF)) ; mm(c) = mr(modulo(j+1-c, NF))
            enddo
            if (k == 1) then
               call seth_dr_detect(mp, mm, lp, lm)
               call seth_dr_override(vp(-1:2), vm(-1:2), lp, lm, SETH_DR_T1, SETH_DR_T2, 1._R_P, sp, sm)
               call seth_teno5_dr(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, SETH_DR_T1, SETH_DR_T2, 1._R_P, rp, rm)
            else
               call seth_thinc_detect(mp, mm, lp, lm)
               call seth_thinc_override(vp(-1:2), vm(-1:2), lp, lm, SETH_THINC_BETA, 1._R_P, sp, sm)
               call seth_teno5_thinc(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, SETH_THINC_BETA, 1._R_P, rp, rm)
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
      print '(A,A,A,I4,A,I4)', '  ', FNAME(f_), ' states overridden (of 400): DR ', nact(1), ', THINC ', nact(2)
   enddo
   print '(A)', 'self-contained kernels equal to the two-pass stored-mask reference, bitwise (N1); smooth fields inert (N2)'
   if (nfail > 0) then ; print '(A,I0)', 'FAIL 12-13: ', nfail ; passed = .false. ; endif
   endsubroutine self_contained

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
