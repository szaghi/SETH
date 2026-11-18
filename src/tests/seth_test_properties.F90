!< SETH test: the invariances and the robustness every scheme must have, on every scheme of the registry.
program seth_test_properties
!< SETH test: the invariances and the robustness every scheme must have, on every scheme of the registry.
!<
!< On random stencils (smooth and rough), for each scheme:
!<
!< 1. Exactness: data of a polynomial of degree < w_min (w_min = the narrowest candidate: S for WENO, 3 for TENO) is
!<    reproduced exactly at the face whatever the nonlinear weights, as every candidate is exact; cell averages for the
!<    reconstructions, point values at the cell centres for the interpolations. Tolerance 64 eps of the data.
!< 2. Constants: a constant stencil gives that constant (16 eps relative: the rounding of sum_j P_j c, sum |P| <= 6).
!< 3. Translation: R(v + c) = R(v) + c, to 1e-11 (|c| + max|v|): the indicators are blind to constants, but the
!<    nonlinear weights amplify the round-off of an indicator by about wexp/(eps + IS) (WENO-JS: up to 1e-13 seen);
!<    an indicator that is not blind to constants errs by O(1)*c, which this bound still catches.
!< 4. Mirror: with vm(m) = vp(1-m) the right state equals the left one: bitwise for WENO-Z and TENO (the same arithmetic
!<    on the mirrored data), to 1e-11 max|v| for WENO-JS (its two biases use mirrored tables, summed in another order,
!<    the difference amplified by the weights as in 3).
!< 5. Power-of-2 scaling: R(2^k v) = 2^k R(v) bitwise for WENO-Z and TENO (eps = 1e-40 is below one ulp of the
!<    indicators of O(1) data) and for WENO-JS with the scale-invariant descaler (rmu = seth_weno_descaling(sigma=1));
!<    plain WENO-JS (absolute eps = 1e-6) is not scale invariant, and the test also asserts that it is not, on data of
!<    magnitude 1e-4, so that the documentation's statement stays true. The adaptive cutoffs the same way (TENO5-AZ:
!<    Tang's min(C_T, tau); TENO5-A, TENO8-A: the sensor's regulariser xi^2): the papers' forms are dimensional
!<    (asserted not invariant), with the descaler (their rmu argument) they are bitwise invariant.
!< 6. Robustness: the zero stencil, and magnitudes from 1e-200 up to the scheme's limit, give finite states inside
!<    the data's range (the debug build traps overflow, invalid, zero divide). The indicators square the data, so 1e154
!<    bounds every scheme; plain WENO-JS also raises (eps + IS) to the power wexp (= S, S-1 for S = 5), which
!<    overflows near 1e77 (WENO-3), 1e51 (WENO-5), 1e38 (WENO-7/9): tested up to 1e30. With the scale-invariant
!<    descaler (eps + IS rmu^2 = O(1)) WENO-JS is tested up to 1e140, as WENO-Z and TENO are. The centred WENO, whose
!<    upwind part is WENO-JS, has WENO-JS's ranges (its blend's indicators take rmu too). These are the R8P
!<    ranges; R4P's (seth_test_parameters' ROBUST_EXP_*) stop where its squares do.
!< 7. A step of BIG_JUMP (1e30 in R8P) at every cell boundary of the stencil, both orientations (TENO family): finite
!<    states inside [0, BIG_JUMP] up to TOL_PROPERTIES of the jump (the composites' overrides included).
!< Every scheme of the registry, those beyond the sweeps (NSCHEMES+1..NALL: the stencil-targeted and Gaussian-RBF
!< kernels, rbf_dx at its default 0.01, so s = 0.03 is fixed while the data scale) included. The RBF rows scale with
!< the data (P + s^2 C, s a length), so 5 holds bitwise; their dependence on the length unit is E5 (seth_vv_nonpoly).
use seth
use seth_test_parameters, only : AI_STEP_TOL, ROBUST_EXP_MAX, ROBUST_EXP_MAX_JS, ROBUST_EXP_MIN, SMALL_STEP, &
                                 TOL_PROPERTIES, NP_GUARD_SCALE, BIG_JUMP
use seth_test_rng, only : seed_rng
use seth_test_schemes
implicit none
integer(I4P), parameter :: NR = 2000_I4P
real(R_P),    parameter :: E = epsilon(1._R_P)
real(R_P)               :: vp(VLO:VHI), vm(VLO:VHI), wp(VLO:VHI), wm(VLO:VHI), rp, rm, sp, sm, c, x, vmax, sc, worst
real(R_P)               :: coef(0:4), rmu
integer(I4P)            :: id, n_, m, wmin, deg, k, fam, s, failures
logical                 :: passed

passed = .true.
call seed_rng(1004_I4P)
do id=1, NALL
   fam = SCHEMES(id)%family ; s = SCHEMES(id)%s
   wmin = merge(3_I4P, s, fam == FAM_TENO)
   failures = 0
   do n_=1, NR
      ! 1. exactness on polynomials of degree < wmin
      deg = wmin - 1
      call random_number(coef) ; coef = 2._R_P*coef - 1._R_P
      do m=VLO, VHI
         if (SCHEMES(id)%interp) then
            vp(m) = poly(coef, deg, m - 0.5_R_P)
         else
            vp(m) = poly_avg(coef, deg, m)
         endif
      enddo
      call reconstruct(id, vp, vp, rp, rm)
      vmax = maxval(abs(vp(1-s:s)))
      if (abs(rp - coef(0)) > 64._R_P*E*vmax .or. abs(rm - coef(0)) > 64._R_P*E*vmax) call fail('exactness', rp, rm)
      ! 2. constants
      call random_number(c) ; c = 10._R_P*(c - 0.5_R_P)
      vp = c
      call reconstruct(id, vp, vp, rp, rm)
      if (abs(rp - c) > 16._R_P*E*abs(c) .or. abs(rm - c) > 16._R_P*E*abs(c)) call fail('constant', rp, rm)
      ! random stencils: half smooth (a random cubic), half rough
      call stencil(vp, vm)
      vmax = max(maxval(abs(vp)), maxval(abs(vm)))
      ! 3. translation
      call random_number(c) ; c = 4._R_P*(c - 0.5_R_P)
      call reconstruct(id, vp, vm, rp, rm)
      call reconstruct(id, vp + c, vm + c, sp, sm)
      if (abs(sp - c - rp) > TOL_PROPERTIES*(abs(c) + vmax) .or. abs(sm - c - rm) > TOL_PROPERTIES*(abs(c) + vmax)) &
         call fail('translation', sp - c - rp, sm - c - rm)
      ! 4. mirror
      do m=VLO, VHI
         wm(m) = 0._R_P
      enddo
      do m=merge(VLO, 1-s, SCHEMES(id)%wide), merge(VHI, s, SCHEMES(id)%wide) ! the composites read all of -5:6
         wm(m) = vp(1-m)
      enddo
      call reconstruct(id, vp, wm, rp, rm)
      if (fam == FAM_Z .or. fam == FAM_TENO) then
         if (rp /= rm) call fail('mirror (bitwise)', rp, rm)
      else
         if (abs(rp - rm) > TOL_PROPERTIES*vmax) call fail('mirror', rp, rm)
      endif
      ! 5. power-of-2 scaling
      call random_number(x) ; k = int(40._R_P*x) - 20
      sc = 2._R_P**k
      select case(fam)
      case(FAM_Z, FAM_TENO)
         if (dimensional(SCHEMES(id)%name)) then ! the papers' adaptive cutoffs are dimensional: the descaled forms
            call descaled(SCHEMES(id)%name, vp, vm, rp, rm)
            call descaled(SCHEMES(id)%name, sc*vp, sc*vm, sp, sm)
            if (sp /= sc*rp .or. sm /= sc*rm) call fail('scaling, descaled adaptive cutoff (bitwise)', sp/sc - rp, sm/sc - rm)
         else
            call reconstruct(id, vp, vm, rp, rm)
            call reconstruct(id, sc*vp, sc*vm, sp, sm)
            if (sp /= sc*rp .or. sm /= sc*rm) call fail('scaling (bitwise)', sp/sc - rp, sm/sc - rm)
         endif
      case default
         rmu = seth_weno_descaling(s, 1._R_P, vp(1-s:s), vm(1-s:s))
         call reconstruct(id, vp, vm, rp, rm, par=rmu)
         rmu = seth_weno_descaling(s, 1._R_P, sc*vp(1-s:s), sc*vm(1-s:s))
         call reconstruct(id, sc*vp, sc*vm, sp, sm, par=rmu)
         if (sp /= sc*rp .or. sm /= sc*rm) call fail('scaling, scale-invariant JS (bitwise)', sp/sc - rp, sm/sc - rm)
      endselect
      ! 6. robustness
      call random_number(x)
      if (fam == FAM_JS .or. fam == FAM_JS_INTERP .or. fam == FAM_C) then ! the centred WENO's upwind part is WENO-JS
         sc = 10._R_P**((ROBUST_EXP_MAX_JS - ROBUST_EXP_MIN)*x + ROBUST_EXP_MIN) ! plain WENO-JS (R8P: 1e-200 .. 1e30)
      else
         sc = 10._R_P**((ROBUST_EXP_MAX - ROBUST_EXP_MIN)*x + ROBUST_EXP_MIN) ! R8P: 1e-200 .. 1e140
      endif
      wp = sc*vp ; wm = sc*vm
      call reconstruct(id, wp, wm, rp, rm)
      vmax = sc*max(maxval(abs(vp)), maxval(abs(vm)))
      if (.not.(abs(rp) <= 4._R_P*vmax .and. abs(rm) <= 4._R_P*vmax)) call fail('robustness (finite, bounded)', rp, rm)
      if (fam == FAM_JS .or. fam == FAM_JS_INTERP .or. fam == FAM_C) then ! scale-invariant (R8P: 1e-200 .. 1e140)
         sc = 10._R_P**((ROBUST_EXP_MAX - ROBUST_EXP_MIN)*x + ROBUST_EXP_MIN)
         wp = sc*vp ; wm = sc*vm
         rmu = seth_weno_descaling(s, 1._R_P, wp(1-s:s), wm(1-s:s))
         call reconstruct(id, wp, wm, rp, rm, par=rmu)
         vmax = sc*max(maxval(abs(vp)), maxval(abs(vm)))
         if (.not.(abs(rp) <= 4._R_P*vmax .and. abs(rm) <= 4._R_P*vmax)) &
            call fail('robustness, scale-invariant JS (finite, bounded)', rp, rm)
      endif
      if (failures > 0) exit
   enddo
   ! the zero stencil
   vp = 0._R_P
   call reconstruct(id, vp, vp, rp, rm)
   if (rp /= 0._R_P .or. rm /= 0._R_P) call fail('zero stencil', rp, rm)
   ! 7. a step of BIG_JUMP at every cell boundary, both orientations
   if (fam == FAM_TENO) then
      do k=VLO+1, VHI
         do n_=0, 1
            do m=VLO, VHI
               vp(m) = merge(BIG_JUMP, 0._R_P, (m >= k) .neqv. (n_ == 1)) ; vm(m) = vp(m)
            enddo
            call reconstruct(id, vp, vm, rp, rm)
            if (.not.(rp >= -TOL_PROPERTIES*BIG_JUMP .and. rp <= (1._R_P + TOL_PROPERTIES)*BIG_JUMP .and. &
                      rm >= -TOL_PROPERTIES*BIG_JUMP .and. rm <= (1._R_P + TOL_PROPERTIES)*BIG_JUMP)) &
               call fail('a step of BIG_JUMP: a state not finite or outside the data', rp, rm)
         enddo
      enddo
   endif
   ! affine-invariant descaler (Wang & Don 2022): bitwise under 2^k scaling; a step of height 1e-7 on an offset of 1 keeps
   ! the ENO property (states within the data's range to 1e-12 of the jump), which the mean-|v| descaler loses
   if (fam == FAM_JS .and. s >= 3) then
      x = 0._R_P ; worst = 0._R_P
      do n_=1, 200
         call stencil(vp, vm, rough=.true.)
         call random_number(sc) ; k = int(40._R_P*sc) - 20 ; sc = 2._R_P**k
         call reconstruct(id, vp, vm, rp, rm, par=seth_weno_descaling_ai(s, vp(1-s:s), vm(1-s:s)))
         call reconstruct(id, sc*vp, sc*vm, sp, sm, par=seth_weno_descaling_ai(s, sc*vp(1-s:s), sc*vm(1-s:s)))
         if (sp /= sc*rp .or. sm /= sc*rm) call fail('scaling, affine-invariant JS (bitwise)', sp/sc - rp, sm/sc - rm)
      enddo
      c = 0._R_P
      do k=1-s, s ! a step of height SMALL_STEP (R8P 1e-7) at every cell boundary, on the offset 1
         do n_=VLO, VHI
            vp(n_) = 1._R_P + merge(SMALL_STEP, 0._R_P, n_ >= k) ; vm(n_) = vp(n_)
         enddo
         call reconstruct(id, vp, vm, rp, rm, par=seth_weno_descaling_ai(s, vp(1-s:s), vm(1-s:s)))
         x = max(x, (max(rp, rm) - (1._R_P + SMALL_STEP))/SMALL_STEP, (1._R_P - min(rp, rm))/SMALL_STEP)
         ! with the stencil's mean subtracted first and added back (Wang & Don, Remark 2): the expanded indicators of
         ! WENO-7/9 (Balsara-Shu integer forms) cancel the offset only to its round-off, which rmu^2 then amplifies
         wp = vp - sum(vp(1-s:s-1))/real(2*s-1, R_P) ; wm = wp
         call reconstruct(id, wp, wm, sp, sm, par=seth_weno_descaling_ai(s, wp(1-s:s), wm(1-s:s)))
         sp = sp + sum(vp(1-s:s-1))/real(2*s-1, R_P) ; sm = sm + sum(vp(1-s:s-1))/real(2*s-1, R_P)
         c = max(c, (max(sp, sm) - (1._R_P + SMALL_STEP))/SMALL_STEP, (1._R_P - min(sp, sm))/SMALL_STEP)
         call reconstruct(id, vp, vm, rp, rm, par=seth_weno_descaling(s, 1._R_P, vp(1-s:s), vm(1-s:s)))
         worst = max(worst, (max(rp, rm) - (1._R_P + SMALL_STEP))/SMALL_STEP, (1._R_P - min(rp, rm))/SMALL_STEP)
      enddo
      print '(A16,A,3ES10.2)', SCHEMES(id)%name, ': small step on 1, overshoot / jump: Ai, Ai on v - mean, mean-|v|: ', &
            x, c, worst
      if (c > 1.e-5_R_P) call fail('affine-invariant JS (mean subtracted): overshoot on a small step on an offset', c, x)
      ! the descaler alone rests on the round-off of the expanded indicators: WENO-5 1.6e-6 (gfortran), 3.1e-4
      ! (nvfortran, fused multiply-adds); asserted below 1e-3, the mean-subtracted form above to 1e-5
      if (s == 3 .and. x > AI_STEP_TOL) call fail('affine-invariant WENO-5: overshoot on a small step on an offset', x, c)
      if (worst < 1.e-2_R_P) call fail('mean-|v| descaler unexpectedly ENO on a small step on an offset', x, worst)
   endif
   ! plain WENO-JS and the adaptive cutoffs with rmu = 1 (the papers' forms) are not scale invariant: data of magnitude
   ! 1e-4 is reconstructed differently from the same data at 1
   ! (the composites' guard compares a product of two differences with 1e-15: it acts on data near 1e-9, NP_GUARD_SCALE)
   if (((fam == FAM_JS .or. fam == FAM_JS_INTERP) .and. s >= 2) .or. dimensional(SCHEMES(id)%name)) then
      x = 0._R_P
      sc = merge(NP_GUARD_SCALE, 1.e-4_R_P, SCHEMES(id)%wide)
      do n_=1, merge(2000, 50, SCHEMES(id)%wide) ! the composites differ only where a detector fires: a few % of stencils
         call stencil(vp, vm, rough=.true.)
         call reconstruct(id, vp, vm, rp, rm)
         call reconstruct(id, sc*vp, sc*vm, sp, sm)
         x = max(x, abs(sp/sc - rp))
      enddo
      if (x < 1.e-6_R_P) call fail('plain JS or an adaptive cutoff unexpectedly scale invariant', x, 0._R_P)
   endif
   print '(A16,A,I0,A)', SCHEMES(id)%name, ': ', failures, ' failures'
   if (failures > 0) passed = .false.
enddo
if (.not.passed) error stop 'property test failed'
print '(A)', 'seth_test_properties: ok'

contains
   subroutine fail(what, a, b)
   character(*), intent(in) :: what
   real(R_P),    intent(in) :: a, b
   failures = failures + 1
   if (failures <= 3) print '(A16,A,A,2ES24.16)', SCHEMES(id)%name, ': ', what//' ', a, b
   endsubroutine fail

   subroutine stencil(vp, vm, rough)
   !< Random stencils: a random cubic (smooth) or independent values (rough), two different fields.
   real(R_P), intent(out)          :: vp(VLO:VHI), vm(VLO:VHI)
   logical,   intent(in), optional :: rough
   real(R_P)                       :: y, cf(0:3)
   integer(I4P)                    :: m
   call random_number(y)
   if (present(rough)) y = 1._R_P
   if (y < 0.5_R_P) then
      call random_number(cf) ; cf = cf - 0.5_R_P
      do m=WLO, WHI
         vp(m) = 1._R_P + cf(1)*m + cf(2)*m*m*0.1_R_P + cf(3)*m**3*0.01_R_P
         vm(m) = 2._R_P - cf(1)*m + cf(0)*m*m*0.1_R_P
      enddo
   else
      call random_number(vp(WLO:WHI)) ; call random_number(vm(WLO:WHI))
      vp(WLO:WHI) = vp(WLO:WHI) - 0.5_R_P ; vm(WLO:WHI) = vm(WLO:WHI) - 0.3_R_P
   endif
   call extend(vp) ; call extend(vm) ! the draws and the maxima of the stencils as before the widening to VLO:VHI
   endsubroutine stencil

   pure function poly(cf, deg, x) result(y)
   real(R_P),    intent(in) :: cf(0:4), x
   integer(I4P), intent(in) :: deg
   real(R_P)                :: y
   integer(I4P)             :: n
   y = 0._R_P
   do n=deg, 0, -1
      y = y*x + cf(n)
   enddo
   endfunction poly

   pure function poly_avg(cf, deg, m) result(y)
   !< Average over cell m = [m-1, m] of sum cf(n) x**n.
   real(R_P),    intent(in) :: cf(0:4)
   integer(I4P), intent(in) :: deg, m
   real(R_P)                :: y
   integer(I4P)             :: n
   y = 0._R_P
   do n=0, deg
      y = y + cf(n)*(real(m, R_P)**(n+1) - real(m-1, R_P)**(n+1))/(n+1)
   enddo
   endfunction poly_avg

   pure function dimensional(name) result(yes)
   !< The schemes whose published form is not scale invariant (an adaptive cutoff with a dimensional quantity; the
   !< first-order guard of the DR and THINC overrides, a product of two differences against 1e-15).
   character(*), intent(in) :: name
   logical                  :: yes
   yes = any(trim(name) == [character(16) :: 'teno5_az', 'teno5_a', 'teno8_a', 'teno5_dr', 'teno5_thinc', &
                                              'teno5_st_dr', 'teno5_rbf_dr'])
   endfunction dimensional

   subroutine descaled(name, vp, vm, rp, rm)
   !< The scheme with the descaler of its stencil as rmu: the scale-invariant form of a dimensional cutoff.
   character(*), intent(in)  :: name
   real(R_P),    intent(in)  :: vp(VLO:VHI), vm(VLO:VHI)
   real(R_P),    intent(out) :: rp, rm
   select case(trim(name))
   case('teno5_az')
      call seth_teno5_az(vp(-2:3), vm(-2:3), SETH_TENO_EPS, SETH_TENO5_AZ_CT, &
                         seth_weno_descaling(3_I4P, 1._R_P, vp(-2:3), vm(-2:3)), rp, rm)
   case('teno5_a')
      call seth_teno5_a(vp(-2:3), vm(-2:3), SETH_TENO_EPS, SETH_TENO5_A_CR, SETH_TENO5_A_XI, SETH_TENO5_A_A1, &
                        SETH_TENO5_A_A2, seth_weno_descaling(3_I4P, 1._R_P, vp(-2:3), vm(-2:3)), rp, rm)
   case('teno8_a')
      call seth_teno8_a(vp(-3:4), vm(-3:4), SETH_TENO_EPS, SETH_TENO8_A_CR, SETH_TENO8_A_XI, SETH_TENO8_A_A1, &
                        SETH_TENO8_A_A2, seth_weno_descaling(4_I4P, 1._R_P, vp(-3:4), vm(-3:4)), rp, rm)
   case('teno5_dr')
      call seth_teno5_dr(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, SETH_DR_T1, SETH_DR_T2, &
                         seth_weno_descaling(3_I4P, 1._R_P, vp(-2:3), vm(-2:3)), rp, rm)
   case('teno5_thinc')
      call seth_teno5_thinc(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, SETH_THINC_BETA, &
                            seth_weno_descaling(3_I4P, 1._R_P, vp(-2:3), vm(-2:3)), rp, rm)
   case('teno5_st_dr')
      call seth_teno5_rbf_dr(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, SETH_RBF_CR, 0._R_P, 0._R_P, SETH_DR_T1, SETH_DR_T2, &
                             seth_weno_descaling(3_I4P, 1._R_P, vp(-2:3), vm(-2:3)), rp, rm)
   case('teno5_rbf_dr')
      call seth_teno5_rbf_dr(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, SETH_RBF_CR, SETH_RBF_ALPHA1*rbf_dx, &
                             SETH_RBF_ALPHA2*rbf_dx, SETH_DR_T1, SETH_DR_T2, &
                             seth_weno_descaling(3_I4P, 1._R_P, vp(-2:3), vm(-2:3)), rp, rm)
   endselect
   endsubroutine descaled
endprogram seth_test_properties
