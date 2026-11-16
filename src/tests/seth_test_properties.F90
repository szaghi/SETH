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
!<    descaler (eps + IS rmu^2 = O(1)) WENO-JS is tested up to 1e140, as WENO-Z and TENO are.
use seth
use seth_test_rng, only : seed_rng
use seth_test_schemes
implicit none
integer(I4P), parameter :: NR = 2000_I4P
real(R8P),    parameter :: E = epsilon(1._R8P)
real(R8P)               :: vp(VLO:VHI), vm(VLO:VHI), wp(VLO:VHI), wm(VLO:VHI), rp, rm, sp, sm, c, x, vmax, sc, worst
real(R8P)               :: coef(0:4), rmu
integer(I4P)            :: id, n_, m, wmin, deg, k, fam, s, failures
logical                 :: passed

passed = .true.
call seed_rng(1004_I4P)
do id=1, NSCHEMES
   fam = SCHEMES(id)%family ; s = SCHEMES(id)%s
   wmin = merge(3_I4P, s, fam == FAM_TENO)
   failures = 0
   do n_=1, NR
      ! 1. exactness on polynomials of degree < wmin
      deg = wmin - 1
      call random_number(coef) ; coef = 2._R8P*coef - 1._R8P
      do m=VLO, VHI
         if (SCHEMES(id)%interp) then
            vp(m) = poly(coef, deg, m - 0.5_R8P)
         else
            vp(m) = poly_avg(coef, deg, m)
         endif
      enddo
      call reconstruct(id, vp, vp, rp, rm)
      vmax = maxval(abs(vp(1-s:s)))
      if (abs(rp - coef(0)) > 64._R8P*E*vmax .or. abs(rm - coef(0)) > 64._R8P*E*vmax) call fail('exactness', rp, rm)
      ! 2. constants
      call random_number(c) ; c = 10._R8P*(c - 0.5_R8P)
      vp = c
      call reconstruct(id, vp, vp, rp, rm)
      if (abs(rp - c) > 16._R8P*E*abs(c) .or. abs(rm - c) > 16._R8P*E*abs(c)) call fail('constant', rp, rm)
      ! random stencils: half smooth (a random cubic), half rough
      call stencil(vp, vm)
      vmax = max(maxval(abs(vp)), maxval(abs(vm)))
      ! 3. translation
      call random_number(c) ; c = 4._R8P*(c - 0.5_R8P)
      call reconstruct(id, vp, vm, rp, rm)
      call reconstruct(id, vp + c, vm + c, sp, sm)
      if (abs(sp - c - rp) > 1.e-11_R8P*(abs(c) + vmax) .or. abs(sm - c - rm) > 1.e-11_R8P*(abs(c) + vmax)) &
         call fail('translation', sp - c - rp, sm - c - rm)
      ! 4. mirror
      do m=VLO, VHI
         wm(m) = 0._R8P
      enddo
      do m=1-s, s
         wm(m) = vp(1-m)
      enddo
      call reconstruct(id, vp, wm, rp, rm)
      if (fam == FAM_Z .or. fam == FAM_TENO) then
         if (rp /= rm) call fail('mirror (bitwise)', rp, rm)
      else
         if (abs(rp - rm) > 1.e-11_R8P*vmax) call fail('mirror', rp, rm)
      endif
      ! 5. power-of-2 scaling
      call random_number(x) ; k = int(40._R8P*x) - 20
      sc = 2._R8P**k
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
         rmu = seth_weno_descaling(s, 1._R8P, vp(1-s:s), vm(1-s:s))
         call reconstruct(id, vp, vm, rp, rm, par=rmu)
         rmu = seth_weno_descaling(s, 1._R8P, sc*vp(1-s:s), sc*vm(1-s:s))
         call reconstruct(id, sc*vp, sc*vm, sp, sm, par=rmu)
         if (sp /= sc*rp .or. sm /= sc*rm) call fail('scaling, scale-invariant JS (bitwise)', sp/sc - rp, sm/sc - rm)
      endselect
      ! 6. robustness
      call random_number(x)
      if (fam == FAM_JS .or. fam == FAM_JS_INTERP) then
         sc = 10._R8P**(230._R8P*x - 200._R8P) ! plain WENO-JS: 1e-200 .. 1e30
      else
         sc = 10._R8P**(340._R8P*x - 200._R8P) ! 1e-200 .. 1e140
      endif
      wp = sc*vp ; wm = sc*vm
      call reconstruct(id, wp, wm, rp, rm)
      vmax = sc*max(maxval(abs(vp)), maxval(abs(vm)))
      if (.not.(abs(rp) <= 4._R8P*vmax .and. abs(rm) <= 4._R8P*vmax)) call fail('robustness (finite, bounded)', rp, rm)
      if (fam == FAM_JS .or. fam == FAM_JS_INTERP) then ! scale-invariant WENO-JS: 1e-200 .. 1e140
         sc = 10._R8P**(340._R8P*x - 200._R8P)
         wp = sc*vp ; wm = sc*vm
         rmu = seth_weno_descaling(s, 1._R8P, wp(1-s:s), wm(1-s:s))
         call reconstruct(id, wp, wm, rp, rm, par=rmu)
         vmax = sc*max(maxval(abs(vp)), maxval(abs(vm)))
         if (.not.(abs(rp) <= 4._R8P*vmax .and. abs(rm) <= 4._R8P*vmax)) &
            call fail('robustness, scale-invariant JS (finite, bounded)', rp, rm)
      endif
      if (failures > 0) exit
   enddo
   ! the zero stencil
   vp = 0._R8P
   call reconstruct(id, vp, vp, rp, rm)
   if (rp /= 0._R8P .or. rm /= 0._R8P) call fail('zero stencil', rp, rm)
   ! affine-invariant descaler (Wang & Don 2022): bitwise under 2^k scaling; a step of height 1e-7 on an offset of 1 keeps
   ! the ENO property (states within the data's range to 1e-12 of the jump), which the mean-|v| descaler loses
   if (fam == FAM_JS .and. s >= 3) then
      x = 0._R8P ; worst = 0._R8P
      do n_=1, 200
         call stencil(vp, vm, rough=.true.)
         call random_number(sc) ; k = int(40._R8P*sc) - 20 ; sc = 2._R8P**k
         call reconstruct(id, vp, vm, rp, rm, par=seth_weno_descaling_ai(s, vp(1-s:s), vm(1-s:s)))
         call reconstruct(id, sc*vp, sc*vm, sp, sm, par=seth_weno_descaling_ai(s, sc*vp(1-s:s), sc*vm(1-s:s)))
         if (sp /= sc*rp .or. sm /= sc*rm) call fail('scaling, affine-invariant JS (bitwise)', sp/sc - rp, sm/sc - rm)
      enddo
      c = 0._R8P
      do k=1-s, s ! a step of height 1e-7 at every cell boundary, on the offset 1
         do n_=VLO, VHI
            vp(n_) = 1._R8P + merge(1.e-7_R8P, 0._R8P, n_ >= k) ; vm(n_) = vp(n_)
         enddo
         call reconstruct(id, vp, vm, rp, rm, par=seth_weno_descaling_ai(s, vp(1-s:s), vm(1-s:s)))
         x = max(x, (max(rp, rm) - (1._R8P + 1.e-7_R8P))/1.e-7_R8P, (1._R8P - min(rp, rm))/1.e-7_R8P)
         ! with the stencil's mean subtracted first and added back (Wang & Don, Remark 2): the expanded indicators of
         ! WENO-7/9 (Balsara-Shu integer forms) cancel the offset only to its round-off, which rmu^2 then amplifies
         wp = vp - sum(vp(1-s:s-1))/real(2*s-1, R8P) ; wm = wp
         call reconstruct(id, wp, wm, sp, sm, par=seth_weno_descaling_ai(s, wp(1-s:s), wm(1-s:s)))
         sp = sp + sum(vp(1-s:s-1))/real(2*s-1, R8P) ; sm = sm + sum(vp(1-s:s-1))/real(2*s-1, R8P)
         c = max(c, (max(sp, sm) - (1._R8P + 1.e-7_R8P))/1.e-7_R8P, (1._R8P - min(sp, sm))/1.e-7_R8P)
         call reconstruct(id, vp, vm, rp, rm, par=seth_weno_descaling(s, 1._R8P, vp(1-s:s), vm(1-s:s)))
         worst = max(worst, (max(rp, rm) - (1._R8P + 1.e-7_R8P))/1.e-7_R8P, (1._R8P - min(rp, rm))/1.e-7_R8P)
      enddo
      print '(A16,A,3ES10.2)', SCHEMES(id)%name, ': 1e-7 step on 1, overshoot / jump: Ai, Ai on v - mean, mean-|v|: ', &
            x, c, worst
      if (c > 1.e-5_R8P) call fail('affine-invariant JS (mean subtracted): overshoot on a small step on an offset', c, x)
      if (s == 3 .and. x > 1.e-5_R8P) call fail('affine-invariant WENO-5: overshoot on a small step on an offset', x, c)
      if (worst < 1.e-2_R8P) call fail('mean-|v| descaler unexpectedly ENO on a small step on an offset', x, worst)
   endif
   ! plain WENO-JS and the adaptive cutoffs with rmu = 1 (the papers' forms) are not scale invariant: data of magnitude
   ! 1e-4 is reconstructed differently from the same data at 1
   if (((fam == FAM_JS .or. fam == FAM_JS_INTERP) .and. s >= 2) .or. dimensional(SCHEMES(id)%name)) then
      x = 0._R8P
      do n_=1, 50
         call stencil(vp, vm, rough=.true.)
         call reconstruct(id, vp, vm, rp, rm)
         call reconstruct(id, 1.e-4_R8P*vp, 1.e-4_R8P*vm, sp, sm)
         x = max(x, abs(sp*1.e4_R8P - rp))
      enddo
      if (x < 1.e-6_R8P) call fail('plain JS or an adaptive cutoff unexpectedly scale invariant', x, 0._R8P)
   endif
   print '(A16,A,I0,A)', SCHEMES(id)%name, ': ', failures, ' failures'
   if (failures > 0) passed = .false.
enddo
if (.not.passed) error stop 'property test failed'
print '(A)', 'seth_test_properties: ok'

contains
   subroutine fail(what, a, b)
   character(*), intent(in) :: what
   real(R8P),    intent(in) :: a, b
   failures = failures + 1
   if (failures <= 3) print '(A16,A,A,2ES24.16)', SCHEMES(id)%name, ': ', what//' ', a, b
   endsubroutine fail

   subroutine stencil(vp, vm, rough)
   !< Random stencils: a random cubic (smooth) or independent values (rough), two different fields.
   real(R8P), intent(out)          :: vp(VLO:VHI), vm(VLO:VHI)
   logical,   intent(in), optional :: rough
   real(R8P)                       :: y, cf(0:3)
   integer(I4P)                    :: m
   call random_number(y)
   if (present(rough)) y = 1._R8P
   if (y < 0.5_R8P) then
      call random_number(cf) ; cf = cf - 0.5_R8P
      do m=VLO, VHI
         vp(m) = 1._R8P + cf(1)*m + cf(2)*m*m*0.1_R8P + cf(3)*m**3*0.01_R8P
         vm(m) = 2._R8P - cf(1)*m + cf(0)*m*m*0.1_R8P
      enddo
   else
      call random_number(vp) ; call random_number(vm)
      vp = vp - 0.5_R8P ; vm = vm - 0.3_R8P
   endif
   endsubroutine stencil

   pure function poly(cf, deg, x) result(y)
   real(R8P),    intent(in) :: cf(0:4), x
   integer(I4P), intent(in) :: deg
   real(R8P)                :: y
   integer(I4P)             :: n
   y = 0._R8P
   do n=deg, 0, -1
      y = y*x + cf(n)
   enddo
   endfunction poly

   pure function poly_avg(cf, deg, m) result(y)
   !< Average over cell m = [m-1, m] of sum cf(n) x**n.
   real(R8P),    intent(in) :: cf(0:4)
   integer(I4P), intent(in) :: deg, m
   real(R8P)                :: y
   integer(I4P)             :: n
   y = 0._R8P
   do n=0, deg
      y = y + cf(n)*(real(m, R8P)**(n+1) - real(m-1, R8P)**(n+1))/(n+1)
   enddo
   endfunction poly_avg

   pure function dimensional(name) result(yes)
   !< The schemes whose published form is not scale invariant (an adaptive cutoff with a dimensional quantity).
   character(*), intent(in) :: name
   logical                  :: yes
   yes = any(trim(name) == [character(16) :: 'teno5_az', 'teno5_a', 'teno8_a'])
   endfunction dimensional

   subroutine descaled(name, vp, vm, rp, rm)
   !< The scheme with the descaler of its stencil as rmu: the scale-invariant form of a dimensional cutoff.
   character(*), intent(in)  :: name
   real(R8P),    intent(in)  :: vp(VLO:VHI), vm(VLO:VHI)
   real(R8P),    intent(out) :: rp, rm
   select case(trim(name))
   case('teno5_az')
      call seth_teno5_az(vp(-2:3), vm(-2:3), SETH_TENO_EPS, SETH_TENO5_AZ_CT, &
                         seth_weno_descaling(3_I4P, 1._R8P, vp(-2:3), vm(-2:3)), rp, rm)
   case('teno5_a')
      call seth_teno5_a(vp(-2:3), vm(-2:3), SETH_TENO_EPS, SETH_TENO5_A_CR, SETH_TENO5_A_XI, SETH_TENO5_A_A1, &
                        SETH_TENO5_A_A2, seth_weno_descaling(3_I4P, 1._R8P, vp(-2:3), vm(-2:3)), rp, rm)
   case('teno8_a')
      call seth_teno8_a(vp(-3:4), vm(-3:4), SETH_TENO_EPS, SETH_TENO8_A_CR, SETH_TENO8_A_XI, SETH_TENO8_A_A1, &
                        SETH_TENO8_A_A2, seth_weno_descaling(4_I4P, 1._R8P, vp(-3:4), vm(-3:4)), rp, rm)
   endselect
   endsubroutine descaled
endprogram seth_test_properties
