!< SETH test: the WENO-JS kernels reproduce ADAM's (seth_adam_reference, ADAM's code verbatim) on random stencils.
program seth_test_weno_js_adam
!< SETH test: the WENO-JS kernels reproduce ADAM's (seth_adam_reference, ADAM's code verbatim) on random stencils.
!<
!< Every scheme (S = 1..5), reconstruction and interpolation tables, Jiang-Shu and scale-invariant weights, smooth and
!< rough stencils at magnitudes 1e-6..1e6: states, weights ratios and descaler compared with ADAM's on the same values
!< (ADAM's rows: v(2,m) = vp(m), v(1,m) = vm(m+1)). A state's difference is measured against the magnitude of its data,
!< max|stencil| (a state may cancel to near zero, so its own value is no scale), a weights ratio's absolutely (it lies
!< in [0,1]), the descaler's relatively. Passes when every difference is within TOL = 64 eps: round-off, including the
!< different fused multiply-add contractions an optimising compiler picks for the constant-folded SETH kernels and
!< ADAM's table-driven ones. Without contraction (the debug builds) the results are bitwise equal; the count is reported.
use seth
use seth_weno_js_tables
use seth_adam_reference, only : weno_reconstruct_upwind, weno_reconstruct_upwind_wratio, weno_descaling
implicit none
integer(I4P), parameter :: N = 20000_I4P !< Random stencils per scheme and table set.
real(R8P),    parameter :: TOL = 64._R8P * epsilon(1._R8P)
real(R8P), allocatable  :: a(:,:,:), p(:,:,:,:), d(:,:,:,:), v(:,:), vp(:), vm(:)
real(R8P)               :: vr(2), wr(2), rp, rm, wrp, wrm, rmu, rmu_ref, sigma, x, scale, diff, vmax
integer(I4P)            :: S, n_, m, t, bitwise, total
logical                 :: interp, passed

passed = .true. ; bitwise = 0 ; total = 0
call random_seed()
do S=1, 5
   do t=1, 2
      interp = t == 2
      call tables(S, interp)
      if (allocated(v)) deallocate(v, vp, vm)
      allocate(v(2,1-S:S-1), vp(1-S:S), vm(1-S:S))
      diff = 0._R8P
      do n_=1, N
         call random_number(x) ; scale = 10._R8P**(12._R8P*x - 6._R8P)
         call random_number(x)
         if (x < 0.5_R8P) then ! smooth: samples of a random quadratic
            call random_number(x)
            do m=1-S, S
               vp(m) = scale * (1._R8P + x*m + 0.1_R8P*x*m*m)
               vm(m) = scale * (2._R8P - x*m + 0.2_R8P*x*m*m)
            enddo
         else                  ! rough: independent random values
            call random_number(vp) ; call random_number(vm)
            vp = scale*(vp - 0.5_R8P) ; vm = scale*(vm - 0.5_R8P)
         endif
         do m=1-S, S-1
            v(1,m) = vm(m+1)
            v(2,m) = vp(m)
         enddo
         call random_number(x) ; sigma = merge(1._R8P, 0._R8P, x < 0.5_R8P)
         rmu     = seth_weno_descaling(S, sigma, vp, vm)
         rmu_ref = weno_descaling(S, sigma, v)
         vmax = max(maxval(abs(vp)), maxval(abs(vm)), tiny(1._R8P))
         call compare(rmu, rmu_ref, abs(rmu_ref))
         call weno_reconstruct_upwind_wratio(S, a, p, d, SETH_WENO_JS_EPS, rmu_ref, v, vr, wr)
         call seth_kernel(S, interp, rmu, rp, rm, wrp, wrm)
         call compare(rp, vr(2), vmax) ; call compare(rm, vr(1), vmax)
         call compare(wrp, wr(2), 1._R8P) ; call compare(wrm, wr(1), 1._R8P)
         call weno_reconstruct_upwind(S, a, p, d, SETH_WENO_JS_EPS, rmu_ref, v, vr)
         call seth_kernel_nowr(S, interp, rmu, rp, rm)
         call compare(rp, vr(2), vmax) ; call compare(rm, vr(1), vmax)
      enddo
      print '(A,I0,A,L1,A,ES9.2)', 'S=', S, ' interp=', interp, ' max scaled difference=', diff
      if (diff > TOL) passed = .false.
   enddo
enddo
print '(A,I0,A,I0,A)', 'bitwise equal: ', bitwise, ' of ', total, ' results'
if (.not.passed) error stop 'SETH WENO-JS differs from ADAM beyond round-off'
print '(A)', 'seth_test_weno_js_adam: ok'

contains
   subroutine compare(val, ref, mag)
   real(R8P), intent(in) :: val, ref !< SETH's and ADAM's result.
   real(R8P), intent(in) :: mag      !< The magnitude the difference is measured against.
   total = total + 1
   if (val == ref .or. (val /= val .and. ref /= ref)) then
      bitwise = bitwise + 1
   else
      diff = max(diff, abs(val - ref)/max(mag, tiny(1._R8P)))
   endif
   endsubroutine compare

   subroutine tables(S, interp)
   !< ADAM-layout tables of level S from SETH's constants (the other levels are never read at S).
   integer(I4P), intent(in) :: S
   logical,      intent(in) :: interp
   if (allocated(a)) deallocate(a, p, d)
   allocate(a(1:2,0:S-1,1:S), p(1:2,0:S-1,0:S-1,1:S), d(0:S-1,0:S-1,0:S-1,1:S))
   a = 0._R8P ; p = 0._R8P ; d = 0._R8P
   select case(S)
   case(1)
      a(:,:,S) = merge(WENO1_JS_INTERP_A, WENO1_JS_A, interp) ; p(:,:,:,S) = merge(WENO1_JS_INTERP_P, WENO1_JS_P, interp)
      d(:,:,:,S) = WENO1_JS_D
   case(2)
      a(:,:,S) = merge(WENO3_JS_INTERP_A, WENO3_JS_A, interp) ; p(:,:,:,S) = merge(WENO3_JS_INTERP_P, WENO3_JS_P, interp)
      d(:,:,:,S) = WENO3_JS_D
   case(3)
      a(:,:,S) = merge(WENO5_JS_INTERP_A, WENO5_JS_A, interp) ; p(:,:,:,S) = merge(WENO5_JS_INTERP_P, WENO5_JS_P, interp)
      d(:,:,:,S) = WENO5_JS_D
   case(4)
      a(:,:,S) = merge(WENO7_JS_INTERP_A, WENO7_JS_A, interp) ; p(:,:,:,S) = merge(WENO7_JS_INTERP_P, WENO7_JS_P, interp)
      d(:,:,:,S) = WENO7_JS_D
   case(5)
      a(:,:,S) = merge(WENO9_JS_INTERP_A, WENO9_JS_A, interp) ; p(:,:,:,S) = merge(WENO9_JS_INTERP_P, WENO9_JS_P, interp)
      d(:,:,:,S) = WENO9_JS_D
   endselect
   endsubroutine tables

   subroutine seth_kernel(S, interp, rmu, rp, rm, wrp, wrm)
   integer(I4P), intent(in)  :: S
   logical,      intent(in)  :: interp
   real(R8P),    intent(in)  :: rmu
   real(R8P),    intent(out) :: rp, rm, wrp, wrm
   real(R8P), parameter      :: E = SETH_WENO_JS_EPS
   select case(S)
   case(1) ; if (interp) then ; call seth_weno1_js_interp_wr(vp, vm, E, rmu, rp, rm, wrp, wrm)
             else ; call seth_weno1_js_wr(vp, vm, E, rmu, rp, rm, wrp, wrm) ; endif
   case(2) ; if (interp) then ; call seth_weno3_js_interp_wr(vp, vm, E, rmu, rp, rm, wrp, wrm)
             else ; call seth_weno3_js_wr(vp, vm, E, rmu, rp, rm, wrp, wrm) ; endif
   case(3) ; if (interp) then ; call seth_weno5_js_interp_wr(vp, vm, E, rmu, rp, rm, wrp, wrm)
             else ; call seth_weno5_js_wr(vp, vm, E, rmu, rp, rm, wrp, wrm) ; endif
   case(4) ; if (interp) then ; call seth_weno7_js_interp_wr(vp, vm, E, rmu, rp, rm, wrp, wrm)
             else ; call seth_weno7_js_wr(vp, vm, E, rmu, rp, rm, wrp, wrm) ; endif
   case(5) ; if (interp) then ; call seth_weno9_js_interp_wr(vp, vm, E, rmu, rp, rm, wrp, wrm)
             else ; call seth_weno9_js_wr(vp, vm, E, rmu, rp, rm, wrp, wrm) ; endif
   endselect
   endsubroutine seth_kernel

   subroutine seth_kernel_nowr(S, interp, rmu, rp, rm)
   integer(I4P), intent(in)  :: S
   logical,      intent(in)  :: interp
   real(R8P),    intent(in)  :: rmu
   real(R8P),    intent(out) :: rp, rm
   real(R8P), parameter      :: E = SETH_WENO_JS_EPS
   select case(S)
   case(1) ; if (interp) then ; call seth_weno1_js_interp(vp, vm, E, rmu, rp, rm)
             else ; call seth_weno1_js(vp, vm, E, rmu, rp, rm) ; endif
   case(2) ; if (interp) then ; call seth_weno3_js_interp(vp, vm, E, rmu, rp, rm)
             else ; call seth_weno3_js(vp, vm, E, rmu, rp, rm) ; endif
   case(3) ; if (interp) then ; call seth_weno5_js_interp(vp, vm, E, rmu, rp, rm)
             else ; call seth_weno5_js(vp, vm, E, rmu, rp, rm) ; endif
   case(4) ; if (interp) then ; call seth_weno7_js_interp(vp, vm, E, rmu, rp, rm)
             else ; call seth_weno7_js(vp, vm, E, rmu, rp, rm) ; endif
   case(5) ; if (interp) then ; call seth_weno9_js_interp(vp, vm, E, rmu, rp, rm)
             else ; call seth_weno9_js(vp, vm, E, rmu, rp, rm) ; endif
   endselect
   endsubroutine seth_kernel_nowr
endprogram seth_test_weno_js_adam
