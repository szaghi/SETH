!< SETH test: the centred WENO kernels reproduce ADAM's weno_reconstruct_centered with its defect D1 fixed.
program seth_test_weno_c_adam
!< SETH test: the centred WENO kernels reproduce ADAM's weno_reconstruct_centered with its defect D1 fixed.
!<
!< S = 2..4 (seth_weno4/6/8_c), random stencils, smooth (samples of a random quadratic) and rough (independent random
!< values), at magnitudes 1e-6..1e6, rmu = 1 (ADAM's blend has no descaler; SETH's with rmu = 1 is the same
!< arithmetic). ADAM's rows: v(1,m) = vm(m), v(2,m) = vp(m) on the big stencil 1-S..S; rp = vr(2), rm = vr(1). Against
!< seth_adam_reference:
!<
!< - weno_reconstruct_centered_d1fix (ADAM's routine with the left weights on the shifted cells): both states equal to
!<   TOL = 64 eps of max|stencil| (round-off: other summation orders and contractions);
!< - weno_reconstruct_centered, ADAM's routine as it is: the right state equal to TOL, the left one not (audit D1, its
!<   weights read the right-biased cells): asserted to differ beyond TOL on some rough stencils, the count printed.
!<
!< R8P only (ADAM's code is R8P): left out of the R4P builds (fobos [template-r4-tests] exclude).
use seth
use seth_test_rng, only : seed_rng
use seth_weno_js_tables
use seth_z_teno_tables, only : CEN4_C, CEN6_C, CEN8_C
use seth_adam_reference, only : weno_reconstruct_centered, weno_reconstruct_centered_d1fix
implicit none
integer(I4P), parameter :: N = 20000_I4P !< Random stencils per scheme.
real(R8P),    parameter :: TOL = 64._R8P*epsilon(1._R8P)
real(R8P), allocatable  :: a(:,:,:), p(:,:,:,:), d(:,:,:,:), c(:,:), v(:,:), vp(:), vm(:)
real(R8P)               :: vr(2), vr0(2), rp, rm, x, scale, diff, diff_right, diff_left_max, vmax
integer(I4P)            :: S, n_, m, nleft
logical                 :: passed

passed = .true.
call seed_rng(1012_I4P)
do S=2, 4
   call tables(S)
   if (allocated(v)) deallocate(v, vp, vm)
   allocate(v(2,1-S:S), vp(1-S:S), vm(1-S:S))
   diff = 0._R8P ; diff_right = 0._R8P ; diff_left_max = 0._R8P ; nleft = 0
   do n_=1, N
      call random_number(x) ; scale = 10._R8P**(12._R8P*x - 6._R8P)
      call random_number(x)
      if (x < 0.5_R8P) then ! smooth: samples of a random quadratic
         call random_number(x)
         do m=1-S, S
            vp(m) = scale*(1._R8P + x*m + 0.1_R8P*x*m*m)
            vm(m) = scale*(2._R8P - x*m + 0.2_R8P*x*m*m)
         enddo
      else                  ! rough: independent random values
         call random_number(vp) ; call random_number(vm)
         vp = scale*(vp - 0.5_R8P) ; vm = scale*(vm - 0.5_R8P)
      endif
      v(1,:) = vm ; v(2,:) = vp
      vmax = max(maxval(abs(vp)), maxval(abs(vm)), tiny(1._R8P))
      call seth_kernel(S, rp, rm)
      call weno_reconstruct_centered_d1fix(S, a, p, d, c, SETH_WENO_JS_EPS, 1._R8P, v, vr)
      diff = max(diff, abs(rp - vr(2))/vmax, abs(rm - vr(1))/vmax)
      call weno_reconstruct_centered(S, a, p, d, c, SETH_WENO_JS_EPS, 1._R8P, v, vr0)
      diff_right = max(diff_right, abs(rp - vr0(2))/vmax)
      if (abs(rm - vr0(1))/vmax > TOL) nleft = nleft + 1
      diff_left_max = max(diff_left_max, abs(rm - vr0(1))/vmax)
   enddo
   print '(A,I0,A,ES9.2,A,ES9.2,A,I0,A,I0,A,ES9.2,A)', 'S=', S, ': vs D1-fixed ', diff, '; vs ADAM as is: right ', &
         diff_right, ', left differs on ', nleft, ' of ', N, ' (max ', diff_left_max, ')'
   if (diff > TOL .or. diff_right > TOL .or. nleft == 0) passed = .false.
enddo
if (.not.passed) error stop 'SETH centred WENO differs from ADAM (D1 fixed) beyond round-off'
print '(A)', 'seth_test_weno_c_adam: ok'

contains
   subroutine tables(S)
   !< ADAM-layout tables of level S from SETH's constants (the other levels are never read at S).
   integer(I4P), intent(in) :: S
   if (allocated(a)) deallocate(a, p, d, c)
   allocate(a(1:2,0:S-1,1:S), p(1:2,0:S-1,0:S-1,1:S), d(0:S-1,0:S-1,0:S-1,1:S), c(1-S:S,1:S))
   a = 0._R8P ; p = 0._R8P ; d = 0._R8P ; c = 0._R8P
   select case(S)
   case(2) ; a(:,:,S) = WENO3_JS_A ; p(:,:,:,S) = WENO3_JS_P ; d(:,:,:,S) = WENO3_JS_D ; c(:,S) = CEN4_C
   case(3) ; a(:,:,S) = WENO5_JS_A ; p(:,:,:,S) = WENO5_JS_P ; d(:,:,:,S) = WENO5_JS_D ; c(:,S) = CEN6_C
   case(4) ; a(:,:,S) = WENO7_JS_A ; p(:,:,:,S) = WENO7_JS_P ; d(:,:,:,S) = WENO7_JS_D ; c(:,S) = CEN8_C
   endselect
   endsubroutine tables

   subroutine seth_kernel(S, rp, rm)
   integer(I4P), intent(in)  :: S
   real(R8P),    intent(out) :: rp, rm
   select case(S)
   case(2) ; call seth_weno4_c(vp, vm, SETH_WENO_JS_EPS, 1._R8P, rp, rm)
   case(3) ; call seth_weno6_c(vp, vm, SETH_WENO_JS_EPS, 1._R8P, rp, rm)
   case(4) ; call seth_weno8_c(vp, vm, SETH_WENO_JS_EPS, 1._R8P, rp, rm)
   endselect
   endsubroutine seth_kernel
endprogram seth_test_weno_c_adam
