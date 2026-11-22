!< SETH test: seth_weno6_point and seth_weno6_jumps reproduce ADAM's seam WENO-6 restriction and fill bit for bit.
program seth_test_weno6_point_adam
!< SETH test: seth_weno6_point and seth_weno6_jumps reproduce ADAM's seam WENO-6 restriction and fill bit for bit.
!<
!< The reference is seth_adam_seam_reference, ADAM's routines copied verbatim from ADAM commit 2ea71524 (issue #7, E2).
!< First the tables: every linear and cubic weight SETH compiles at a midpoint t4 = 4 code + 6 equals ADAM's
!< restriction table (SEAM_D_WENO6, SEAM_W_CUBIC_MID) and at a quarter t4 = 4 p + 2 sub - 3 ADAM's fill table
!< (SEAM_D_FILL_WENO6, SEAM_W_CUBIC_Q), bit for bit, and so do the regulariser and the jump ratio. Then the results, on
!< N random lines of six values of seven kinds (smooth samples, independent random values, steps, constants, a flat
!< sub-stencil next to a ramp, ADAM's sharpened rnd**8, small variations on a large offset) at magnitudes
!< 1e-12..1e12: seth_weno6_point(u, 4 code + 6) against seam_restriction_weno6_1d(u, code) for code = 0..4,
!< seth_weno6_point(u, 4 p + 2 sub - 3) against seam_fill_weno6_1d(u, sub, p) for p = 1..6, sub = 1, 2, and
!< seth_weno6_jumps against seam_weno6_jumps. Built by gfortran (debug or optimised) every result must be bitwise
!< equal; another compiler may contract multiply-adds differently in the two copies, and a result then passes within
!< TOL = 64 eps of the magnitude of its data (the count of bitwise results is printed either way). R8P only.
use seth
use seth_point_weno6_tables, only : WENO6_POINT_D, WENO6_POINT_W
use seth_test_rng, only : seed_rng
use seth_adam_seam_reference, only : SEAM_D_WENO6, SEAM_W_CUBIC_MID, SEAM_D_FILL_WENO6, SEAM_W_CUBIC_Q, &
                                     SEAM_WENO6_EPS, SEAM_WENO6_JUMP_RATIO, seam_weno6_jumps, &
                                     seam_restriction_weno6_1d, seam_fill_weno6_1d
implicit none
integer(I4P), parameter :: N = 100000_I4P !< Random lines.
integer(I4P), parameter :: NKIND = 7_I4P  !< Kinds of data.
real(R8P),    parameter :: TOL = 64._R8P * epsilon(1._R8P)
real(R8P)               :: u(1:6), got, want, diff, mag
integer(I4P)            :: n_, code, p, sb, t4, kind_, bitwise, total, nbad, ntab
logical                 :: jump(1:5), jump_ref(1:5), passed
#ifdef __GFORTRAN__
logical, parameter      :: BITWISE_REQUIRED = .true.
#else
logical, parameter      :: BITWISE_REQUIRED = .false.
#endif

passed = .true.
! the tables, the regulariser and the jump ratio
ntab = 0 ; nbad = 0
do code=0, 4
   t4 = 4_I4P * code + 6_I4P
   call table(WENO6_POINT_D(:,t4), SEAM_D_WENO6(:,code))
   do sb=1, 3
      if (SEAM_D_WENO6(sb,code) > 0._R8P) call table(WENO6_POINT_W(:,sb,t4), SEAM_W_CUBIC_MID(:,code-sb+2))
   enddo
enddo
do p=1, 6
   do sb=1, 2
      t4 = 4_I4P * p + 2_I4P * sb - 3_I4P
      call table(WENO6_POINT_D(:,t4), SEAM_D_FILL_WENO6(:,(t4-1)/2))
      do code=1, 3
         if (SEAM_D_FILL_WENO6(code,(t4-1)/2) > 0._R8P) &
            call table(WENO6_POINT_W(:,code,t4), SEAM_W_CUBIC_Q(:,(t4 - 4_I4P * code + 3_I4P) / 2_I4P))
      enddo
   enddo
enddo
call table([SETH_WENO6_EPS, SETH_WENO6_JUMP_RATIO], [SEAM_WENO6_EPS, SEAM_WENO6_JUMP_RATIO])
print '(A,I0,A,I0,A)', 'tables: ', ntab - nbad, ' of ', ntab, ' constants bitwise equal to ADAM''s'
if (nbad > 0) passed = .false.

! the results
call seed_rng(2026112_I4P)
bitwise = 0 ; total = 0 ; diff = 0._R8P ; nbad = 0
do n_=1, N
   kind_ = mod(n_, NKIND)
   call line(kind_, u)
   mag = max(maxval(abs(u)), tiny(1._R8P))
   do code=0, 4
      got  = seth_weno6_point(u, 4_I4P * code + 6_I4P)
      want = seam_restriction_weno6_1d(u=u, code=code)
      call compare(got, want, mag)
   enddo
   do p=1, 6
      do sb=1, 2
         got  = seth_weno6_point(u, 4_I4P * p + 2_I4P * sb - 3_I4P)
         want = seam_fill_weno6_1d(u=u, sub=sb, p=p)
         call compare(got, want, mag)
      enddo
   enddo
   call seth_weno6_jumps(u, jump)
   call seam_weno6_jumps(u=u, jump=jump_ref)
   if (any(jump .neqv. jump_ref)) nbad = nbad + 1
enddo
print '(A,I0,A,I0,A,ES9.2)', 'results: ', bitwise, ' of ', total, ' bitwise equal to ADAM''s; max scaled difference ', &
      diff
print '(A,I0,A,I0)', 'jump flags: lines differing from ADAM''s ', nbad, ' of ', N
if (nbad > 0) passed = .false.
if (BITWISE_REQUIRED .and. bitwise /= total) passed = .false.
if (diff > TOL) passed = .false.
if (.not.passed) error stop 'seth_weno6_point differs from ADAM''s seam WENO-6'
print '(A)', 'seth_test_weno6_point_adam: ok'

contains
   subroutine table(seth_, adam_)
   !< Count the constants of SETH's table that are not bitwise ADAM's.
   real(R8P), intent(in) :: seth_(:), adam_(:)
   ntab = ntab + size(seth_)
   nbad = nbad + count(seth_ /= adam_)
   endsubroutine table

   subroutine compare(val, ref, mag)
   real(R8P), intent(in) :: val, ref !< SETH's and ADAM's result.
   real(R8P), intent(in) :: mag      !< The magnitude of the data.
   total = total + 1
   if (val == ref) then
      bitwise = bitwise + 1
   else
      diff = max(diff, abs(val - ref) / mag)
   endif
   endsubroutine compare

   subroutine line(kind_, u)
   !< Six random values of one kind, at a random magnitude 10^(-12..12).
   integer(I4P), intent(in)  :: kind_
   real(R8P),    intent(out) :: u(1:6)
   real(R8P)                 :: r(1:6), x, scale, h, c(1:4)
   integer(I4P)              :: k, jmp
   call random_number(x) ; scale = 10._R8P**(24._R8P * x - 12._R8P)
   select case(kind_)
   case(0) ! smooth: samples of a cubic plus a wave, spacing 1e-3..1
      call random_number(c) ; call random_number(x) ; h = 10._R8P**(-3._R8P * x)
      do k=1, 6
         u(k) = scale * (c(1) - 0.5_R8P + c(2) * (h * k) + c(3) * (h * k)**2 + sin(6._R8P * c(4) * h * k + c(1)))
      enddo
   case(1) ! independent random values
      call random_number(r) ; u = scale * (r - 0.5_R8P)
   case(2) ! a step at a random interval (or none), two random levels, a little noise
      call random_number(r) ; call random_number(x) ; jmp = int(7._R8P * x, I4P)
      do k=1, 6
         u(k) = scale * (merge(r(1), r(2) - 0.5_R8P, k > jmp) + 1.e-3_R8P * r(k))
      enddo
   case(3) ! constant
      call random_number(x) ; u = scale * (x - 0.5_R8P)
   case(4) ! a flat sub-stencil next to a ramp (beta = 0 for it), either way round
      call random_number(x)
      u = scale * [1._R8P, 1._R8P, 1._R8P, 1._R8P, 1._R8P + x * 1.e-3_R8P, 1._R8P + x * 2.e-3_R8P]
      call random_number(x)
      if (x < 0.5_R8P) u = u(6:1:-1)
   case(5) ! ADAM's sharpened, jump-like data
      call random_number(r) ; u = scale * r**8
   case default ! small variations on a large offset
      call random_number(r) ; u = scale * (1._R8P + 1.e-6_R8P * r)
   endselect
   endsubroutine line
endprogram seth_test_weno6_point_adam
