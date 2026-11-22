!< SETH test: the WENO-6 point interpolation, seth_weno6_point and seth_weno6_jumps, against its definition.
program seth_test_weno6_point
!< SETH test: the WENO-6 point interpolation, seth_weno6_point and seth_weno6_jumps, against its definition.
!<
!< In the working kind R_P (issue #7, E2; ADAM's test_seam_interpolation checks the same properties of its copy):
!< 1. the compiled tables: the linear weights of every target sum to 1, the cubic weights of every candidate too, and
!<    the linear combination of the candidates' cubics reproduces x^n up to the size of its union minus 1;
!< 2. a node (t4 = 0 mod 4) returns its value, bitwise;
!< 3. the kernel is exact on a cubic at every target, under amplitudes and offsets of 1e-3..1e3;
!< 4. design order on smooth data, h = 0.05 and 0.025: 6 at x = 13/4..15/4, 5 one interval further out, 4 toward the
!<    ends (R8P; in R4P the errors of the finer grid are below its round-off and no order is observable);
!< 5. the jump detector flags the interval of a step and nothing else, and nothing on resolved smooth data;
!< 6. steps: a target takes the state of its side of a jump outside its interval; a jump inside its interval gives the
!<    troubled-cell fallback, the mean of the two values at a midpoint (bitwise), the clamped linear interpolation at a
!<    quarter;
!< 7. mirror symmetry on random data: t4 on u equals 28 - t4 on u reversed, to round-off;
!< 8. the regulariser in R4P: a flat sub-stencil next to a ramp (1, 1, 1, 1, 1.001, 1.002 at magnitudes 1e-3..1e3,
!<    either way round) gives a finite value, within the ramp's extent of the data, at every target (in R4P
!<    (tau/eps)^2 alone overflows); constant data return the constant;
!< 9. on smooth data (spacing 0.01, where the nonlinear weights are the linear ones to far below round-off) the value
!<    agrees with the R8P Lagrange interpolant of the candidates' union to the round-off of R_P: in R4P, R4P against
!<    R8P.
use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
use seth
use seth_point_weno6_tables, only : WENO6_POINT_D, WENO6_POINT_W
use seth_test_rng, only : seed_rng
implicit none
real(R_P), parameter :: EPS = epsilon(1._R_P)
real(R_P)            :: u(1:6), x(1:6), got, want, err(1:2), h, x0, amp, off, r(1:6), lo_, hi_, worst
real(R8P)            :: u8(1:6), ref8, l8
integer(I4P)         :: t4, k, l, n_, m, s, jmp, i, s_lo, s_hi
logical              :: jump(1:5), ok, passed

passed = .true.
call seed_rng(2026113_I4P)

! 1. the compiled tables
ok = .true.
do t4=3, 25
   if (mod(t4, 4) == 0) cycle
   if (abs(sum(WENO6_POINT_D(:,t4)) - 1._R_P) > 2._R_P * EPS) ok = .false.
   if (any(WENO6_POINT_D(:,t4) < 0._R_P)) ok = .false.
   call candidates(t4, s_lo, s_hi)
   do s=s_lo, s_hi
      if (abs(sum(WENO6_POINT_W(:,s,t4)) - 1._R_P) > 4._R_P * EPS) ok = .false.
   enddo
   do n_=0, s_hi - s_lo + 3 ! the union's size minus 1
      got = 0._R_P ; amp = 0._R_P
      do s=s_lo, s_hi
         do k=1, 4
            off = WENO6_POINT_D(s,t4) * WENO6_POINT_W(k,s,t4) * (real(s + k - 1, R_P) - 0.25_R_P * t4)**n_
            got = got + off ; amp = amp + abs(off)
         enddo
      enddo
      want = merge(1._R_P, 0._R_P, n_ == 0) ! x^n at the target, the origin
      if (abs(got - want) > 16._R_P * EPS * amp) ok = .false.
   enddo
enddo
call verdict(ok, 'tables: linear and cubic weights sum to 1, the combination reproduces the union''s polynomials')

! 2. a node returns its value
ok = .true.
call random_number(u)
do k=1, 6
   if (seth_weno6_point(u, 4_I4P * k) /= u(k)) ok = .false.
enddo
call verdict(ok, 'a node returns its value, bitwise')

! 3. exact on a cubic, f = x^3 + 6 x + 1 (no jump), every target, scale and offset invariant
ok = .true. ; worst = 0._R_P
do t4=3, 25
   if (mod(t4, 4) == 0) cycle
   do m=1, 6
      x(m) = real(m, R_P) - 0.25_R_P * t4 ! node coordinates, the target at 0
   enddo
   do l=-3, 3
      amp = 10._R_P**l ; off = 10._R_P**(-l)
      u = amp * (x**3 + 6._R_P * x + 1._R_P) + off
      got = seth_weno6_point(u, t4)
      want = amp + off
      worst = max(worst, abs(got - want) / (amp * maxval(abs(x**3 + 6._R_P * x + 1._R_P)) + off))
   enddo
enddo
if (worst > 64._R_P * EPS) ok = .false.
call verdict(ok, 'exact on cubics at every target, scale and offset invariant', worst)

! 4. design order on smooth data
#ifdef _R_P_IS_R4P
print '(A)', 'SKIP: design order (R4P: the errors at h = 0.025, 1e-9..1e-12, are below its round-off)'
#else
ok = .true.
do t4=3, 25
   if (mod(t4, 4) == 0) cycle
   do l=1, 2
      h = 0.05_R_P / real(l, R_P)
      err(l) = 0._R_P
      do n_=0, 100
         x0 = real(n_, R_P) / 100._R_P
         do m=1, 6
            x(m) = x0 + h * (real(m, R_P) - 0.25_R_P * t4)
         enddo
         u = sin(2._R_P * acos(-1._R_P) * x) + 0.3_R_P * cos(6._R_P * x)
         got = seth_weno6_point(u, t4)
         err(l) = max(err(l), abs(got - (sin(2._R_P * acos(-1._R_P) * x0) + 0.3_R_P * cos(6._R_P * x0))))
      enddo
   enddo
   want = 3.7_R_P ! the design order minus 0.3
   if (t4 >= 9 .and. t4 <= 19) want = 4.7_R_P
   if (t4 >= 13 .and. t4 <= 15) want = 5.7_R_P
   if (log(err(1) / err(2)) / log(2._R_P) < want) then
      ok = .false.
      print '(A,I3,A,F6.2)', '   t4 =', t4, ': order ', log(err(1) / err(2)) / log(2._R_P)
   endif
enddo
call verdict(ok, 'design order on smooth data: 6 at x = 13/4..15/4, 5 further out, 4 toward the ends')
#endif

! 5. the jump detector: a step flags its interval alone; resolved smooth data flag nothing
ok = .true.
do jmp=1, 5
   do m=1, 6
      u(m) = merge(1._R_P, 1.e-3_R_P, m > jmp)
   enddo
   call seth_weno6_jumps(u, jump)
   do i=1, 5
      if (jump(i) .neqv. (i == jmp)) ok = .false.
   enddo
enddo
do n_=0, 100
   x0 = real(n_, R_P) / 100._R_P
   do m=1, 6
      u(m) = sin(2._R_P * acos(-1._R_P) * (x0 + 0.05_R_P * m))
   enddo
   call seth_weno6_jumps(u, jump)
   if (any(jump)) ok = .false.
enddo
call verdict(ok, 'jump detector: the interval of a step alone, nothing on resolved smooth data')

! 6. steps (1e-3 / 1)
ok = .true.
do t4=3, 25
   if (mod(t4, 4) == 0) cycle
   m = min(max(t4 / 4, 1), 5) ! the target's interval
   do jmp=1, 5
      do i=1, 6
         u(i) = merge(1._R_P, 1.e-3_R_P, i > jmp)
      enddo
      got = seth_weno6_point(u, t4)
      if (jmp == m) then
         if (mod(t4, 4) == 2) then
            if (got /= 0.5_R_P * (u(m) + u(m+1))) ok = .false.
         else
            want = min(max(u(m) + 0.25_R_P * real(t4 - 4 * m, R_P) * (u(m+1) - u(m)), 1.e-3_R_P), 1._R_P)
            if (got /= want) ok = .false.
         endif
      elseif (jmp < m) then
         if (abs(got - 1._R_P) > 8._R_P * EPS) ok = .false.
      else
         if (abs(got - 1.e-3_R_P) > 8._R_P * EPS) ok = .false.
      endif
   enddo
enddo
call verdict(ok, 'steps: the side of the target, the troubled-cell fallback for a jump in its interval')

! 7. mirror symmetry on random data (rnd**8: jump-like)
ok = .true. ; worst = 0._R_P
do n_=1, 2000
   call random_number(r)
   if (mod(n_, 2) == 0) r = r**8
   do t4=3, 25
      got = seth_weno6_point(r, t4)
      u = r(6:1:-1)
      want = seth_weno6_point(u, 28_I4P - t4)
      worst = max(worst, abs(got - want) / maxval(abs(r)))
   enddo
enddo
if (worst > 64._R_P * EPS) ok = .false.
call verdict(ok, 'mirror symmetry: t4 on u equals 28 - t4 on u reversed', worst)

! 8. a flat sub-stencil next to a ramp, and constants
ok = .true.
do l=-3, 3
   do i=1, 2
      u = 10._R_P**l * [1._R_P, 1._R_P, 1._R_P, 1._R_P, 1.001_R_P, 1.002_R_P]
      if (i == 2) u = u(6:1:-1)
      lo_ = minval(u) ; hi_ = maxval(u)
      do t4=3, 25 ! finite, and within the ramp's extent of the data (the end targets extrapolate a quarter)
         got = seth_weno6_point(u, t4)
         if (.not.ieee_is_finite(got)) ok = .false.
         if (.not.(got >= 2._R_P * lo_ - hi_ .and. got <= 2._R_P * hi_ - lo_)) ok = .false.
      enddo
   enddo
   u = -0.7_R_P * 10._R_P**l
   do t4=3, 25
      if (abs(seth_weno6_point(u, t4) - u(1)) > 8._R_P * EPS * abs(u(1))) ok = .false.
   enddo
enddo
call verdict(ok, 'flat sub-stencil next to a ramp: finite, bounded; constants returned')

! 9. smooth data against the R8P Lagrange interpolant of the candidates' union
ok = .true. ; worst = 0._R_P
do n_=1, 2000
   call random_number(r)
   do t4=3, 25
      if (mod(t4, 4) == 0) cycle
      do m=1, 6
         u(m) = r(1) + sin(2._R_P * (r(2) + 0.01_R_P * m)) + r(3) * (0.01_R_P * m)**2
      enddo
      call candidates(t4, s_lo, s_hi)
      u8 = real(u, R8P)
      ref8 = 0._R8P
      do k=s_lo, s_hi + 3
         l8 = 1._R8P
         do m=s_lo, s_hi + 3
            if (m /= k) l8 = l8 * (0.25_R8P * t4 - m) / real(k - m, R8P)
         enddo
         ref8 = ref8 + l8 * u8(k)
      enddo
      got = seth_weno6_point(u, t4)
      worst = max(worst, real(abs(real(got, R8P) - ref8) / maxval(abs(u8)), R_P))
   enddo
enddo
if (worst > 16._R_P * EPS) ok = .false.
call verdict(ok, 'smooth data: the R8P Lagrange interpolant of the candidates'' union, to the round-off of R_P', worst)

if (.not.passed) error stop 'seth_weno6_point violates its definition'
print '(A)', 'seth_test_weno6_point: ok'

contains
   subroutine candidates(t4, s_lo, s_hi)
   !< The first and last sub-stencil with a positive linear weight at t4.
   integer(I4P), intent(in)  :: t4
   integer(I4P), intent(out) :: s_lo, s_hi
   integer(I4P)              :: s
   s_lo = 0 ; s_hi = 0
   do s=1, 3
      if (WENO6_POINT_D(s,t4) > 0._R_P) then
         if (s_lo == 0) s_lo = s
         s_hi = s
      endif
   enddo
   endsubroutine candidates

   subroutine verdict(ok, what, measured)
   logical,      intent(in)           :: ok
   character(*), intent(in)           :: what
   real(R_P),    intent(in), optional :: measured
   if (present(measured)) then
      print '(A,A,A,ES9.2,A)', merge('PASS: ', 'FAIL: ', ok), what, ' (worst ', measured / EPS, ' eps)'
   else
      print '(A,A)', merge('PASS: ', 'FAIL: ', ok), what
   endif
   if (.not.ok) passed = .false.
   endsubroutine verdict
endprogram seth_test_weno6_point
