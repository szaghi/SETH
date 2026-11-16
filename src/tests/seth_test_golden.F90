!< SETH test: regression against stored golden values of every kernel on deterministic stencils.
program seth_test_golden
!< SETH test: regression against stored golden values of every kernel on deterministic stencils.
!<
!< The verification tests assert orders, properties and bounds: a change of a few ulp, or a 30% change of an error
!< that stays inside a bound, passes them. This test pins the values themselves. Every scheme of the registry (WENO-Z
!< also with p = 2) is evaluated on NST deterministic stencils (smooth polynomials, jumps at every cell boundary,
!< data scaled by powers of 2, rough data from an integer LCG): exact arithmetic only (no libm, no random-number
!< generator), so the inputs are bit-identical on every compiler and any difference is the kernels'. The two states
!< (and the weights ratios of WENO-JS) are compared with the golden file:
!<
!<   seth_test_golden                    compare with src/tests/golden/seth_golden_v6.txt (run from the repository root)
!<   seth_test_golden write <file>       write a new golden file (refuses to overwrite an existing one)
!<
!< A state passes when it differs from the golden one by at most 1e-13 of the stencil's magnitude: other compilers
!< and flags contract multiply-adds differently from the build that wrote the file (its compiler and options are in
!< the file's header). A weights ratio passes within 1e-9: it reads the weights themselves, which move by about
!< wexp dIS/(eps + IS) when the indicators round differently (dIS: the round-off of IS, large against a small IS with
!< the Balsara-Shu coefficients of WENO-7/9), while the states, a combination of nearly equal candidates on smooth
!< data, barely move (measured: gfortran -O3 2.6e-13, nvfortran -fast 3.5e-10). The numbers of bitwise-equal values
!< and the worst differences are printed. A golden file is never edited: a deliberate change of results gets a new
!< file and this test's default path. v2 (2026-11-15) adds the WENO-Z and TENO interpolations; its lines for the schemes
!< of v1 are identical to v1's, character for character (checked when v2 was written, same build), and v1 is kept.
!< v3 (2026-11-15) adds TENO5-opt, v4 TENO5-AZ, v5 TENO8, v6 TENO5-A and TENO8-A; each keeps the lines of the previous
!< file identical (same build), and the previous files are kept.
use, intrinsic :: iso_fortran_env, only : compiler_version, compiler_options
use seth
use seth_test_schemes
implicit none
character(*), parameter :: GOLDEN = 'src/tests/golden/seth_golden_v6.txt'
integer(I4P), parameter :: NST = 64_I4P
real(R8P),    parameter :: TOL = 1.e-13_R8P, TOL_WR = 1.e-9_R8P
real(R8P)               :: vp(VLO:VHI,NST), vm(VLO:VHI,NST), r(4), g(4), mag, d, worst, worst_wr
integer(I4P)            :: id, p, k, u, n_values, n_bitwise, n_fail, ios, kk
character(16)           :: name, gname
character(256)          :: arg, path, line
logical                 :: write_mode, exists

call stencils(vp, vm)
write_mode = .false.
if (command_argument_count() >= 1) then
   call get_command_argument(1, arg)
   write_mode = trim(arg) == 'write'
endif
if (write_mode) then
   call get_command_argument(2, path)
   inquire(file=trim(path), exist=exists)
   if (exists) error stop 'seth_test_golden: refusing to overwrite an existing golden file'
   open(newunit=u, file=trim(path), status='new', action='write')
   write(u, '(A)') '# SETH golden values: scheme p stencil rp rm wrp wrm (wr = -1: no weights ratio)'
   write(u, '(A)') '# compiler: '//compiler_version()
   write(u, '(A)') '# options: '//compiler_options()
   do id=1, NSCHEMES
      do p=1, merge(2, 1, SCHEMES(id)%family == FAM_Z)
         do k=1, NST
            call evaluate(id, p, k, r)
            write(u, '(A16,1X,I1,1X,I3,4(1X,ES25.17E3))') SCHEMES(id)%name, p, k, r
         enddo
      enddo
   enddo
   close(u)
   print '(A)', 'golden values written to '//trim(path)
   stop
endif
open(newunit=u, file=GOLDEN, status='old', action='read', iostat=ios)
if (ios /= 0) error stop 'seth_test_golden: cannot open '//GOLDEN//' (run from the repository root)'
n_values = 0 ; n_bitwise = 0 ; n_fail = 0 ; worst = 0._R8P ; worst_wr = 0._R8P
do
   read(u, '(A)', iostat=ios) line
   if (ios /= 0) exit
   if (line(1:1) == '#') cycle
   read(line, *) gname, p, k, g
   id = scheme_index(gname)
   if (id < 1) then ; print '(A)', 'unknown scheme in the golden file: '//trim(gname) ; n_fail = n_fail + 1 ; cycle ; endif
   call evaluate(id, p, k, r)
   mag = max(maxval(abs(vp(:,k))), maxval(abs(vm(:,k))), tiny(1._R8P))
   do kk=1, 4
      n_values = n_values + 1
      if (r(kk) == g(kk)) then
         n_bitwise = n_bitwise + 1
      else
         if (kk <= 2) then
            d = abs(r(kk) - g(kk))/mag ; worst = max(worst, d)
         else
            d = abs(r(kk) - g(kk)) ; worst_wr = max(worst_wr, d)
         endif
         if (.not.(d <= merge(TOL, TOL_WR, kk <= 2))) then
            n_fail = n_fail + 1
            if (n_fail <= 10) print '(A16,I2,I4,I2,2ES25.16)', gname, p, k, kk, r(kk), g(kk)
         endif
      endif
   enddo
enddo
close(u)
print '(A,I0,A,I0,A,ES9.2,A,ES9.2)', 'values: ', n_values, ', bitwise equal: ', n_bitwise, &
      ', worst state difference (scaled): ', worst, ', worst weights-ratio difference: ', worst_wr
if (n_values /= expected_values()) then
   print '(A,I0)', 'the golden file does not cover every scheme and stencil; expected values: ', expected_values()
   n_fail = n_fail + 1
endif
if (n_fail > 0) error stop 'golden regression failed'
print '(A)', 'seth_test_golden: ok'

contains
   subroutine evaluate(id, p, k, r)
   !< rp, rm, wrp, wrm of scheme id (WENO-Z: exponent p) on stencil k.
   integer(I4P), intent(in)  :: id, p, k
   real(R8P),    intent(out) :: r(4)
   if (SCHEMES(id)%family == FAM_Z) then
      call reconstruct(id, vp(:,k), vm(:,k), r(1), r(2), r(3), r(4), par=real(p, R8P))
   else
      call reconstruct(id, vp(:,k), vm(:,k), r(1), r(2), r(3), r(4))
   endif
   endsubroutine evaluate

   pure function scheme_index(n) result(id)
   character(*), intent(in) :: n
   integer(I4P)             :: id
   do id=1, NSCHEMES
      if (trim(SCHEMES(id)%name) == trim(n)) return
   enddo
   id = -1
   endfunction scheme_index

   pure function expected_values() result(n)
   integer(I4P) :: n, id
   n = 0
   do id=1, NSCHEMES
      n = n + 4*NST*merge(2, 1, SCHEMES(id)%family == FAM_Z)
   enddo
   endfunction expected_values

   subroutine stencils(vp, vm)
   !< NST deterministic stencil pairs, exact in binary: 16 smooth cubics with dyadic coefficients, 20 jumps (each cell
   !< boundary, two shapes; VLO: constants), 12 scaled (a cubic with a jump times 2^-498 .. 2^96, about 1e-150 .. 1e29),
   !< 16 rough (an integer LCG).
   real(R8P), intent(out) :: vp(VLO:VHI,NST), vm(VLO:VHI,NST)
   integer(I4P)           :: k, m, j
   integer(I8P)           :: state
   real(R8P)              :: sc
   k = 0
   do j=1, 16 ! smooth: cubics with dyadic coefficients (exact)
      k = k + 1
      do m=VLO, VHI
         vp(m,k) = 1._R8P + (j/8._R8P)*m - (j/64._R8P)*m*m + (1._R8P/512._R8P)*m**3
         vm(m,k) = -0.5_R8P + (j/16._R8P)*m + (j/128._R8P)*m*m - (j/1024._R8P)*m**3
      enddo
   enddo
   do j=VLO, VHI ! jumps between cells j-1 and j, both orientations (j = VLO: the constant stencils)
      k = k + 1
      do m=VLO, VHI
         vp(m,k) = merge(1._R8P, 0._R8P, m >= j)
         vm(m,k) = merge(0._R8P, 2._R8P, m >= j)
      enddo
      k = k + 1
      do m=VLO, VHI
         vp(m,k) = merge(-3._R8P, 0.5_R8P, m >= j)
         vm(m,k) = merge(1._R8P, 1.25_R8P, m >= j)
      enddo
   enddo
   do j=1, 12 ! scaled: a cubic with a jump, times 2^-498 .. 2^96 (exact)
      k = k + 1
      sc = 2._R8P**(-498 + (j-1)*54)
      do m=VLO, VHI
         vp(m,k) = sc*(0.25_R8P*m - 0.03125_R8P*m*m + merge(1._R8P, 0._R8P, m >= 1))
         vm(m,k) = sc*(1._R8P - 0.125_R8P*m + 0.0078125_R8P*m**3 - merge(0.5_R8P, 0._R8P, m >= 0))
      enddo
   enddo
   state = 20261115_I8P
   do j=1, 16 ! rough: a 31-bit LCG (Park-Miller), the same numbers on every compiler
      k = k + 1
      do m=VLO, VHI
         state = mod(16807_I8P*state, 2147483647_I8P) ; vp(m,k) = real(state, R8P)/2147483647._R8P - 0.5_R8P
         state = mod(16807_I8P*state, 2147483647_I8P) ; vm(m,k) = real(state, R8P)/2147483647._R8P - 0.3_R8P
      enddo
   enddo
   if (k /= NST) error stop 'stencils: count mismatch'
   endsubroutine stencils
endprogram seth_test_golden
