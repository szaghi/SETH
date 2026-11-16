!< SETH test: seth_keep_pe6 and seth_ducros against ADAM's compute_kep_flux and compute_ducros (R8P only).
program seth_test_central_adam
!< SETH test: seth_keep_pe6 and seth_ducros against ADAM's compute_kep_flux and compute_ducros (R8P only).
!<
!< On 20000 random states: seth_keep_pe6 equals ADAM's compute_kep_flux bitwise (gfortran; other compilers within 2e-15
!< of the flux scale, fused multiply-adds), and seth_ducros equals compute_ducros with its override off.
!<
!< ADAM's code (seth_adam_reference) is R8P: the comparison is R8P only, and the R4P test builds leave this program out
!< (the `exclude` of fobos' [template-r4-tests]; scripts/check_literals.py skips it for the same reason). It was item 5
!< of seth_test_central until the tests were made per kind (roadmap issue #4, 6.2, phase 3).
use seth
use seth_test_rng, only : seed_rng
use seth_adam_reference, only : compute_kep_flux, compute_ducros
implicit none
real(R8P), parameter :: GAM = 1.4_R8P
real(R8P)            :: w(5,-2:3), f(5), fa(5), phi, psi, g(3,3), worst
integer(I4P)         :: n, i, d, nbit, ntot
logical              :: passed

passed = .true.
call seed_rng(1013_I4P)
nbit = 0 ; ntot = 0 ; worst = 0._R8P
do n=1, 20000
   do i=-2, 3
      call random_number(w(:,i))
      w(1,i) = 0.1_R8P + w(1,i) ; w(5,i) = 0.1_R8P + w(5,i) ; w(2:4,i) = 2._R8P*w(2:4,i) - 1._R8P
   enddo
   d = 1 + mod(n, 3)
   call seth_keep_pe6(d, GAM, w, f)
   call compute_kep_flux(d, GAM, w, fa)
   ntot = ntot + 5 ; nbit = nbit + count(f == fa)
   worst = max(worst, maxval(abs(f - fa))/maxval(abs(fa)))
   call random_number(g) ; g = 2._R8P*g - 1._R8P
   call seth_ducros(1.e-3_R8P, g, phi)
   call compute_ducros(1.e-3_R8P, 0._R8P, 1._R8P, g, psi)
   if (phi /= psi) then ; print '(A)', 'seth_ducros differs from ADAM''s compute_ducros' ; passed = .false. ; exit ; endif
enddo
print '(A,I0,A,I0,A,ES10.2)', 'seth_keep_pe6 against ADAM''s compute_kep_flux: ', nbit, ' of ', ntot, &
      ' bitwise, worst scaled difference ', worst
if (worst > 2.e-15_R8P) passed = .false.
if (.not.passed) error stop 'central flux test against ADAM failed'
print '(A)', 'seth_test_central_adam: ok'
endprogram seth_test_central_adam
