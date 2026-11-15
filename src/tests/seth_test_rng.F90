!< SETH tests: a fixed, printed seed for every test that draws random stencils, so that a failure reproduces.
module seth_test_rng
!< SETH tests: a fixed, printed seed for every test that draws random stencils, so that a failure reproduces.
!<
!< `call seed_rng(base)` puts the seed base, base+1, ... into random_seed and prints `base`. The default base of each
!< test is fixed in its source; the environment variable SETH_TEST_SEED, when set, overrides it (to explore other
!< streams without editing the tests).
use penf, only : I4P
implicit none
private
public :: seed_rng

contains
   subroutine seed_rng(base)
   !< Seed the random number generator with base, base+1, ... (or SETH_TEST_SEED) and print the base used.
   integer(I4P), intent(in) :: base !< Default base of the seed.
   integer, allocatable     :: seed(:)
   integer                  :: n, i, b, status
   character(32)            :: env
   b = base
   call get_environment_variable('SETH_TEST_SEED', env, status=status)
   if (status == 0 .and. len_trim(env) > 0) read(env, *) b
   call random_seed(size=n)
   allocate(seed(n))
   seed = [(b + i, i=0, n-1)]
   call random_seed(put=seed)
   print '(A,I0)', 'random seed base: ', b
   endsubroutine seed_rng
endmodule seth_test_rng
