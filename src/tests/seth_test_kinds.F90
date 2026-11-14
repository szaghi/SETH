!< SETH test: the working kinds reach a client through `use seth` (the PENF dependency is wired).
program seth_test_kinds
!< SETH test: the working kinds reach a client through `use seth` (the PENF dependency is wired).
use seth, only : I4P, I8P, R4P, R8P
implicit none

if (storage_size(1_I4P) /= 32)  error stop 'I4P is not a 32-bit integer'
if (storage_size(1_I8P) /= 64)  error stop 'I8P is not a 64-bit integer'
if (storage_size(1._R4P) /= 32) error stop 'R4P is not a 32-bit real'
if (storage_size(1._R8P) /= 64) error stop 'R8P is not a 64-bit real'
print '(A)', 'seth_test_kinds: ok'
endprogram seth_test_kinds
