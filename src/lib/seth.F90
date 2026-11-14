!< SETH, Smoothness-weighted ENO, TENO and Hybrids: the library's entry point.
module seth
!< SETH, Smoothness-weighted ENO, TENO and Hybrids: the library's entry point.
!<
!< Seed of the repository: it re-exports the working kinds, so a client needs `use seth` only. The reconstruction and
!< interpolation machinery (extracted from ADAM) is added behind this module.
use penf, only : I4P, I8P, R4P, R8P
implicit none
private
public :: I4P, I8P, R4P, R8P
endmodule seth
