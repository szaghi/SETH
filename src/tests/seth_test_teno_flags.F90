!< SETH tests: a twin of seth_teno5 that also returns the masks of its two states (the reference of seth_teno5_mask).
module seth_test_teno_flags
!< SETH tests: a twin of seth_teno5 that also returns the masks of its two states (the reference of seth_teno5_mask).
!<
!< The library's body (seth_z_teno_body.INC) compiled with SETH_TENO_FLAGS: the kernel's code, plus the outputs mp, mm =
!< sum_k delta_k 2**k of the left and right states. seth_test_nonpoly checks that its states are seth_teno5's bitwise
!< (so it is the kernel) and that seth_teno5_mask returns its masks (so the flags are the kernel's cutoff).
use penf, only : I4P, R_P
use seth_z_teno_tables
implicit none
private
public :: teno5_flags

contains
   pure subroutine teno5_flags(vp, vm, eps, ct, rp, rm, mp, mm)
   !< seth_teno5 with the masks of its two states.
   integer(I4P), parameter :: NC = TENO5_NC, LO = TENO5_LO, UMAX = TENO5_UMAX, MAXW = TENO5_MAXW
   integer(I4P), parameter :: CLO(0:NC-1) = TENO5_CLO, CW(0:NC-1) = TENO5_CW
   integer(I4P), parameter :: GLO = TENO5_GLO, GW = TENO5_GW
   real(R_P),    parameter :: CP(0:MAXW-1,0:NC-1) = TENO5_P
   real(R_P),    parameter :: CL(0:MAXW-1,0:MAXW-1,0:NC-1) = TENO5_L, CDQ(0:MAXW-1,0:NC-1) = TENO5_DQ
   real(R_P),    parameter :: D(0:NC-1) = TENO5_D, TC(0:NC-1) = TENO5_TC, TG = TENO5_TG
   real(R_P),    parameter :: CLG(0:GW-1,0:GW-1) = TENO5_LG, CDG(0:GW-1) = TENO5_DG
#define SETH_TENO
#define SETH_TENO_FLAGS
#include "../lib/seth_z_teno_body.INC"
#undef SETH_TENO_FLAGS
#undef SETH_TENO
   endsubroutine teno5_flags
endmodule seth_test_teno_flags
