!< SETH, WENO-Z kernels: WENO-Z3/5/7/9 reconstruction at one face.
module seth_wenoz
!< SETH, WENO-Z kernels: WENO-Z3/5/7/9 reconstruction at one face.
!<
!< One kernel per scheme, its coefficients compile-time constants (seth_z_teno_tables), the body shared with the other
!< family (seth_z_teno_body.INC). Every kernel is pure, `acc routine seq` and `omp declare target`:
!<
!<    call seth_wenoz3(vp, vm, eps, p, rp, rm)
!<
!< `vp`/`vm` the values reconstructed with a left/right bias, cells by offset from the face (cell 0 left of it, bounds
!< LO:1-LO of the scheme), `eps` the regulariser, `p` the exponent of tau/beta in the weights (1 or 2),
!< `rp`/`rm` the left/right states at the face. Sources and derivations: notes/m2-sources.md.
use penf, only : I4P, R8P
use seth_z_teno_tables
implicit none
private
public :: seth_wenoz3
public :: seth_wenoz5
public :: seth_wenoz7
public :: seth_wenoz9

contains
   pure subroutine seth_wenoz3(vp, vm, eps, p, rp, rm)
   !< WENO-Z3 (Borges et al. weights), tau_3 = |beta_0 - beta_1|.
   integer(I4P), parameter :: NC = WENOZ3_NC, LO = WENOZ3_LO, UMAX = WENOZ3_UMAX, MAXW = WENOZ3_MAXW
   integer(I4P), parameter :: CLO(0:NC-1) = WENOZ3_CLO, CW(0:NC-1) = WENOZ3_CW
   integer(I4P), parameter :: GLO = WENOZ3_GLO, GW = WENOZ3_GW
   real(R8P),    parameter :: CP(0:MAXW-1,0:NC-1) = WENOZ3_P
   real(R8P),    parameter :: CL(0:MAXW-1,0:MAXW-1,0:NC-1) = WENOZ3_L, CDQ(0:MAXW-1,0:NC-1) = WENOZ3_DQ
   real(R8P),    parameter :: D(0:NC-1) = WENOZ3_D, TC(0:NC-1) = WENOZ3_TC, TG = WENOZ3_TG
   real(R8P),    parameter :: CLG(0:GW-1,0:GW-1) = WENOZ3_LG, CDG(0:GW-1) = WENOZ3_DG
#include "seth_z_teno_body.INC"
   endsubroutine seth_wenoz3

   pure subroutine seth_wenoz5(vp, vm, eps, p, rp, rm)
   !< WENO-Z5 (Borges et al. 2008), tau_5 = |beta_0 - beta_2|.
   integer(I4P), parameter :: NC = WENOZ5_NC, LO = WENOZ5_LO, UMAX = WENOZ5_UMAX, MAXW = WENOZ5_MAXW
   integer(I4P), parameter :: CLO(0:NC-1) = WENOZ5_CLO, CW(0:NC-1) = WENOZ5_CW
   integer(I4P), parameter :: GLO = WENOZ5_GLO, GW = WENOZ5_GW
   real(R8P),    parameter :: CP(0:MAXW-1,0:NC-1) = WENOZ5_P
   real(R8P),    parameter :: CL(0:MAXW-1,0:MAXW-1,0:NC-1) = WENOZ5_L, CDQ(0:MAXW-1,0:NC-1) = WENOZ5_DQ
   real(R8P),    parameter :: D(0:NC-1) = WENOZ5_D, TC(0:NC-1) = WENOZ5_TC, TG = WENOZ5_TG
   real(R8P),    parameter :: CLG(0:GW-1,0:GW-1) = WENOZ5_LG, CDG(0:GW-1) = WENOZ5_DG
#include "seth_z_teno_body.INC"
   endsubroutine seth_wenoz5

   pure subroutine seth_wenoz7(vp, vm, eps, p, rp, rm)
   !< WENO-Z7, tau_7 = |beta_0 + 3 beta_1 - 3 beta_2 - beta_3| (Don & Borges).
   integer(I4P), parameter :: NC = WENOZ7_NC, LO = WENOZ7_LO, UMAX = WENOZ7_UMAX, MAXW = WENOZ7_MAXW
   integer(I4P), parameter :: CLO(0:NC-1) = WENOZ7_CLO, CW(0:NC-1) = WENOZ7_CW
   integer(I4P), parameter :: GLO = WENOZ7_GLO, GW = WENOZ7_GW
   real(R8P),    parameter :: CP(0:MAXW-1,0:NC-1) = WENOZ7_P
   real(R8P),    parameter :: CL(0:MAXW-1,0:MAXW-1,0:NC-1) = WENOZ7_L, CDQ(0:MAXW-1,0:NC-1) = WENOZ7_DQ
   real(R8P),    parameter :: D(0:NC-1) = WENOZ7_D, TC(0:NC-1) = WENOZ7_TC, TG = WENOZ7_TG
   real(R8P),    parameter :: CLG(0:GW-1,0:GW-1) = WENOZ7_LG, CDG(0:GW-1) = WENOZ7_DG
#include "seth_z_teno_body.INC"
   endsubroutine seth_wenoz7

   pure subroutine seth_wenoz9(vp, vm, eps, p, rp, rm)
   !< WENO-Z9, tau_9 = |beta_0 - beta_4| (Castro et al.'s rule for odd r, as restated by Li et al. 2022b).
   integer(I4P), parameter :: NC = WENOZ9_NC, LO = WENOZ9_LO, UMAX = WENOZ9_UMAX, MAXW = WENOZ9_MAXW
   integer(I4P), parameter :: CLO(0:NC-1) = WENOZ9_CLO, CW(0:NC-1) = WENOZ9_CW
   integer(I4P), parameter :: GLO = WENOZ9_GLO, GW = WENOZ9_GW
   real(R8P),    parameter :: CP(0:MAXW-1,0:NC-1) = WENOZ9_P
   real(R8P),    parameter :: CL(0:MAXW-1,0:MAXW-1,0:NC-1) = WENOZ9_L, CDQ(0:MAXW-1,0:NC-1) = WENOZ9_DQ
   real(R8P),    parameter :: D(0:NC-1) = WENOZ9_D, TC(0:NC-1) = WENOZ9_TC, TG = WENOZ9_TG
   real(R8P),    parameter :: CLG(0:GW-1,0:GW-1) = WENOZ9_LG, CDG(0:GW-1) = WENOZ9_DG
#include "seth_z_teno_body.INC"
   endsubroutine seth_wenoz9

endmodule seth_wenoz
