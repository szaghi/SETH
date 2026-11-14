!< SETH, TENO kernels: TENO5, TENO6 and TENO6-opt reconstruction at one face.
module seth_teno
!< SETH, TENO kernels: TENO5, TENO6 and TENO6-opt reconstruction at one face.
!<
!< Targeted ENO (Fu, Hu & Adams 2016, as restated by Li, Fu & Adams 2020): candidate stencils of incremental width, the
!< smoothness measures gamma_k = (1 + tau/(beta_k + eps))**6 normalised and cut at C_T, the linear weights of the kept
!< candidates renormalised. One kernel per scheme, its coefficients compile-time constants (seth_z_teno_tables), the
!< body shared with WENO-Z (seth_z_teno_body.INC). Every kernel is pure, `acc routine seq` and `omp declare target`:
!<
!<    call seth_teno6(vp, vm, eps, ct, rp, rm)
!<
!< `vp`/`vm` the values reconstructed with a left/right bias, cells by offset from the face (cell 0 left of it, bounds
!< -2:3), `eps` the regulariser (SETH_TENO_EPS), `ct` the cutoff, `rp`/`rm` the left/right states at the face. Only the
!< sourced C_T is published (TENO6-opt); for TENO5 and TENO6 C_T is the caller's choice (notes/m2-sources.md).
use penf, only : I4P, R8P
use seth_z_teno_tables
implicit none
private
public :: SETH_TENO_EPS, SETH_TENO6_OPT_CT
public :: seth_teno5
public :: seth_teno6
public :: seth_teno6_opt

real(R8P), parameter :: SETH_TENO_EPS     = 1.e-40_R8P !< eps of the TENO papers (Li, Fu & Adams 2020, eq. 9).
real(R8P), parameter :: SETH_TENO6_OPT_CT = 1.e-6_R8P  !< C_T of TENO6-opt (NENO6 paper, section 2.2).

contains
   pure subroutine seth_teno5(vp, vm, eps, ct, rp, rm)
   !< TENO5: candidates S0..S2 (3-point), tau_5 = |beta_1 - beta_2|, linear weights 6/10, 3/10, 1/10.
   integer(I4P), parameter :: NC = TENO5_NC, LO = TENO5_LO, UMAX = TENO5_UMAX, MAXW = TENO5_MAXW
   integer(I4P), parameter :: CLO(0:NC-1) = TENO5_CLO, CW(0:NC-1) = TENO5_CW
   integer(I4P), parameter :: GLO = TENO5_GLO, GW = TENO5_GW
   real(R8P),    parameter :: CP(0:MAXW-1,0:NC-1) = TENO5_P
   real(R8P),    parameter :: CU(0:MAXW-1,0:MAXW-1,0:NC-1) = TENO5_U
   real(R8P),    parameter :: D(0:NC-1) = TENO5_D, TC(0:NC-1) = TENO5_TC, TG = TENO5_TG
   real(R8P),    parameter :: CUG(0:GW-1,0:GW-1) = TENO5_UG
#define SETH_TENO
#include "seth_z_teno_body.INC"
#undef SETH_TENO
   endsubroutine seth_teno5

   pure subroutine seth_teno6(vp, vm, eps, ct, rp, rm)
   !< TENO6: candidates S0..S3 (3,3,3,4 points), tau_6 = beta_6 - (beta_1 + beta_2 + 4 beta_0)/6, max-order weights.
   integer(I4P), parameter :: NC = TENO6_NC, LO = TENO6_LO, UMAX = TENO6_UMAX, MAXW = TENO6_MAXW
   integer(I4P), parameter :: CLO(0:NC-1) = TENO6_CLO, CW(0:NC-1) = TENO6_CW
   integer(I4P), parameter :: GLO = TENO6_GLO, GW = TENO6_GW
   real(R8P),    parameter :: CP(0:MAXW-1,0:NC-1) = TENO6_P
   real(R8P),    parameter :: CU(0:MAXW-1,0:MAXW-1,0:NC-1) = TENO6_U
   real(R8P),    parameter :: D(0:NC-1) = TENO6_D, TC(0:NC-1) = TENO6_TC, TG = TENO6_TG
   real(R8P),    parameter :: CUG(0:GW-1,0:GW-1) = TENO6_UG
#define SETH_TENO
#include "seth_z_teno_body.INC"
#undef SETH_TENO
   endsubroutine seth_teno6

   pure subroutine seth_teno6_opt(vp, vm, eps, ct, rp, rm)
   !< TENO6-opt: TENO6 with the spectrally optimised linear weights (NENO6 paper); use SETH_TENO6_OPT_CT.
   integer(I4P), parameter :: NC = TENO6_OPT_NC, LO = TENO6_OPT_LO, UMAX = TENO6_OPT_UMAX, MAXW = TENO6_OPT_MAXW
   integer(I4P), parameter :: CLO(0:NC-1) = TENO6_OPT_CLO, CW(0:NC-1) = TENO6_OPT_CW
   integer(I4P), parameter :: GLO = TENO6_OPT_GLO, GW = TENO6_OPT_GW
   real(R8P),    parameter :: CP(0:MAXW-1,0:NC-1) = TENO6_OPT_P
   real(R8P),    parameter :: CU(0:MAXW-1,0:MAXW-1,0:NC-1) = TENO6_OPT_U
   real(R8P),    parameter :: D(0:NC-1) = TENO6_OPT_D, TC(0:NC-1) = TENO6_OPT_TC, TG = TENO6_OPT_TG
   real(R8P),    parameter :: CUG(0:GW-1,0:GW-1) = TENO6_OPT_UG
#define SETH_TENO
#include "seth_z_teno_body.INC"
#undef SETH_TENO
   endsubroutine seth_teno6_opt

endmodule seth_teno
