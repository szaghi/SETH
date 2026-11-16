!< SETH, TENO kernels: TENO5, TENO5-opt, TENO5-AZ, TENO6, TENO6-opt, TENO6-NENO6; TENO5/6 interpolation.
module seth_teno
!< SETH, TENO kernels: TENO5, TENO5-opt, TENO5-AZ, TENO6, TENO6-opt, TENO6-NENO6; TENO5/6 interpolation.
!<
!< Targeted ENO (Fu, Hu & Adams 2016, as restated by Li, Fu & Adams 2020): candidate stencils of incremental width, the
!< smoothness measures gamma_k = (1 + tau/(beta_k + eps))**6 normalised and cut at C_T, the linear weights of the kept
!< candidates renormalised. One kernel per scheme, its coefficients compile-time constants (seth_z_teno_tables), the
!< body shared with WENO-Z (seth_z_teno_body.INC). Every kernel is pure, `acc routine seq` and `omp declare target`:
!<
!<    call seth_teno6(vp, vm, eps, ct, rp, rm)
!<
!< `vp`/`vm` the values reconstructed with a left/right bias, cells by offset from the face (cell 0 left of it, bounds
!< -2:3), `eps` the regulariser (SETH_TENO_EPS), `ct` the cutoff (SETH_TENO*_CT: Fu, Hu & Adams 2016, Table 5, and
!< the NENO6 paper), `rp`/`rm` the left/right states at the face. Sources: notes/m2-sources.md. The `_interp`
!< kernels take point values at the cell centres: other face coefficients and linear weights, the same indicators,
!< tau and cutoff (notes/interpolation-plan.md); TENO6-opt and NENO6 have none (their weights are tuned for the
!< reconstruction).
use penf, only : I4P, R8P
use seth_z_teno_tables
implicit none
private
public :: SETH_TENO_EPS, SETH_TENO5_CT, SETH_TENO6_CT, SETH_TENO6_OPT_CT, SETH_TENO6_NENO6_CT
public :: SETH_TENO5_OPT_CT, SETH_TENO5_AZ_CT
public :: seth_teno5
public :: seth_teno5_interp
public :: seth_teno5_opt
public :: seth_teno5_az
public :: seth_teno6
public :: seth_teno6_interp
public :: seth_teno6_opt
public :: seth_teno6_neno6

real(R8P), parameter :: SETH_TENO_EPS     = 1.e-40_R8P !< eps of the TENO papers (Li, Fu & Adams 2020, eq. 9).
real(R8P), parameter :: SETH_TENO5_CT       = 1.e-5_R8P !< C_T of TENO5 (Fu, Hu & Adams 2016, Table 5).
real(R8P), parameter :: SETH_TENO6_CT       = 1.e-7_R8P !< C_T of TENO6 (Fu, Hu & Adams 2016, Table 5).
real(R8P), parameter :: SETH_TENO5_OPT_CT   = 1.e-5_R8P !< C_T of TENO5-opt (Fu, Hu & Adams 2016, Table 5).
real(R8P), parameter :: SETH_TENO5_AZ_CT    = 1.e-5_R8P !< lambda_d of TENO5-AZ, the cutoff's ceiling (Tang 2024, Eq. 20).
real(R8P), parameter :: SETH_TENO6_OPT_CT   = 1.e-7_R8P !< C_T of TENO6-opt (Fu, Hu & Adams 2016, Table 5).
real(R8P), parameter :: SETH_TENO6_NENO6_CT = 1.e-6_R8P !< C_T used with the NENO6 weights (NENO6 paper, section 2.2).

contains
   pure subroutine seth_teno5(vp, vm, eps, ct, rp, rm)
   !< TENO5: candidates S0..S2 (3-point), tau_5 = |beta_1 - beta_2|, linear weights 6/10, 3/10, 1/10.
   integer(I4P), parameter :: NC = TENO5_NC, LO = TENO5_LO, UMAX = TENO5_UMAX, MAXW = TENO5_MAXW
   integer(I4P), parameter :: CLO(0:NC-1) = TENO5_CLO, CW(0:NC-1) = TENO5_CW
   integer(I4P), parameter :: GLO = TENO5_GLO, GW = TENO5_GW
   real(R8P),    parameter :: CP(0:MAXW-1,0:NC-1) = TENO5_P
   real(R8P),    parameter :: CL(0:MAXW-1,0:MAXW-1,0:NC-1) = TENO5_L, CDQ(0:MAXW-1,0:NC-1) = TENO5_DQ
   real(R8P),    parameter :: D(0:NC-1) = TENO5_D, TC(0:NC-1) = TENO5_TC, TG = TENO5_TG
   real(R8P),    parameter :: CLG(0:GW-1,0:GW-1) = TENO5_LG, CDG(0:GW-1) = TENO5_DG
#define SETH_TENO
#include "seth_z_teno_body.INC"
#undef SETH_TENO
   endsubroutine seth_teno5

   pure subroutine seth_teno5_interp(vp, vm, eps, ct, rp, rm)
   !< TENO5 interpolation from point values: weights 5/8, 5/16, 1/16 (same indicators, tau and cutoff).
   integer(I4P), parameter :: NC = TENO5_NC, LO = TENO5_LO, UMAX = TENO5_UMAX, MAXW = TENO5_MAXW
   integer(I4P), parameter :: CLO(0:NC-1) = TENO5_CLO, CW(0:NC-1) = TENO5_CW
   integer(I4P), parameter :: GLO = TENO5_GLO, GW = TENO5_GW
   real(R8P),    parameter :: CP(0:MAXW-1,0:NC-1) = TENO5_INTERP_P
   real(R8P),    parameter :: CL(0:MAXW-1,0:MAXW-1,0:NC-1) = TENO5_L, CDQ(0:MAXW-1,0:NC-1) = TENO5_DQ
   real(R8P),    parameter :: D(0:NC-1) = TENO5_INTERP_D, TC(0:NC-1) = TENO5_TC, TG = TENO5_TG
   real(R8P),    parameter :: CLG(0:GW-1,0:GW-1) = TENO5_LG, CDG(0:GW-1) = TENO5_DG
#define SETH_TENO
#include "seth_z_teno_body.INC"
#undef SETH_TENO
   endsubroutine seth_teno5_interp

   pure subroutine seth_teno5_opt(vp, vm, eps, ct, rp, rm)
   !< TENO5-opt (Fu, Hu & Adams 2016, Table 5): optimised linear weights 11/20, 2/5, 1/20 (eta = 0.2), 4th order.
   integer(I4P), parameter :: NC = TENO5_OPT_NC, LO = TENO5_OPT_LO, UMAX = TENO5_OPT_UMAX, MAXW = TENO5_OPT_MAXW
   integer(I4P), parameter :: CLO(0:NC-1) = TENO5_OPT_CLO, CW(0:NC-1) = TENO5_OPT_CW
   integer(I4P), parameter :: GLO = TENO5_OPT_GLO, GW = TENO5_OPT_GW
   real(R8P),    parameter :: CP(0:MAXW-1,0:NC-1) = TENO5_OPT_P
   real(R8P),    parameter :: CL(0:MAXW-1,0:MAXW-1,0:NC-1) = TENO5_OPT_L, CDQ(0:MAXW-1,0:NC-1) = TENO5_OPT_DQ
   real(R8P),    parameter :: D(0:NC-1) = TENO5_OPT_D, TC(0:NC-1) = TENO5_OPT_TC, TG = TENO5_OPT_TG
   real(R8P),    parameter :: CLG(0:GW-1,0:GW-1) = TENO5_OPT_LG, CDG(0:GW-1) = TENO5_OPT_DG
#define SETH_TENO
#include "seth_z_teno_body.INC"
#undef SETH_TENO
   endsubroutine seth_teno5_opt

   pure subroutine seth_teno5_az(vp, vm, eps, ct, rmu, rp, rm)
   !< TENO5-AZ (Tang 2024): TENO5 with tau = |beta_1 + beta_2 - 2 beta_0| and the cutoff min(ct, tau rmu^2), ct = lambda_d;
   !< rmu = 1 is Tang's form, rmu = seth_weno_descaling(3, 1, vp, vm) makes the cutoff scale invariant.
   integer(I4P), parameter :: NC = TENO5_AZ_NC, LO = TENO5_AZ_LO, UMAX = TENO5_AZ_UMAX, MAXW = TENO5_AZ_MAXW
   integer(I4P), parameter :: CLO(0:NC-1) = TENO5_AZ_CLO, CW(0:NC-1) = TENO5_AZ_CW
   integer(I4P), parameter :: GLO = TENO5_AZ_GLO, GW = TENO5_AZ_GW
   real(R8P),    parameter :: CP(0:MAXW-1,0:NC-1) = TENO5_AZ_P
   real(R8P),    parameter :: CL(0:MAXW-1,0:MAXW-1,0:NC-1) = TENO5_AZ_L, CDQ(0:MAXW-1,0:NC-1) = TENO5_AZ_DQ
   real(R8P),    parameter :: D(0:NC-1) = TENO5_AZ_D, TC(0:NC-1) = TENO5_AZ_TC, TG = TENO5_AZ_TG
   real(R8P),    parameter :: CLG(0:GW-1,0:GW-1) = TENO5_AZ_LG, CDG(0:GW-1) = TENO5_AZ_DG
#define SETH_TENO
#define SETH_TENO_AZ
#include "seth_z_teno_body.INC"
#undef SETH_TENO_AZ
#undef SETH_TENO
   endsubroutine seth_teno5_az

   pure subroutine seth_teno6(vp, vm, eps, ct, rp, rm)
   !< TENO6: candidates S0..S3 (3,3,3,4 points), tau_6 = beta_6 - (beta_1 + beta_2 + 4 beta_0)/6, max-order weights.
   integer(I4P), parameter :: NC = TENO6_NC, LO = TENO6_LO, UMAX = TENO6_UMAX, MAXW = TENO6_MAXW
   integer(I4P), parameter :: CLO(0:NC-1) = TENO6_CLO, CW(0:NC-1) = TENO6_CW
   integer(I4P), parameter :: GLO = TENO6_GLO, GW = TENO6_GW
   real(R8P),    parameter :: CP(0:MAXW-1,0:NC-1) = TENO6_P
   real(R8P),    parameter :: CL(0:MAXW-1,0:MAXW-1,0:NC-1) = TENO6_L, CDQ(0:MAXW-1,0:NC-1) = TENO6_DQ
   real(R8P),    parameter :: D(0:NC-1) = TENO6_D, TC(0:NC-1) = TENO6_TC, TG = TENO6_TG
   real(R8P),    parameter :: CLG(0:GW-1,0:GW-1) = TENO6_LG, CDG(0:GW-1) = TENO6_DG
#define SETH_TENO
#include "seth_z_teno_body.INC"
#undef SETH_TENO
   endsubroutine seth_teno6

   pure subroutine seth_teno6_interp(vp, vm, eps, ct, rp, rm)
   !< TENO6 interpolation from point values: weights 15/32, 5/16, 1/32, 3/16 (derived; same indicators, tau, cutoff).
   integer(I4P), parameter :: NC = TENO6_NC, LO = TENO6_LO, UMAX = TENO6_UMAX, MAXW = TENO6_MAXW
   integer(I4P), parameter :: CLO(0:NC-1) = TENO6_CLO, CW(0:NC-1) = TENO6_CW
   integer(I4P), parameter :: GLO = TENO6_GLO, GW = TENO6_GW
   real(R8P),    parameter :: CP(0:MAXW-1,0:NC-1) = TENO6_INTERP_P
   real(R8P),    parameter :: CL(0:MAXW-1,0:MAXW-1,0:NC-1) = TENO6_L, CDQ(0:MAXW-1,0:NC-1) = TENO6_DQ
   real(R8P),    parameter :: D(0:NC-1) = TENO6_INTERP_D, TC(0:NC-1) = TENO6_TC, TG = TENO6_TG
   real(R8P),    parameter :: CLG(0:GW-1,0:GW-1) = TENO6_LG, CDG(0:GW-1) = TENO6_DG
#define SETH_TENO
#include "seth_z_teno_body.INC"
#undef SETH_TENO
   endsubroutine seth_teno6_interp

   pure subroutine seth_teno6_opt(vp, vm, eps, ct, rp, rm)
   !< TENO6-opt (Fu, Hu & Adams 2016): optimised linear weights 0.462, 0.3, 0.054, 0.184, 5th order; SETH_TENO6_OPT_CT.
   integer(I4P), parameter :: NC = TENO6_OPT_NC, LO = TENO6_OPT_LO, UMAX = TENO6_OPT_UMAX, MAXW = TENO6_OPT_MAXW
   integer(I4P), parameter :: CLO(0:NC-1) = TENO6_OPT_CLO, CW(0:NC-1) = TENO6_OPT_CW
   integer(I4P), parameter :: GLO = TENO6_OPT_GLO, GW = TENO6_OPT_GW
   real(R8P),    parameter :: CP(0:MAXW-1,0:NC-1) = TENO6_OPT_P
   real(R8P),    parameter :: CL(0:MAXW-1,0:MAXW-1,0:NC-1) = TENO6_OPT_L, CDQ(0:MAXW-1,0:NC-1) = TENO6_OPT_DQ
   real(R8P),    parameter :: D(0:NC-1) = TENO6_OPT_D, TC(0:NC-1) = TENO6_OPT_TC, TG = TENO6_OPT_TG
   real(R8P),    parameter :: CLG(0:GW-1,0:GW-1) = TENO6_OPT_LG, CDG(0:GW-1) = TENO6_OPT_DG
#define SETH_TENO
#include "seth_z_teno_body.INC"
#undef SETH_TENO
   endsubroutine seth_teno6_opt

   pure subroutine seth_teno6_neno6(vp, vm, eps, ct, rp, rm)
   !< TENO6 with the optimised weights of the NENO6 paper (4th order); use SETH_TENO6_NENO6_CT.
   integer(I4P), parameter :: NC = TENO6_NENO6_NC, LO = TENO6_NENO6_LO, UMAX = TENO6_NENO6_UMAX, MAXW = TENO6_NENO6_MAXW
   integer(I4P), parameter :: CLO(0:NC-1) = TENO6_NENO6_CLO, CW(0:NC-1) = TENO6_NENO6_CW
   integer(I4P), parameter :: GLO = TENO6_NENO6_GLO, GW = TENO6_NENO6_GW
   real(R8P),    parameter :: CP(0:MAXW-1,0:NC-1) = TENO6_NENO6_P
   real(R8P),    parameter :: CL(0:MAXW-1,0:MAXW-1,0:NC-1) = TENO6_NENO6_L, CDQ(0:MAXW-1,0:NC-1) = TENO6_NENO6_DQ
   real(R8P),    parameter :: D(0:NC-1) = TENO6_NENO6_D, TC(0:NC-1) = TENO6_NENO6_TC, TG = TENO6_NENO6_TG
   real(R8P),    parameter :: CLG(0:GW-1,0:GW-1) = TENO6_NENO6_LG, CDG(0:GW-1) = TENO6_NENO6_DG
#define SETH_TENO
#include "seth_z_teno_body.INC"
#undef SETH_TENO
   endsubroutine seth_teno6_neno6

endmodule seth_teno
