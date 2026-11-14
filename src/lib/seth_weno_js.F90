!< SETH, WENO-JS kernels: WENO-1/3/5/7/9 reconstruction and interpolation at one face.
module seth_weno_js
!< SETH, WENO-JS kernels: WENO-1/3/5/7/9 reconstruction and interpolation at one face.
!<
!< One kernel per scheme, its coefficients compile-time constants (seth_weno_js_tables): no tables to pass, nothing to
!< mirror on a device. Every kernel is pure, `acc routine seq` and `omp declare target`, and has the same dummies:
!<
!<    call seth_weno5_js(vp, vm, eps, rmu, rp, rm)
!<
!< `vp`/`vm` the values reconstructed with a left/right bias, cells by offset from the face (cell 0 left of it,
!< bounds 1-S:S), `eps` the regulariser (SETH_WENO_JS_EPS: ADAM's value), `rmu` the inverse descaler (1 for the
!< Jiang-Shu weights, see seth_weno_descaling), `rp`/`rm` the left/right states at the face. The `_wr` kernels also
!< return the weights ratios `wrp`, `wrm`; the `_interp` kernels take point values instead of cell averages.
use penf, only : I4P, R8P
use seth_weno_js_tables
implicit none
private
public :: SETH_WENO_JS_EPS
public :: seth_weno1_js
public :: seth_weno1_js_wr
public :: seth_weno1_js_interp
public :: seth_weno1_js_interp_wr
public :: seth_weno3_js
public :: seth_weno3_js_wr
public :: seth_weno3_js_interp
public :: seth_weno3_js_interp_wr
public :: seth_weno5_js
public :: seth_weno5_js_wr
public :: seth_weno5_js_interp
public :: seth_weno5_js_interp_wr
public :: seth_weno7_js
public :: seth_weno7_js_wr
public :: seth_weno7_js_interp
public :: seth_weno7_js_interp_wr
public :: seth_weno9_js
public :: seth_weno9_js_wr
public :: seth_weno9_js_interp
public :: seth_weno9_js_interp_wr

real(R8P), parameter :: SETH_WENO_JS_EPS = 1.e-6_R8P !< ADAM's regulariser (absolute, on the IS as the tables scale them).

contains
   pure subroutine seth_weno1_js(vp, vm, eps, rmu, rp, rm)
   !< WENO-1 (S = 1) reconstruction from cell averages.
   integer(I4P), parameter :: S = 1_I4P
   integer(I4P), parameter :: WEXP = 1_I4P
   real(R8P),    parameter :: WA(1:2,0:S-1) = WENO1_JS_A
   real(R8P),    parameter :: WP(1:2,0:S-1,0:S-1) = WENO1_JS_P
   real(R8P),    parameter :: WD(0:S-1,0:S-1,0:S-1) = WENO1_JS_D
#include "seth_weno_js_body.INC"
   endsubroutine seth_weno1_js

   pure subroutine seth_weno1_js_wr(vp, vm, eps, rmu, rp, rm, wrp, wrm)
   !< WENO-1 (S = 1) reconstruction from cell averages, with the weights ratios.
   integer(I4P), parameter :: S = 1_I4P
   integer(I4P), parameter :: WEXP = 1_I4P
   real(R8P),    parameter :: WA(1:2,0:S-1) = WENO1_JS_A
   real(R8P),    parameter :: WP(1:2,0:S-1,0:S-1) = WENO1_JS_P
   real(R8P),    parameter :: WD(0:S-1,0:S-1,0:S-1) = WENO1_JS_D
#define SETH_WR
#include "seth_weno_js_body.INC"
#undef SETH_WR
   endsubroutine seth_weno1_js_wr

   pure subroutine seth_weno1_js_interp(vp, vm, eps, rmu, rp, rm)
   !< WENO-1 (S = 1) interpolation from point values.
   integer(I4P), parameter :: S = 1_I4P
   integer(I4P), parameter :: WEXP = 1_I4P
   real(R8P),    parameter :: WA(1:2,0:S-1) = WENO1_JS_INTERP_A
   real(R8P),    parameter :: WP(1:2,0:S-1,0:S-1) = WENO1_JS_INTERP_P
   real(R8P),    parameter :: WD(0:S-1,0:S-1,0:S-1) = WENO1_JS_D
#include "seth_weno_js_body.INC"
   endsubroutine seth_weno1_js_interp

   pure subroutine seth_weno1_js_interp_wr(vp, vm, eps, rmu, rp, rm, wrp, wrm)
   !< WENO-1 (S = 1) interpolation from point values, with the weights ratios.
   integer(I4P), parameter :: S = 1_I4P
   integer(I4P), parameter :: WEXP = 1_I4P
   real(R8P),    parameter :: WA(1:2,0:S-1) = WENO1_JS_INTERP_A
   real(R8P),    parameter :: WP(1:2,0:S-1,0:S-1) = WENO1_JS_INTERP_P
   real(R8P),    parameter :: WD(0:S-1,0:S-1,0:S-1) = WENO1_JS_D
#define SETH_WR
#include "seth_weno_js_body.INC"
#undef SETH_WR
   endsubroutine seth_weno1_js_interp_wr

   pure subroutine seth_weno3_js(vp, vm, eps, rmu, rp, rm)
   !< WENO-3 (S = 2) reconstruction from cell averages.
   integer(I4P), parameter :: S = 2_I4P
   integer(I4P), parameter :: WEXP = 2_I4P
   real(R8P),    parameter :: WA(1:2,0:S-1) = WENO3_JS_A
   real(R8P),    parameter :: WP(1:2,0:S-1,0:S-1) = WENO3_JS_P
   real(R8P),    parameter :: WD(0:S-1,0:S-1,0:S-1) = WENO3_JS_D
#include "seth_weno_js_body.INC"
   endsubroutine seth_weno3_js

   pure subroutine seth_weno3_js_wr(vp, vm, eps, rmu, rp, rm, wrp, wrm)
   !< WENO-3 (S = 2) reconstruction from cell averages, with the weights ratios.
   integer(I4P), parameter :: S = 2_I4P
   integer(I4P), parameter :: WEXP = 2_I4P
   real(R8P),    parameter :: WA(1:2,0:S-1) = WENO3_JS_A
   real(R8P),    parameter :: WP(1:2,0:S-1,0:S-1) = WENO3_JS_P
   real(R8P),    parameter :: WD(0:S-1,0:S-1,0:S-1) = WENO3_JS_D
#define SETH_WR
#include "seth_weno_js_body.INC"
#undef SETH_WR
   endsubroutine seth_weno3_js_wr

   pure subroutine seth_weno3_js_interp(vp, vm, eps, rmu, rp, rm)
   !< WENO-3 (S = 2) interpolation from point values.
   integer(I4P), parameter :: S = 2_I4P
   integer(I4P), parameter :: WEXP = 2_I4P
   real(R8P),    parameter :: WA(1:2,0:S-1) = WENO3_JS_INTERP_A
   real(R8P),    parameter :: WP(1:2,0:S-1,0:S-1) = WENO3_JS_INTERP_P
   real(R8P),    parameter :: WD(0:S-1,0:S-1,0:S-1) = WENO3_JS_D
#include "seth_weno_js_body.INC"
   endsubroutine seth_weno3_js_interp

   pure subroutine seth_weno3_js_interp_wr(vp, vm, eps, rmu, rp, rm, wrp, wrm)
   !< WENO-3 (S = 2) interpolation from point values, with the weights ratios.
   integer(I4P), parameter :: S = 2_I4P
   integer(I4P), parameter :: WEXP = 2_I4P
   real(R8P),    parameter :: WA(1:2,0:S-1) = WENO3_JS_INTERP_A
   real(R8P),    parameter :: WP(1:2,0:S-1,0:S-1) = WENO3_JS_INTERP_P
   real(R8P),    parameter :: WD(0:S-1,0:S-1,0:S-1) = WENO3_JS_D
#define SETH_WR
#include "seth_weno_js_body.INC"
#undef SETH_WR
   endsubroutine seth_weno3_js_interp_wr

   pure subroutine seth_weno5_js(vp, vm, eps, rmu, rp, rm)
   !< WENO-5 (S = 3) reconstruction from cell averages.
   integer(I4P), parameter :: S = 3_I4P
   integer(I4P), parameter :: WEXP = 3_I4P
   real(R8P),    parameter :: WA(1:2,0:S-1) = WENO5_JS_A
   real(R8P),    parameter :: WP(1:2,0:S-1,0:S-1) = WENO5_JS_P
   real(R8P),    parameter :: WD(0:S-1,0:S-1,0:S-1) = WENO5_JS_D
#include "seth_weno_js_body.INC"
   endsubroutine seth_weno5_js

   pure subroutine seth_weno5_js_wr(vp, vm, eps, rmu, rp, rm, wrp, wrm)
   !< WENO-5 (S = 3) reconstruction from cell averages, with the weights ratios.
   integer(I4P), parameter :: S = 3_I4P
   integer(I4P), parameter :: WEXP = 3_I4P
   real(R8P),    parameter :: WA(1:2,0:S-1) = WENO5_JS_A
   real(R8P),    parameter :: WP(1:2,0:S-1,0:S-1) = WENO5_JS_P
   real(R8P),    parameter :: WD(0:S-1,0:S-1,0:S-1) = WENO5_JS_D
#define SETH_WR
#include "seth_weno_js_body.INC"
#undef SETH_WR
   endsubroutine seth_weno5_js_wr

   pure subroutine seth_weno5_js_interp(vp, vm, eps, rmu, rp, rm)
   !< WENO-5 (S = 3) interpolation from point values.
   integer(I4P), parameter :: S = 3_I4P
   integer(I4P), parameter :: WEXP = 3_I4P
   real(R8P),    parameter :: WA(1:2,0:S-1) = WENO5_JS_INTERP_A
   real(R8P),    parameter :: WP(1:2,0:S-1,0:S-1) = WENO5_JS_INTERP_P
   real(R8P),    parameter :: WD(0:S-1,0:S-1,0:S-1) = WENO5_JS_D
#include "seth_weno_js_body.INC"
   endsubroutine seth_weno5_js_interp

   pure subroutine seth_weno5_js_interp_wr(vp, vm, eps, rmu, rp, rm, wrp, wrm)
   !< WENO-5 (S = 3) interpolation from point values, with the weights ratios.
   integer(I4P), parameter :: S = 3_I4P
   integer(I4P), parameter :: WEXP = 3_I4P
   real(R8P),    parameter :: WA(1:2,0:S-1) = WENO5_JS_INTERP_A
   real(R8P),    parameter :: WP(1:2,0:S-1,0:S-1) = WENO5_JS_INTERP_P
   real(R8P),    parameter :: WD(0:S-1,0:S-1,0:S-1) = WENO5_JS_D
#define SETH_WR
#include "seth_weno_js_body.INC"
#undef SETH_WR
   endsubroutine seth_weno5_js_interp_wr

   pure subroutine seth_weno7_js(vp, vm, eps, rmu, rp, rm)
   !< WENO-7 (S = 4) reconstruction from cell averages.
   integer(I4P), parameter :: S = 4_I4P
   integer(I4P), parameter :: WEXP = 4_I4P
   real(R8P),    parameter :: WA(1:2,0:S-1) = WENO7_JS_A
   real(R8P),    parameter :: WP(1:2,0:S-1,0:S-1) = WENO7_JS_P
   real(R8P),    parameter :: WD(0:S-1,0:S-1,0:S-1) = WENO7_JS_D
#include "seth_weno_js_body.INC"
   endsubroutine seth_weno7_js

   pure subroutine seth_weno7_js_wr(vp, vm, eps, rmu, rp, rm, wrp, wrm)
   !< WENO-7 (S = 4) reconstruction from cell averages, with the weights ratios.
   integer(I4P), parameter :: S = 4_I4P
   integer(I4P), parameter :: WEXP = 4_I4P
   real(R8P),    parameter :: WA(1:2,0:S-1) = WENO7_JS_A
   real(R8P),    parameter :: WP(1:2,0:S-1,0:S-1) = WENO7_JS_P
   real(R8P),    parameter :: WD(0:S-1,0:S-1,0:S-1) = WENO7_JS_D
#define SETH_WR
#include "seth_weno_js_body.INC"
#undef SETH_WR
   endsubroutine seth_weno7_js_wr

   pure subroutine seth_weno7_js_interp(vp, vm, eps, rmu, rp, rm)
   !< WENO-7 (S = 4) interpolation from point values.
   integer(I4P), parameter :: S = 4_I4P
   integer(I4P), parameter :: WEXP = 4_I4P
   real(R8P),    parameter :: WA(1:2,0:S-1) = WENO7_JS_INTERP_A
   real(R8P),    parameter :: WP(1:2,0:S-1,0:S-1) = WENO7_JS_INTERP_P
   real(R8P),    parameter :: WD(0:S-1,0:S-1,0:S-1) = WENO7_JS_D
#include "seth_weno_js_body.INC"
   endsubroutine seth_weno7_js_interp

   pure subroutine seth_weno7_js_interp_wr(vp, vm, eps, rmu, rp, rm, wrp, wrm)
   !< WENO-7 (S = 4) interpolation from point values, with the weights ratios.
   integer(I4P), parameter :: S = 4_I4P
   integer(I4P), parameter :: WEXP = 4_I4P
   real(R8P),    parameter :: WA(1:2,0:S-1) = WENO7_JS_INTERP_A
   real(R8P),    parameter :: WP(1:2,0:S-1,0:S-1) = WENO7_JS_INTERP_P
   real(R8P),    parameter :: WD(0:S-1,0:S-1,0:S-1) = WENO7_JS_D
#define SETH_WR
#include "seth_weno_js_body.INC"
#undef SETH_WR
   endsubroutine seth_weno7_js_interp_wr

   pure subroutine seth_weno9_js(vp, vm, eps, rmu, rp, rm)
   !< WENO-9 (S = 5) reconstruction from cell averages.
   integer(I4P), parameter :: S = 5_I4P
   integer(I4P), parameter :: WEXP = 4_I4P
   real(R8P),    parameter :: WA(1:2,0:S-1) = WENO9_JS_A
   real(R8P),    parameter :: WP(1:2,0:S-1,0:S-1) = WENO9_JS_P
   real(R8P),    parameter :: WD(0:S-1,0:S-1,0:S-1) = WENO9_JS_D
#include "seth_weno_js_body.INC"
   endsubroutine seth_weno9_js

   pure subroutine seth_weno9_js_wr(vp, vm, eps, rmu, rp, rm, wrp, wrm)
   !< WENO-9 (S = 5) reconstruction from cell averages, with the weights ratios.
   integer(I4P), parameter :: S = 5_I4P
   integer(I4P), parameter :: WEXP = 4_I4P
   real(R8P),    parameter :: WA(1:2,0:S-1) = WENO9_JS_A
   real(R8P),    parameter :: WP(1:2,0:S-1,0:S-1) = WENO9_JS_P
   real(R8P),    parameter :: WD(0:S-1,0:S-1,0:S-1) = WENO9_JS_D
#define SETH_WR
#include "seth_weno_js_body.INC"
#undef SETH_WR
   endsubroutine seth_weno9_js_wr

   pure subroutine seth_weno9_js_interp(vp, vm, eps, rmu, rp, rm)
   !< WENO-9 (S = 5) interpolation from point values.
   integer(I4P), parameter :: S = 5_I4P
   integer(I4P), parameter :: WEXP = 4_I4P
   real(R8P),    parameter :: WA(1:2,0:S-1) = WENO9_JS_INTERP_A
   real(R8P),    parameter :: WP(1:2,0:S-1,0:S-1) = WENO9_JS_INTERP_P
   real(R8P),    parameter :: WD(0:S-1,0:S-1,0:S-1) = WENO9_JS_D
#include "seth_weno_js_body.INC"
   endsubroutine seth_weno9_js_interp

   pure subroutine seth_weno9_js_interp_wr(vp, vm, eps, rmu, rp, rm, wrp, wrm)
   !< WENO-9 (S = 5) interpolation from point values, with the weights ratios.
   integer(I4P), parameter :: S = 5_I4P
   integer(I4P), parameter :: WEXP = 4_I4P
   real(R8P),    parameter :: WA(1:2,0:S-1) = WENO9_JS_INTERP_A
   real(R8P),    parameter :: WP(1:2,0:S-1,0:S-1) = WENO9_JS_INTERP_P
   real(R8P),    parameter :: WD(0:S-1,0:S-1,0:S-1) = WENO9_JS_D
#define SETH_WR
#include "seth_weno_js_body.INC"
#undef SETH_WR
   endsubroutine seth_weno9_js_interp_wr

endmodule seth_weno_js
