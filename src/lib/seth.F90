!< SETH, Smoothness-weighted ENO, TENO and Hybrids: the library's entry point.
module seth
!< SETH, Smoothness-weighted ENO, TENO and Hybrids: the library's entry point.
!<
!< `use seth` gives the working kinds and every scheme family. Each kernel reconstructs (or interpolates) the two
!< states at one face from two stencils, its coefficients compile-time constants; see seth_weno_js.
use penf, only : I4P, I8P, R4P, R8P
use seth_descaling
use seth_weno_js
use seth_weno_js_tables
use seth_wenoz
use seth_teno
use seth_positivity
use seth_mp
use seth_linear
implicit none
private
public :: I4P, I8P, R4P, R8P
public :: seth_weno_descaling, seth_weno_descaling_ai
public :: SETH_WENO_JS_EPS
public :: seth_weno1_js, seth_weno1_js_wr, seth_weno1_js_interp, seth_weno1_js_interp_wr
public :: seth_weno3_js, seth_weno3_js_wr, seth_weno3_js_interp, seth_weno3_js_interp_wr
public :: seth_weno5_js, seth_weno5_js_wr, seth_weno5_js_interp, seth_weno5_js_interp_wr
public :: seth_weno7_js, seth_weno7_js_wr, seth_weno7_js_interp, seth_weno7_js_interp_wr
public :: seth_weno9_js, seth_weno9_js_wr, seth_weno9_js_interp, seth_weno9_js_interp_wr
public :: seth_wenoz3, seth_wenoz5, seth_wenoz7, seth_wenoz9
public :: seth_wenoz3_interp, seth_wenoz5_interp, seth_wenoz7_interp, seth_wenoz9_interp
public :: seth_teno5, seth_teno5_opt, seth_teno5_az, seth_teno6, seth_teno6_opt, seth_teno6_neno6, seth_teno8
public :: seth_teno5_a, seth_teno8_a
public :: seth_teno5_interp, seth_teno6_interp
public :: SETH_TENO_EPS, SETH_TENO5_CT, SETH_TENO5_OPT_CT, SETH_TENO5_AZ_CT, SETH_TENO6_CT, SETH_TENO6_OPT_CT
public :: SETH_TENO6_NENO6_CT, SETH_TENO8_CT
public :: SETH_TENO5_A_CR, SETH_TENO5_A_XI, SETH_TENO5_A_A1, SETH_TENO5_A_A2
public :: SETH_TENO8_A_CR, SETH_TENO8_A_XI, SETH_TENO8_A_A1, SETH_TENO8_A_A2
public :: SETH_MP_ALPHA, seth_mp_bound
public :: seth_linear3, seth_linear5, seth_linear7, seth_linear9, seth_linear6, seth_linear8
public :: seth_linear3_interp, seth_linear5_interp, seth_linear7_interp, seth_linear9_interp, seth_linear6_interp
public :: SETH_POSITIVITY_EPS, SETH_POSITIVITY_NV_MAX, seth_positivity_euler, seth_euler_admissible
public :: WENO1_JS_D_SCALE, WENO3_JS_D_SCALE, WENO5_JS_D_SCALE, WENO7_JS_D_SCALE, WENO9_JS_D_SCALE
endmodule seth
