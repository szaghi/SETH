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
implicit none
private
public :: I4P, I8P, R4P, R8P
public :: seth_weno_descaling
public :: SETH_WENO_JS_EPS
public :: seth_weno1_js, seth_weno1_js_wr, seth_weno1_js_interp, seth_weno1_js_interp_wr
public :: seth_weno3_js, seth_weno3_js_wr, seth_weno3_js_interp, seth_weno3_js_interp_wr
public :: seth_weno5_js, seth_weno5_js_wr, seth_weno5_js_interp, seth_weno5_js_interp_wr
public :: seth_weno7_js, seth_weno7_js_wr, seth_weno7_js_interp, seth_weno7_js_interp_wr
public :: seth_weno9_js, seth_weno9_js_wr, seth_weno9_js_interp, seth_weno9_js_interp_wr
public :: seth_wenoz3, seth_wenoz5, seth_wenoz7
public :: WENO1_JS_D_SCALE, WENO3_JS_D_SCALE, WENO5_JS_D_SCALE, WENO7_JS_D_SCALE, WENO9_JS_D_SCALE
endmodule seth
