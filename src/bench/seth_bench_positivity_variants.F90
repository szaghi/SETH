!< SETH benchmark: variants of the positivity limiter with constant array bounds and/or loop extents.
module seth_bench_positivity_max
!< Local arrays over-dimensioned to 5 (constant bounds), loops to the dummy nv: ADAM's strategy 1 alone.
use penf, only : I4P, R8P
implicit none
private
public :: positivity_euler
contains
#define NVD 5
#define NVL nv
#include "seth_bench_positivity_body.INC"
#undef NVD
#undef NVL
endmodule seth_bench_positivity_max

module seth_bench_positivity_nv3
!< Constant bounds and constant loop extents, nv = 3 (the dummy nv is ignored).
use penf, only : I4P, R8P
implicit none
private
public :: positivity_euler
contains
#define NVD 3
#define NVL 3
#include "seth_bench_positivity_body.INC"
#undef NVD
#undef NVL
endmodule seth_bench_positivity_nv3

module seth_bench_positivity_nv5
!< Constant bounds and constant loop extents, nv = 5 (the dummy nv is ignored).
use penf, only : I4P, R8P
implicit none
private
public :: positivity_euler
contains
#define NVD 5
#define NVL 5
#include "seth_bench_positivity_body.INC"
#undef NVD
#undef NVL
endmodule seth_bench_positivity_nv5
