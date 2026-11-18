!< SETH tests: the generic smooth test function of the order tests, its derivative and its exact cell averages.
module seth_test_functions
!< SETH tests: the generic smooth test function of the order tests, its derivative and its exact cell averages.
!<
!< g(x) = sin(0.9 x + 0.4) + 0.5 exp(0.6 x) + 0.3 cos(1.7 x - 1.1): three modes of no common pattern. The order tests
!< used exp(x) before, which is degenerate for them: the cell averages of exp are a constant multiple of its point
!< values (sinh(h/2)/(h/2)), so a test cannot tell a reconstruction from an interpolation, and exp satisfies
!< f' f'''' = f'' f''' and similar identities, which cancel the leading error terms of the nonlinear weights and can
!< make an order look higher than it is. g has no such identity, and g' > 1 near the test point x = 0.3 (no critical
!< point: the order tests measure the design order). Its wavenumbers stay below 2: with a mode of wavenumber 2.9 the
!< 9th-order schemes, whose stencils span 10 cells, were still pre-asymptotic when their error reached round-off
!< (observed ratios swinging between 7.5 and 12.7).
!<
!< The cell average is a 10-point Gauss-Legendre quadrature: exact to degree 19, an error far below round-off for the
!< cell widths of the tests (h <= 0.8), and no cancellation, unlike (G(b) - G(a))/h with the primitive G.
use penf, only : I4P, R_P
implicit none
private
public :: gfun, gder, gder2, gavg, GL_X, GL_W

real(R_P), parameter :: GL_X(5) = [0.1488743389816312_R_P, 0.4333953941292472_R_P, 0.6794095682990244_R_P, &
                                   0.8650633666889845_R_P, 0.9739065285171717_R_P] !< Gauss-Legendre nodes (> 0).
real(R_P), parameter :: GL_W(5) = [0.2955242247147529_R_P, 0.2692667193099963_R_P, 0.2190863625159820_R_P, &
                                   0.1494513491505806_R_P, 0.0666713443086881_R_P] !< Their weights.

contains
   elemental function gfun(x) result(y)
   !< g(x).
   real(R_P), intent(in) :: x
   real(R_P)             :: y
   y = sin(0.9_R_P*x + 0.4_R_P) + 0.5_R_P*exp(0.6_R_P*x) + 0.3_R_P*cos(1.7_R_P*x - 1.1_R_P)
   endfunction gfun

   elemental function gder(x) result(y)
   !< g'(x).
   real(R_P), intent(in) :: x
   real(R_P)             :: y
   y = 0.9_R_P*cos(0.9_R_P*x + 0.4_R_P) + 0.3_R_P*exp(0.6_R_P*x) - 0.51_R_P*sin(1.7_R_P*x - 1.1_R_P)
   endfunction gder

   elemental function gder2(x) result(y)
   !< g''(x).
   real(R_P), intent(in) :: x
   real(R_P)             :: y
   y = -0.81_R_P*sin(0.9_R_P*x + 0.4_R_P) + 0.18_R_P*exp(0.6_R_P*x) - 0.867_R_P*cos(1.7_R_P*x - 1.1_R_P)
   endfunction gder2

   elemental function gavg(a, b) result(y)
   !< The average of g over [a, b] (10-point Gauss-Legendre).
   real(R_P), intent(in) :: a, b
   real(R_P)             :: y, c, r
   integer(I4P)          :: k
   c = 0.5_R_P*(a + b) ; r = 0.5_R_P*(b - a)
   y = 0._R_P
   do k=1, 5
      y = y + GL_W(k)*(gfun(c - r*GL_X(k)) + gfun(c + r*GL_X(k)))
   enddo
   y = 0.5_R_P*y
   endfunction gavg
endmodule seth_test_functions
