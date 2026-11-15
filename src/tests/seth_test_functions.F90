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
use penf, only : I4P, R8P
implicit none
private
public :: gfun, gder, gder2, gavg

real(R8P), parameter :: GL_X(5) = [0.1488743389816312_R8P, 0.4333953941292472_R8P, 0.6794095682990244_R8P, &
                                   0.8650633666889845_R8P, 0.9739065285171717_R8P] !< Gauss-Legendre nodes (> 0).
real(R8P), parameter :: GL_W(5) = [0.2955242247147529_R8P, 0.2692667193099963_R8P, 0.2190863625159820_R8P, &
                                   0.1494513491505806_R8P, 0.0666713443086881_R8P] !< Their weights.

contains
   elemental function gfun(x) result(y)
   !< g(x).
   real(R8P), intent(in) :: x
   real(R8P)             :: y
   y = sin(0.9_R8P*x + 0.4_R8P) + 0.5_R8P*exp(0.6_R8P*x) + 0.3_R8P*cos(1.7_R8P*x - 1.1_R8P)
   endfunction gfun

   elemental function gder(x) result(y)
   !< g'(x).
   real(R8P), intent(in) :: x
   real(R8P)             :: y
   y = 0.9_R8P*cos(0.9_R8P*x + 0.4_R8P) + 0.3_R8P*exp(0.6_R8P*x) - 0.51_R8P*sin(1.7_R8P*x - 1.1_R8P)
   endfunction gder

   elemental function gder2(x) result(y)
   !< g''(x).
   real(R8P), intent(in) :: x
   real(R8P)             :: y
   y = -0.81_R8P*sin(0.9_R8P*x + 0.4_R8P) + 0.18_R8P*exp(0.6_R8P*x) - 0.867_R8P*cos(1.7_R8P*x - 1.1_R8P)
   endfunction gder2

   elemental function gavg(a, b) result(y)
   !< The average of g over [a, b] (10-point Gauss-Legendre).
   real(R8P), intent(in) :: a, b
   real(R8P)             :: y, c, r
   integer(I4P)          :: k
   c = 0.5_R8P*(a + b) ; r = 0.5_R8P*(b - a)
   y = 0._R8P
   do k=1, 5
      y = y + GL_W(k)*(gfun(c - r*GL_X(k)) + gfun(c + r*GL_X(k)))
   enddo
   y = 0.5_R8P*y
   endfunction gavg
endmodule seth_test_functions
