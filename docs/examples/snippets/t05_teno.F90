!< Tutorial 5: TENO and its cutoff. Exactly linear where smooth, and what a critical point does to the cutoff.
program t05_teno
use seth
implicit none
real(R_P), parameter :: PI = acos(-1._R_P)
real(R_P)            :: dmax, u(-3:67), h, rp, rm, lp, lm, dev(3), ct(3)
real(R_P)            :: e5(6), eaz(6), hk
integer              :: i, k, j
character(28)        :: label(4) = [character(28) :: 'sine, 16 cells per wave', 'sine, 4 cells per wave', &
                                       'Gaussian, sigma 2.6 cells', 'two steps']

! four fields on 64 periodic cells: the faces where TENO5 is not its linear scheme, for three cutoffs
ct = [1.e-3_R_P, SETH_TENO5_CT, 1.e-7_R_P]
print '(A)', 'faces (of 64) where TENO5 is not its linear scheme'
print '(A,3ES10.0,A)', '  data                  C_T =', ct, '  max |TENO5 - linear|'
h = 1._R_P/64
do j=1, 4
   do i=1, 64
      select case(j)
      case(1) ; u(i) = (cos(2*PI*4*(i - 1)*h) - cos(2*PI*4*i*h))/(2*PI*4*h)        ! 16 cells per wave
      case(2) ; u(i) = (cos(2*PI*16*(i - 1)*h) - cos(2*PI*16*i*h))/(2*PI*16*h)     ! 4 cells per wave
      case(3) ; u(i) = exp(-300*((i - 0.5_R_P)*h - 0.5_R_P)**2)                    ! sigma = 2.6 cells
      case(4) ; u(i) = merge(1._R_P, 0._R_P, i > 32)                                ! steps at x = 0.5 and 1
      endselect
   enddo
   u(-2:0) = u(62:64) ; u(65:67) = u(1:3)
   dev = 0 ; dmax = 0
   do k=1, 3
      do i=0, 63
         call seth_teno5(u(i-2:i+3), u(i-2:i+3), SETH_TENO_EPS, ct(k), rp, rm)
         call seth_linear5(u(i-2:i+3), u(i-2:i+3), lp, lm)
         if (abs(rp - lp) > 4*epsilon(1._R_P) .or. abs(rm - lm) > 4*epsilon(1._R_P)) dev(k) = dev(k) + 1
         if (k == 2) dmax = max(dmax, abs(rp - lp), abs(rm - lm))
      enddo
   enddo
   if (dmax < 1.e-14_R_P) then                  ! round-off: printed as such, its digits depend on the machine
      print '(2X,A28,3F10.0,A14)', label(j), dev, 'round-off'
   else
      print '(2X,A28,3F10.0,ES14.2)', label(j), dev, dmax
   endif
enddo

! f = x^3/3! + x^4/4! + ... + x^7/7! (e^x - 1 - x - x^2/2, truncated) has f' = f'' = 0 at x = 0. Finite differences
! at the node 0 (Henrick et al. 2005): the point values f(m h) as the kernel's data, the derivative
! (F_{1/2} - F_{-1/2})/h, exactly f'(0) = 0
hk = 0.04_R_P
do k=1, 6
   do i=-3, 3
      u(i) = f(i*hk)
   enddo
   call seth_teno5(u(-2:3), u(-2:3), SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)
   call seth_teno5(u(-3:2), u(-3:2), SETH_TENO_EPS, SETH_TENO5_CT, lp, lm)
   e5(k) = abs(rp - lp)/hk
   call seth_teno5_az(u(-2:3), u(-2:3), SETH_TENO_EPS, SETH_TENO5_AZ_CT, 1._R_P, rp, rm)
   call seth_teno5_az(u(-3:2), u(-3:2), SETH_TENO_EPS, SETH_TENO5_AZ_CT, 1._R_P, lp, lm)
   eaz(k) = abs(rp - lp)/hk
   hk = hk/2
enddo
print '(A)', 'error of the derivative at a second-order critical point'
print '(A)', '        h       TENO5  order    TENO5-AZ  order'
hk = 0.04_R_P
do k=1, 6
   if (k == 1) then
      print '(F9.5,ES12.3,7X,ES12.3)', hk, e5(k), eaz(k)
   else
      print '(F9.5,ES12.3,F7.2,ES12.3,F7.2)', hk, e5(k), log(e5(k-1)/e5(k))/log(2._R_P), &
                                              eaz(k), log(eaz(k-1)/eaz(k))/log(2._R_P)
   endif
   hk = hk/2
enddo

contains
   pure function f(x) result(v)
   real(R_P), intent(in) :: x
   real(R_P)             :: v
   v = x**3*(1._R_P/6 + x*(1._R_P/24 + x*(1._R_P/120 + x*(1._R_P/720 + x/5040))))
   endfunction f
endprogram t05_teno
