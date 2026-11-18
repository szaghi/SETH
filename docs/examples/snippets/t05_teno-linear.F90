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
