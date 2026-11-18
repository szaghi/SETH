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
