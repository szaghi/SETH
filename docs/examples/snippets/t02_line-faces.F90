do i=0, n-1                                   ! the face x_{i+1/2} = i h, between cells i and i+1
   call seth_weno5_js(u(i-2:i+3), u(i-2:i+3), SETH_WENO_JS_EPS, 1._R_P, rp, rm)
   err(1) = err(1) + abs(rp - sin(2*PI*i*h))*h
   call seth_wenoz5(u(i-2:i+3), u(i-2:i+3), SETH_TENO_EPS, 1, rp, rm)
   err(2) = err(2) + abs(rp - sin(2*PI*i*h))*h
   call seth_teno5(u(i-2:i+3), u(i-2:i+3), SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)
   err(3) = err(3) + abs(rp - sin(2*PI*i*h))*h
enddo
