do i=0, N                                    ! a face is flagged where either weights ratio is below 1/2
   call seth_weno5_js_wr(u(i-2:i+3), u(i-2:i+3), SETH_WENO_JS_EPS, 1._R_P, rp, rm, wrp, wrm)
   line(i+1:i+1) = merge('|', '.', min(wrp, wrm) < 0.5_R_P)
enddo
