! cell averages of sin(2 pi x) on the six cells around the face: cell m is [XF + (m-1) H, XF + m H]
do m=-2, 3
   u(m) = (cos(2*PI*(XF + (m - 1)*H)) - cos(2*PI*(XF + m*H)))/(2*PI*H)
enddo
s = u                                           ! kept for the mirror, below
call seth_weno5_js_wr(u, u, SETH_WENO_JS_EPS, 1._R_P, rp, rm, wrp, wrm)
call seth_linear5(u, u, lp, lm)
