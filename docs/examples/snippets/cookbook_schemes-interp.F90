call seth_weno5_js(v, v, SETH_WENO_JS_EPS, 1._R_P, rp, rm)          ! wrong data: it expects cell averages
e_rec(k) = abs(rp - sin(2*PI*xf))
call seth_weno5_js_interp(v, v, SETH_WENO_JS_EPS, 1._R_P, rp, rm)   ! point values: the interpolation kernel
e_int(k) = abs(rp - sin(2*PI*xf))
