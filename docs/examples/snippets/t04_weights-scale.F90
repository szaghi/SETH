call seth_weno5_js(u, u, SETH_WENO_JS_EPS, 1._R_P, rp, rm)        ! Jiang and Shu's weights
o1 = over(rm, scale, scale)
rmu = seth_weno_descaling(3, 1._R_P, u, u)                         ! 1/mean|v| of the stencils
call seth_weno5_js(u, u, SETH_WENO_JS_EPS, rmu, rp, rm)
