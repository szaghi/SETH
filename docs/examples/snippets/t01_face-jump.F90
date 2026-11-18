! a step between cells 1 and 2: the left state's stencil (cells -2..2) crosses it, the right state's (-1..3) too
u(-2:1) = 1._R_P
u(2:3)  = 0._R_P
call seth_weno5_js_wr(u, u, SETH_WENO_JS_EPS, 1._R_P, rp, rm, wrp, wrm)
call seth_linear5(u, u, lp, lm)
