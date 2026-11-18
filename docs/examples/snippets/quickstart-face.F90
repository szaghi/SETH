! vp and vm are the six cells around the face, cell 0 the one left of it: u(i-2:i+3)
call seth_weno5_js(u(i-2:i+3), u(i-2:i+3), SETH_WENO_JS_EPS, 1._R_P, rp, rm)
call seth_linear5(u(i-2:i+3), u(i-2:i+3), lp, lm)
