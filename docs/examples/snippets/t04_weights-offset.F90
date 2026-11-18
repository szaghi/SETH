mean = sum(u)/size(u)                           ! subtract the stencil's mean, descale, add the mean back
rmu = seth_weno_descaling_ai(3, u - mean, u - mean)
call seth_weno5_js(u - mean, u - mean, SETH_WENO_JS_EPS, rmu, rp, rm)
rp = rp + mean ; rm = rm + mean
