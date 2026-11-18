use seth
implicit none
print '(A,I0,A,I0,A)', 'R_P: ', storage_size(1._R_P), ' bits, ', precision(1._R_P), ' decimal digits'
print '(A,ES10.3)', '  epsilon(1._R_P)   ', epsilon(1._R_P)
print '(A,ES10.3)', '  tiny(1._R_P)      ', tiny(1._R_P)
print '(A,ES10.3)', '  SETH_TENO_EPS     ', SETH_TENO_EPS
print '(A,ES10.3)', '  SETH_WENO_JS_EPS  ', SETH_WENO_JS_EPS
