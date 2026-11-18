call seth_teno5_dr(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, SETH_DR_T1, SETH_DR_T2, 1._R_P, rp, rm)
print '(A16,2F10.6)', 'TENO5-DR', rp, rm
call seth_teno5_thinc(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, SETH_THINC_BETA, 1._R_P, rp, rm)
print '(A16,2F10.6)', 'TENO5-THINC', rp, rm
call seth_teno5_st(vp(-2:3), vm(-2:3), SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)
print '(A16,2F10.6)', 'TENO5-ST', rp, rm
call seth_teno5_rbf(vp(-2:3), vm(-2:3), SETH_TENO_EPS, SETH_TENO5_CT, SETH_RBF_CR, SETH_RBF_ALPHA1*h, &
                    SETH_RBF_ALPHA2*h, rp, rm)
print '(A16,2F10.6)', 'TENO5-RBF', rp, rm
call seth_teno5_rbf_dr(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, SETH_RBF_CR, SETH_RBF_ALPHA1*h, SETH_RBF_ALPHA2*h, &
                       SETH_DR_T1, SETH_DR_T2, 1._R_P, rp, rm)
print '(A16,2F10.6)', 'TENO5-RBF+DR', rp, rm
