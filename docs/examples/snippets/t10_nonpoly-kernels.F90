select case(s)
case(1)
   call seth_teno5(fp(i-2:i+3), fm(i-2:i+3), SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)
case(2)                                   ! twelve cells, -5..6
   call seth_teno5_dr(fp(i-5:i+6), fm(i-5:i+6), SETH_TENO_EPS, SETH_TENO5_CT, SETH_DR_T1, SETH_DR_T2, &
                      1._R_P, rp, rm)
case(3)
   call seth_teno5_thinc(fp(i-5:i+6), fm(i-5:i+6), SETH_TENO_EPS, SETH_TENO5_CT, SETH_THINC_BETA, 1._R_P, &
                         rp, rm)
case(4)                                   ! the RBF shape parameters in the grid's units: s = 3 h
   call seth_teno5_rbf_dr(fp(i-5:i+6), fm(i-5:i+6), SETH_TENO_EPS, SETH_TENO5_CT, SETH_RBF_CR, &
                          SETH_RBF_ALPHA1*h, SETH_RBF_ALPHA2*h, SETH_DR_T1, SETH_DR_T2, 1._R_P, rp, rm)
endselect
