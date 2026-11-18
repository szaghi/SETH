do i=-3, NC+3                                        ! pass 1, at each face i+1/2: msl(i) of cell i, msr(i+1) of cell i+1
   call seth_teno5_face_masks(u(i-2:i+3), u(i-2:i+3), SETH_TENO_EPS, SETH_TENO5_CT, msl(i), msr(i+1))
enddo
do i=0, NC                                           ! pass 2, the faces 0..NC: mp(c) = msl(i+c), mm(c) = msr(i+1-c)
   do m=-3, 3
      mp(m) = msl(i+m) ; mm(m) = msr(i+1-m)
   enddo
   call seth_teno5_dr_stored(u(i-2:i+3), u(i-2:i+3), mp, mm, SETH_TENO_EPS, SETH_TENO5_CT, SETH_DR_T1, SETH_DR_T2, &
                             1._R_P, fp(i), fm(i))
enddo
