call seth_teno5(vp(-2:3), vm(-2:3), SETH_TENO_EPS, SETH_TENO5_CT, sp, sm)  ! TENO5's two states
do k=-3, 3                                           ! the masks of the cells -3..3, both states
   do m=-2, 2
      w(m) = vp(k+m)                                 ! the left state: cell k, upwind to the left
   enddo
   call seth_teno5_mask(w, SETH_TENO_EPS, SETH_TENO5_CT, mp(k))
   do m=-2, 2
      w(m) = vm(1-k-m)                               ! the right state: the mirrored cell 1 - k
   enddo
   call seth_teno5_mask(w, SETH_TENO_EPS, SETH_TENO5_CT, mm(k))
enddo
call seth_dr_detect(mp, mm, lp, lm)                  ! the levels: 0 (none), 1 (Eq. 24), 2 (Eq. 32)
do m=-1, 2
   c(m) = vp(m) ; d(m) = vm(m)
enddo
call seth_dr_override(c, d, lp, lm, SETH_DR_T1, SETH_DR_T2, 1._R_P, sp, sm)  ! sp, sm: inout
