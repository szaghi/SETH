!< Reference, non-polynomial reconstructions: every kernel called on one stencil, and TENO5-DR built from its pieces.
program reference_nonpoly
use seth
implicit none
real(R_P) :: vp(-5:6), vm(-5:6), rp, rm, sp, sm, h
real(R_P) :: w(-2:2), c(-1:2), d(-1:2)
integer   :: m, mk, mp(-3:3), mm(-3:3), lp, lm, k

h = 0.01_R_P                                         ! the grid spacing: the RBF kernels' s = alpha h
do m=-5, 6                                           ! a unit jump between the cells 0 and 1, on a slope
   vp(m) = 0.02_R_P*m + merge(1._R_P, 0._R_P, m >= 1)
   vm(m) = vp(m)                                     ! one scalar field: both biases read the same values
enddo
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
mk = mp(0)
print '(A,I2,A,I3,A,2I2)', 'cell 0: mask', mk, ', zeta', seth_teno5_zeta(mk), '; DR levels of the face', lp, lm
call seth_teno5_dr(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, SETH_DR_T1, SETH_DR_T2, 1._R_P, rp, rm)
print '(A,L2)', 'pieces equal to seth_teno5_dr, bitwise:', sp == rp .and. sm == rm
endprogram reference_nonpoly
