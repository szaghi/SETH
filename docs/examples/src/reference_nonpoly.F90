!< Reference, non-polynomial reconstructions: every kernel called on one stencil, TENO5-DR built from its pieces, and
!< the stored-mask pair on a field.
program reference_nonpoly
!run reference-nonpoly reference_nonpoly
use seth
implicit none
real(R_P) :: vp(-5:6), vm(-5:6), rp, rm, sp, sm, h
real(R_P) :: w(-2:2), c(-1:2), d(-1:2)
integer   :: m, mk, mp(-3:3), mm(-3:3), lp, lm, k
integer, parameter :: NC = 16                        ! the stored-mask pair: a field of NC cells and its halo
real(R_P) :: u(-5:NC+6)
real(R_P) :: fp(0:NC), fm(0:NC)
integer   :: msl(-3:NC+3), msr(-2:NC+4), i, nequal

h = 0.01_R_P                                         ! the grid spacing: the RBF kernels' s = alpha h
do m=-5, 6                                           ! a unit jump between the cells 0 and 1, on a slope
   vp(m) = 0.02_R_P*m + merge(1._R_P, 0._R_P, m >= 1)
   vm(m) = vp(m)                                     ! one scalar field: both biases read the same values
enddo
!region composites
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
!endregion composites

!region pieces
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
!endregion pieces
mk = mp(0)
print '(A,I2,A,I3,A,2I2)', 'cell 0: mask', mk, ', zeta', seth_teno5_zeta(mk), '; DR levels of the face', lp, lm
call seth_teno5_dr(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, SETH_DR_T1, SETH_DR_T2, 1._R_P, rp, rm)
print '(A,L2)', 'pieces equal to seth_teno5_dr, bitwise:', sp == rp .and. sm == rm

do i=-5, NC+6                                        ! the same jump on a field, between the cells 8 and 9
   u(i) = 0.02_R_P*i + merge(1._R_P, 0._R_P, i >= 9)
enddo
!region stored
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
!endregion stored
nequal = 0
do i=0, NC
   call seth_teno5_dr(u(i-5:i+6), u(i-5:i+6), SETH_TENO_EPS, SETH_TENO5_CT, SETH_DR_T1, SETH_DR_T2, 1._R_P, rp, rm)
   if (fp(i) == rp .and. fm(i) == rm) nequal = nequal + 1
enddo
print '(A,I3,A,I3,A)', 'stored-mask pair equal to seth_teno5_dr, bitwise, on', nequal, ' of', NC + 1, ' faces'
endprogram reference_nonpoly
