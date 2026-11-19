subroutine residual(u, res)
!< res = -(F_{i+1/2} - F_{i-1/2})/h, every loop on the device; the kernel called inside the face loop. u and res are
!< device memory (FUNDAL's and SELAH's: deviceptr), the work arrays live in the data region (present).
real(R_P), intent(in)  :: u(N)
real(R_P), intent(out) :: res(N)
real(R_P)              :: vp(-2:3), vm(-2:3), rp, rm
integer                :: i, m
!$acc parallel loop deviceptr(u) present(ug)
do i=1-NG, N+NG                              ! periodic ghost cells
   ug(i) = u(modulo(i - 1, N) + 1)
enddo
!$acc parallel loop present(ug, fp, fm)
do i=1-NG, N+NG
   fp(i) = 0.5_R_P*(A*ug(i) + abs(A)*ug(i))
   fm(i) = 0.5_R_P*(A*ug(i) - abs(A)*ug(i))
enddo
!$acc parallel loop present(fp, fm, f) private(vp, vm, rp, rm)
do i=0, N
   do m=-2, 3                                ! the face's stencils, in arrays of constant bounds
      vp(m) = fp(i+m) ; vm(m) = fm(i+m)
   enddo
   call seth_teno5(vp, vm, SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)   ! acc routine seq: one thread per face
   f(i) = rp + rm
enddo
!$acc parallel loop present(f) deviceptr(res)
do i=1, N
   res(i) = -(f(i) - f(i-1))/h
enddo
endsubroutine residual
