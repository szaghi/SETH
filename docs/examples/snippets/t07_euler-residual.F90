subroutine residual(q, res, characteristic)
!< res = -(F_{i+1/2} - F_{i-1/2})/h: Lax-Friedrichs splitting per field, TENO5, transmissive ends.
real(R_P), intent(in)  :: q(3,N)
real(R_P), intent(out) :: res(3,N)
logical,   intent(in)  :: characteristic
real(R_P)              :: qg(3,1-NG:N+NG), fg(3,1-NG:N+NG), wc(3,-2:3), gc(3,-2:3), l(3,3), r(3,3)
real(R_P)              :: alpha(3), fc(3), fhat(3,0:N), wp(3), rp, rm, c
integer                :: i, k, m
qg(:,1:N) = q
do i=1, NG
   qg(:,1-i) = q(:,1) ; qg(:,N+i) = q(:,N)   ! transmissive: copies of the end cells
enddo
do i=1-NG, N+NG
   fg(:,i) = flux(qg(:,i))
enddo
alpha = 0                                    ! alpha_k = 1.1 max |lambda_k| over the domain
do i=1, N
   wp = cons2prim(q(:,i)) ; c = sqrt(GAMMA*wp(3)/wp(1))
   alpha = max(alpha, abs([wp(2) - c, wp(2), wp(2) + c]))
enddo
alpha = 1.1_R_P*alpha
if (.not.characteristic) alpha = maxval(alpha)
do i=0, N
   if (characteristic) then                  ! project the stencil on the face's characteristic fields
      call eigenvectors(qg(:,i), qg(:,i+1), l, r)
      do m=-2, 3
         wc(:,m) = matmul(l, qg(:,i+m)) ; gc(:,m) = matmul(l, fg(:,i+m))
      enddo
   else                                      ! the conservative variables themselves
      wc = qg(:,i-2:i+3) ; gc = fg(:,i-2:i+3)
   endif
   do k=1, 3                                 ! each field: split, reconstruct, sum
      call seth_teno5(0.5_R_P*(gc(k,:) + alpha(k)*wc(k,:)), 0.5_R_P*(gc(k,:) - alpha(k)*wc(k,:)), &
                      SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)
      fc(k) = rp + rm
   enddo
   if (characteristic) then
      fhat(:,i) = matmul(r, fc)              ! back to the conservative variables
   else
      fhat(:,i) = fc
   endif
enddo
res = -(fhat(:,1:N) - fhat(:,0:N-1))/h
endsubroutine residual
