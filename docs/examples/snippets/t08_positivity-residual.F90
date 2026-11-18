subroutine residual(q, dt, limit, res, nlim)
!< res = -(F_{i+1/2} - F_{i-1/2})/h; with limit, every face flux limited by seth_positivity_euler.
real(R_P), intent(in)    :: q(3,N), dt
logical,   intent(in)    :: limit
real(R_P), intent(out)   :: res(3,N)
integer,   intent(inout) :: nlim
real(R_P)                :: qg(3,1-NG:N+NG), fg(3,1-NG:N+NG), wc(3,-2:3), gc(3,-2:3), l(3,3), r(3,3)
real(R_P)                :: alpha(3), fc(3), fhat(3,0:N), flow(3), wp(3), rp, rm, c, theta, a
integer                  :: i, k, m
qg(:,1:N) = q
do i=1, NG                                   ! reflecting walls: mirrored cells, the momentum reversed
   qg(:,1-i) = q(:,i) ; qg(2,1-i) = -q(2,i)
   qg(:,N+i) = q(:,N+1-i) ; qg(2,N+i) = -q(2,N+1-i)
enddo
do i=1-NG, N+NG
   fg(:,i) = flux(qg(:,i))
enddo
alpha = 0
do i=1, N
   wp = cons2prim(q(:,i)) ; c = sqrt(GAMMA*wp(3)/wp(1))
   alpha = max(alpha, abs([wp(2) - c, wp(2), wp(2) + c]))
enddo
a = maxval(alpha)                            ! the first-order flux's alpha: max(|u| + c)
alpha = 1.1_R_P*alpha
do i=0, N
   call eigenvectors(qg(:,i), qg(:,i+1), l, r)
   do m=-2, 3
      wc(:,m) = matmul(l, qg(:,i+m)) ; gc(:,m) = matmul(l, fg(:,i+m))
   enddo
   do k=1, 3
      call seth_teno5(0.5_R_P*(gc(k,:) + alpha(k)*wc(k,:)), 0.5_R_P*(gc(k,:) - alpha(k)*wc(k,:)), &
                      SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)
      fc(k) = rp + rm
   enddo
   fhat(:,i) = matmul(r, fc)
   if (limit) then
      ! the first-order Lax-Friedrichs flux, positive for 2 dt/h a <= 1; theta blends the two
      flow = 0.5_R_P*(fg(:,i) + fg(:,i+1)) - 0.5_R_P*a*(qg(:,i+1) - qg(:,i))
      call seth_positivity_euler(3, GAMMA, 2*dt/h, SETH_POSITIVITY_EPS, SETH_POSITIVITY_EPS, qg(:,i), qg(:,i+1), &
                                 fhat(:,i), flow, fc, theta)
      fhat(:,i) = fc
      if (theta < 1) nlim = nlim + 1
   endif
enddo
res = -(fhat(:,1:N) - fhat(:,0:N-1))/h
endsubroutine residual
