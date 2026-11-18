pure subroutine eigenvectors(ql, qr, l, r)
!< The left (rows of l) and right (columns of r) eigenvectors of the flux Jacobian at the Roe average of ql, qr.
real(R_P), intent(in)  :: ql(3), qr(3)
real(R_P), intent(out) :: l(3,3), r(3,3)
real(R_P)              :: wl(3), wr(3), sl, sr, u, hh, c, b1, b2
wl = cons2prim(ql) ; wr = cons2prim(qr)
sl = sqrt(wl(1)) ; sr = sqrt(wr(1))
u  = (sl*wl(2) + sr*wr(2))/(sl + sr)
hh = (sl*(ql(3) + wl(3))/wl(1) + sr*(qr(3) + wr(3))/wr(1))/(sl + sr)   ! enthalpy
c  = sqrt((GAMMA - 1)*(hh - 0.5_R_P*u*u))
r(:,1) = [1._R_P, u - c, hh - u*c]
r(:,2) = [1._R_P, u, 0.5_R_P*u*u]
r(:,3) = [1._R_P, u + c, hh + u*c]
b1 = (GAMMA - 1)/(c*c) ; b2 = 0.5_R_P*b1*u*u
l(1,:) = 0.5_R_P*[b2 + u/c, -b1*u - 1/c, b1]
l(2,:) = [1 - b2, b1*u, -b1]
l(3,:) = 0.5_R_P*[b2 - u/c, -b1*u + 1/c, b1]
endsubroutine eigenvectors
