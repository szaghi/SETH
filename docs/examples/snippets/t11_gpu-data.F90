! the fields live on the device for the whole run: one copy in, one copy out
!$acc data copy(u) create(u1, u2, res, ug, fp, fm, f)
t = 0 ; steps = 0
do while (t < 1._R_P)
   dt = min(CFL*h/abs(A), 1._R_P - t)
   call residual(u, res)
   !$acc parallel loop present(u, u1, res)
   do i=1, N
      u1(i) = u(i) + dt*res(i)
   enddo
   call residual(u1, res)
   !$acc parallel loop present(u, u1, u2, res)
   do i=1, N
      u2(i) = 0.75_R_P*u(i) + 0.25_R_P*(u1(i) + dt*res(i))
   enddo
   call residual(u2, res)
   !$acc parallel loop present(u, u2, res)
   do i=1, N
      u(i) = u(i)/3 + 2._R_P/3*(u2(i) + dt*res(i))
   enddo
   t = t + dt ; steps = steps + 1
enddo
!$acc end data
