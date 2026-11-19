! the field in device memory (FUNDAL: on the GPU with OpenACC, on the host without), and SELAH's registers with it:
! one copy in, one copy out
call dev_alloc(fptr_dev=u, ubounds=[N], lbounds=[1], ierr=ierr)
call dev_memcpy_to_device(dst=u, src=u0)
bx = box_interior(lb=[1], ub=[N], ngc=0, components=.false.)
call st%initialize(SCHEME_SSP_33, FORM_SHU_OSHER_2R, lb=[1], ub=[N], components=.false.)
!$acc data create(ug, fp, fm, f)
t = 0 ; steps = 0
do while (t < 1._R_P)
   dt = min(CFL*h/abs(A), 1._R_P - t)
   call st%begin_step(bx, u)
   do k=1, st%stages()                           ! SELAH's SSP(3,3): the stage state is u itself
      call st%begin_stage(k, bx, dt, u)
      call st%rhs_register_view(r)
      call residual(u, r)
      call st%end_stage(k, bx, dt, u)
   enddo
   call st%finish_step(bx, dt, u)
   t = t + dt ; steps = steps + 1
enddo
!$acc end data
call dev_memcpy_from_device(dst=uh, src=u)
