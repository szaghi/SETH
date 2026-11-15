!< SETH verification: long-time advection of the Jiang-Shu composite profile (smooth, discontinuous and kinked waves).
program seth_vv_advection
!< SETH verification: long-time advection of the Jiang-Shu composite profile (smooth, discontinuous and kinked waves).
!<
!< u_t + u_x = 0 on [-1, 1], periodic, from the composite profile of Jiang & Shu (1996), example 2: a Gaussian pulse,
!< a square wave, a triangle and a half ellipse (constants as llm-wiki high-order-capturing-schemes
!< `jiang-shu-composite-advection.md`), 200 cells, CFL 0.4, SSP-RK3, to t = 8 (four periods): dissipation, dispersion
!< and the non-oscillatory property at once. Every reconstruction (upwind splitting, seth_vv_scalar).
!<
!< Asserted: the solution stays in [-0.05, 1.05] (no growing oscillation: WENO is not bounded exactly, but an overshoot
!< of 5% after four periods would be an instability); the
!< low-dissipation claims of the literature on this problem: WENO-Z5 (Borges et al. 2008) and TENO5 (Fu, Hu & Adams
!< 2016) have a smaller L1 error than WENO5-JS.
!< Report mode: advection.dat (scheme L1 max min), the final profiles, 81 frames over the four periods (WENO5-JS,
!< WENO-Z5, TENO5, TENO6).
use seth
use seth_test_schemes
use seth_vv_scalar
use seth_vv_cases, only : vv_report, vv_open, vv_schemes, vv_label, l1_norm
use seth_vv_euler, only : PATH_SPLIT
implicit none
integer(I4P), parameter   :: N = 200_I4P, NF = 80_I4P
real(R8P), parameter      :: TEND = 8._R8P, CFL = 0.4_R8P
integer(I4P), allocatable :: ids(:)
real(R8P)                 :: u0(0:N-1), u(0:N-1), h, l1, js5, fr(0:N-1,0:NF)
real(R8P), allocatable    :: l1s(:)
integer(I4P)              :: k, i, ua, uf, f
logical                   :: passed, frames

passed = .true.
h = 2._R8P/N
do i=0, N-1
   u0(i) = jiang_shu(-1._R8P + (i + 0.5_R8P)*h)
enddo
call vv_schemes(PATH_SPLIT, .not.vv_report(), ids)
allocate(l1s(size(ids)))
if (vv_report()) call vv_open('advection.dat', ua, 'scheme L1 max min (composite profile, N = 200, t = 8)')
js5 = -1._R8P
do k=1, size(ids)
   frames = vv_report() .and. any(vv_label(ids(k)) == [character(8) :: 'weno5_js', 'wenoz5', 'teno5', 'teno6'])
   u = u0
   call advect(ids(k), u, frames)
   l1 = l1_norm(u - u0, h) ; l1s(k) = l1
   if (vv_label(ids(k)) == 'weno5_js') js5 = l1
   print '(A16,3ES12.4)', vv_label(ids(k)), l1, maxval(u), minval(u)
   if (.not.(maxval(u) <= 1.05_R8P .and. minval(u) >= -0.05_R8P)) call fail(vv_label(ids(k))//': out of [-0.05, 1.05]')
   if (vv_report()) then
      write(ua, '(A,1X,3ES24.15)') vv_label(ids(k)), l1, maxval(u), minval(u)
      call vv_open('advection-'//vv_label(ids(k))//'.dat', uf, 'x u u_exact (t = 8)')
      do i=0, N-1
         write(uf, '(3ES24.15)') -1._R8P + (i + 0.5_R8P)*h, u(i), u0(i)
      enddo
      close(uf)
      if (frames) then
         call vv_open('advection-frames-'//vv_label(ids(k))//'.dat', uf, &
                      'frames: line 1 x(1..N); then per frame f = 0..80 (t = f/10): t u(1..N) u_exact(1..N)')
         write(uf, '(*(ES16.8))') (-1._R8P + (i + 0.5_R8P)*h, i=0, N-1)
         do f=0, NF
            write(uf, '(*(ES16.8))') TEND*f/NF, fr(:,f), &
               (jiang_shu(modulo(-1._R8P + (i + 0.5_R8P)*h - TEND*f/NF + 1._R8P, 2._R8P) - 1._R8P), i=0, N-1)
         enddo
         close(uf)
      endif
   endif
enddo
do k=1, size(ids)
   if (any(vv_label(ids(k)) == [character(8) :: 'wenoz5', 'teno5']) .and. .not.(l1s(k) < js5)) &
      call fail(vv_label(ids(k))//': not less dissipative than WENO5-JS')
enddo
if (vv_report()) close(ua)
if (.not.passed) error stop 'advection verification failed'
print '(A)', 'seth_vv_advection: ok'

contains
   subroutine advect(id, u, frames)
   integer(I4P), intent(in)    :: id
   real(R8P),    intent(inout) :: u(0:N-1)
   logical,      intent(in)    :: frames
   real(R8P)                   :: t, dt, tf
   integer(I4P)                :: kf
   t = 0._R8P ; kf = 1
   if (frames) fr(:,0) = u
   do while (t < TEND*(1._R8P - 1.e-14_R8P))
      dt = min(CFL*h, TEND - t)
      if (frames) then
         tf = TEND*kf/NF ; dt = min(dt, tf - t)
      endif
      call scalar_step(id, FLUX_ADVECTION, u, N, h, dt)
      t = t + dt
      if (frames .and. kf <= NF) then
         if (t >= (TEND*kf/NF)*(1._R8P - 1.e-14_R8P)) then
            fr(:,kf) = u ; kf = kf + 1
         endif
      endif
   enddo
   endsubroutine advect

   subroutine fail(what)
   character(*), intent(in) :: what
   print '(A)', '   FAILED: '//what
   passed = .false.
   endsubroutine fail

   pure function jiang_shu(x) result(u)
   !< The composite profile, constants as llm-wiki `jiang-shu-composite-advection.md` prints them.
   real(R8P), intent(in) :: x
   real(R8P)             :: u
   real(R8P), parameter  :: A = 0.5_R8P, Z = -0.7_R8P, DG = 0.005_R8P, AF = 10._R8P
   real(R8P)             :: bg
   bg = log(2._R8P)/(36._R8P*DG*DG)
   if (x >= -0.8_R8P .and. x <= -0.6_R8P) then
      u = (g(x, bg, Z - DG) + g(x, bg, Z + DG) + 4._R8P*g(x, bg, Z))/6._R8P
   elseif (x >= -0.4_R8P .and. x <= -0.2_R8P) then
      u = 1._R8P
   elseif (x >= 0._R8P .and. x <= 0.2_R8P) then
      u = 1._R8P - abs(10._R8P*(x - 0.1_R8P))
   elseif (x >= 0.4_R8P .and. x <= 0.6_R8P) then
      u = (ef(x, AF, A - DG) + ef(x, AF, A + DG) + 4._R8P*ef(x, AF, A))/6._R8P
   else
      u = 0._R8P
   endif
   endfunction jiang_shu

   pure function g(x, b, c) result(y)
   real(R8P), intent(in) :: x, b, c
   real(R8P)             :: y
   y = exp(-b*(x - c)**2)
   endfunction g

   pure function ef(x, al, c) result(y)
   real(R8P), intent(in) :: x, al, c
   real(R8P)             :: y
   y = sqrt(max(1._R8P - al*al*(x - c)**2, 0._R8P))
   endfunction ef
endprogram seth_vv_advection
