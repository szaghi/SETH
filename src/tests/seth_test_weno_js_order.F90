!< SETH test: the WENO-JS kernels converge at their design order 2S-1 on smooth data.
program seth_test_weno_js_order
!< SETH test: the WENO-JS kernels converge at their design order 2S-1 on smooth data.
!<
!< u = g(x), the generic function of seth_test_functions (not exp, whose averages are a constant multiple of its point
!< values), no critical point near the face (where Jiang-Shu weights lose order); the face at x = 0.3, cells of width h
!< halved from h0. Reconstruction from exact cell averages and interpolation from point values, both states
!< (vp = vm = u).
!< The observed order is taken on the finest pair of h whose errors stay above round-off (200 eps) and must reach
!< 2S-1 - 0.5.
use seth
use seth_test_functions, only : gfun, gavg
use seth_vv_cases, only : vv_report, vv_open
implicit none
integer(I4P), parameter :: NH = 6_I4P            !< Number of grid sizes.
real(R_P),    parameter :: XF = 0.3_R_P          !< Face position.
real(R_P),    parameter :: H0(5) = [0.4_R_P, 0.4_R_P, 0.4_R_P, 0.4_R_P, 0.8_R_P] !< Coarsest h per S.
real(R_P), allocatable  :: vp(:)
real(R_P)               :: h, err(2,NH,2), order, floor_
integer(I4P)            :: S, j, m, t, side
logical                 :: passed, interp
integer                 :: ur

passed = .true.
if (vv_report()) call vv_open('orders-js.dat', ur, 'S interp state observed-order design')
floor_ = 2.e2_R_P * epsilon(1._R_P) ! 100 eps |g|, |g| < 2 near XF
do S=1, 5
   if (allocated(vp)) deallocate(vp)
   allocate(vp(1-S:S))
   do t=1, 2
      interp = t == 2
      do j=1, NH
         h = H0(S) / 2._R_P**(j-1)
         do m=1-S, S ! cell m spans [XF + (m-1)h, XF + m h]
            if (interp) then
               vp(m) = gfun(XF + (m - 0.5_R_P)*h)
            else
               vp(m) = gavg(XF + (m-1)*h, XF + m*h)
            endif
         enddo
         call kernel(S, interp, err(1,j,t), err(2,j,t))
         err(:,j,t) = abs(err(:,j,t) - gfun(XF))
      enddo
      do side=1, 2
         order = -1._R_P
         do j=NH-1, 1, -1
            if (err(side,j+1,t) > floor_) then
               order = log(err(side,j,t)/err(side,j+1,t)) / log(2._R_P)
               exit
            endif
         enddo
         print '(A,I0,A,L1,A,A,A,F6.2,A,I0)', 'S=', S, ' interp=', interp, ' state=', merge('+','-',side==1), &
               ' observed order=', order, ' design=', 2*S-1
         if (vv_report()) write(ur, '(I2,1X,L1,1X,A1,F8.3,I3)') S, interp, merge('+','-',side==1), order, 2*S-1
         if (order < 2*S - 1 - 0.5_R_P) passed = .false.
      enddo
   enddo
enddo
if (vv_report()) close(ur)
if (.not.passed) error stop 'a WENO-JS kernel misses its design order'
print '(A)', 'seth_test_weno_js_order: ok'

contains
   subroutine kernel(S, interp, rp, rm)
   integer(I4P), intent(in)  :: S
   logical,      intent(in)  :: interp
   real(R_P),    intent(out) :: rp, rm
   real(R_P), parameter      :: E = SETH_WENO_JS_EPS
   select case(S)
   case(1) ; if (interp) then ; call seth_weno1_js_interp(vp, vp, E, 1._R_P, rp, rm)
             else ; call seth_weno1_js(vp, vp, E, 1._R_P, rp, rm) ; endif
   case(2) ; if (interp) then ; call seth_weno3_js_interp(vp, vp, E, 1._R_P, rp, rm)
             else ; call seth_weno3_js(vp, vp, E, 1._R_P, rp, rm) ; endif
   case(3) ; if (interp) then ; call seth_weno5_js_interp(vp, vp, E, 1._R_P, rp, rm)
             else ; call seth_weno5_js(vp, vp, E, 1._R_P, rp, rm) ; endif
   case(4) ; if (interp) then ; call seth_weno7_js_interp(vp, vp, E, 1._R_P, rp, rm)
             else ; call seth_weno7_js(vp, vp, E, 1._R_P, rp, rm) ; endif
   case(5) ; if (interp) then ; call seth_weno9_js_interp(vp, vp, E, 1._R_P, rp, rm)
             else ; call seth_weno9_js(vp, vp, E, 1._R_P, rp, rm) ; endif
   endselect
   endsubroutine kernel
endprogram seth_test_weno_js_order
