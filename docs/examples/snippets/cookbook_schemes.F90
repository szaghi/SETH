!< Cookbook: recipes on the reconstruction kernels, the recipe chosen on the command line.
program cookbook_schemes
use seth
implicit none
real(R_P), parameter :: PI = acos(-1._R_P)
character(16)        :: recipe

call get_command_argument(1, recipe)
select case(trim(recipe))
case('interp') ; call point_values
case('flags')  ; call shock_flags
case('mp')     ; call mp_bound
case default   ; print '(A)', 'usage: cookbook_schemes interp|flags|mp'
endselect

contains
   subroutine point_values
   !< The face value of sin(2 pi x) from point values at the cell centres: interpolation, not reconstruction.
   real(R_P) :: v(-2:3), h, xf, rp, rm, e_rec(5), e_int(5)
   integer   :: k, m
   print '(A)', 'the face value of sin(2 pi x) from its point values at the cell centres'
   print '(A)', '        h  reconstruction  order  interpolation  order'
   h = 0.05_R_P ; xf = 0.1_R_P
   do k=1, 5
      do m=-2, 3                                ! the point values at the centres of the cells around the face
         v(m) = sin(2*PI*(xf + (m - 0.5_R_P)*h))
      enddo
      call seth_weno5_js(v, v, SETH_WENO_JS_EPS, 1._R_P, rp, rm)          ! wrong data: it expects cell averages
      e_rec(k) = abs(rp - sin(2*PI*xf))
      call seth_weno5_js_interp(v, v, SETH_WENO_JS_EPS, 1._R_P, rp, rm)   ! point values: the interpolation kernel
      e_int(k) = abs(rp - sin(2*PI*xf))
      if (k == 1) then
         print '(F9.5,ES16.3,7X,ES15.3)', h, e_rec(k), e_int(k)
      else
         print '(F9.5,ES16.3,F7.2,ES15.3,F7.2)', h, e_rec(k), log(e_rec(k-1)/e_rec(k))/log(2._R_P), &
                                                 e_int(k), log(e_int(k-1)/e_int(k))/log(2._R_P)
      endif
      h = h/2
   enddo
   endsubroutine point_values

   subroutine shock_flags
   !< The weights ratios of the _wr kernels as a free shock sensor: the faces they flag on a sine with a step.
   integer, parameter :: N = 40
   real(R_P)          :: u(-2:N+3), h, rp, rm, wrp, wrm
   character(N+1)     :: line
   integer            :: i
   h = 1._R_P/N
   do i=1, N
      u(i) = (cos(2*PI*(i - 1)*h) - cos(2*PI*i*h))/(2*PI*h) + merge(2._R_P, 0._R_P, i > N/2)
   enddo
   u(-2:0) = u(N-2:N) ; u(N+1:N+3) = u(1:3)
   do i=0, N                                    ! a face is flagged where either weights ratio is below 1/2
      call seth_weno5_js_wr(u(i-2:i+3), u(i-2:i+3), SETH_WENO_JS_EPS, 1._R_P, rp, rm, wrp, wrm)
      line(i+1:i+1) = merge('|', '.', min(wrp, wrm) < 0.5_R_P)
   enddo
   print '(A)', 'faces 0..40 of a sine plus a step at x = 0.5 (and its periodic image at x = 0): | flagged'
   print '(2X,A)', line
   endsubroutine shock_flags

   subroutine mp_bound
   !< The ninth-order linear scheme across a step, then the monotonicity-preserving bound on its states; and the bound
   !< on smooth data, where it must not act.
   real(R_P) :: u(-4:5), rp, rm, bp, bm, h
   integer   :: j, m
   print '(A)', 'linear9 on a step, the step moved across the stencil'
   print '(A)', '  step after cell    linear9: rp       rm    bounded: rp       rm'
   do j=-3, 4
      u = merge(1._R_P, 0._R_P, [(m, m=-4, 5)] <= j)
      call seth_linear9(u, u, rp, rm)                                   ! any reconstruction kernel
      bp = rp ; bm = rm
      call seth_mp_bound(u(-2:3), u(-2:3), SETH_MP_ALPHA, bp, bm)       ! the five cells around each upwind cell
      print '(I17,2F11.5,4X,2F11.5)', j, shown(rp), shown(rm), shown(bp), shown(bm)
   enddo
   h = 0.02_R_P                                 ! the cell averages of sin(2 pi x) around its maximum at x = 0.25
   do m=-4, 5
      u(m) = (cos(2*PI*(0.25_R_P + (m - 1)*h)) - cos(2*PI*(0.25_R_P + m*h)))/(2*PI*h)
   enddo
   call seth_linear9(u, u, rp, rm)
   bp = rp ; bm = rm
   call seth_mp_bound(u(-2:3), u(-2:3), SETH_MP_ALPHA, bp, bm)
   print '(A,L2)', 'at the smooth maximum of a sine, the bound leaves the states unchanged:', bp == rp .and. bm == rm
   endsubroutine mp_bound

   pure real(R_P) function shown(v)
   !< v to 5 decimals for printing, a negative zero printed as 0.
   real(R_P), intent(in) :: v
   shown = anint(v*1e5_R_P)/1e5_R_P + 0._R_P
   endfunction shown
endprogram cookbook_schemes
