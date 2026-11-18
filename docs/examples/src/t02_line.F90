!< Tutorial 2: a line of faces. Every face of a periodic grid, and the order of accuracy it gives.
program t02_line
!run t02_line t02_line
use seth
implicit none
real(R_P), parameter :: PI = acos(-1._R_P)
integer,   parameter :: NG = 3                  ! ghost cells: WENO5 reads cells -2..3 around a face
real(R_P)            :: err(3), err_old(3)
integer              :: n, k, j

print '(A)', '    N     WENO5-JS  order     WENO-Z5  order       TENO5  order'
err_old = 0
n = 10
do k=1, 6
   n = 2*n
   err = l1_error(n)
   if (k == 1) then
      print '(I5,3(ES12.3,7X))', n, err
   else
      print '(I5,3(ES12.3,F7.2))', n, (err(j), log(err_old(j)/err(j))/log(2._R_P), j=1, 3)
   endif
   err_old = err
enddo

contains
   function l1_error(n) result(err)
   !< The L1 error of the left states at the faces of n cells on [0, 1], for each scheme.
   integer, intent(in) :: n
   real(R_P)           :: err(3)
   real(R_P)           :: u(1-NG:n+NG), h, rp, rm
   integer             :: i
   !region fill
   h = 1._R_P/n
   do i=1, n                                     ! cell i is [(i-1) h, i h]: its average of sin(2 pi x)
      u(i) = (cos(2*PI*(i - 1)*h) - cos(2*PI*i*h))/(2*PI*h)
   enddo
   u(1-NG:0)   = u(n+1-NG:n)                     ! periodic ghost cells
   u(n+1:n+NG) = u(1:NG)
   !endregion fill
   err = 0
   !region faces
   do i=0, n-1                                   ! the face x_{i+1/2} = i h, between cells i and i+1
      call seth_weno5_js(u(i-2:i+3), u(i-2:i+3), SETH_WENO_JS_EPS, 1._R_P, rp, rm)
      err(1) = err(1) + abs(rp - sin(2*PI*i*h))*h
      call seth_wenoz5(u(i-2:i+3), u(i-2:i+3), SETH_TENO_EPS, 1, rp, rm)
      err(2) = err(2) + abs(rp - sin(2*PI*i*h))*h
      call seth_teno5(u(i-2:i+3), u(i-2:i+3), SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)
      err(3) = err(3) + abs(rp - sin(2*PI*i*h))*h
   enddo
   !endregion faces
   endfunction l1_error
endprogram t02_line
