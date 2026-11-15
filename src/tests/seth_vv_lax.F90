!< SETH verification: the Lax shock tube, against the exact solution.
program seth_vv_lax
!< SETH verification: the Lax shock tube, against the exact solution.
!<
!< Lax (1954): (rho, u, p) = (0.445, 0.698, 3.528) left of x = 0.5, (0.5, 0, 0.571) right of it, gamma = 1.4, t = 0.13
!< (Jiang & Shu's t = 1.3 on [-5, 5]): a strong contact behind a shock, where low-dissipation schemes oscillate first;
!< the exact solution by the exact Riemann solver (seth_vv_riemann, Toro 2009).
!<
!< Both paths (the splitting path with every reconstruction, the interpolation path with every interpolation), N = 200
!< cells on [0, 1], transmissive ends, CFL 0.5. Every scheme runs without any positivity safeguard (recording whether
!< rho and p stay positive) and with the positivity-preserving flux limiter of Hu, Adams & Shu (2013) (counting the
!< limited faces). Asserted: with and without the limiter positivity throughout and no limited face; L1 of the density
!< below the recorded bounds of the quick schemes (measured + 5%); WENO-Z5 and TENO5 closer to the exact solution than
!< WENO5-JS on each path. Report mode: every scheme, the profiles, the exact solution, L1 on 100 .. 1600 cells, and
!< lax-summary.dat (per scheme: positive without the limiter, limited faces, L1).
use seth
use seth_vv_euler
use seth_vv_cases
implicit none
integer(I4P), allocatable :: ids(:), nfb(:), nfbu(:), nlim(:)
real(R8P), allocatable    :: l1(:), l1u(:)
logical, allocatable      :: pos(:), posu(:)
integer(I4P)              :: path, k, us
real(R8P)                 :: bound, js5
logical                   :: passed
real(R8P), parameter      :: WL(3) = [0.445_R8P, 0.698_R8P, 3.528_R8P]
real(R8P), parameter      :: WR(3) = [0.5_R8P, 0._R8P, 0.571_R8P]
real(R8P), parameter      :: X0 = 0.5_R8P, TEND = 0.13_R8P
passed = .true.
if (vv_report()) call vv_open('lax-summary.dat', us, 'scheme positive-without-limiter limited-faces L1(rho)')
do path=PATH_SPLIT, PATH_INTERP
   call vv_schemes(path, .not.vv_report(), ids)
   allocate(l1(size(ids)), pos(size(ids)), l1u(size(ids)), posu(size(ids)), nfb(size(ids)), nfbu(size(ids)), &
            nlim(size(ids)))
   limiter = .false.
   call riemann_problem('lax', WL, WR, X0, TEND, 0.5_R8P, 200, ids, path, l1u, posu, nfbu, files=.false.)
   limiter = .true.
   call riemann_problem('lax', WL, WR, X0, TEND, 0.5_R8P, 200, ids, path, l1, pos, nfb, nlim)
   limiter = .false.
   js5 = -1._R8P
   do k=1, size(ids)
      if (vv_label(ids(k)) == 'weno5_js' .or. vv_label(ids(k)) == 'weno5_js_interp') js5 = l1(k)
   enddo
   do k=1, size(ids)
      print '(A18,2ES12.4,A,L2,A,I6)', vv_label(ids(k)), l1(k), l1u(k), '   positive without limiter:', posu(k), &
            '   limited faces:', nlim(k)
      if (vv_report()) write(us, '(A,1X,L1,1X,I8,ES24.15)') vv_label(ids(k)), posu(k), nlim(k), l1(k)
      if (.not.pos(k)) call fail(vv_label(ids(k))//': not positive with the limiter')
      if (.not.posu(k)) call fail(vv_label(ids(k))//': not positive without the limiter')
      if (nlim(k) > 0) call fail(vv_label(ids(k))//': the limiter acted')
      select case(vv_label(ids(k)))
      case('wenoz5', 'teno5', 'wenoz5_interp', 'teno5_interp')
         if (.not.(l1(k) < js5)) call fail(vv_label(ids(k))//': not closer to the exact solution than WENO5-JS')
      end select
      bound = recorded(vv_label(ids(k)))
      if (.not.vv_report() .and. bound > 0._R8P .and. l1(k) > bound) call fail(vv_label(ids(k))//': L1 above bound')
   enddo
   deallocate(l1, pos, l1u, posu, nfb, nfbu, nlim)
enddo
if (vv_report()) close(us)
if (.not.passed) error stop 'lax verification failed'
print '(A)', 'seth_vv_lax: ok'

contains
   pure function recorded(name) result(b)
   !< L1 bounds of the test run (N = 200, with the limiter): measured (gfortran -O2) + 5%; -1: none.
   character(*), intent(in) :: name
   real(R8P)                :: b
   select case(name)
   case('weno5_js') ; b = 1.3929e-02_R8P ! 1.3266e-02
   case('wenoz5') ; b = 9.7912e-03_R8P ! 9.3250e-03
   case('teno5') ; b = 9.6596e-03_R8P ! 9.1996e-03
   case('teno6') ; b = 8.3436e-03_R8P ! 7.9463e-03
   case('weno5_js_interp') ; b = 1.1004e-02_R8P ! 1.0480e-02
   case('wenoz5_interp') ; b = 8.2395e-03_R8P ! 7.8471e-03
   case('teno5_interp') ; b = 8.2674e-03_R8P ! 7.8737e-03
   case('teno6_interp') ; b = 7.6446e-03_R8P ! 7.2806e-03
   case default ; b = -1._R8P
   end select
   endfunction recorded

   subroutine fail(what)
   character(*), intent(in) :: what
   print '(A)', '   FAILED: '//what
   passed = .false.
   endsubroutine fail
endprogram seth_vv_lax
