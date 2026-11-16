!< SETH verification: the Sod shock tube, against the exact solution.
program seth_vv_sod
!< SETH verification: the Sod shock tube, against the exact solution.
!<
!< Sod (1978): (rho, u, p) = (1, 0, 1) left of x = 0.5, (0.125, 0, 0.1) right of it, gamma = 1.4, t = 0.2: a
!< rarefaction, a contact and a shock; the exact solution by the exact Riemann solver (seth_vv_riemann, Toro 2009).
!<
!< Both paths (the splitting path with every reconstruction, the interpolation path with every interpolation), N = 200
!< cells on [0, 1], transmissive ends, CFL 0.5. Every scheme runs without any positivity safeguard (recording whether
!< rho and p stay positive) and with the positivity-preserving flux limiter of Hu, Adams & Shu (2013) (counting the
!< limited faces). Asserted: with and without the limiter positivity throughout and no limited face (the limiter is
!< inactive on this problem); L1 of the density below the recorded bounds of the quick schemes (measured + 5%); the
!< low-dissipation claims of WENO-Z5 (Borges et al. 2008) and TENO5 (Fu, Hu & Adams 2016): closer to the exact solution
!< than WENO5-JS on each path. Report mode: every scheme, the profiles, the exact solution, L1 on 100 .. 1600 cells, and
!< sod-summary.dat (per scheme: positive without the limiter, limited faces, L1).
use seth
use seth_vv_euler
use seth_vv_cases
implicit none
integer(I4P), allocatable :: ids(:), nfb(:), nfbu(:), nlim(:)
real(R_P), allocatable    :: l1(:), l1u(:)
logical, allocatable      :: pos(:), posu(:)
integer(I4P)              :: path, k, us
real(R_P)                 :: bound, js5
logical                   :: passed
real(R_P), parameter      :: WL(3) = [1._R_P, 0._R_P, 1._R_P]
real(R_P), parameter      :: WR(3) = [0.125_R_P, 0._R_P, 0.1_R_P]
real(R_P), parameter      :: X0 = 0.5_R_P, TEND = 0.2_R_P
passed = .true.
if (vv_report()) call vv_open('sod-summary.dat', us, 'scheme positive-without-limiter limited-faces L1(rho)')
do path=PATH_SPLIT, PATH_INTERP
   call vv_schemes(path, .not.vv_report(), ids)
   allocate(l1(size(ids)), pos(size(ids)), l1u(size(ids)), posu(size(ids)), nfb(size(ids)), nfbu(size(ids)), &
            nlim(size(ids)))
   limiter = .false.
   call riemann_problem('sod', WL, WR, X0, TEND, 0.5_R_P, 200, ids, path, l1u, posu, nfbu, files=.false.)
   limiter = .true.
   call riemann_problem('sod', WL, WR, X0, TEND, 0.5_R_P, 200, ids, path, l1, pos, nfb, nlim)
   limiter = .false.
   js5 = -1._R_P
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
      if (.not.vv_report() .and. bound > 0._R_P .and. l1(k) > bound) call fail(vv_label(ids(k))//': L1 above bound')
   enddo
   deallocate(l1, pos, l1u, posu, nfb, nfbu, nlim)
enddo
if (vv_report()) close(us)
if (.not.passed) error stop 'sod verification failed'
print '(A)', 'seth_vv_sod: ok'

contains
   pure function recorded(name) result(b)
   !< L1 bounds of the test run (N = 200, with the limiter): measured (gfortran -O2) + 5%; -1: none.
   character(*), intent(in) :: name
   real(R_P)                :: b
   select case(name)
   case('weno5_js') ; b = 4.0102e-03_R_P ! 3.8192e-03
   case('wenoz5') ; b = 2.9696e-03_R_P ! 2.8282e-03
   case('teno5') ; b = 2.8670e-03_R_P ! 2.7305e-03
   case('teno6') ; b = 2.4527e-03_R_P ! 2.3359e-03
   case('weno5_js_interp') ; b = 2.9953e-03_R_P ! 2.8527e-03
   case('wenoz5_interp') ; b = 2.4393e-03_R_P ! 2.3231e-03
   case('teno5_interp') ; b = 2.3849e-03_R_P ! 2.2713e-03
   case('teno6_interp') ; b = 2.3867e-03_R_P ! 2.2730e-03
   case default ; b = -1._R_P
   end select
   endfunction recorded

   subroutine fail(what)
   character(*), intent(in) :: what
   print '(A)', '   FAILED: '//what
   passed = .false.
   endsubroutine fail
endprogram seth_vv_sod
