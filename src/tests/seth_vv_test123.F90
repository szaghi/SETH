!< SETH verification: the 123 problem (two strong rarefactions, near vacuum), against the exact solution.
program seth_vv_test123
!< SETH verification: the 123 problem (two strong rarefactions, near vacuum), against the exact solution.
!<
!< Einfeldt et al. (1991), Toro's test 2: (rho, u, p) = (1, -2, 0.4) left of x = 0.5, (1, 2, 0.4) right of it, gamma =
!< 1.4, t = 0.15: two rarefactions leave a near-vacuum region (rho, p ~ 0.02) in the middle, the standard positivity
!< stress test; the exact solution by the exact Riemann solver (seth_vv_riemann, Toro 2009).
!<
!< Both paths (the splitting path with every reconstruction, the interpolation path with every interpolation), N = 200
!< cells on [0, 1], transmissive ends, CFL 0.5. Every scheme runs without any positivity safeguard (recording whether
!< rho and p stay positive) and with the positivity-preserving flux limiter of Hu, Adams & Shu (2013) (counting the
!< limited faces). Asserted: with the limiter positivity throughout; without it, every reconstruction up to 7th order
!< on the splitting path and the WENO3-JS, WENO5-JS and WENO-Z3 interpolations stay positive (measured: the
!< WENO7/9-JS interpolations and the higher-order WENO-Z and TENO interpolations lose positivity without the limiter,
!< the finding documented on the case page; WENO-Z9 on the splitting path stays positive with p = 2); L1 of the
!< density below the recorded bounds of the
!< quick schemes (measured + 5%).
!< Report mode: every scheme, the profiles, the exact solution, L1 on 100 .. 1600 cells, and 123-summary.dat (per
!< scheme: positive without the limiter, limited faces, L1).
use seth
use seth_vv_euler
use seth_vv_cases
use seth_test_schemes, only : SCHEMES
implicit none
integer(I4P), allocatable :: ids(:), nfb(:), nfbu(:), nlim(:)
real(R_P), allocatable    :: l1(:), l1u(:)
logical, allocatable      :: pos(:), posu(:)
integer(I4P)              :: path, k, us
real(R_P)                 :: bound, js5
logical                   :: passed
real(R_P), parameter      :: WL(3) = [1._R_P, -2._R_P, 0.4_R_P]
real(R_P), parameter      :: WR(3) = [1._R_P, 2._R_P, 0.4_R_P]
real(R_P), parameter      :: X0 = 0.5_R_P, TEND = 0.15_R_P
passed = .true.
if (vv_report()) call vv_open('123-summary.dat', us, 'scheme positive-without-limiter limited-faces L1(rho)')
do path=PATH_SPLIT, PATH_INTERP
   call vv_schemes(path, .not.vv_report(), ids)
   allocate(l1(size(ids)), pos(size(ids)), l1u(size(ids)), posu(size(ids)), nfb(size(ids)), nfbu(size(ids)), &
            nlim(size(ids)))
   limiter = .false.
   call riemann_problem('123', WL, WR, X0, TEND, 0.5_R_P, 200, ids, path, l1u, posu, nfbu, files=.false.)
   limiter = .true.
   call riemann_problem('123', WL, WR, X0, TEND, 0.5_R_P, 200, ids, path, l1, pos, nfb, nlim)
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
      if (((path == PATH_SPLIT .and. SCHEMES(ids(k))%order <= 7._R_P) .or. &
           any(vv_label(ids(k)) == [character(16) :: 'weno3_js_interp', 'weno5_js_interp', 'wenoz3_interp'])) &
          .and. .not.posu(k)) call fail(vv_label(ids(k))//': not positive without the limiter')
      bound = recorded(vv_label(ids(k)))
      if (.not.vv_report() .and. bound > 0._R_P .and. l1(k) > bound) call fail(vv_label(ids(k))//': L1 above bound')
   enddo
   deallocate(l1, pos, l1u, posu, nfb, nfbu, nlim)
enddo
if (vv_report()) close(us)
if (.not.passed) error stop '123 verification failed'
print '(A)', 'seth_vv_test123: ok'

contains
   pure function recorded(name) result(b)
   !< L1 bounds of the test run (N = 200, with the limiter): measured (gfortran -O2) + 5%; -1: none.
   character(*), intent(in) :: name
   real(R_P)                :: b
   select case(name)
   case('weno5_js') ; b = 5.6569e-03_R_P ! 5.3875e-03
   case('wenoz5') ; b = 5.2971e-03_R_P ! 5.0449e-03
   case('teno5') ; b = 5.2417e-03_R_P ! 4.9921e-03
   case('teno6') ; b = 5.4952e-03_R_P ! 5.2335e-03
   case('weno5_js_interp') ; b = 5.8725e-03_R_P ! 5.5929e-03
   case('wenoz5_interp') ; b = 6.0822e-03_R_P ! 5.7926e-03
   case('teno5_interp') ; b = 6.0884e-03_R_P ! 5.7985e-03
   case('teno6_interp') ; b = 5.9806e-03_R_P ! 5.6958e-03
   case default ; b = -1._R_P
   end select
   endfunction recorded

   subroutine fail(what)
   character(*), intent(in) :: what
   print '(A)', '   FAILED: '//what
   passed = .false.
   endsubroutine fail
endprogram seth_vv_test123
