!< SETH verification: the Shu-Osher shock / entropy-wave interaction against a 16 times finer reference.
program seth_vv_shu_osher
!< SETH verification: the Shu-Osher shock / entropy-wave interaction against a 16 times finer reference.
!<
!< Shu & Osher (1989), example 8: a Mach 3 shock at x = -4 runs into the density waves rho = 1 + 0.2 sin(5 x) on
!< [-5, 5], transmissive ends, gamma = 1.4, t = 1.8; behind the shock the waves are compressed into a band the scheme
!< must resolve. No exact solution: the reference is WENO5-JS (splitting) on 16 times the cells (the test run: 4 times),
!< and L1 of the density is measured against it. Both paths: the splitting path with every reconstruction, the
!< interpolation path with every interpolation (Riemann solver HLLC with the high-order flux corrections).
!<
!< Asserted: positivity throughout; in the test run (N = 200, reference 800) L1 below the recorded bounds of the quick
!< schemes (measured + 5%), and the low-dissipation claims of the literature: WENO-Z5 (Borges et al. 2008) and TENO5
!< (Fu, Hu & Adams 2016) closer to the reference than WENO5-JS, on each path.
!< Report mode: N = 400 against N = 6400, every scheme, the profiles, 60 frames of every scheme and
!< shu-osher-summary.dat (per scheme L1).
use seth
use seth_vv_euler
use seth_vv_cases
implicit none
integer(I4P), allocatable :: ids(:)
real(R_P), allocatable    :: l1(:), tf(:)
logical, allocatable      :: pos(:)
integer(I4P)              :: path, k, n, nref, us
real(R_P)                 :: js5, bound
real(R_P), allocatable    :: qref(:,:)
logical                   :: passed
passed = .true.
if (vv_report()) call vv_open('shu-osher-summary.dat', us, 'scheme L1(rho) (against the reference)')
n = merge(400, 200, vv_report()) ; nref = merge(16, 4, vv_report())*n
tf = [(1.8_R_P*k/60._R_P, k=1, 60)]
do path=PATH_SPLIT, PATH_INTERP
   call vv_schemes(path, .not.vv_report(), ids)
   allocate(l1(size(ids)), pos(size(ids)))
   call reference_problem('shu-osher', init_shu_osher, -5._R_P, 5._R_P, BC_TRANSMISSIVE, 1.8_R_P, 0.5_R_P, n, nref, &
                          3_I4P, ids, path, l1, pos, tf, qref_in=qref)
   js5 = -1._R_P
   do k=1, size(ids)
      if (vv_label(ids(k)) == 'weno5_js' .or. vv_label(ids(k)) == 'weno5_js_interp') js5 = l1(k)
   enddo
   do k=1, size(ids)
      print '(A18,ES12.4,L3)', vv_label(ids(k)), l1(k), pos(k)
      if (vv_report()) write(us, '(A,1X,ES24.15)') vv_label(ids(k)), l1(k)
      if (.not.pos(k)) call fail(vv_label(ids(k))//': not positive')
      bound = recorded(vv_label(ids(k)))
      if (.not.vv_report() .and. bound > 0._R_P .and. l1(k) > bound) call fail(vv_label(ids(k))//': L1 above bound')
      select case(vv_label(ids(k)))
      case('wenoz5', 'teno5', 'wenoz5_interp', 'teno5_interp')
         if (.not.(l1(k) < js5)) call fail(vv_label(ids(k))//': not closer to the reference than WENO5-JS')
      end select
   enddo
   deallocate(l1, pos)
enddo
if (vv_report()) close(us)
! G7 of issue #6: the non-polynomial schemes on the splitting path against TENO5's recorded bound (the reference reused)
call vv_nonpoly_ids(ids)
allocate(l1(size(ids)), pos(size(ids)))
call reference_problem('shu-osher', init_shu_osher, -5._R_P, 5._R_P, BC_TRANSMISSIVE, 1.8_R_P, 0.5_R_P, n, nref, &
                       3_I4P, ids, PATH_SPLIT, l1, pos, files=.false., qref_in=qref)
do k=1, size(ids)
   print '(A18,ES12.4,L3,A,ES12.4)', vv_label(ids(k)), l1(k), pos(k), '   TENO5''s bound', recorded('teno5')
   if (.not.pos(k)) call fail(vv_label(ids(k))//': not positive')
   if (.not.vv_report() .and. l1(k) > recorded('teno5')) call fail(vv_label(ids(k))//': L1 above TENO5''s bound')
enddo
deallocate(l1, pos)
if (.not.passed) error stop 'Shu-Osher verification failed'
print '(A)', 'seth_vv_shu_osher: ok'

contains
   pure function recorded(name) result(b)
   !< L1 bounds of the test run (N = 200 against N = 800): measured (gfortran -O2) + 5%; -1: none.
   character(*), intent(in) :: name
   real(R_P)                :: b
   select case(name)
   case('weno5_js') ; b = 8.0728e-1_R_P ! 7.6884e-1
   case('wenoz5') ; b = 5.8472e-1_R_P ! 5.5688e-1
   case('teno5') ; b = 5.2884e-1_R_P ! 5.0366e-1
   case('teno6') ; b = 2.7250e-1_R_P ! 2.5952e-1
   case('weno5_js_interp') ; b = 6.8476e-1_R_P ! 6.5215e-1
   case('wenoz5_interp') ; b = 3.6405e-1_R_P ! 3.4671e-1
   case('teno5_interp') ; b = 3.1265e-1_R_P ! 2.9776e-1
   case('teno6_interp') ; b = 2.2110e-1_R_P ! 2.1057e-1
   case default ; b = -1._R_P
   end select
   endfunction recorded

   subroutine fail(what)
   character(*), intent(in) :: what
   print '(A)', '   FAILED: '//what
   passed = .false.
   endsubroutine fail
endprogram seth_vv_shu_osher
