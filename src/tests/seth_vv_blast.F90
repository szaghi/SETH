!< SETH verification: the Woodward-Colella interacting blast waves against a 16 times finer reference.
program seth_vv_blast
!< SETH verification: the Woodward-Colella interacting blast waves against a 16 times finer reference.
!<
!< Woodward & Colella (1984): rho = 1, u = 0 on [0, 1] between reflecting walls, p = 1000 for x < 0.1, 0.01 in the
!< middle, 100 for x > 0.9, gamma = 1.4, t = 0.038: two strong blast waves, their reflections and their collision
!< leave narrow density peaks and contacts. No exact solution: the reference is WENO5-JS (splitting) on 16 times the
!< cells (the test run: 4 times), and L1 of the density is measured against it. Both paths.
!<
!< Every scheme runs twice: without any positivity safeguard, recording whether rho and p stay positive (just before
!< the collision the low-pressure gap shrinks to one or two cells, where any polynomial reconstruction can go negative:
!< Harten et al. 1987, section 7), then with the positivity-preserving flux limiter of Hu, Adams & Shu (2013), counting
!< the limited faces. Asserted: with the limiter, positivity throughout and, in the test run (N = 200, reference 800),
!< L1 below the recorded bounds of the quick schemes (measured + 5%); without it, WENO5-JS (splitting) stays positive.
!< Report mode: N = 400 against N = 6400, every scheme, profiles, 60 frames, and blast-summary.dat (per scheme:
!< positive without the limiter, limited faces, L1).
use seth
use seth_vv_euler
use seth_vv_cases
implicit none
integer(I4P), allocatable :: ids(:)
real(R8P), allocatable    :: l1(:), l1u(:), tf(:), qref(:,:)
logical, allocatable      :: pos(:), posu(:)
integer(I4P), allocatable :: nlim(:)
integer(I4P)              :: path, k, n, nref, us
real(R8P)                 :: bound
logical                   :: passed
passed = .true.
if (vv_report()) call vv_open('blast-summary.dat', us, 'scheme positive-without-limiter limited-faces L1(rho)')
n = merge(400, 200, vv_report()) ; nref = merge(16, 4, vv_report())*n
tf = [(0.038_R8P*k/60._R8P, k=1, 60)]
do path=PATH_SPLIT, PATH_INTERP
   call vv_schemes(path, .not.vv_report(), ids)
   allocate(l1(size(ids)), pos(size(ids)), l1u(size(ids)), posu(size(ids)), nlim(size(ids)))
   limiter = .false.
   call reference_problem('blast', init_blast, 0._R8P, 1._R8P, BC_REFLECTIVE, 0.038_R8P, 0.5_R8P, n, nref, 3_I4P, &
                          ids, path, l1u, posu, files=.false., qref_in=qref)
   limiter = .true.
   call reference_problem('blast', init_blast, 0._R8P, 1._R8P, BC_REFLECTIVE, 0.038_R8P, 0.5_R8P, n, nref, 3_I4P, &
                          ids, path, l1, pos, tf, nlim, qref_in=qref)
   do k=1, size(ids)
      print '(A18,ES12.4,L3,A,L2,A,I7)', vv_label(ids(k)), l1(k), pos(k), '   positive without limiter:', posu(k), &
            '   limited faces:', nlim(k)
      if (vv_report()) write(us, '(A,1X,L1,1X,I8,ES24.15)') vv_label(ids(k)), posu(k), nlim(k), l1(k)
      if (.not.pos(k)) call fail(vv_label(ids(k))//': not positive with the limiter')
      if (vv_label(ids(k)) == 'weno5_js' .and. .not.posu(k)) call fail('weno5_js not positive without the limiter')
      bound = recorded(vv_label(ids(k)))
      if (.not.vv_report() .and. bound > 0._R8P .and. l1(k) > bound) call fail(vv_label(ids(k))//': L1 above bound')
   enddo
   deallocate(l1, pos, l1u, posu, nlim)
enddo
if (vv_report()) close(us)
if (.not.passed) error stop 'blast-wave verification failed'
print '(A)', 'seth_vv_blast: ok'

contains
   pure function recorded(name) result(b)
   !< L1 bounds of the test run (N = 200 against N = 800, with the limiter): measured (gfortran -O2) + 5%; -1: none.
   character(*), intent(in) :: name
   real(R8P)                :: b
   select case(name)
   case('weno5_js') ; b = 1.4765e-01_R8P ! 1.4062e-01
   case('wenoz5') ; b = 1.0963e-01_R8P ! 1.0441e-01
   case('teno5') ; b = 9.8958e-02_R8P ! 9.4246e-02
   case('teno6') ; b = 6.6681e-02_R8P ! 6.3506e-02
   case('weno5_js_interp') ; b = 7.8476e-02_R8P ! 7.4739e-02
   case('wenoz5_interp') ; b = 4.8645e-02_R8P ! 4.6329e-02
   case('teno5_interp') ; b = 4.9426e-02_R8P ! 4.7072e-02
   case('teno6_interp') ; b = 4.5190e-02_R8P ! 4.3038e-02
   case default ; b = -1._R8P
   end select
   endfunction recorded

   subroutine fail(what)
   character(*), intent(in) :: what
   print '(A)', '   FAILED: '//what
   passed = .false.
   endsubroutine fail
endprogram seth_vv_blast
