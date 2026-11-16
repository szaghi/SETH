!< SETH test: the test registry is complete and consistent (every scheme dispatched, names unique, _wr kernels give
!< the states of the plain ones bitwise).
program seth_test_registry
!< SETH test: the test registry is complete and consistent (every scheme dispatched, names unique, _wr kernels give
!< the states of the plain ones bitwise).
use seth
use seth_test_rng, only : seed_rng
use seth_test_schemes
implicit none
real(R_P)    :: vp(VLO:VHI), vm(VLO:VHI), rp, rm, rp2, rm2, wrp, wrm
integer(I4P) :: id, j, n_
logical      :: passed

passed = .true.
call seed_rng(1002_I4P)
do id=1, NSCHEMES
   do j=id+1, NSCHEMES
      if (SCHEMES(id)%name == SCHEMES(j)%name) then ; print '(A)', 'duplicate '//SCHEMES(id)%name ; passed = .false. ; endif
   enddo
   do n_=1, 200
      call random_number(vp) ; call random_number(vm)
      call reconstruct(id, vp, vm, rp, rm)
      call reconstruct(id, vp, vm, rp2, rm2, wrp, wrm)
      if (rp /= rp2 .or. rm /= rm2) then
         print '(A)', trim(SCHEMES(id)%name)//': _wr kernel states differ from the plain kernel' ; passed = .false. ; exit
      endif
      if (SCHEMES(id)%family <= FAM_JS_INTERP .and. .not.(wrp >= 0._R_P .and. wrp <= 1._R_P + 4._R_P*epsilon(1._R_P) &
                                                      .and. wrm >= 0._R_P .and. wrm <= 1._R_P + 4._R_P*epsilon(1._R_P))) then
         print '(A,2ES12.4)', trim(SCHEMES(id)%name)//': weights ratio outside [0,1]: ', wrp, wrm ; passed = .false. ; exit
      endif
   enddo
enddo
print '(A,I0,A)', 'registry: ', NSCHEMES, ' schemes dispatched'
if (.not.passed) error stop 'registry test failed'
print '(A)', 'seth_test_registry: ok'
endprogram seth_test_registry
