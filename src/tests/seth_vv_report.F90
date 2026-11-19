!< SETH verification: the report utilities and the scheme lists of the verification cases.
module seth_vv_report
!< SETH verification: the report utilities and the scheme lists of the verification cases.
!<
!< The report mode (SETH_VV_DIR, vv_report, vv_open), the scheme lists of a path (vv_schemes, vv_nonpoly_ids,
!< vv_label), the L1 norm, and the path constants. They depend on the scheme registry only, not on the Euler solver
!< and its time integrator (SELAH's, seth_vv_euler), so that a unit test writing report data needs neither:
!< seth_vv_cases re-exports them for the verification cases.
use seth
use seth_test_schemes, only : SCHEMES, NSCHEMES, NALL
implicit none
private
public :: PATH_SPLIT, PATH_INTERP, PATH_ROE
public :: vv_dir, vv_report, vv_open, vv_schemes, vv_label, vv_nonpoly_ids, l1_norm

integer(I4P), parameter :: PATH_SPLIT = 1, PATH_INTERP = 2, PATH_ROE = 3 !< Paths of the Euler solver.

contains
   function vv_dir() result(d)
   !< The report directory (SETH_VV_DIR), '' outside the report mode.
   character(:), allocatable :: d
   character(512)            :: buf
   integer                   :: l, st
   call get_environment_variable('SETH_VV_DIR', buf, l, st)
   d = ''
   if (st == 0 .and. l > 0) d = trim(buf)
   endfunction vv_dir

   function vv_report() result(r)
   logical :: r
   r = len(vv_dir()) > 0
   endfunction vv_report

   subroutine vv_open(name, u, header)
   !< Open <SETH_VV_DIR>/<name> for writing, with a '#' header line.
   character(*), intent(in)  :: name, header
   integer,      intent(out) :: u
   call execute_command_line('mkdir -p '//vv_dir())
   open(newunit=u, file=vv_dir()//'/'//name, status='replace', action='write')
   write(u, '(A)') '# '//header
   endsubroutine vv_open

   subroutine vv_schemes(path, quick, ids)
   !< The schemes of a path: the reconstructions (PATH_SPLIT) or the interpolations (PATH_INTERP), first order
   !< excluded; `quick`: the subset of the test-suite run (WENO5-JS, WENO-Z5, TENO5, TENO6 and their interpolations).
   integer(I4P),              intent(in)  :: path
   logical,                   intent(in)  :: quick
   integer(I4P), allocatable, intent(out) :: ids(:)
   integer(I4P)                           :: id
   character(16)                          :: n
   allocate(ids(0))
   do id=1, NSCHEMES
      if (SCHEMES(id)%s < 2) cycle
      if (SCHEMES(id)%interp .neqv. (path == PATH_INTERP)) cycle
      n = SCHEMES(id)%name
      if (quick .and. .not.any(n == [character(16) :: 'weno5_js', 'wenoz5', 'teno5', 'teno6', 'weno5_js_interp', &
                                     'wenoz5_interp', 'teno5_interp', 'teno6_interp'])) cycle
      ids = [ids, id]
   enddo
   endsubroutine vv_schemes

   subroutine vv_nonpoly_ids(ids)
   !< The non-polynomial schemes of issue #6 (splitting path): TENO5-DR, -THINC, -ST, -RBF, -ST-DR, -RBF-DR, the last
   !< four beyond the sweeps (NSCHEMES+1..NALL). G7 runs them against TENO5's recorded bounds.
   integer(I4P), allocatable, intent(out) :: ids(:)
   character(16), parameter :: NP(6) = [character(16) :: 'teno5_dr', 'teno5_thinc', 'teno5_st', 'teno5_rbf', &
                                        'teno5_st_dr', 'teno5_rbf_dr']
   integer(I4P)             :: k, j
   allocate(ids(size(NP)))
   do k=1, size(NP)
      ids(k) = 0
      do j=1, NALL
         if (trim(SCHEMES(j)%name) == trim(NP(k))) ids(k) = j
      enddo
      if (ids(k) == 0) error stop 'vv_nonpoly_ids: scheme not in the registry'
   enddo
   endsubroutine vv_nonpoly_ids

   function vv_label(id) result(s)
   character(:), allocatable :: s
   integer(I4P), intent(in)  :: id
   s = trim(SCHEMES(id)%name)
   endfunction vv_label

   pure function l1_norm(e, h) result(s)
   !< h sum |e| (pairwise, so that the sum does not lose the small terms).
   real(R_P), intent(in) :: e(:), h
   real(R_P)             :: s
   s = h*pairwise(abs(e))
   endfunction l1_norm

   pure recursive function pairwise(a) result(s)
   real(R_P), intent(in) :: a(:)
   real(R_P)             :: s
   integer               :: m
   if (size(a) <= 8) then
      s = sum(a)
   else
      m = size(a)/2
      s = pairwise(a(1:m)) + pairwise(a(m+1:))
   endif
   endfunction pairwise
endmodule seth_vv_report
