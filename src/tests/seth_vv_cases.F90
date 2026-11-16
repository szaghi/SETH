!< SETH verification: what the verification cases share (scheme lists, runs, norms, report files).
module seth_vv_cases
!< SETH verification: what the verification cases share (scheme lists, runs, norms, report files).
!<
!< Every verification case is a program that asserts its bounds and exits non-zero on a violation, so that it runs in
!< the test suite (scripts/run_tests.sh) on a small set of schemes and grids. With the environment variable
!< SETH_VV_DIR set to a directory, the same program runs every scheme on the report grids and writes its data there
!< (profiles, error tables, frames) for the documentation's tables, figures and animations (scripts/run_vv.sh,
!< docs/scripts/make_doc_figures.py): the report mode. The asserted bounds are the same in both modes.
use seth
use seth_test_parameters, only : T_SLACK
use seth_test_schemes
use seth_vv_euler
use seth_vv_riemann, only : riemann_sample
implicit none
private
public :: vv_dir, vv_report, vv_open, vv_schemes, vv_label
public :: euler_run, l1_norm, riemann_problem, reference_problem, write_frames
public :: init_shu_osher, init_blast

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

   function vv_label(id) result(s)
   character(:), allocatable :: s
   integer(I4P), intent(in)  :: id
   s = trim(SCHEMES(id)%name)
   endfunction vv_label

   subroutine euler_run(id, path, w0, n, h, bc, tend, cfl, q, positive, frames, tframes, nframes)
   !< Advance the Euler state from the primitive w0 to tend; `positive`: rho and p stayed > 0 at every step. With
   !< `frames`, the primitive state at the times tframes(1:nframes) (the step is cut to land on them).
   integer(I4P), intent(in)            :: id, path, n, bc
   real(R_P),    intent(in)            :: w0(3,n), h, tend, cfl
   real(R_P),    intent(out)           :: q(3,1-NG:n+NG)
   logical,      intent(out)           :: positive
   real(R_P),    intent(out), optional :: frames(:,:,:) !< (3, n, nframes)
   real(R_P),    intent(in),  optional :: tframes(:)
   integer(I4P), intent(in),  optional :: nframes
   real(R_P)                           :: t, dt, tnext
   integer(I4P)                        :: i, kf
   do i=1, n
      q(:,i) = prim2cons(w0(:,i))
   enddo
   if (bc == BC_INFLOW) then ! the inflow state, held by the left ghost cells for the whole run
      do i=1-NG, 0
         q(:,i) = q(:,1)
      enddo
   endif
   t = 0._R_P ; positive = .true. ; kf = 1
   do while (t < tend*(1._R_P - T_SLACK))
      dt = euler_dt(q, n, h, cfl)
      tnext = tend
      if (present(frames)) then
         if (kf <= nframes) tnext = min(tend, tframes(kf))
      endif
      dt = min(dt, tnext - t)
      call euler_step(id, path, q, n, h, dt, bc, positive)
      t = t + dt
      if (.not.positive) return
      if (present(frames)) then
         if (kf <= nframes) then
            if (t >= tframes(kf)*(1._R_P - T_SLACK)) then
               do i=1, n
                  frames(:,i,kf) = cons2prim(q(:,i))
               enddo
               kf = kf + 1
            endif
         endif
      endif
   enddo
   endsubroutine euler_run

   subroutine riemann_problem(name, wl, wr, x0, tend, cfl, nmain, ids, path, l1, positive, nfb, nlim, files)
   !< A shock tube on [0, 1] against its exact solution: every scheme of ids on nmain cells, transmissive ends; L1 of
   !< the density over the domain (point values at the cell centres against the exact ones); `positive`: rho, p > 0
   !< throughout; `nfb`: interpolated face states that fell back to first order; `nlim`: faces limited by the
   !< positivity limiter (when the caller switched it on). Report mode (unless files = .false.): the profiles
   !< (<name>-<scheme>.dat), the exact solution (<name>-exact.dat), and L1 on 100 .. 1600 cells (<name>-l1.dat).
   character(*), intent(in)  :: name
   real(R_P),    intent(in)  :: wl(3), wr(3), x0, tend, cfl
   integer(I4P), intent(in)  :: nmain, ids(:), path
   real(R_P),    intent(out) :: l1(size(ids))
   logical,      intent(out) :: positive(size(ids))
   integer(I4P), intent(out) :: nfb(size(ids))
   integer(I4P), intent(out), optional :: nlim(size(ids))
   logical,      intent(in),  optional :: files
   real(R_P), allocatable    :: w0(:,:), q(:,:), ex(:,:), e(:)
   real(R_P)                 :: h, x, w(3), l1n
   integer(I4P)              :: k, i, n, u, uc, j
   logical                   :: pos, write_files
   write_files = vv_report() ; if (present(files)) write_files = write_files .and. files
   do k=1, size(ids)
      n = nmain ; h = 1._R_P/n
      call setup(n)
      nfallback = 0 ; nlimited = 0
      call euler_run(ids(k), path, w0, n, h, BC_TRANSMISSIVE, tend, cfl, q, positive(k))
      nfb(k) = nfallback
      if (present(nlim)) nlim(k) = nlimited
      if (.not.positive(k)) then
         l1(k) = huge(1._R_P) ; cycle
      endif
      do i=1, n
         e(i) = q(1,i) - ex(1,i)
      enddo
      l1(k) = l1_norm(e, h)
      if (write_files) then
         call vv_open(name//'-'//vv_label(ids(k))//'.dat', u, 'x rho u p (scheme '//vv_label(ids(k))//', '// &
                      'N = '//trim(str_i(n))//', t = '//trim(str_r(tend))//')')
         do i=1, n
            w = cons2prim(q(:,i))
            write(u, '(4ES24.15)') (i - 0.5_R_P)*h, w
         enddo
         close(u)
      endif
   enddo
   if (write_files) then
      call vv_open(name//'-exact.dat', u, 'x rho u p (exact, t = '//trim(str_r(tend))//')')
      do i=1, 4000
         x = (i - 0.5_R_P)/4000._R_P
         write(u, '(4ES24.15)') x, riemann_sample(wl, wr, GAMMA, (x - x0)/tend)
      enddo
      close(u)
      inquire(file=vv_dir()//'/'//name//'-l1.dat', exist=pos)
      if (pos) then
         open(newunit=uc, file=vv_dir()//'/'//name//'-l1.dat', position='append', action='write')
      else
         call vv_open(name//'-l1.dat', uc, 'scheme N L1(rho) (point values against the exact solution)')
      endif
      do k=1, size(ids)
         do j=0, 4
            n = 100*2**j ; h = 1._R_P/n
            call setup(n)
            call euler_run(ids(k), path, w0, n, h, BC_TRANSMISSIVE, tend, cfl, q, pos)
            l1n = -1._R_P
            if (pos) then
               do i=1, n
                  e(i) = q(1,i) - ex(1,i)
               enddo
               l1n = l1_norm(e, h)
            endif
            write(uc, '(A,1X,I6,ES24.15)') vv_label(ids(k)), n, l1n
         enddo
      enddo
      close(uc)
   endif
   contains
      subroutine setup(nn)
      integer(I4P), intent(in) :: nn
      integer(I4P)             :: ii
      real(R_P)                :: xx
      if (allocated(w0)) deallocate(w0, q, ex, e)
      allocate(w0(3,nn), q(3,1-NG:nn+NG), ex(3,nn), e(nn))
      do ii=1, nn
         xx = (ii - 0.5_R_P)/nn
         w0(:,ii) = merge(wl, wr, xx < x0)
         ex(:,ii) = riemann_sample(wl, wr, GAMMA, (xx - x0)/tend)
      enddo
      endsubroutine setup
   endsubroutine riemann_problem

   subroutine reference_problem(name, init, xl, xr, bc, tend, cfl, nmain, nref, refid, ids, path, l1, positive, &
                                tframes, nlim, files, qref_in)
   !< A problem without exact solution against a reference computed on nref = r nmain cells (r even) by scheme refid
   !< (splitting path; with the positivity limiter if the caller switched it on): L1 of the density at the nmain cell
   !< centres against the mean of the two reference cells that share that centre; `nlim`: faces limited per scheme.
   !< Report mode (unless files = .false.): the profiles, the reference, with `tframes` the frames of every scheme.
   character(*), intent(in)           :: name
   interface
      pure function init(x) result(w)
      import :: R_P
      real(R_P), intent(in) :: x
      real(R_P)             :: w(3)
      endfunction init
   endinterface
   real(R_P),    intent(in)           :: xl, xr, tend, cfl
   integer(I4P), intent(in)           :: bc, nmain, nref, refid, ids(:), path
   real(R_P),    intent(out)          :: l1(size(ids))
   logical,      intent(out)          :: positive(size(ids))
   real(R_P),    intent(in), optional :: tframes(:)
   integer(I4P), intent(out), optional :: nlim(size(ids))
   logical,      intent(in),  optional :: files
   real(R_P),    intent(inout), optional, allocatable :: qref_in(:,:) !< The reference: computed once, then reused.
   real(R_P), allocatable             :: w0(:,:), q(:,:), qref(:,:), e(:), fr(:,:,:)
   real(R_P)                          :: h, href, w(3)
   integer(I4P)                       :: k, i, r, u
   logical                            :: pos, write_files
   write_files = vv_report() ; if (present(files)) write_files = write_files .and. files
   r = nref/nmain
   href = (xr - xl)/nref
   allocate(w0(3,nref), qref(3,1-NG:nref+NG))
   reuse: block
   if (present(qref_in)) then
      if (allocated(qref_in)) then
         qref = qref_in ; exit reuse
      endif
   endif
   do i=1, nref
      w0(:,i) = init(xl + (i - 0.5_R_P)*href)
   enddo
   call euler_run(refid, PATH_SPLIT, w0, nref, href, bc, tend, cfl, qref, pos)
   if (.not.pos) error stop 'reference_problem: the reference run lost positivity'
   if (present(qref_in)) qref_in = qref
   endblock reuse
   if (write_files) then
      call vv_open(name//'-reference.dat', u, 'x rho u p (reference: '//vv_label(refid)//', N = '// &
                   trim(str_i(nref))//')')
      do i=1, nref
         write(u, '(4ES24.15)') xl + (i - 0.5_R_P)*href, cons2prim(qref(:,i))
      enddo
      close(u)
   endif
   deallocate(w0)
   h = (xr - xl)/nmain
   allocate(w0(3,nmain), q(3,1-NG:nmain+NG), e(nmain))
   do i=1, nmain
      w0(:,i) = init(xl + (i - 0.5_R_P)*h)
   enddo
   do k=1, size(ids)
      nlimited = 0
      if (present(tframes) .and. write_files) then
         if (allocated(fr)) deallocate(fr)
         allocate(fr(3,nmain,size(tframes)))
         call euler_run(ids(k), path, w0, nmain, h, bc, tend, cfl, q, positive(k), fr, tframes, size(tframes))
         if (positive(k)) call write_frames(name//'-frames-'//vv_label(ids(k))//'.dat', xl, h, nmain, tframes, fr)
      else
         call euler_run(ids(k), path, w0, nmain, h, bc, tend, cfl, q, positive(k))
      endif
      if (present(nlim)) nlim(k) = nlimited
      l1(k) = huge(1._R_P)
      if (.not.positive(k)) cycle
      do i=1, nmain
         e(i) = q(1,i) - 0.5_R_P*(qref(1,(i-1)*r + r/2) + qref(1,(i-1)*r + r/2 + 1))
      enddo
      l1(k) = l1_norm(e, h)
      if (write_files) then
         call vv_open(name//'-'//vv_label(ids(k))//'.dat', u, 'x rho u p (scheme '//vv_label(ids(k))//', N = '// &
                      trim(str_i(nmain))//', t = '//trim(str_r(tend))//')')
         do i=1, nmain
            w = cons2prim(q(:,i))
            write(u, '(4ES24.15)') xl + (i - 0.5_R_P)*h, w
         enddo
         close(u)
      endif
   enddo
   endsubroutine reference_problem

   subroutine write_frames(name, xl, h, n, tframes, fr)
   !< Frames file: line 1 the cell centres, then one line per frame: t rho(1..n) u(1..n) p(1..n).
   character(*), intent(in) :: name
   real(R_P),    intent(in) :: xl, h, tframes(:), fr(:,:,:)
   integer(I4P), intent(in) :: n
   integer(I4P)             :: u, i, f
   call vv_open(name, u, 'frames: line 1 x(1..N); then per frame: t rho(1..N) u(1..N) p(1..N)')
   write(u, '(*(ES16.8))') (xl + (i - 0.5_R_P)*h, i=1, n)
   do f=1, size(tframes)
      write(u, '(*(ES16.8))') tframes(f), fr(1,:,f), fr(2,:,f), fr(3,:,f)
   enddo
   close(u)
   endsubroutine write_frames

   pure function init_shu_osher(x) result(w)
   !< Shu & Osher (1989), example 8: a Mach 3 shock at x = -4 running into rho = 1 + 0.2 sin(5 x), on [-5, 5].
   real(R_P), intent(in) :: x
   real(R_P)             :: w(3)
   if (x < -4._R_P) then
      w = [3.857143_R_P, 2.629369_R_P, 10.33333_R_P]
   else
      w = [1._R_P + 0.2_R_P*sin(5._R_P*x), 0._R_P, 1._R_P]
   endif
   endfunction init_shu_osher

   pure function init_blast(x) result(w)
   !< Woodward & Colella (1984): p = 1000, 0.01, 100 on [0, 0.1), [0.1, 0.9], (0.9, 1]; rho = 1, u = 0.
   real(R_P), intent(in) :: x
   real(R_P)             :: w(3)
   w = [1._R_P, 0._R_P, merge(1000._R_P, merge(100._R_P, 0.01_R_P, x > 0.9_R_P), x < 0.1_R_P)]
   endfunction init_blast

   pure function str_i(i) result(s)
   integer(I4P), intent(in) :: i
   character(16)            :: s
   write(s, '(I0)') i
   endfunction str_i

   pure function str_r(r) result(s)
   real(R_P), intent(in) :: r
   character(24)         :: s
   write(s, '(G0.6)') r
   endfunction str_r

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
endmodule seth_vv_cases
