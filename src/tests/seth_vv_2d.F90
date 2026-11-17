!< SETH verification: the non-polynomial schemes on three 2-D cases of the minimal 2-D harness (issue #6, step 4.3).
program seth_vv_2d
!< SETH verification: the non-polynomial schemes on three 2-D cases of the minimal 2-D harness (issue #6, step 4.3).
!<
!< The harness is seth_vv_euler2d (fenced by owner decision D8). Three cases on [0, 1]^2, CFL 0.4:
!< - an oblique contact: rho = 1.5 on the diagonal band 0.25 <= frac(x + y) < 0.75, 1 elsewhere, u = v = 1, p = 1,
!<   periodic, t = 1 (one period: the exact solution is the initial one), symmetric about the diagonal;
!< - Kelvin-Helmholtz (Liang & Fu 2024, Eq. 46): (2, -0.5, 0.01 sin 2 pi x, 2.5) for 0.25 <= y <= 0.75, (1, 0.5,
!<   0.01 sin 2 pi x, 2.5) elsewhere, periodic, t = 0.5;
!< - the 2-D Riemann problem of configuration 3 (Liang & Fu 2024, Eq. 47), transmissive, t = 0.3, symmetric about the
!<   diagonal.
!< Asserted: the harness's eigenvectors (L R = I to 64 ulps at a random state); positivity; on the periodic cases the
!< sums of the conservative variables conserved to NP_TOL of their magnitude; on the contact the solution symmetric
!< about the diagonal (q(i,j) and q(j,i) with the momenta swapped) to 1e5 ulps of the data. Configuration 3's symmetry
!< is reported: bitwise with gfortran and ifx, but where the compiler breaks the identity of the x and y passes
!< (nvfortran -fast, -acc) its slip lines amplify the round-off asymmetry, to 2.6e-8 by t = 0.3 with RBF+DR at N = 48
!< (measured 2026-11-17): the reason the papers add a symmetry-preserving technique.
!< The suite runs N = 48 (the contact and configuration 3 to their final time, KH to t = 0.1) with TENO5, TENO5-DR and
!< the complete RBF+DR scheme; report mode (SETH_VV_DIR) N = 400 and the seven schemes of the ablation to the final
!< times, writing the density fields (nonpoly2d-<case>-<scheme>.dat) and a summary (nonpoly2d-summary.dat). The
!< papers' 2-D runs also use a symmetry-preserving technique (Fleischmann, Adami & Adams 2019) that this harness does
!< not: its symmetry rests on the identical x and y passes.
use seth
use seth_vv_euler, only : NG
use seth_vv_euler2d
use seth_vv_cases, only : vv_report, vv_open
use seth_test_schemes, only : SCHEMES, NALL
use seth_test_parameters, only : NP_TOL
implicit none
integer(I4P), parameter  :: NA = 7_I4P
character(16), parameter :: AN(NA) = [character(16) :: 'teno5', 'teno5_st', 'teno5_rbf', 'teno5_thinc', 'teno5_dr', &
                                      'teno5_st_dr', 'teno5_rbf_dr']
logical                  :: passed, rep
integer(I4P)             :: k, ns, n, us
real(R_P)                :: l(4,4), r(4,4), lam(4), e
passed = .true.
rep = vv_report()
! the eigenvectors
call eigen2(prim2cons2([1.3_R_P, 0.4_R_P, -0.7_R_P, 0.9_R_P]), prim2cons2([0.8_R_P, -0.2_R_P, 0.3_R_P, 1.4_R_P]), l, r, lam)
e = maxval(abs(matmul(l, r) - reshape([1._R_P, 0._R_P, 0._R_P, 0._R_P, 0._R_P, 1._R_P, 0._R_P, 0._R_P, &
                                       0._R_P, 0._R_P, 1._R_P, 0._R_P, 0._R_P, 0._R_P, 0._R_P, 1._R_P], [4, 4])))
print '(A,ES10.3)', 'eigenvectors: max |L R - I| = ', e
if (e > 64._R_P*epsilon(1._R_P)) then ; print '(A)', 'FAIL: eigenvectors' ; passed = .false. ; endif
n = merge(400_I4P, 48_I4P, rep)
if (rep) call vv_open('nonpoly2d-summary.dat', us, 'case scheme N positive conservation-error symmetry-error')
print '(A,I0)', '2-D cases, N = ', n
do k=1, NA
   if (.not.rep .and. .not.any(trim(AN(k)) == [character(16) :: 'teno5', 'teno5_dr', 'teno5_rbf_dr'])) cycle
   call run_case('contact', AN(k))
   call run_case('kh', AN(k))
   call run_case('config3', AN(k))
enddo
if (rep) close(us)
if (.not.passed) error stop 'seth_vv_2d failed'
print '(A)', 'seth_vv_2d: ok'

contains
   function scheme_id(name) result(id)
   character(*), intent(in) :: name
   integer(I4P)             :: id, j
   id = 0
   do j=1, NALL
      if (trim(SCHEMES(j)%name) == trim(name)) id = j
   enddo
   if (id == 0) error stop 'seth_vv_2d: scheme not in the registry'
   endfunction scheme_id

   subroutine run_case(cname, sname)
   character(*), intent(in) :: cname, sname
   real(R_P), allocatable   :: q(:,:,:)
   real(R_P)                :: h, x, y, t, tend, dt, tot0(4), tot(4), cons, sym
   integer(I4P)             :: i, j, bc, id, uf
   logical                  :: ok, periodic, symmetric
   id = scheme_id(sname)
   allocate(q(4,1-NG:n+NG,1-NG:n+NG))
   h = 1._R_P/n
   periodic = cname /= 'config3' ; symmetric = cname /= 'kh'
   bc = merge(BC2_PERIODIC, BC2_TRANSMISSIVE, periodic)
   do j=1, n
      do i=1, n
         x = (i - 0.5_R_P)*h ; y = (j - 0.5_R_P)*h
         select case(cname)
         case('contact')
            q(:,i,j) = prim2cons2([merge(1.5_R_P, 1._R_P, modulo(x + y, 1._R_P) >= 0.25_R_P .and. &
                                         modulo(x + y, 1._R_P) < 0.75_R_P), 1._R_P, 1._R_P, 1._R_P])
         case('kh')
            if (y >= 0.25_R_P .and. y <= 0.75_R_P) then
               q(:,i,j) = prim2cons2([2._R_P, -0.5_R_P, 0.01_R_P*sin(2._R_P*acos(-1._R_P)*x), 2.5_R_P])
            else
               q(:,i,j) = prim2cons2([1._R_P, 0.5_R_P, 0.01_R_P*sin(2._R_P*acos(-1._R_P)*x), 2.5_R_P])
            endif
         case default
            if (x > 0.5_R_P .and. y > 0.5_R_P) then
               q(:,i,j) = prim2cons2([1.5_R_P, 0._R_P, 0._R_P, 1.5_R_P])
            elseif (x < 0.5_R_P .and. y > 0.5_R_P) then
               q(:,i,j) = prim2cons2([0.5323_R_P, 1.206_R_P, 0._R_P, 0.3_R_P])
            elseif (x < 0.5_R_P .and. y < 0.5_R_P) then
               q(:,i,j) = prim2cons2([0.138_R_P, 1.206_R_P, 1.206_R_P, 0.029_R_P])
            else
               q(:,i,j) = prim2cons2([0.5323_R_P, 0._R_P, 1.206_R_P, 0.3_R_P])
            endif
         endselect
      enddo
   enddo
   select case(cname)
   case('contact') ; tend = 1._R_P
   case('kh') ; tend = merge(0.5_R_P, 0.1_R_P, rep)
   case default ; tend = 0.3_R_P
   endselect
   tot0 = sum(sum(q(:,1:n,1:n), dim=3), dim=2)
   t = 0._R_P ; ok = .true.
   do while (t < tend .and. ok)
      dt = euler2_dt(q, n, h, 0.4_R_P)
      if (t + dt > tend) dt = tend - t
      call euler2_step(id, q, n, h, dt, bc, ok)
      t = t + dt
   enddo
   tot = sum(sum(q(:,1:n,1:n), dim=3), dim=2)
   cons = maxval(abs(tot - tot0)/max(abs(tot0), 1._R_P))
   sym = 0._R_P
   if (symmetric) then
      do j=1, n
         do i=1, n
            sym = max(sym, abs(q(1,i,j) - q(1,j,i)), abs(q(2,i,j) - q(3,j,i)), abs(q(4,i,j) - q(4,j,i)))
         enddo
      enddo
   endif
   print '(A10,A14,A,L2,A,ES10.3,A,ES10.3)', cname, trim(sname), ': positive', ok, ', conservation ', cons, &
         ', symmetry ', sym
   if (.not.ok) then ; print '(A)', 'FAIL: not positive' ; passed = .false. ; endif
   ! the sums gather n^2 cells over ~n steps: NP_TOL at the suite's N = 48, scaled by (n/48)^2 beyond (N = 400: 7e-12;
   ! measured 3e-13 on Kelvin-Helmholtz)
   if (periodic .and. cons > NP_TOL*max(1._R_P, (n/48._R_P)**2)) then
      print '(A)', 'FAIL: not conservative' ; passed = .false.
   endif
   if (cname == 'contact' .and. sym > 1.e5_R_P*epsilon(1._R_P)*maxval(abs(q(:,1:n,1:n)))) then
      print '(A)', 'FAIL: not symmetric' ; passed = .false.
   endif
   if (rep) then
      write(us, '(A,1X,A,1X,I0,1X,L1,2ES14.5)') cname, trim(sname), n, ok, cons, sym
      call vv_open('nonpoly2d-'//cname//'-'//trim(sname)//'.dat', uf, 'density on the N x N cells, one row per y')
      do j=1, n
         write(uf, '(*(ES12.5,1X))') q(1,1:n,j)
      enddo
      close(uf)
   endif
   endsubroutine run_case
endprogram seth_vv_2d
