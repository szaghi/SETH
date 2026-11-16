!< SETH test: accuracy at critical points, against the orders measured in the literature.
program seth_test_critical_points
!< SETH test: accuracy at critical points, against the orders measured in the literature.
!<
!< The setting of Henrick et al. 2005 and Borges et al. 2008 (finite-difference WENO): point values f_i = f(x_i) on
!< nodes x_i = i h, the critical point at node 0, the fluxes at the faces 1/2 and -1/2 reconstructed from them (the
!< kernels' coefficients are those of the finite-difference flux reconstruction), and the derivative
!< D = (F_{1/2} - F_{-1/2})/h compared with f'(0) = 0. The order is that of D: where the nonlinear weights vary
!< non-smoothly from face to face, the face errors stop cancelling in D and the order drops (Henrick et al. 2005).
!< Both biases (the left state: a > 0; the right state: a < 0). eps = 1e-40 for every scheme, as Borges et al. use to
!< see the weights' own behaviour.
!<
!< - first-order critical point (f'(0) = 0, f''(0) /= 0): f = e^x - 1 - x;
!< - second-order critical point (f'(0) = f''(0) = 0, f'''(0) /= 0): f = e^x - 1 - x - x^2/2
!< (Tang et al. 2024's family e^x - sum_{k<=n} x^k/k!; the constant removed keeps the differences out of the
!< cancellation of data near 1). h is halved from h0 = 0.1, 0.2, 0.4 for stencils of half-width <= 3, 4, 5 (the
!< wide stencils reach round-off sooner), 8 levels; the order is taken on the finest pair of h whose errors stay above
!< the round-off of D, 10 eps max|f|/h (the data are small and computed without cancellation, so errors far below
!< 1e-13 are still exact: WENO9-JS converges to Henrick's 7 only below h = 0.025, after a pre-asymptotic dip).
!<
!< Asserted, at the first-order point (sources via llm-wiki, weno-teno `properties/critical-point-accuracy.md`):
!<
!< | scheme            | expected                  | source                                        |
!< |-------------------|---------------------------|-----------------------------------------------|
!< | WENO5-JS          | 3 (2.5 .. 3.6)            | Henrick 2005; Borges 2008: 2.95-2.99          |
!< | WENO5-Z, p = 1    | 4 (3.5 .. 4.6)            | Borges 2008: analysis 4, measured 3.93-3.99   |
!< | WENO5-Z, p = 2    | 5 (>= 4.5)                | Borges 2008: analysis 5, measured 4.80-4.98   |
!< | TENO5             | 5 (>= 4.5)                | Tang et al. 2024                              |
!< | TENO6 / TENO6-opt | 6 / 5 (>= 5.5 / >= 4.5)   | Zhang et al. 2024b                            |
!< | TENO5-opt         | 4 (>= 3.5)                | its linear scheme's order (Fu et al. 2016)    |
!< | TENO5-AZ          | 5 at both points (>= 4.5) | Tang 2024, Tables 1-3                         |
!< | TENO8             | 6 (>= 5.5)                | its sixth-order background (Liang & Fu 2024b) |
!<
!< WENO-JS of every order, at both points: Henrick et al.'s rate max(2r - 2 - n_cp, n_cp) (eq. 60: r <= 6, eps = 0;
!< llm-wiki high-order-capturing-schemes `properties/accuracy-at-critical-points.md`), within 0.5.
!<
!< At the second-order point the literature measures 2 for WENO5-Z and disagrees on TENO (Fu 2023: no degeneration;
!< Zhang et al. 2024b: 2; SETH measures 2): those schemes must reach at least 1.7, a floor that catches a collapse
!< without deciding the contradiction. The other schemes are printed, not asserted. The interpolation kernels are not
!< tested here: the literature has no critical-point test for them.
use seth
use seth_test_parameters, only : EPS_ZT
use seth_test_schemes
use seth_vv_cases, only : vv_report, vv_open
implicit none
integer(I4P), parameter :: NH = 8_I4P
real(R_P),    parameter :: FLOOR_ = 10._R_P*epsilon(1._R_P)
real(R_P)               :: o1(2), o2(2), lo, hi
integer(I4P)            :: id, p, r
real(R_P)               :: rate1, rate2
logical                 :: passed
integer                 :: ur

passed = .true.
if (vv_report()) call vv_open('critical-points.dat', ur, 'scheme p cp1(a>0) cp1(a<0) cp2(a>0) cp2(a<0)')
print '(A16,A4,4A12)', 'scheme', 'p', '1st-cp a>0', '1st-cp a<0', '2nd-cp a>0', '2nd-cp a<0'
do id=1, NSCHEMES
   if (SCHEMES(id)%interp) cycle ! the finite-difference flux setting reconstructs: interpolations elsewhere
   do p=1, merge(2, 1, SCHEMES(id)%family == FAM_Z)
      call orders(id, 1, p, o1)
      call orders(id, 2, p, o2)
      print '(A16,I4,4F12.2)', SCHEMES(id)%name, p, o1, o2
      if (vv_report()) write(ur, '(A,1X,I1,4F8.3)') trim(SCHEMES(id)%name), p, o1, o2
      lo = -1._R_P ; hi = 1.e9_R_P
      if (SCHEMES(id)%family == FAM_JS) then ! Henrick et al.'s rate
         r = SCHEMES(id)%s
         if (r >= 2) then
            rate1 = real(max(2*r - 3, 1), R_P) ; rate2 = real(max(2*r - 4, 2), R_P)
            if (any(abs(o1 - rate1) > 0.5_R_P) .or. any(abs(o2 - rate2) > 0.5_R_P)) then
               print '(A,2F6.1)', '   ^ not Henrick et al.''s rates ', rate1, rate2 ; passed = .false.
            endif
         endif
      endif
      select case(trim(SCHEMES(id)%name))
      case('weno5_js')  ; lo = 2.5_R_P ; hi = 3.6_R_P
      case('wenoz5')    ; if (p == 1) then ; lo = 3.5_R_P ; hi = 4.6_R_P ; else ; lo = 4.5_R_P ; endif
      case('teno5')     ; lo = 4.5_R_P
      case('teno6')     ; lo = 5.5_R_P
      case('teno6_opt') ; lo = 4.5_R_P
      case('teno5_opt') ; lo = 3.5_R_P
      case('teno5_az')  ; lo = 4.5_R_P
      case('teno8')     ; lo = 5.5_R_P
      case('teno5_a')   ; lo = 4.5_R_P
      case('teno8_a')   ; lo = 5.5_R_P
      endselect
      if (any(o1 < lo) .or. any(o1 > hi)) then
         print '(A,2ES10.2)', '   ^ first-order critical point outside the literature range: ', lo, hi
         passed = .false.
      endif
      if (SCHEMES(id)%family /= FAM_JS .and. SCHEMES(id)%order >= 3._R_P .and. any(o2 < 1.7_R_P)) then
         print '(A)', '   ^ second-order critical point below the floor 1.7' ; passed = .false.
      endif
      if (trim(SCHEMES(id)%name) == 'teno5_az' .and. any(o2 < 4.5_R_P)) then ! Tang 2024, Tables 2-3
         print '(A)', '   ^ TENO5-AZ below 4.5 at the second-order critical point (Tang 2024: 5)' ; passed = .false.
      endif
   enddo
enddo
if (.not.passed) error stop 'critical-point test failed'
if (vv_report()) close(ur)
print '(A)', 'seth_test_critical_points: ok'

contains
   subroutine orders(id, kind, p, o)
   !< Observed orders of D for both biases at the critical point `kind` (1 or 2).
   integer(I4P), intent(in)  :: id, kind, p
   real(R_P),    intent(out) :: o(2)
   real(R_P)                 :: up(VLO:VHI), um(VLO:VHI), h, err(2,NH), flo(NH), rp1, rm1, rp0, rm0, par, h0
   integer(I4P)              :: j, m, side
   par = SCHEMES(id)%par
   if (SCHEMES(id)%family == FAM_Z) par = real(p, R_P)
   h0 = 0.1_R_P*2._R_P**max(0, SCHEMES(id)%s - 3)
   do j=1, NH
      h = h0/2._R_P**(j-1)
      do m=VLO, VHI ! face 1/2: cell m is node m; face -1/2: cell m is node m-1
         up(m) = f(kind, m*h)
         um(m) = f(kind, (m-1)*h)
      enddo
      call reconstruct(id, up, up, rp1, rm1, eps=EPS_ZT, par=par)
      call reconstruct(id, um, um, rp0, rm0, eps=EPS_ZT, par=par)
      err(1,j) = abs((rp1 - rp0)/h) ! f'(0) = 0
      err(2,j) = abs((rm1 - rm0)/h)
      flo(j) = FLOOR_*max(maxval(abs(up)), maxval(abs(um)))/h
   enddo
   do side=1, 2
      o(side) = 99._R_P ! exact to round-off on every h
      do j=NH-1, 1, -1
         if (err(side,j+1) > flo(j+1)) then
            o(side) = log(err(side,j)/err(side,j+1))/log(2._R_P) ; exit
         endif
      enddo
   enddo
   endsubroutine orders

   pure function f(kind, x) result(y)
   !< e^x - sum_{k=0}^{kind} x^k/k!, by its Taylor series beyond the subtracted terms (no cancellation).
   integer(I4P), intent(in) :: kind
   real(R_P),    intent(in) :: x
   real(R_P)                :: y, t
   integer(I4P)             :: k
   y = 0._R_P ; t = 1._R_P
   do k=1, 30
      t = t*x/k
      if (k > kind) y = y + t
   enddo
   endfunction f
endprogram seth_test_critical_points
