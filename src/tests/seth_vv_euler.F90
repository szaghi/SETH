!< SETH verification: a minimal 1-D Euler solver around the SETH kernels (test harness, not library API).
module seth_vv_euler
!< SETH verification: a minimal 1-D Euler solver around the SETH kernels (test harness, not library API).
!<
!< Ideal gas, gamma = 1.4, conservative variables q = (rho, rho u, E) on N cells of width h with NG = 6 ghost cells
!< (the -5:6 footprint of the composite kernels of issue #6; WENO-9 needs 5), point values at the cell centres
!< (conservative finite differences: the fluxes at the faces, du/dt = -(F_{i+1/2} - F_{i-1/2})/h). Three space
!< discretizations, the two ways SETH's kernels are used and a variant of the first:
!<
!< - PATH_SPLIT, the reconstruction kernels: flux splitting (Shu & Osher 1989; Jiang & Shu 1996) on characteristic
!<   variables. At each face the Roe-averaged left and right eigenvectors L, R project the stencil's states and fluxes,
!<   w = L q, g = L f; global Lax-Friedrichs splitting per field, g+- = (g +- alpha_k w)/2, alpha_k = max |lambda_k|
!<   over the domain times 1.1; g+ is reconstructed with the left bias (vp), g- with the right bias (vm), and
!<   F = R (rp + rm). The 10% margin keeps the split fluxes free of degenerate critical points: with alpha = max|lambda|
!<   exactly, d g-/dx = (lambda - alpha) dw/dx vanishes to second order where |lambda| peaks, and TENO drops to second
!<   order at second-order critical points (seth_test_critical_points);
!< - PATH_INTERP, the interpolation kernels: the alternative finite-difference formulation (Shu & Osher 1988; Jiang,
!<   Shu & Zhang 2013). The characteristic variables w = L q are interpolated to the face (left and right states), a
!<   Riemann solver (HLLC, Toro et al. 1994 with Davis-Einfeldt speeds) gives F0, and the flux is corrected to high
!<   order with the even derivatives of the flux, F = F0 - h^2/24 f_xx + 7 h^4/5760 f_xxxx, f_xx h^2 and f_xxxx h^4
!<   from the six node fluxes around the face (exact coefficients: O(h^4) and O(h^2), so the correction errs by
!<   O(h^6)). A face state with rho <= 0 or p <= 0 falls back to the first-order states (counted in `nfallback`).
!< - PATH_ROE, the reconstruction kernels with Roe's flux splitting with entropy fix (issue #5; the "RF" of Jiang & Shu
!<   1996, Shu & Osher 1989 Algorithm 2.3 in characteristic fields, notes/issue5-roe-sources.md): w and g as on the
!<   splitting path; per field k, when lambda_k of cell i, of cell i+1 and of the Roe average all have one sign, g is
!<   reconstructed unsplit with the upwind bias of that sign; otherwise (a sign change or a zero: the sonic case) local
!<   Lax-Friedrichs, g+- = (g +- alpha w)/2, alpha = max(|lambda_k(q_i)|, |lambda_k(q_{i+1})|) (Shu & Osher Eq. (4.4)),
!<   reconstructed as on the splitting path; such fields are counted in `nsonic`. With `roe_fix` off, every field is
!<   upwinded by the sign of the Roe eigenvalue (Algorithm 2.1, no entropy fix: the negative control).
!<
!< Positivity: with `limiter` on, the flux of every face is blended with the first-order Lax-Friedrichs flux,
!< F = theta F_high + (1 - theta) F_LF, theta in [0, 1] as large as keeps rho and p of the two one-sided states
!< q_i - 2 lambda F and q_{i+1} + 2 lambda F above 1e-13 (lambda = dt/h): the positivity-preserving flux limiter of
!< Hu, Adams & Shu (2013) for finite-difference WENO, SETH's kernel seth_positivity_euler (theta in closed form). Each stage is then a convex combination of
!< positive states when dt max(|u| + c) <= h/2 (CFL <= 0.5), and SSP-RK3 keeps that. Every face with theta < 1 is
!< counted in `nlimited`: the cases assert that the limiter never acts on smooth flows, and report how often each
!< scheme needs it on the stress tests.
!<
!< Boundaries: periodic, transmissive (zeroth-order extrapolation), reflecting walls, or inflow (BC_INFLOW: the left
!< ghost cells keep the state euler_run gave them, the first cell's initial state; the right end transmissive). A
!< supersonic inflow needs BC_INFLOW under a central flux: the central stencil reads cells downstream, and a
!< transmissive end copies what reaches the first cell into the inflow state, which then drifts (Shu-Osher's cost
!< hybrid of TENO6: from round-off to 0.1-0.15 in density by t = 1.8). Time: SSP-RK3 (Shu & Osher 1988).
!<
!< Hybrids (PATH_SPLIT only; roadmap issue #4, steps 3.4 and 4.3), set by the module variable `hybrid`: the density of
!< cells i-2..i+3 feeds the sensor `hyb_sensor` (default seth_sensor_lcr, Li, Chen & Ren 2022; or seth_sensor_teno,
!< _ren, _tau, at their published thresholds), normalised by max(|rho|, 1); a face whose sensor is
!< below 1 is flagged, and the flags are dilated by `hyb_buffer` faces on each side. Flagged faces take the
!< characteristic flux-split nonlinear scheme as above; the others
!< - HYB_COST: the linear background of the scheme on component-wise Lax-Friedrichs-split fluxes (alpha the largest of
!<   the characteristic alphas): no projection, no weights;
!< - HYB_KEP: the sixth-order KEEP-PE central flux (seth_keep_pe6) of the primitive state;
!< - HYB_LEFT_KEP: KEEP-PE on the faces left of `hyb_split` (a fraction of the domain), the nonlinear scheme right of
!<   it, regardless of the sensor (the junction experiment J1);
!< - HYB_RANDOM: KEEP-PE or the nonlinear scheme at random faces (J3, conservation).
!< `nflagged` and `nhybfaces` count the flagged and the evaluated faces.
use seth
use seth_test_schemes, only : reconstruct, VLO, VHI, SCHEMES, rbf_dx
implicit none
private
public :: GAMMA, NG, PATH_SPLIT, PATH_INTERP, PATH_ROE, BC_PERIODIC, BC_TRANSMISSIVE, BC_REFLECTIVE, BC_INFLOW
public :: prim2cons, cons2prim, euler_dt, euler_step, nfallback, corrections, limiter, nlimited, positive_state
public :: roe_fix, nsonic, stored_masks, count_activations, nactivations
public :: SENS_TENO, SENS_REN, SENS_LCR, SENS_TAU, hyb_sensor
public :: HYB_NONE, HYB_COST, HYB_KEP, HYB_LEFT_KEP, HYB_RANDOM, hybrid, hyb_buffer, hyb_split, nflagged, nhybfaces

real(R_P),    parameter :: GAMMA = 1.4_R_P
integer(I4P), parameter :: NG = 6_I4P
integer(I4P), parameter :: PATH_SPLIT = 1, PATH_INTERP = 2, PATH_ROE = 3
integer(I4P), parameter :: BC_PERIODIC = 1, BC_TRANSMISSIVE = 2, BC_REFLECTIVE = 3, BC_INFLOW = 4
! f'' h^2 and f'''' h^4 at the face from the node fluxes of cells -2..3 (exact rationals)
real(R_P), parameter :: C2(-2:3) = [-5._R_P/48._R_P, 13._R_P/16._R_P, -17._R_P/24._R_P, -17._R_P/24._R_P, &
                                    13._R_P/16._R_P, -5._R_P/48._R_P]
real(R_P), parameter :: C4(-2:3) = [0.5_R_P, -1.5_R_P, 1._R_P, 1._R_P, -1.5_R_P, 0.5_R_P]
integer(I4P) :: nfallback = 0 !< Interpolated face states replaced by first-order ones (non-positive rho or p).
logical      :: corrections = .true. !< PATH_INTERP: add the flux corrections (off: the second-order flux F0).
logical      :: limiter = .false.    !< Positivity-preserving flux limiter (Hu, Adams & Shu 2013).
integer(I4P) :: nlimited = 0         !< Faces where the limiter blended in the first-order flux (theta < 1).
logical      :: count_activations = .false. !< Count, on the splitting path, the states of a composite (`wide`) scheme
                                            !< whose detector fires, per characteristic field (issue #6, N5's table).
integer(I8P) :: nactivations = 0_I8P  !< Those states, summed over the run.
logical      :: roe_fix = .true.     !< PATH_ROE: the entropy fix (local Lax-Friedrichs at sign changes).
integer(I8P) :: nsonic = 0           !< PATH_ROE: characteristic fields evaluated with local Lax-Friedrichs.
logical      :: stored_masks = .false. !< PATH_SPLIT, the composite kernels (registry flag `wide`): the stored-mask
                                       !< variant of notes/kernel-api-design.md, section 11.4, instead of the
                                       !< self-contained kernel. Every cell's masks are computed once, in the basis of
                                       !< that cell's own face (the left-biased mask of cell i in face i+1/2's, the
                                       !< right-biased mask of cell i+1 in face i+1/2's), and read by the faces around
                                       !< it; the default recomputes them in the basis of the face being reconstructed.
integer(I4P), parameter :: HYB_NONE = 0, HYB_COST = 1, HYB_KEP = 2, HYB_LEFT_KEP = 3, HYB_RANDOM = 4
integer(I4P) :: hybrid = HYB_NONE    !< The hybrid of PATH_SPLIT (HYB_NONE: the nonlinear scheme everywhere).
integer(I4P) :: hyb_buffer = 2       !< Faces each side of a flagged one that are flagged too.
real(R_P)    :: hyb_split = 0.5_R_P  !< HYB_LEFT_KEP: the fraction of the faces that take KEEP-PE.
integer(I8P) :: nflagged = 0         !< Faces the hybrid's sensor flagged (nonlinear scheme), summed over residuals.
integer(I8P) :: nhybfaces = 0        !< Faces a hybrid evaluated, summed over residuals.
integer(I8P) :: hyb_seed = 12345     !< State of the LCG of HYB_RANDOM.
integer(I4P), parameter :: SENS_TENO = 1, SENS_REN = 2, SENS_LCR = 3, SENS_TAU = 4
integer(I4P) :: hyb_sensor = SENS_LCR !< The hybrid's sensor (SETH's four kernels; Phase 5 compares them).

contains
   pure function prim2cons(w) result(q)
   !< (rho, u, p) -> (rho, rho u, E).
   real(R_P), intent(in) :: w(3)
   real(R_P)             :: q(3)
   q = [w(1), w(1)*w(2), w(3)/(GAMMA - 1._R_P) + 0.5_R_P*w(1)*w(2)**2]
   endfunction prim2cons

   pure function cons2prim(q) result(w)
   !< (rho, rho u, E) -> (rho, u, p).
   real(R_P), intent(in) :: q(3)
   real(R_P)             :: w(3)
   w(1) = q(1) ; w(2) = q(2)/q(1) ; w(3) = (GAMMA - 1._R_P)*(q(3) - 0.5_R_P*q(2)*w(2))
   endfunction cons2prim

   pure function flux(q) result(f)
   real(R_P), intent(in) :: q(3)
   real(R_P)             :: f(3), w(3)
   w = cons2prim(q)
   f = [q(2), q(2)*w(2) + w(3), (q(3) + w(3))*w(2)]
   endfunction flux

   pure subroutine eigen(ql, qr, l, r, lam)
   !< Left and right eigenvectors of the flux Jacobian at the Roe average of ql and qr, and its eigenvalues.
   real(R_P), intent(in)  :: ql(3), qr(3)
   real(R_P), intent(out) :: l(3,3), r(3,3), lam(3)
   real(R_P)              :: wl(3), wr(3), sl, sr, u, hh, c, b1, b2
   wl = cons2prim(ql) ; wr = cons2prim(qr)
   sl = sqrt(wl(1)) ; sr = sqrt(wr(1))
   u = (sl*wl(2) + sr*wr(2))/(sl + sr)
   hh = (sl*(ql(3) + wl(3))/wl(1) + sr*(qr(3) + wr(3))/wr(1))/(sl + sr)
   c = sqrt(max((GAMMA - 1._R_P)*(hh - 0.5_R_P*u*u), tiny(1._R_P)))
   lam = [u - c, u, u + c]
   r(:,1) = [1._R_P, u - c, hh - u*c]
   r(:,2) = [1._R_P, u,     0.5_R_P*u*u]
   r(:,3) = [1._R_P, u + c, hh + u*c]
   b2 = (GAMMA - 1._R_P)/(c*c) ; b1 = 0.5_R_P*b2*u*u
   l(1,:) = 0.5_R_P*[b1 + u/c, -(b2*u + 1._R_P/c), b2]
   l(2,:) = [1._R_P - b1, b2*u, -b2]
   l(3,:) = 0.5_R_P*[b1 - u/c, -(b2*u - 1._R_P/c), b2]
   endsubroutine eigen

   pure function hllc(ql, qr) result(f)
   !< HLLC flux (Toro, Spruce & Speares 1994), Davis-Einfeldt wave speed estimates.
   real(R_P), intent(in) :: ql(3), qr(3)
   real(R_P)             :: f(3), wl(3), wr(3), cl, cr, sl, sr, sm, qs(3)
   wl = cons2prim(ql) ; wr = cons2prim(qr)
   cl = sqrt(GAMMA*wl(3)/wl(1)) ; cr = sqrt(GAMMA*wr(3)/wr(1))
   sl = min(wl(2) - cl, wr(2) - cr) ; sr = max(wl(2) + cl, wr(2) + cr)
   sm = (wr(3) - wl(3) + ql(2)*(sl - wl(2)) - qr(2)*(sr - wr(2)))/(wl(1)*(sl - wl(2)) - wr(1)*(sr - wr(2)))
   if (sl >= 0._R_P) then
      f = flux(ql)
   elseif (sr <= 0._R_P) then
      f = flux(qr)
   elseif (sm >= 0._R_P) then
      qs = wl(1)*(sl - wl(2))/(sl - sm)*[1._R_P, sm, ql(3)/wl(1) + (sm - wl(2))*(sm + wl(3)/(wl(1)*(sl - wl(2))))]
      f = flux(ql) + sl*(qs - ql)
   else
      qs = wr(1)*(sr - wr(2))/(sr - sm)*[1._R_P, sm, qr(3)/wr(1) + (sm - wr(2))*(sm + wr(3)/(wr(1)*(sr - wr(2))))]
      f = flux(qr) + sr*(qs - qr)
   endif
   endfunction hllc

   subroutine ghosts(q, n, bc)
   !< Fill the NG ghost cells of each end.
   integer(I4P), intent(in)    :: n, bc
   real(R_P),    intent(inout) :: q(3,1-NG:n+NG)
   integer(I4P)                :: g
   do g=1, NG
      select case(bc)
      case(BC_PERIODIC)
         q(:,1-g) = q(:,n+1-g) ; q(:,n+g) = q(:,g)
      case(BC_TRANSMISSIVE)
         q(:,1-g) = q(:,1) ; q(:,n+g) = q(:,n)
      case(BC_REFLECTIVE)
         q(:,1-g) = q(:,g) ; q(2,1-g) = -q(2,g)
         q(:,n+g) = q(:,n+1-g) ; q(2,n+g) = -q(2,n+1-g)
      case(BC_INFLOW) ! the left ghosts are left as they are
         q(:,n+g) = q(:,n)
      endselect
   enddo
   endsubroutine ghosts

   function euler_dt(q, n, h, cfl) result(dt)
   !< dt = cfl h / max(|u| + c).
   integer(I4P), intent(in) :: n
   real(R_P),    intent(in) :: q(3,1-NG:n+NG), h, cfl
   real(R_P)                :: dt, w(3), smax
   integer(I4P)             :: i
   smax = 0._R_P
   do i=1, n
      w = cons2prim(q(:,i))
      smax = max(smax, abs(w(2)) + sqrt(GAMMA*w(3)/w(1)))
   enddo
   dt = cfl*h/smax
   endfunction euler_dt

   subroutine residual(id, path, q, n, h, dt, bc, res, par)
   !< res = -(F_{i+1/2} - F_{i-1/2})/h for i = 1..n (dt: for the positivity limiter only).
   integer(I4P), intent(in)           :: id, path, n, bc
   real(R_P),    intent(inout)        :: q(3,1-NG:n+NG)
   real(R_P),    intent(in)           :: h, dt
   real(R_P),    intent(out)          :: res(3,n)
   real(R_P),    intent(in), optional :: par
   real(R_P)                          :: fhat(3,0:n), fn(3,1-NG:n+NG), l(3,3), r(3,3), lam(3), alpha(3)
   real(R_P)                          :: w(3,VLO:VHI), g(3,VLO:VHI), vp(VLO:VHI), vm(VLO:VHI), rp, rm, fc(3)
   real(R_P)                          :: wl(3), wr(3), ql(3), qr(3), pl(3), pr(3)
   real(R_P)                          :: ws(5,-2:3), fk(5), sv, laml(3), lamr(3), a
   logical                            :: flagged(0:n), raw(0:n)
   real(R_P)                          :: u5(-2:2), v5(-2:2)
   integer(I4P)                       :: msl(3,-3:n+3), msr(3,-2:n+4), mpk(-3:3), mmk(-3:3), lpk, lmk
   integer(I4P)                       :: i, m, k
   rbf_dx = h ! the RBF kernels' s = alpha dx (seth_test_schemes)
   call ghosts(q, n, bc)
   do i=1-NG, n+NG
      fn(:,i) = flux(q(:,i))
   enddo
   alpha = 0._R_P
   if (path == PATH_SPLIT) then
      do i=0, n
         call eigen(q(:,i), q(:,i+1), l, r, lam)
         alpha = max(alpha, abs(lam))
      enddo
      do i=1, n ! the cell eigenvalues too: the splitting must dominate every |lambda| on the stencils
         pl = cons2prim(q(:,i))
         alpha = max(alpha, abs(pl(2)) + sqrt(GAMMA*pl(3)/pl(1)))
      enddo
      alpha = 1.1_R_P*alpha ! margin: alpha = max|lambda| exactly makes g- = (lambda - alpha) w/2 a second-order
                            ! critical point where |lambda| peaks, where TENO drops to second order
   endif
   if (hybrid /= HYB_NONE) then
      do i=0, n
         select case(hybrid)
         case(HYB_LEFT_KEP)
            raw(i) = real(i, R_P) > hyb_split*n
         case(HYB_RANDOM)
            hyb_seed = mod(48271_I8P*hyb_seed, 2147483647_I8P) ! Park & Miller's minimal standard generator
            raw(i) = hyb_seed > 1073741823_I8P
         case default
            vp(-2:3) = q(1,i-2:i+3)
            select case(hyb_sensor)
            case(SENS_TENO)
               call seth_sensor_teno(vp(-2:3), 1._R_P/max(maxval(abs(vp(-2:3))), 1._R_P), SETH_SENSOR_TENO_CT, sv)
            case(SENS_REN)
               call seth_sensor_ren(vp(-2:3), 1._R_P/max(maxval(abs(vp(-2:3))), 1._R_P), SETH_SENSOR_REN_PSIC, sv)
            case(SENS_TAU)
               call seth_sensor_tau(vp(-2:3), 1._R_P/max(maxval(abs(vp(-2:3))), 1._R_P), SETH_SENSOR_TAU_ETA, sv)
            case default
               call seth_sensor_lcr(vp(-2:3), 1._R_P/max(maxval(abs(vp(-2:3))), 1._R_P), SETH_SENSOR_LCR_PSIC, sv)
            endselect
            raw(i) = sv < 1._R_P
         endselect
      enddo
      if (bc == BC_PERIODIC) raw(n) = raw(0) ! faces 0 and n are the same face: one flux, or mass leaks
      if (hybrid == HYB_COST .or. hybrid == HYB_KEP) then
         do i=0, n
            flagged(i) = any(raw(max(0, i-hyb_buffer):min(n, i+hyb_buffer)))
         enddo
      else
         flagged = raw
      endif
      nhybfaces = nhybfaces + n + 1 ; nflagged = nflagged + count(flagged)
   endif
   if (stored_masks .and. path == PATH_SPLIT .and. SCHEMES(id)%wide) then ! pass 1 of the stored-mask variant
      do i=-3, n+3
         call eigen(q(:,i), q(:,i+1), l, r, lam)
         do m=-2, 3
            w(:,m) = matmul(l, q(:,i+m))
            g(:,m) = matmul(l, fn(:,i+m))
         enddo
         do k=1, 3
            do m=-2, 2
               u5(m) = 0.5_R_P*(g(k,m) + alpha(k)*w(k,m))          ! the left-biased mask of cell i
               v5(m) = 0.5_R_P*(g(k,1-m) - alpha(k)*w(k,1-m))      ! the right-biased mask of cell i+1, mirrored
            enddo
            call seth_teno5_mask(u5, SCHEMES(id)%eps, SCHEMES(id)%par, msl(k,i))
            call seth_teno5_mask(v5, SCHEMES(id)%eps, SCHEMES(id)%par, msr(k,i+1))
         enddo
      enddo
   endif
   do i=0, n ! face i+1/2: stencil cell m is cell i+m
      if (hybrid /= HYB_NONE) then
         if (.not.flagged(i)) then
            if (hybrid == HYB_COST) then ! component-wise split fluxes, the linear background
               do k=1, 3
                  vp = 0.5_R_P*(fn(k,i+VLO:i+VHI) + maxval(alpha)*q(k,i+VLO:i+VHI))
                  vm = 0.5_R_P*(fn(k,i+VLO:i+VHI) - maxval(alpha)*q(k,i+VLO:i+VHI))
                  call linear_background(id, vp, vm, rp, rm)
                  fhat(k,i) = rp + rm
               enddo
            else ! KEEP-PE of the primitive state
               do m=-2, 3
                  pl = cons2prim(q(:,i+m))
                  ws(:,m) = [pl(1), pl(2), 0._R_P, 0._R_P, pl(3)]
               enddo
               call seth_keep_pe6(1_I4P, GAMMA, ws, fk)
               fhat(:,i) = [fk(1), fk(2), fk(5)]
            endif
            cycle
         endif
      endif
      call eigen(q(:,i), q(:,i+1), l, r, lam)
      do m=VLO, VHI
         w(:,m) = matmul(l, q(:,i+m))
         g(:,m) = matmul(l, fn(:,i+m))
      enddo
      if (path == PATH_SPLIT) then
         do k=1, 3
            vp = 0.5_R_P*(g(k,:) + alpha(k)*w(k,:))
            vm = 0.5_R_P*(g(k,:) - alpha(k)*w(k,:))
            if (stored_masks .and. SCHEMES(id)%wide) then ! pass 2: TENO5, the stored masks, detector, override
               call seth_teno5(vp(-2:3), vm(-2:3), SCHEMES(id)%eps, SCHEMES(id)%par, rp, rm)
               do m=-3, 3
                  mpk(m) = msl(k,i+m) ; mmk(m) = msr(k,i+1-m)
               enddo
               if (trim(SCHEMES(id)%name) == 'teno5_thinc') then
                  call seth_thinc_detect(mpk, mmk, lpk, lmk)
                  call seth_thinc_override(vp(-1:2), vm(-1:2), lpk, lmk, SETH_THINC_BETA, 1._R_P, rp, rm)
               else
                  call seth_dr_detect(mpk, mmk, lpk, lmk)
                  call seth_dr_override(vp(-1:2), vm(-1:2), lpk, lmk, SETH_DR_T1, SETH_DR_T2, 1._R_P, rp, rm)
               endif
            else
               call reconstruct(id, vp, vm, rp, rm, par=par)
            endif
            if (count_activations .and. SCHEMES(id)%wide) call activations(trim(SCHEMES(id)%name) == 'teno5_thinc', &
                                                                          vp, vm)
            fc(k) = rp + rm
         enddo
         fhat(:,i) = matmul(r, fc)
      elseif (path == PATH_ROE) then
         pl = cons2prim(q(:,i)) ; pr = cons2prim(q(:,i+1))
         a = sqrt(GAMMA*pl(3)/pl(1)) ; laml = [pl(2) - a, pl(2), pl(2) + a]
         a = sqrt(GAMMA*pr(3)/pr(1)) ; lamr = [pr(2) - a, pr(2), pr(2) + a]
         do k=1, 3
            if (roe_fix .and. .not.((laml(k) > 0._R_P .and. lamr(k) > 0._R_P .and. lam(k) > 0._R_P) .or. &
                                    (laml(k) < 0._R_P .and. lamr(k) < 0._R_P .and. lam(k) < 0._R_P))) then
               a = max(abs(laml(k)), abs(lamr(k)))
               vp = 0.5_R_P*(g(k,:) + a*w(k,:))
               vm = 0.5_R_P*(g(k,:) - a*w(k,:))
               call reconstruct(id, vp, vm, rp, rm, par=par)
               fc(k) = rp + rm ; nsonic = nsonic + 1
            else
               vp = g(k,:)
               call reconstruct(id, vp, vp, rp, rm, par=par)
               fc(k) = merge(rp, rm, lam(k) >= 0._R_P)
            endif
         enddo
         fhat(:,i) = matmul(r, fc)
      else
         do k=1, 3
            vp = w(k,:)
            call reconstruct(id, vp, vp, wl(k), wr(k), par=par)
         enddo
         ql = matmul(r, wl) ; qr = matmul(r, wr)
         pl = cons2prim(ql) ; pr = cons2prim(qr)
         if (.not.(pl(1) > 0._R_P .and. pl(3) > 0._R_P .and. pr(1) > 0._R_P .and. pr(3) > 0._R_P)) then
            ql = q(:,i) ; qr = q(:,i+1) ; nfallback = nfallback + 1
         endif
         fhat(:,i) = hllc(ql, qr)
         if (corrections) then
         do m=-2, 3
            fhat(:,i) = fhat(:,i) - (C2(m)/24._R_P - 7._R_P*C4(m)/5760._R_P)*fn(:,i+m)
         enddo
         endif
      endif
   enddo
   if (limiter) call limit(q, n, h, dt, fhat)
   do i=1, n
      res(:,i) = -(fhat(:,i) - fhat(:,i-1))/h
   enddo
   endsubroutine residual

   subroutine linear_background(id, vp, vm, rp, rm)
   !< The linear background of scheme id (its order and family), the smooth branch of HYB_COST.
   integer(I4P), intent(in)  :: id
   real(R_P),    intent(in)  :: vp(VLO:VHI), vm(VLO:VHI)
   real(R_P),    intent(out) :: rp, rm
   select case(trim(SCHEMES(id)%name))
   case('weno3_js', 'wenoz3')                       ; call seth_linear3(vp(-1:2), vm(-1:2), rp, rm)
   case('teno6', 'teno6_opt', 'teno6_neno6')        ; call seth_linear6(vp(-2:3), vm(-2:3), rp, rm)
   case('weno7_js', 'wenoz7')                       ; call seth_linear7(vp(-3:4), vm(-3:4), rp, rm)
   case('teno8', 'teno8_a')                         ; call seth_linear8(vp(-3:4), vm(-3:4), rp, rm)
   case('weno9_js', 'wenoz9')                       ; call seth_linear9(vp(-4:5), vm(-4:5), rp, rm)
   case default                                     ; call seth_linear5(vp(-2:3), vm(-2:3), rp, rm)
   endselect
   endsubroutine linear_background

   subroutine limit(q, n, h, dt, fhat)
   !< The positivity-preserving flux limiter of SETH (seth_positivity_euler, Hu, Adams & Shu 2013) at every face: the
   !< low-order flux is Lax-Friedrichs with alpha = max(|u| + c) over the domain, lam = 2 dt/h, floors 1e-13.
   integer(I4P), intent(in)    :: n
   real(R_P),    intent(in)    :: q(3,1-NG:n+NG), h, dt
   real(R_P),    intent(inout) :: fhat(3,0:n)
   real(R_P)                   :: lam, a, w(3), flf(3), f(3), theta
   integer(I4P)                :: i
   lam = 2._R_P*dt/h
   a = 0._R_P
   do i=0, n+1
      w = cons2prim(q(:,i))
      a = max(a, abs(w(2)) + sqrt(max(GAMMA*w(3)/w(1), 0._R_P)))
   enddo
   do i=0, n
      flf = 0.5_R_P*(flux(q(:,i)) + flux(q(:,i+1))) - 0.5_R_P*a*(q(:,i+1) - q(:,i))
      call seth_positivity_euler(3_I4P, GAMMA, lam, SETH_POSITIVITY_EPS, SETH_POSITIVITY_EPS, q(:,i), q(:,i+1), &
                                 fhat(:,i), flf, f, theta)
      if (theta < 1._R_P) then
         fhat(:,i) = f ; nlimited = nlimited + 1
      endif
   enddo
   endsubroutine limit

   subroutine euler_step(id, path, q, n, h, dt, bc, ok, par)
   !< One SSP-RK3 step; `ok` false (and q left at the failed stage) when a stage has a cell with rho <= 0 or p <= 0:
   !< the next residual would take the square root of a negative pressure (a floating-point trap in the debug builds).
   integer(I4P), intent(in)           :: id, path, n, bc
   real(R_P),    intent(inout)        :: q(3,1-NG:n+NG)
   real(R_P),    intent(in)           :: h, dt
   logical,      intent(out)          :: ok
   real(R_P),    intent(in), optional :: par
   real(R_P)                          :: q1(3,1-NG:n+NG), q2(3,1-NG:n+NG), res(3,n)
   call residual(id, path, q, n, h, dt, bc, res, par)
   q1 = q ; q1(:,1:n) = q(:,1:n) + dt*res
   ok = positive_state(q1, n) ; if (.not.ok) then ; q = q1 ; return ; endif
   call residual(id, path, q1, n, h, dt, bc, res, par)
   q2 = q ; q2(:,1:n) = 0.75_R_P*q(:,1:n) + 0.25_R_P*(q1(:,1:n) + dt*res)
   ok = positive_state(q2, n) ; if (.not.ok) then ; q = q2 ; return ; endif
   call residual(id, path, q2, n, h, dt, bc, res, par)
   q(:,1:n) = q(:,1:n)/3._R_P + 2._R_P/3._R_P*(q2(:,1:n) + dt*res)
   ok = positive_state(q, n)
   endsubroutine euler_step

   pure function positive_state(q, n) result(ok)
   !< rho > 0 and p > 0 in every interior cell (no division: rho is tested first, p as rho E - (rho u)^2/2 > 0).
   integer(I4P), intent(in) :: n
   real(R_P),    intent(in) :: q(3,1-NG:n+NG)
   logical                  :: ok
   integer(I4P)             :: i
   ok = .true.
   do i=1, n
      if (.not.(q(1,i) > 0._R_P)) then ; ok = .false. ; return ; endif
      if (.not.(q(1,i)*q(3,i) - 0.5_R_P*q(2,i)**2 > 0._R_P)) then ; ok = .false. ; return ; endif
   enddo
   endfunction positive_state
   subroutine activations(thinc, vp, vm)
   !< Add to nactivations the states (0, 1 or 2) of the face whose detector fires on vp, vm (TENO5's masks of the cells
   !< -3..3, the right state's mirrored), the detector of the scheme: THINC's or DR's (DR, ST+DR, RBF+DR).
   logical,   intent(in) :: thinc
   real(R_P), intent(in) :: vp(VLO:VHI), vm(VLO:VHI)
   real(R_P)             :: u(-2:2)
   integer(I4P)          :: mp(-3:3), mm(-3:3), c, m, lp, lm
   do c=-3, 3
      do m=-2, 2
         u(m) = vp(c+m)
      enddo
      call seth_teno5_mask(u, SETH_TENO_EPS, SETH_TENO5_CT, mp(c))
      do m=-2, 2
         u(m) = vm(1-c-m)
      enddo
      call seth_teno5_mask(u, SETH_TENO_EPS, SETH_TENO5_CT, mm(c))
   enddo
   if (thinc) then
      call seth_thinc_detect(mp, mm, lp, lm)
   else
      call seth_dr_detect(mp, mm, lp, lm)
   endif
   nactivations = nactivations + merge(1_I8P, 0_I8P, lp > 0) + merge(1_I8P, 0_I8P, lm > 0)
   endsubroutine activations

endmodule seth_vv_euler
