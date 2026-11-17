!< SETH, non-polynomial reconstructions of the TENO line (issue #6): the flags, detectors and overrides of TENO5-THINC,
!< TENO5-DR and the Gaussian-RBF TENO.
module seth_nonpoly
!< SETH, non-polynomial reconstructions of the TENO line (issue #6): the flags, detectors and overrides of TENO5-THINC,
!< TENO5-DR and the Gaussian-RBF TENO.
!<
!< The three schemes share one structure (notes/kernel-api-design.md, section 11): a flag per cell (TENO5's cutoff
!< outcome), a detector reading the flags of the neighbouring cells, an override of the face state where the detector
!< fires. This module holds the pieces; every kernel is pure, `acc routine seq` and `omp declare target`.
!<
!<    call seth_teno5_mask(v, eps, ct, m)
!<    call seth_dr_detect(mp, mm, lp, lm)
!<    call seth_thinc_detect(mp, mm, lp, lm)
!<    call seth_thinc_override(vp, vm, lp, lm, beta, rmu, rp, rm)
!<    call seth_dr_override(vp, vm, lp, lm, t1, t2, rmu, rp, rm)
!<    call seth_teno5_dr(vp, vm, eps, ct, t1, t2, rmu, rp, rm)       ! vp, vm(-5:6)
!<    call seth_teno5_thinc(vp, vm, eps, ct, beta, rmu, rp, rm)      ! vp, vm(-5:6)
!<
!< `v(-2:2)` the five values centred on one cell, upwind to the left (the caller mirrors the data for the other bias),
!< `eps` and `ct` TENO5's regulariser and cutoff (the papers: SETH_TENO_EPS, C_T = 1e-5, SETH_TENO5_CT), `m` the
!< mask delta_0 + 2 delta_1 + 4 delta_2: delta_k = 1 where candidate k of seth_teno5 passes the cutoff (S0 = {-1,0,1}
!< central, S1 = {0,1,2} downwind, S2 = {-2,-1,0} upwind, the papers' order). The flag of a cell is the cutoff of the
!< left-biased TENO5 reconstruction at its right face. The detectors read the masks of the cells around the face i+1/2
!< by offset from it, mp(-3:3) for the left state (cell 0 left of the face) and mm(-3:3) for the right state (the
!< masks of the mirrored data w(m) = vm(1-m), so mm(c) is the mask of cell 1-c computed on mirrored values), and return
!< the levels lp, lm: DR 0 (none), 1 (Liang, Xiao, Shyy & Fu 2022, Eq. 24, steepness 0.91), 2 (Eq. 32 too, 1.23);
!< THINC 0 or 1 (Takagi et al. 2022, Eq. 20). zeta (seth_teno5_zeta) is Eq. 17 / Takagi's Eq. 19. DR reads the
!< masks -3..2, THINC -3..3. The overrides take the four cells vp, vm(-1:2) and the levels, and replace the states
!< rp, rm (inout) where the level is positive: by the face value of Takagi's THINC function (Eq. 16) or by Liang et
!< al.'s two-steepness interpolation (Eq. 31, tanh(0.91) at level 1, tanh(1.23) at level 2), both with the first-order
!< guard: where (f_{i+1}-f_i)(f_i-f_{i-1}) rmu^2 <= SETH_NONPOLY_GUARD the state is the cell value f_i (Takagi Eq. 17,
!< Liang Eq. 33; rmu = 1 the papers' form). The self-contained kernels seth_teno5_dr and seth_teno5_thinc are the
!< composite schemes in one call per face (notes/kernel-api-design.md, section 11.3): TENO5's two states, the masks
!< of the neighbouring cells recomputed from the twelve values vp, vm(-5:6) (in the basis the caller projected them
!< into: the face's own), the detector, the override. The face's own cell is tested first and the downwind cells next;
!< detection is ruled out as soon as it cannot hold, without changing the result (an early exit). Sources:
!< notes/m5-sources.md.
use penf, only : I4P, R_P
use seth_z_teno_tables
use seth_nonpoly_tables, only : SETH_DR_T1, SETH_DR_T2
use seth_kind_parameters, only : SETH_NONPOLY_GUARD
use seth_teno, only : seth_teno5
implicit none
private
public :: seth_teno5_mask, seth_teno5_zeta, seth_dr_detect, seth_thinc_detect
public :: seth_thinc_override, seth_dr_override
public :: seth_teno5_dr, seth_teno5_thinc
public :: SETH_DR_T1, SETH_DR_T2, SETH_THINC_BETA, SETH_NONPOLY_GUARD

real(R_P), parameter :: SETH_THINC_BETA = 1.6_R_P !< Steepness of the THINC function (Takagi et al. 2022, Eq. 15).

contains
   pure subroutine seth_teno5_mask(v, eps, ct, mk)
   !< The mask of one cell: TENO5's indicators, tau_5 = |beta_1 - beta_2| and cutoff, no candidate values (the body of
   !< seth_teno5, seth_z_teno_body.INC, with SETH_TENO_MASK: the flags are the kernel's own).
   integer(I4P), parameter :: NC = TENO5_NC, LO = TENO5_LO, UMAX = TENO5_UMAX, MAXW = TENO5_MAXW
   integer(I4P), parameter :: CLO(0:NC-1) = TENO5_CLO, CW(0:NC-1) = TENO5_CW
   integer(I4P), parameter :: GLO = TENO5_GLO, GW = TENO5_GW
   real(R_P),    parameter :: CL(0:MAXW-1,0:MAXW-1,0:NC-1) = TENO5_L, CDQ(0:MAXW-1,0:NC-1) = TENO5_DQ
   real(R_P),    parameter :: D(0:NC-1) = TENO5_D, TC(0:NC-1) = TENO5_TC, TG = TENO5_TG
   real(R_P),    parameter :: CLG(0:GW-1,0:GW-1) = TENO5_LG, CDG(0:GW-1) = TENO5_DG
#define SETH_TENO
#define SETH_TENO_MASK
#include "seth_z_teno_body.INC"
#undef SETH_TENO_MASK
#undef SETH_TENO
   endsubroutine seth_teno5_mask

   pure function seth_teno5_zeta(mk) result(z)
   !< zeta of a mask (Liang et al. 2022, Eq. 17; Takagi et al. 2022, Eq. 19): 2 where the downwind and the upwind
   !< candidates are both cut, 1 where only the downwind one is, -1 where only the upwind one is, 0 otherwise.
   integer(I4P), intent(in) :: mk !< Mask delta_0 + 2 delta_1 + 4 delta_2.
   integer(I4P)             :: z  !< zeta.
   integer(I4P)             :: d1, d2
   !$acc routine seq
   !$omp declare target
   d1 = ibits(mk, 1, 1) ; d2 = ibits(mk, 2, 1)
   if (d2 == 0 .and. d1 == 0) then
      z = 2
   elseif (d1 == 0 .and. d2 == 1) then
      z = 1
   elseif (d2 == 0 .and. d1 == 1) then
      z = -1
   else
      z = 0
   endif
   endfunction seth_teno5_zeta

   pure subroutine seth_dr_detect(mp, mm, lp, lm)
   !< DR detection levels of the two states (Liang et al. 2022, Eqs. 24 and 32).
   integer(I4P), intent(in)  :: mp(-3:3) !< Masks of the cells -3..3 for the left state (-3..2 read).
   integer(I4P), intent(in)  :: mm(-3:3) !< Masks of the mirrored cells for the right state (-3..2 read).
   integer(I4P), intent(out) :: lp       !< Level of the left state: 0, 1 or 2.
   integer(I4P), intent(out) :: lm       !< Level of the right state.
   !$acc routine seq
   !$omp declare target
   lp = dr_level(mp)
   lm = dr_level(mm)
   endsubroutine seth_dr_detect

   pure subroutine seth_thinc_detect(mp, mm, lp, lm)
   !< THINC detection levels of the two states (Takagi et al. 2022, Eq. 20; Liang et al. 2022, Eq. 18).
   integer(I4P), intent(in)  :: mp(-3:3) !< Masks of the cells -3..3 for the left state.
   integer(I4P), intent(in)  :: mm(-3:3) !< Masks of the mirrored cells for the right state.
   integer(I4P), intent(out) :: lp       !< Level of the left state: 0 or 1.
   integer(I4P), intent(out) :: lm       !< Level of the right state.
   !$acc routine seq
   !$omp declare target
   lp = thinc_level(mp)
   lm = thinc_level(mm)
   endsubroutine seth_thinc_detect

   pure subroutine seth_thinc_override(vp, vm, lp, lm, beta, rmu, rp, rm)
   !< Takagi et al. 2022: where the level is 1, the face value of the THINC function h = f+ + f- tanh(beta (X - d))
   !< whose cell average is f_i (Eqs. 15-16), in the form A = tanh(alpha beta/2)/tanh(beta/2) (odd in alpha: mirror
   !< symmetric), with the first-order guard (Eq. 17). Level 0: the states unchanged.
   real(R_P),    intent(in)    :: vp(-1:2) !< Left-biased values, cells -1..1 read (f_{i-1}, f_i, f_{i+1}).
   real(R_P),    intent(in)    :: vm(-1:2) !< Right-biased values, cells 0..2 read, mirrored.
   integer(I4P), intent(in)    :: lp, lm   !< Levels of the two states (seth_thinc_detect).
   real(R_P),    intent(in)    :: beta     !< Steepness (SETH_THINC_BETA).
   real(R_P),    intent(in)    :: rmu      !< Inverse descaler of the guard (1: the paper's form).
   real(R_P),    intent(inout) :: rp, rm   !< States: replaced where the level is 1.
   !$acc routine seq
   !$omp declare target
   if (lp > 0) rp = thinc_face(vp(-1), vp(0), vp(1), beta, rmu)
   if (lm > 0) rm = thinc_face(vm(2), vm(1), vm(0), beta, rmu)
   endsubroutine seth_thinc_override

   pure subroutine seth_dr_override(vp, vm, lp, lm, t1, t2, rmu, rp, rm)
   !< Liang et al. 2022, Eq. 31: the face value m + h (t + B)/(1 + t B), m and h the mean and half difference of
   !< f_{i+1}, f_{i-1}, B = (2 f_i - f_{i+1} - f_{i-1})/(f_{i+1} - f_{i-1}), t = t1 at level 1, t2 at level 2 (Eq. 32),
   !< with the first-order guard (Eq. 33). Level 0: the states unchanged.
   real(R_P),    intent(in)    :: vp(-1:2) !< Left-biased values, cells -1..1 read (f_{i-1}, f_i, f_{i+1}).
   real(R_P),    intent(in)    :: vm(-1:2) !< Right-biased values, cells 0..2 read, mirrored.
   integer(I4P), intent(in)    :: lp, lm   !< Levels of the two states (seth_dr_detect).
   real(R_P),    intent(in)    :: t1, t2   !< tanh of the two steepness constants (SETH_DR_T1, SETH_DR_T2).
   real(R_P),    intent(in)    :: rmu      !< Inverse descaler of the guard (1: the paper's form).
   real(R_P),    intent(inout) :: rp, rm   !< States: replaced where the level is 1 or 2.
   !$acc routine seq
   !$omp declare target
   if (lp > 0) rp = dr_face(vp(-1), vp(0), vp(1), merge(t2, t1, lp > 1), rmu)
   if (lm > 0) rm = dr_face(vm(2), vm(1), vm(0), merge(t2, t1, lm > 1), rmu)
   endsubroutine seth_dr_override

   pure subroutine seth_teno5_dr(vp, vm, eps, ct, t1, t2, rmu, rp, rm)
   !< TENO5-DR (Liang, Xiao, Shyy & Fu 2022), self-contained: seth_teno5, the DR detector on the masks of the cells -3..2
   !< of each state, the DR override.
   real(R_P), intent(in)  :: vp(-5:6) !< Left-biased values, cells by offset from the face (-5..4 read).
   real(R_P), intent(in)  :: vm(-5:6) !< Right-biased values (-3..6 read).
   real(R_P), intent(in)  :: eps, ct  !< TENO5's regulariser and cutoff (the papers: SETH_TENO_EPS, 1e-5).
   real(R_P), intent(in)  :: t1, t2   !< tanh of the two steepness constants (SETH_DR_T1, SETH_DR_T2).
   real(R_P), intent(in)  :: rmu      !< Inverse descaler of the guard (1: the paper's form).
   real(R_P), intent(out) :: rp, rm   !< Left and right states at the face.
   real(R_P)              :: a(-2:3), b(-2:3), w(-5:5), c(-1:2), d(-1:2)
   integer(I4P)           :: m, lp, lm
   !$acc routine seq
   !$omp declare target
   do m=-2, 3
      a(m) = vp(m) ; b(m) = vm(m)
   enddo
   call seth_teno5(a, b, eps, ct, rp, rm)
   do m=-5, 5
      w(m) = vp(m)
   enddo
   lp = side_level(w, eps, ct, .false.)
   do m=-5, 5
      w(m) = vm(1-m)
   enddo
   lm = side_level(w, eps, ct, .false.)
   if (lp > 0 .or. lm > 0) then
      do m=-1, 2
         c(m) = vp(m) ; d(m) = vm(m)
      enddo
      call seth_dr_override(c, d, lp, lm, t1, t2, rmu, rp, rm)
   endif
   endsubroutine seth_teno5_dr

   pure subroutine seth_teno5_thinc(vp, vm, eps, ct, beta, rmu, rp, rm)
   !< TENO5-THINC (Takagi, Fu, Wakimura & Xiao 2022), self-contained: seth_teno5, the THINC detector on the masks of the
   !< cells -3..3 of each state, the THINC override.
   real(R_P), intent(in)  :: vp(-5:6) !< Left-biased values, cells by offset from the face (-5..5 read).
   real(R_P), intent(in)  :: vm(-5:6) !< Right-biased values (-4..6 read).
   real(R_P), intent(in)  :: eps, ct  !< TENO5's regulariser and cutoff (the papers: SETH_TENO_EPS, 1e-5).
   real(R_P), intent(in)  :: beta     !< Steepness of the THINC function (SETH_THINC_BETA).
   real(R_P), intent(in)  :: rmu      !< Inverse descaler of the guard (1: the paper's form).
   real(R_P), intent(out) :: rp, rm   !< Left and right states at the face.
   real(R_P)              :: a(-2:3), b(-2:3), w(-5:5), c(-1:2), d(-1:2)
   integer(I4P)           :: m, lp, lm
   !$acc routine seq
   !$omp declare target
   do m=-2, 3
      a(m) = vp(m) ; b(m) = vm(m)
   enddo
   call seth_teno5(a, b, eps, ct, rp, rm)
   do m=-5, 5
      w(m) = vp(m)
   enddo
   lp = side_level(w, eps, ct, .true.)
   do m=-5, 5
      w(m) = vm(1-m)
   enddo
   lm = side_level(w, eps, ct, .true.)
   if (lp > 0 .or. lm > 0) then
      do m=-1, 2
         c(m) = vp(m) ; d(m) = vm(m)
      enddo
      call seth_thinc_override(c, d, lp, lm, beta, rmu, rp, rm)
   endif
   endsubroutine seth_teno5_thinc

   pure function side_level(w, eps, ct, thinc) result(l)
   !< The detection level of one state from its left-biased data w(-5:5) (the right state's mirrored): the face's own
   !< cell first (zeta 2 or -1 rules both detectors out), then the downwind cells (one zeta in {-1, 2} among 1..2, DR,
   !< or 1..3, THINC, is necessary), then the upwind ones. Every mask not needed stays 7 (zeta 0), unread by the
   !< conditions that the early exits have already decided.
   real(R_P), intent(in) :: w(-5:5)  !< Left-biased values, cells by offset from the face.
   real(R_P), intent(in) :: eps, ct  !< TENO5's regulariser and cutoff.
   logical,   intent(in) :: thinc    !< THINC's detector (cells -3..3) instead of DR's (-3..2).
   integer(I4P)          :: l        !< Level.
   integer(I4P)          :: mk(-3:3), c, z, chi
   logical               :: dn
   !$acc routine seq
   !$omp declare target
   l = 0
   mk = 7
   mk(0) = cell_mask(w, 0_I4P, eps, ct)
   z = seth_teno5_zeta(mk(0))
   if (z == 2 .or. z == -1) return
   chi = merge(3_I4P, 2_I4P, thinc)
   dn = .false.
   do c=1, chi
      mk(c) = cell_mask(w, c, eps, ct)
      z = seth_teno5_zeta(mk(c))
      dn = dn .or. z == -1 .or. z == 2
   enddo
   if (.not.dn) return
   do c=-3, -1
      mk(c) = cell_mask(w, c, eps, ct)
   enddo
   if (thinc) then
      l = thinc_level(mk)
   else
      l = dr_level(mk)
   endif
   endfunction side_level

   pure function cell_mask(w, c, eps, ct) result(mk)
   !< The mask of the cell c of w: seth_teno5_mask on w(c-2:c+2).
   real(R_P),    intent(in) :: w(-5:5) !< Left-biased values.
   integer(I4P), intent(in) :: c       !< Cell, -3..3.
   real(R_P),    intent(in) :: eps, ct !< TENO5's regulariser and cutoff.
   integer(I4P)             :: mk      !< Mask.
   real(R_P)                :: u(-2:2)
   integer(I4P)             :: m
   !$acc routine seq
   !$omp declare target
   do m=-2, 2
      u(m) = w(c+m)
   enddo
   call seth_teno5_mask(u, eps, ct, mk)
   endfunction cell_mask

   pure function thinc_face(fm1, fi, fp1, beta, rmu) result(f)
   !< Takagi's Eq. 16 with the guard of Eq. 17. Past the guard the triple is strictly monotone, so |alpha| < 1.
   real(R_P), intent(in) :: fm1, fi, fp1 !< f_{i-1}, f_i, f_{i+1}.
   real(R_P), intent(in) :: beta, rmu    !< Steepness, inverse descaler of the guard.
   real(R_P)             :: f            !< Face value at i+1/2.
   real(R_P)             :: fs, fd, a, t
   !$acc routine seq
   !$omp declare target
   if ((rmu*(fp1 - fi))*(rmu*(fi - fm1)) <= SETH_NONPOLY_GUARD) then
      f = fi
   else
      fs = 0.5_R_P*(fm1 + fp1) ; fd = 0.5_R_P*(fp1 - fm1)
      t = tanh(0.5_R_P*beta)
      a = tanh(0.5_R_P*beta*((fi - fs)/fd))/t
      f = fs + fd*((a + t)/(1._R_P + a*t))
   endif
   endfunction thinc_face

   pure function dr_face(fm1, fi, fp1, t, rmu) result(f)
   !< Liang's Eq. 31 with the guard of Eq. 33.
   real(R_P), intent(in) :: fm1, fi, fp1 !< f_{i-1}, f_i, f_{i+1}.
   real(R_P), intent(in) :: t, rmu       !< tanh of the steepness, inverse descaler of the guard.
   real(R_P)             :: f            !< Face value at i+1/2.
   real(R_P)             :: b
   !$acc routine seq
   !$omp declare target
   if ((rmu*(fp1 - fi))*(rmu*(fi - fm1)) <= SETH_NONPOLY_GUARD) then
      f = fi
   else
      b = (2._R_P*fi - (fp1 + fm1))/(fp1 - fm1)
      f = 0.5_R_P*(fp1 + fm1) + 0.5_R_P*(fp1 - fm1)*((t + b)/(1._R_P + t*b))
   endif
   endfunction dr_face

   pure function dr_level(mk) result(l)
   !< Eq. 24 (condition I or II) gives 1; Eq. 32 as well gives 2 (each of its conditions implies Eq. 24's).
   integer(I4P), intent(in) :: mk(-3:3) !< Masks by offset from the face.
   integer(I4P)             :: l        !< Level.
   integer(I4P)             :: z(-3:2), c
   logical                  :: up(-3:2), dn(-3:2) !< zeta in {1, 2} (a jump on the upwind side), in {-1, 2}.
   logical                  :: c1, c2
   !$acc routine seq
   !$omp declare target
   do c=-3, 2
      z(c) = seth_teno5_zeta(mk(c))
      up(c) = z(c) == 1 .or. z(c) == 2
      dn(c) = z(c) == -1 .or. z(c) == 2
   enddo
   c1 = (up(-3) .or. up(-2) .or. up(-1)) .and. z(0) == 0 .and. (dn(1) .or. dn(2))
   c2 = (up(-2) .or. up(-1)) .and. (z(0) == 1 .and. z(1) == -1) .and. dn(2)
   l = 0
   if (c1 .or. c2) then
      l = 1
      if ((up(-2) .and. up(-1)) .and. ((z(0) == 0 .and. dn(1) .and. dn(2)) .or. &
                                       (z(0) == 1 .and. z(1) == -1 .and. dn(2)))) l = 2
   endif
   endfunction dr_level

   pure function thinc_level(mk) result(l)
   !< Takagi et al. 2022, Eq. 20.
   integer(I4P), intent(in) :: mk(-3:3) !< Masks by offset from the face.
   integer(I4P)             :: l        !< Level.
   integer(I4P)             :: z(-3:3), c
   logical                  :: up, dn
   !$acc routine seq
   !$omp declare target
   do c=-3, 3
      z(c) = seth_teno5_zeta(mk(c))
   enddo
   up = .false. ; dn = .false.
   do c=-3, -1
      up = up .or. z(c) == 1 .or. z(c) == 2
   enddo
   do c=1, 3
      dn = dn .or. z(c) == -1 .or. z(c) == 2
   enddo
   l = 0
   if (up .and. z(0) == 0 .and. dn) l = 1
   endfunction thinc_level
endmodule seth_nonpoly
