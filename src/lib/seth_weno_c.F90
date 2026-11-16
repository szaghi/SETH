!< SETH, centred WENO: the upwind WENO-JS blended with the central reconstruction of one order more.
module seth_weno_c
!< SETH, centred WENO: the upwind WENO-JS blended with the central reconstruction of one order more.
!<
!< ADAM's weno_reconstruct_centered (adam_weno_object.F90, WENO_C_4/6/8), with its defect D1 fixed (the left state's
!< weights were computed on the unshifted cells, one cell off: notes/adam-weno-audit.md). With S stencils, on the big
!< stencil of cells 1-S..S around the face:
!<
!<    r_up = WENO-JS of order 2S-1 (seth_weno<2S-1>_js, on its own left-biased cells: D1 fixed by construction)
!<    r_c  = sum_m C(m) v(m), the central reconstruction of order 2S from cell averages (CEN<2S>_C)
!<    I_k  = (rmu (v(k) - v(k-1)))^2 + eps, k = 0, 1, 2 (cells -1..2)
!<    phi  = max(|I_1 - I_0|/(I_1 + I_0), |I_2 - I_1|/(I_2 + I_1))
!<    r    = r_c + phi (r_up - r_c)
!<
!< On smooth data phi = O(h) and the error is O(h^2S); at a jump inside cells -1..2 phi tends to 1 and r to the upwind
!< state. The right state is the same on the mirrored vm (the central row is symmetric, and mirroring swaps I_0 and
!< I_2, which phi does not see). `rmu` is the inverse descaler (1: ADAM's form, bitwise in the blend), passed to the
!< upwind weights and to I_k alike; `eps` > 0 (SETH_WENO_JS_EPS), the regulariser of both. ADAM's S = 1 (WENO_C_2)
!< never reaches the blend and is WENO1: use seth_weno1_js. Same convention as the other kernels: vp, vm by offset from
!< the face. Pure, `acc routine seq`, `omp declare target`:
!<
!<    call seth_weno6_c(vp, vm, eps, rmu, rp, rm)          ! vp, vm(-2:3)
use penf, only : I4P, R_P
use seth_weno_js, only : seth_weno3_js, seth_weno5_js, seth_weno7_js
use seth_z_teno_tables, only : CEN4_LO, CEN4_HI, CEN4_C, CEN6_LO, CEN6_HI, CEN6_C, CEN8_LO, CEN8_HI, CEN8_C
implicit none
private
public :: seth_weno4_c
public :: seth_weno6_c
public :: seth_weno8_c

contains
   pure subroutine seth_weno4_c(vp, vm, eps, rmu, rp, rm)
   !< Centred WENO of order 4 (S = 2): WENO3-JS blended with the fourth-order central reconstruction, cells -1..2.
   real(R_P), intent(in)  :: vp(-1:2) !< Values reconstructed with a left bias.
   real(R_P), intent(in)  :: vm(-1:2) !< Values reconstructed with a right bias.
   real(R_P), intent(in)  :: eps      !< Regulariser of the indicators (SETH_WENO_JS_EPS).
   real(R_P), intent(in)  :: rmu      !< Inverse descaler (1: ADAM's form).
   real(R_P), intent(out) :: rp, rm   !< Left and right states at the face.
   real(R_P), parameter   :: C(CEN4_LO:CEN4_HI) = CEN4_C
   real(R_P)              :: cp, cm
   integer(I4P)           :: m
   !$acc routine seq
   !$omp declare target

   call seth_weno3_js(vp, vm, eps, rmu, rp, rm)
   cp = 0._R_P ; cm = 0._R_P
   do m=CEN4_LO, CEN4_HI
      cp = cp + C(m)*vp(m)
      cm = cm + C(m)*vm(1-m)
   enddo
   rp = cp + blend(vp(-1), vp(0), vp(1), vp(2), eps, rmu)*(rp - cp)
   rm = cm + blend(vm(2), vm(1), vm(0), vm(-1), eps, rmu)*(rm - cm)
   endsubroutine seth_weno4_c

   pure subroutine seth_weno6_c(vp, vm, eps, rmu, rp, rm)
   !< Centred WENO of order 6 (S = 3): WENO5-JS blended with the sixth-order central reconstruction, cells -2..3.
   real(R_P), intent(in)  :: vp(-2:3) !< Values reconstructed with a left bias.
   real(R_P), intent(in)  :: vm(-2:3) !< Values reconstructed with a right bias.
   real(R_P), intent(in)  :: eps      !< Regulariser of the indicators (SETH_WENO_JS_EPS).
   real(R_P), intent(in)  :: rmu      !< Inverse descaler (1: ADAM's form).
   real(R_P), intent(out) :: rp, rm   !< Left and right states at the face.
   real(R_P), parameter   :: C(CEN6_LO:CEN6_HI) = CEN6_C
   real(R_P)              :: cp, cm
   integer(I4P)           :: m
   !$acc routine seq
   !$omp declare target

   call seth_weno5_js(vp, vm, eps, rmu, rp, rm)
   cp = 0._R_P ; cm = 0._R_P
   do m=CEN6_LO, CEN6_HI
      cp = cp + C(m)*vp(m)
      cm = cm + C(m)*vm(1-m)
   enddo
   rp = cp + blend(vp(-1), vp(0), vp(1), vp(2), eps, rmu)*(rp - cp)
   rm = cm + blend(vm(2), vm(1), vm(0), vm(-1), eps, rmu)*(rm - cm)
   endsubroutine seth_weno6_c

   pure subroutine seth_weno8_c(vp, vm, eps, rmu, rp, rm)
   !< Centred WENO of order 8 (S = 4): WENO7-JS blended with the eighth-order central reconstruction, cells -3..4.
   real(R_P), intent(in)  :: vp(-3:4) !< Values reconstructed with a left bias.
   real(R_P), intent(in)  :: vm(-3:4) !< Values reconstructed with a right bias.
   real(R_P), intent(in)  :: eps      !< Regulariser of the indicators (SETH_WENO_JS_EPS).
   real(R_P), intent(in)  :: rmu      !< Inverse descaler (1: ADAM's form).
   real(R_P), intent(out) :: rp, rm   !< Left and right states at the face.
   real(R_P), parameter   :: C(CEN8_LO:CEN8_HI) = CEN8_C
   real(R_P)              :: cp, cm
   integer(I4P)           :: m
   !$acc routine seq
   !$omp declare target

   call seth_weno7_js(vp, vm, eps, rmu, rp, rm)
   cp = 0._R_P ; cm = 0._R_P
   do m=CEN8_LO, CEN8_HI
      cp = cp + C(m)*vp(m)
      cm = cm + C(m)*vm(1-m)
   enddo
   rp = cp + blend(vp(-1), vp(0), vp(1), vp(2), eps, rmu)*(rp - cp)
   rm = cm + blend(vm(2), vm(1), vm(0), vm(-1), eps, rmu)*(rm - cm)
   endsubroutine seth_weno8_c

   pure function blend(a, b, c, d, eps, rmu) result(phi)
   !< ADAM's centred-upwind blending coefficient on four consecutive cells a, b, c, d: in [0, 1), O(h) on smooth data.
   real(R_P), intent(in) :: a, b, c, d !< Consecutive cell values (cells -1..2, or 2..-1 for the mirrored state).
   real(R_P), intent(in) :: eps        !< Regulariser.
   real(R_P), intent(in) :: rmu        !< Inverse descaler.
   real(R_P)             :: phi        !< Weight of the upwind state.
   real(R_P)             :: i0, i1, i2
   !$acc routine seq
   !$omp declare target

   i0 = (rmu*(b - a))**2 + eps
   i1 = (rmu*(c - b))**2 + eps
   i2 = (rmu*(d - c))**2 + eps
   phi = max(abs(i1 - i0)/(i1 + i0), abs(i2 - i1)/(i2 + i1))
   endfunction blend
endmodule seth_weno_c
