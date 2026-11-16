!< SETH, linear kernels: the backgrounds of the WENO and TENO schemes, the smooth branch of a hybrid.
module seth_linear
!< SETH, linear kernels: the backgrounds of the WENO and TENO schemes, the smooth branch of a hybrid.
!<
!< Each kernel is the linear scheme sum_k d_k q_k of a nonlinear one, collapsed by the table generator to one row of
!< face coefficients (seth_z_teno_tables, LIN*_C): no indicators, no weights, one dot product per state. It is what
!< TENO returns on smooth data (to round-off) and what WENO-Z approaches, and what a hybrid uses where its sensor says
!< smooth (roadmap issue #4, Phase 3). Same convention as the nonlinear kernels: vp, vm by offset from the face, the
!< right state the left-biased row applied to the mirrored vm. Pure, `acc routine seq`, `omp declare target`:
!<
!<    call seth_linear5(vp, vm, rp, rm)          ! vp, vm(-2:3)
!<
!< seth_linear<2S-1>[_interp] (S = 2..5) are the upwind backgrounds of WENO-(2S-1) (reconstruction and interpolation);
!< seth_linear6[_interp] TENO6's sixth-order one, seth_linear8 TENO8's (Liang & Fu 2024b).
use penf, only : I4P, R8P
use seth_z_teno_tables
implicit none
private
public :: seth_linear3
public :: seth_linear5
public :: seth_linear7
public :: seth_linear9
public :: seth_linear3_interp
public :: seth_linear5_interp
public :: seth_linear7_interp
public :: seth_linear9_interp
public :: seth_linear6
public :: seth_linear6_interp
public :: seth_linear8

contains
   pure subroutine seth_linear3(vp, vm, rp, rm)
   !< Linear background of WENO-3 (3rd order, cells -1..1).
   real(R8P), intent(in)  :: vp(-1:2) !< Values reconstructed with a left bias.
   real(R8P), intent(in)  :: vm(-1:2) !< Values reconstructed with a right bias.
   real(R8P), intent(out) :: rp, rm   !< Left and right states at the face.
   real(R8P), parameter   :: C(LIN3_LO:LIN3_HI) = LIN3_C
   integer(I4P)           :: m
   !$acc routine seq
   !$omp declare target

   rp = 0._R8P ; rm = 0._R8P
   do m=LIN3_LO, LIN3_HI
      rp = rp + C(m)*vp(m)
      rm = rm + C(m)*vm(1-m)
   enddo
   endsubroutine seth_linear3

   pure subroutine seth_linear5(vp, vm, rp, rm)
   !< Linear background of WENO-5 (5th order, cells -2..2).
   real(R8P), intent(in)  :: vp(-2:3) !< Values reconstructed with a left bias.
   real(R8P), intent(in)  :: vm(-2:3) !< Values reconstructed with a right bias.
   real(R8P), intent(out) :: rp, rm   !< Left and right states at the face.
   real(R8P), parameter   :: C(LIN5_LO:LIN5_HI) = LIN5_C
   integer(I4P)           :: m
   !$acc routine seq
   !$omp declare target

   rp = 0._R8P ; rm = 0._R8P
   do m=LIN5_LO, LIN5_HI
      rp = rp + C(m)*vp(m)
      rm = rm + C(m)*vm(1-m)
   enddo
   endsubroutine seth_linear5

   pure subroutine seth_linear7(vp, vm, rp, rm)
   !< Linear background of WENO-7 (7th order, cells -3..3).
   real(R8P), intent(in)  :: vp(-3:4) !< Values reconstructed with a left bias.
   real(R8P), intent(in)  :: vm(-3:4) !< Values reconstructed with a right bias.
   real(R8P), intent(out) :: rp, rm   !< Left and right states at the face.
   real(R8P), parameter   :: C(LIN7_LO:LIN7_HI) = LIN7_C
   integer(I4P)           :: m
   !$acc routine seq
   !$omp declare target

   rp = 0._R8P ; rm = 0._R8P
   do m=LIN7_LO, LIN7_HI
      rp = rp + C(m)*vp(m)
      rm = rm + C(m)*vm(1-m)
   enddo
   endsubroutine seth_linear7

   pure subroutine seth_linear9(vp, vm, rp, rm)
   !< Linear background of WENO-9 (9th order, cells -4..4).
   real(R8P), intent(in)  :: vp(-4:5) !< Values reconstructed with a left bias.
   real(R8P), intent(in)  :: vm(-4:5) !< Values reconstructed with a right bias.
   real(R8P), intent(out) :: rp, rm   !< Left and right states at the face.
   real(R8P), parameter   :: C(LIN9_LO:LIN9_HI) = LIN9_C
   integer(I4P)           :: m
   !$acc routine seq
   !$omp declare target

   rp = 0._R8P ; rm = 0._R8P
   do m=LIN9_LO, LIN9_HI
      rp = rp + C(m)*vp(m)
      rm = rm + C(m)*vm(1-m)
   enddo
   endsubroutine seth_linear9

   pure subroutine seth_linear3_interp(vp, vm, rp, rm)
   !< Linear background of WENO-3 interpolation.
   real(R8P), intent(in)  :: vp(-1:2) !< Values reconstructed with a left bias.
   real(R8P), intent(in)  :: vm(-1:2) !< Values reconstructed with a right bias.
   real(R8P), intent(out) :: rp, rm   !< Left and right states at the face.
   real(R8P), parameter   :: C(LIN3_INTERP_LO:LIN3_INTERP_HI) = LIN3_INTERP_C
   integer(I4P)           :: m
   !$acc routine seq
   !$omp declare target

   rp = 0._R8P ; rm = 0._R8P
   do m=LIN3_INTERP_LO, LIN3_INTERP_HI
      rp = rp + C(m)*vp(m)
      rm = rm + C(m)*vm(1-m)
   enddo
   endsubroutine seth_linear3_interp

   pure subroutine seth_linear5_interp(vp, vm, rp, rm)
   !< Linear background of WENO-5 interpolation.
   real(R8P), intent(in)  :: vp(-2:3) !< Values reconstructed with a left bias.
   real(R8P), intent(in)  :: vm(-2:3) !< Values reconstructed with a right bias.
   real(R8P), intent(out) :: rp, rm   !< Left and right states at the face.
   real(R8P), parameter   :: C(LIN5_INTERP_LO:LIN5_INTERP_HI) = LIN5_INTERP_C
   integer(I4P)           :: m
   !$acc routine seq
   !$omp declare target

   rp = 0._R8P ; rm = 0._R8P
   do m=LIN5_INTERP_LO, LIN5_INTERP_HI
      rp = rp + C(m)*vp(m)
      rm = rm + C(m)*vm(1-m)
   enddo
   endsubroutine seth_linear5_interp

   pure subroutine seth_linear7_interp(vp, vm, rp, rm)
   !< Linear background of WENO-7 interpolation.
   real(R8P), intent(in)  :: vp(-3:4) !< Values reconstructed with a left bias.
   real(R8P), intent(in)  :: vm(-3:4) !< Values reconstructed with a right bias.
   real(R8P), intent(out) :: rp, rm   !< Left and right states at the face.
   real(R8P), parameter   :: C(LIN7_INTERP_LO:LIN7_INTERP_HI) = LIN7_INTERP_C
   integer(I4P)           :: m
   !$acc routine seq
   !$omp declare target

   rp = 0._R8P ; rm = 0._R8P
   do m=LIN7_INTERP_LO, LIN7_INTERP_HI
      rp = rp + C(m)*vp(m)
      rm = rm + C(m)*vm(1-m)
   enddo
   endsubroutine seth_linear7_interp

   pure subroutine seth_linear9_interp(vp, vm, rp, rm)
   !< Linear background of WENO-9 interpolation.
   real(R8P), intent(in)  :: vp(-4:5) !< Values reconstructed with a left bias.
   real(R8P), intent(in)  :: vm(-4:5) !< Values reconstructed with a right bias.
   real(R8P), intent(out) :: rp, rm   !< Left and right states at the face.
   real(R8P), parameter   :: C(LIN9_INTERP_LO:LIN9_INTERP_HI) = LIN9_INTERP_C
   integer(I4P)           :: m
   !$acc routine seq
   !$omp declare target

   rp = 0._R8P ; rm = 0._R8P
   do m=LIN9_INTERP_LO, LIN9_INTERP_HI
      rp = rp + C(m)*vp(m)
      rm = rm + C(m)*vm(1-m)
   enddo
   endsubroutine seth_linear9_interp

   pure subroutine seth_linear6(vp, vm, rp, rm)
   !< Linear background of TENO6 (6th order, cells -2..3).
   real(R8P), intent(in)  :: vp(-2:3) !< Values reconstructed with a left bias.
   real(R8P), intent(in)  :: vm(-2:3) !< Values reconstructed with a right bias.
   real(R8P), intent(out) :: rp, rm   !< Left and right states at the face.
   real(R8P), parameter   :: C(LIN6_LO:LIN6_HI) = LIN6_C
   integer(I4P)           :: m
   !$acc routine seq
   !$omp declare target

   rp = 0._R8P ; rm = 0._R8P
   do m=LIN6_LO, LIN6_HI
      rp = rp + C(m)*vp(m)
      rm = rm + C(m)*vm(1-m)
   enddo
   endsubroutine seth_linear6

   pure subroutine seth_linear6_interp(vp, vm, rp, rm)
   !< Linear background of TENO6 interpolation.
   real(R8P), intent(in)  :: vp(-2:3) !< Values reconstructed with a left bias.
   real(R8P), intent(in)  :: vm(-2:3) !< Values reconstructed with a right bias.
   real(R8P), intent(out) :: rp, rm   !< Left and right states at the face.
   real(R8P), parameter   :: C(LIN6_INTERP_LO:LIN6_INTERP_HI) = LIN6_INTERP_C
   integer(I4P)           :: m
   !$acc routine seq
   !$omp declare target

   rp = 0._R8P ; rm = 0._R8P
   do m=LIN6_INTERP_LO, LIN6_INTERP_HI
      rp = rp + C(m)*vp(m)
      rm = rm + C(m)*vm(1-m)
   enddo
   endsubroutine seth_linear6_interp

   pure subroutine seth_linear8(vp, vm, rp, rm)
   !< Linear background of TENO8, Liang & Fu (6th order, cells -3..4).
   real(R8P), intent(in)  :: vp(-3:4) !< Values reconstructed with a left bias.
   real(R8P), intent(in)  :: vm(-3:4) !< Values reconstructed with a right bias.
   real(R8P), intent(out) :: rp, rm   !< Left and right states at the face.
   real(R8P), parameter   :: C(LIN8_LO:LIN8_HI) = LIN8_C
   integer(I4P)           :: m
   !$acc routine seq
   !$omp declare target

   rp = 0._R8P ; rm = 0._R8P
   do m=LIN8_LO, LIN8_HI
      rp = rp + C(m)*vp(m)
      rm = rm + C(m)*vm(1-m)
   enddo
   endsubroutine seth_linear8

endmodule seth_linear
