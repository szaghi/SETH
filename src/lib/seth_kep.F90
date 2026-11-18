!< SETH, kinetic-energy- and pressure-equilibrium-preserving central fluxes (KEEP-PE) for the Euler equations.
module seth_kep
!< SETH, kinetic-energy- and pressure-equilibrium-preserving central fluxes (KEEP-PE) for the Euler equations.
!<
!< The non-dissipative branch of a hybrid scheme. A kinetic-energy-preserving face flux above second order cannot be a
!< reconstruction followed by a flux function: it is a sum of two-point fluxes of the full state (Coppola & Veldman;
!< llm-wiki weno-teno `properties/kinetic-energy-preservation.md`), so these kernels take the primitive state of the
!< cells around the face and return the flux vector. Pirozzoli's (2010) locally conservative form of order 2L,
!<
!<    F_{i+1/2} = 2 sum_{l=1}^{L} a_l sum_{m=0}^{l-1} f(i-m, i-m+l),
!<
!< a_l the coefficients of the central first derivative of order 2L, with the two-point flux of Shima et al. (2021,
!< KEEP-PE), `<a> = (a_c + a_r)/2`:
!<
!<    mass      f_rho  = <rho> <u_n>
!<    momentum  f_rhou = f_rho <u> + <p> n
!<    energy    f_E    = f_rho (u_c . u_r)/2 + <p> <u_n>/(gamma - 1) + (u_n,c p_r + u_n,r p_c)/2
!<
!< the kinetic energy flux of a kinetic-energy-preserving split form, the internal energy flux of the pressure-
!< equilibrium-preserving scheme (De Michele et al. 2023, Eq. (15): KEEP-PE) and the split pressure work. A uniform
!< velocity and pressure stay so across any density profile; on a periodic domain of uniform pressure the discrete
!< kinetic energy is conserved by the convective terms. `seth_keep_pe6` reproduces ADAM's `compute_kep_flux` bitwise
!< (the same arithmetic, in the same order). `seth_keep_pe8_opt` uses the dispersion-optimised sixth-order a_l of Liang
!< & Fu (2024b), the non-dissipative branch of their TENO-S.
!<
!<    call seth_keep_pe6(d, gamma, w, f)
!<
!< d: the normal direction (1..3); gamma: the ratio of specific heats; w(5,lo:hi): rho, u, v, w, p of the cells around
!< the face, cell 0 on its left (bounds 0:1, -1:2, -2:3, -3:4 for orders 2, 4, 6, 8); f(5): the flux of rho, rho u,
!< rho v, rho w, E. Pure, `acc routine seq`, `omp declare target`.
use penf, only : I4P, R_P
implicit none
private
public :: seth_keep_pe2, seth_keep_pe4, seth_keep_pe6, seth_keep_pe8, seth_keep_pe8_opt

contains
   pure subroutine seth_keep_pe2(d, gamma, w, f)
   !< KEEP-PE flux of order 2 (the two-point flux itself).
   integer(I4P), parameter :: NL = 1_I4P
   real(R_P),    parameter :: A(NL) = [0.5_R_P]
#include "seth_kep_body.INC"
   endsubroutine seth_keep_pe2

   pure subroutine seth_keep_pe4(d, gamma, w, f)
   !< KEEP-PE flux of order 4.
   integer(I4P), parameter :: NL = 2_I4P
   real(R_P),    parameter :: A(NL) = [2._R_P/3._R_P, -1._R_P/12._R_P]
#include "seth_kep_body.INC"
   endsubroutine seth_keep_pe4

   pure subroutine seth_keep_pe6(d, gamma, w, f)
   !< KEEP-PE flux of order 6 (ADAM's compute_kep_flux).
   integer(I4P), parameter :: NL = 3_I4P
   real(R_P),    parameter :: A(NL) = [0.75_R_P, -0.15_R_P, 1._R_P/60._R_P]
#include "seth_kep_body.INC"
   endsubroutine seth_keep_pe6

   pure subroutine seth_keep_pe8(d, gamma, w, f)
   !< KEEP-PE flux of order 8.
   integer(I4P), parameter :: NL = 4_I4P
   real(R_P),    parameter :: A(NL) = [0.8_R_P, -0.2_R_P, 4._R_P/105._R_P, -1._R_P/280._R_P]
#include "seth_kep_body.INC"
   endsubroutine seth_keep_pe8

   pure subroutine seth_keep_pe8_opt(d, gamma, w, f)
   !< KEEP-PE flux on eight points with the dispersion-optimised sixth-order a_l of Liang & Fu (2024b, section 3.2).
   integer(I4P), parameter :: NL = 4_I4P
   real(R_P),    parameter :: A(NL) = [0.8441336457279522_R_P, -0.2441336457279522_R_P, 0.0570096576929319_R_P, &
                                      -0.006723831837710873_R_P]
#include "seth_kep_body.INC"
   endsubroutine seth_keep_pe8_opt
endmodule seth_kep
