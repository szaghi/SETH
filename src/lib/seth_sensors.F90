!< SETH, smoothness sensors of one face stencil: where a hybrid may use its linear branch.
module seth_sensors
!< SETH, smoothness sensors of one face stencil: where a hybrid may use its linear branch.
!<
!< Each sensor reads the six values v(-2:3) around the face i+1/2 (cell 0 = i, SETH's convention), multiplied by rmu
!< (an inverse magnitude: 1, or the normalisation the caller chooses, e.g. Li, Chen & Ren's 1/max(|v|_max, 1)), and
!< returns s in [0, 1] with one meaning for all of them:
!<
!<    s = 1   the published binary test says smooth: the linear branch may be used;
!<    s < 1   not smooth; s itself is the published continuous form where there is one, a blending weight.
!<
!<    call seth_sensor_teno(v, rmu, thr, s)   Fu's TENO-based detector (Fu 2019, as written by Li, Chen & Ren 2022,
!<                                            Eqs. (55)-(57)): four three-point stencils i-2..i, i-1..i+1, i..i+2,
!<                                            i+1..i+3, Jiang-Shu indicators on cell i, gamma_k = (beta_k + 1e-4)^-6,
!<                                            normalised; smooth iff every chi_k > thr = 5e-4. s = min(1, min chi/thr).
!<    call seth_sensor_ren(v, rmu, thr, s)    Ren et al.'s ratio (Li, Chen & Ren, Eqs. (48)-(53)): psi_j =
!<                                            (2 a b + e)/(a^2 + b^2 + e), a = |v_j - v_{j-1}|, b = |v_j - v_{j+1}|,
!<                                            the minimum over j = i-1 .. i+2, e = 0.9 thr/(1 - 0.9 thr) xi^2,
!<                                            xi = 1e-3, thr = psi_c = 0.3. s = min(1, psi/thr) (Eq. (52)).
!<    call seth_sensor_lcr(v, rmu, thr, s)    Li, Chen & Ren's detector (Eq. (60)): Ren's ratio with second
!<                                            differences added to a and b, the minimum over j = i, i+1, xi = 1e-2,
!<                                            thr = psi_c = 0.4. s = min(1, psi/thr).
!<    call seth_sensor_tau(v, rmu, thr, s)    Wan & Xia's (2022, Eqs. (3.13)-(3.15)) eta = tau_5/(e + mean beta_k) on
!<                                            WENO5's three candidates, e = 1e-6: smooth for eta <= thr = 1/2,
!<                                            non-smooth for eta >= 2 thr, the linear blend s = 2 - eta/thr between.
!<
!< The regularisers (1e-4, xi^2, 1e-6) are absolute: the sensors are not scale invariant, which rmu lets the caller
!< correct. Pure, `acc routine seq`, `omp declare target`.
use penf, only : I4P, R_P
implicit none
private
public :: SETH_SENSOR_TENO_CT, SETH_SENSOR_REN_PSIC, SETH_SENSOR_LCR_PSIC, SETH_SENSOR_TAU_ETA
public :: seth_sensor_teno, seth_sensor_ren, seth_sensor_lcr, seth_sensor_tau

real(R_P), parameter :: SETH_SENSOR_TENO_CT  = 5.e-4_R_P !< C_T* of Fu's detector (Li, Chen & Ren 2022, Eq. 57).
real(R_P), parameter :: SETH_SENSOR_REN_PSIC = 0.3_R_P   !< psi_c of Ren's detector (Li, Chen & Ren 2022, Eq. 52).
real(R_P), parameter :: SETH_SENSOR_LCR_PSIC = 0.4_R_P   !< psi_c of Li, Chen & Ren's detector (Eq. 60).
real(R_P), parameter :: SETH_SENSOR_TAU_ETA  = 0.5_R_P   !< Smooth limit of Wan & Xia's eta (non-smooth at twice it).

contains
   pure subroutine seth_sensor_teno(v, rmu, thr, s)
   !< Fu's TENO-based detector: s = min(1, min_k chi_k / thr).
   real(R_P), intent(in)  :: v(-2:3) !< Values around the face.
   real(R_P), intent(in)  :: rmu     !< Inverse magnitude the values are multiplied by.
   real(R_P), intent(in)  :: thr     !< Threshold C_T* (SETH_SENSOR_TENO_CT).
   real(R_P), intent(out) :: s       !< Smoothness, 1: smooth.
   real(R_P), parameter   :: EPS = 1.e-4_R_P !< Regulariser of the indicators (Li, Chen & Ren, after Eq. 55).
   real(R_P)              :: w(-2:3), beta(4), g(4), gmin, gsum
   integer(I4P)           :: k
   !$acc routine seq
   !$omp declare target

   do k=-2, 3
      w(k) = rmu*v(k)
   enddo
   beta(1) = 13._R_P/12._R_P*(w(-2) - 2._R_P*w(-1) + w(0))**2 + 0.25_R_P*(w(-2) - 4._R_P*w(-1) + 3._R_P*w(0))**2
   beta(2) = 13._R_P/12._R_P*(w(-1) - 2._R_P*w(0) + w(1))**2 + 0.25_R_P*(w(-1) - w(1))**2
   beta(3) = 13._R_P/12._R_P*(w(0) - 2._R_P*w(1) + w(2))**2 + 0.25_R_P*(3._R_P*w(0) - 4._R_P*w(1) + w(2))**2
   beta(4) = 13._R_P/12._R_P*(w(1) - 2._R_P*w(2) + w(3))**2 + 0.25_R_P*(5._R_P*w(1) - 8._R_P*w(2) + 3._R_P*w(3))**2
   ! chi_k = gamma_k / sum gamma, gamma_k = (beta_k + eps)^-6 = ((beta_min + eps)/(beta_k + eps))^6 / (beta_min + eps)^6
   gmin = minval(beta) + EPS
   gsum = 0._R_P
   do k=1, 4
      g(k) = gmin/(beta(k) + EPS)
      g(k) = (g(k)*g(k))*(g(k)*g(k))*(g(k)*g(k))
      gsum = gsum + g(k)
   enddo
   s = min(1._R_P, (minval(g)/gsum)/thr)
   endsubroutine seth_sensor_teno

   pure subroutine seth_sensor_ren(v, rmu, thr, s)
   !< Ren et al.'s ratio of first differences, the minimum over the nodes i-1 .. i+2: s = min(1, psi/thr).
   real(R_P), intent(in)  :: v(-2:3) !< Values around the face.
   real(R_P), intent(in)  :: rmu     !< Inverse magnitude the values are multiplied by.
   real(R_P), intent(in)  :: thr     !< Threshold psi_c (SETH_SENSOR_REN_PSIC).
   real(R_P), intent(out) :: s       !< Smoothness, 1: smooth.
   real(R_P), parameter   :: XI = 1.e-3_R_P !< xi of Eq. (51).
   real(R_P)              :: a, b, e, psi
   integer(I4P)           :: j
   !$acc routine seq
   !$omp declare target

   e = 0.9_R_P*thr/(1._R_P - 0.9_R_P*thr)*XI*XI
   psi = 1._R_P
   do j=-1, 2
      a = abs(rmu*(v(j) - v(j-1)))
      b = abs(rmu*(v(j) - v(j+1)))
      psi = min(psi, (2._R_P*a*b + e)/(a*a + b*b + e))
   enddo
   s = min(1._R_P, psi/thr)
   endsubroutine seth_sensor_ren

   pure subroutine seth_sensor_lcr(v, rmu, thr, s)
   !< Li, Chen & Ren's detector: Ren's ratio with second differences, the minimum over the nodes i, i+1.
   real(R_P), intent(in)  :: v(-2:3) !< Values around the face.
   real(R_P), intent(in)  :: rmu     !< Inverse magnitude the values are multiplied by.
   real(R_P), intent(in)  :: thr     !< Threshold psi_c (SETH_SENSOR_LCR_PSIC).
   real(R_P), intent(out) :: s       !< Smoothness, 1: smooth.
   real(R_P), parameter   :: XI = 1.e-2_R_P !< xi of Eq. (60).
   real(R_P)              :: a, b, e, psi
   integer(I4P)           :: j
   !$acc routine seq
   !$omp declare target

   e = 0.9_R_P*thr/(1._R_P - 0.9_R_P*thr)*XI*XI
   psi = 1._R_P
   do j=0, 1
      a = abs(rmu*(v(j) - v(j-1))) + abs(rmu*(v(j) - 2._R_P*v(j-1) + v(j-2)))
      b = abs(rmu*(v(j) - v(j+1))) + abs(rmu*(v(j) - 2._R_P*v(j+1) + v(j+2)))
      psi = min(psi, (2._R_P*a*b + e)/(a*a + b*b + e))
   enddo
   s = min(1._R_P, psi/thr)
   endsubroutine seth_sensor_lcr

   pure subroutine seth_sensor_tau(v, rmu, thr, s)
   !< Wan & Xia's eta = tau_5/(e + mean beta_k) on WENO5's left-biased candidates: s = 1 for eta <= thr, 0 for
   !< eta >= 2 thr, 2 - eta/thr between.
   real(R_P), intent(in)  :: v(-2:3) !< Values around the face (the left-biased cells -2..2 are read).
   real(R_P), intent(in)  :: rmu     !< Inverse magnitude the values are multiplied by.
   real(R_P), intent(in)  :: thr     !< Smooth limit of eta (SETH_SENSOR_TAU_ETA).
   real(R_P), intent(out) :: s       !< Smoothness, 1: smooth.
   real(R_P), parameter   :: EPS = 1.e-6_R_P !< Regulariser (Wan & Xia, after Eq. (3.13)).
   real(R_P)              :: w(-2:2), b0, b1, b2, eta
   integer(I4P)           :: k
   !$acc routine seq
   !$omp declare target

   do k=-2, 2
      w(k) = rmu*v(k)
   enddo
   b0 = 13._R_P/12._R_P*(w(-2) - 2._R_P*w(-1) + w(0))**2 + 0.25_R_P*(w(-2) - 4._R_P*w(-1) + 3._R_P*w(0))**2
   b1 = 13._R_P/12._R_P*(w(-1) - 2._R_P*w(0) + w(1))**2 + 0.25_R_P*(w(-1) - w(1))**2
   b2 = 13._R_P/12._R_P*(w(0) - 2._R_P*w(1) + w(2))**2 + 0.25_R_P*(3._R_P*w(0) - 4._R_P*w(1) + w(2))**2
   eta = abs(b0 - b2)/(EPS + (b0 + b1 + b2)/3._R_P)
   s = max(0._R_P, min(1._R_P, 2._R_P - eta/thr))
   endsubroutine seth_sensor_tau
endmodule seth_sensors
