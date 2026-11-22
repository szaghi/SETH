!< SETH tests: ADAM's former WENO kernels, verbatim, as the reference of the extraction (frozen history).
module seth_adam_reference
!< SETH tests: ADAM's former WENO kernels, verbatim, as the reference of the extraction (frozen history).
!<
!< The non-TBP procedures of ADAM's adam_weno_object.F90 (commit 302f093e), copied by a script without any edit;
!< only the module around them is new (S_MAX_M1 as ADAM defines it). weno_reconstruct_centered (copied at ADAM commit
!< f46b56fb, where the file equals 302f093e's) and weno_reconstruct_centered_d1fix, the same with the weights computed on the
!< shifted cells (audit D1: the one line changed), the references of seth_weno4/6/8_c.
!< compute_kep_flux and compute_ducros: ADAM's src/app/flume/common/adam_flume_hybrid_library.F90 (last changed in
!< commit 5e4fe007), the reference of seth_keep_pe6 and seth_ducros, copied the same way.
!<
!< Frozen history: ADAM no longer holds this code. Since ADAM issue 139 ADAM's applications call SETH: FLUME (CPU and
!< OpenACC), PRISM's CPU backend and RIPPLE take SETH's WENO-JS kernels (PRISM also seth_weno_descaling), FLUME's
!< hybrid seth_keep_pe6 and seth_ducros (ADAM commits dede5114, 73ac233a, 0e085764, b4398d8c, b1683349), every ADAM
!< regression case bitwise the same on CPU and OpenACC; ADAM commit ba37c4cb then deleted ADAM's WENO kernels and
!< tables and adam_flume_hybrid_library.F90. The centred WENO was not migrated: ADAM refuses the weno-c-* names since
!< commit b1683349 and deleted weno_reconstruct_centered in ba37c4cb, so seth_weno4/6/8_c have no counterpart in
!< ADAM. The copies stay as they were copied, the fixed reference of SETH's equivalence tests.
use penf, only : I4P, R8P
implicit none
private
public :: weno_reconstruct_upwind, weno_reconstruct_upwind_wratio, weno_descaling, weno_weights_exponent
public :: compute_kep_flux, compute_ducros
public :: weno_reconstruct_centered, weno_reconstruct_centered_d1fix
integer(I4P), parameter :: S_MAX_M1=4 !< Maximum number/dimensions of stencils minus 1 (ADAM).

contains
   pure subroutine weno_compute_convolution(S, vp, w, vr)
   !< Compute WENO convulution, non TBP.
   integer(I4P), intent(in)  :: S             !< Number of stencils used.
   real(R8P),    intent(in)  :: vp(1:2,0:S-1) !< Polynomial reconstructions.
   real(R8P),    intent(in)  :: w (1:2,0:S-1) !< Weights of the stencils.
   real(R8P),    intent(out) :: vr(1:2      ) !< Left and right (1,2) interface value of reconstructed v.
   integer(I4P)              :: k,f           !< Counter.

#ifdef DEV_OMP
   !$omp declare target(weno_compute_convolution)
#endif
   vr = 0._R8P
   do k=0, S-1
      do f=1, 2 ! 1 => left interface, 2 => right interface
         vr(f) = vr(f) + w(f,k)*vp(f,k)
      enddo
   enddo
   endsubroutine weno_compute_convolution

   pure subroutine weno_compute_polynomials_upwind(S, weno_p, v, vp)
   !< Compute WENO polynomials upwind, non TBP.
   integer(I4P), intent(in)  :: S                   !< Number of stencils used.
   real(R8P),    intent(in)  :: weno_p(1:,0:,0:,1:) !< Polinomials coefficients.
   real(R8P),    intent(in)  :: v (1:2,1-S:-1+S)    !< Variable to be reconstructed.
   real(R8P),    intent(out) :: vp(1:2,0:S-1   )    !< Polynomial reconstructions.
   integer(I4P)              :: s1,s2,f             !< Counter.

#ifdef DEV_OMP
   !$omp declare target(weno_compute_polynomials_upwind)
#endif
   vp = 0._R8P
   do s1=0, S-1 ! stencil counter
      do s2=0, S-1 ! cell counter counter
         do f=1, 2 ! 1 => left interface, 2 => right interface
            vp(f,s1) = vp(f,s1) + weno_p(f,s2,s1,S)*v(f,-s2+s1)
         enddo
      enddo
   enddo
   endsubroutine weno_compute_polynomials_upwind

   pure subroutine weno_compute_weights_upwind(S, weno_a, weno_d, weno_zeps, weno_rmu, v, w)
   !< Compute WENO weights, non TBP: `a = d / (zeps + (IS rmu) rmu)**wexp`, `rmu` the inverse of the descaler of the
   !< smoothness indicators (1 for the Jiang-Shu weights; see [[weno_descaling]]).
   integer(I4P), intent(in)  :: S                   !< Number of stencils used.
   real(R8P),    intent(in)  :: weno_a(1:,0:,1:)    !< Optimal weights.
   real(R8P),    intent(in)  :: weno_d(0:,0:,0:,1:) !< Smoothness indicators coefficients.
   real(R8P),    intent(in)  :: weno_zeps           !< Parameter for avoiding division by zero in computing IS.
   real(R8P),    intent(in)  :: weno_rmu            !< Inverse of the descaler (1: Jiang-Shu weights).
   real(R8P),    intent(in)  :: v    (1:2,1-S:-1+S) !< Variable to be reconstructed.
   real(R8P),    intent(out) :: w    (1:2,0:S-1)    !< Weights of the stencils.
   ! locals of constant bounds (S_MAX): bounds that depend on S make automatic arrays, which gfortran places on the heap
   ! without -fstack-arrays, a malloc and a free per call (10% of a CPU run, FLUME profile-couette3d)
   real(R8P)                 :: IS   (1:2,0:S_MAX_M1) !< Smoothness indicators of the stencils (0:S-1 used).
   real(R8P)                 :: a    (1:2,0:S_MAX_M1) !< Alpha coifficients for the weights (0:S-1 used).
   real(R8P)                 :: a_tot(1:2)          !< Summ of the alpha coefficients.
   real(R8P)                 :: d                   !< Denominator base, zeps + IS rmu^2.
   real(R8P)                 :: dp                  !< The base to the power wexp.
   integer(I4P)              :: wexp                !< Exponent of the smoothness indicators.
   integer(I4P)              :: s1,s2,s3,f          !< Counter.

#ifdef DEV_OMP
   !$omp declare target(weno_compute_weights_upwind)
#endif
   wexp = weno_weights_exponent(S=S)
   ! computing smoothness indicators
   do s1=0,S-1 ! stencil counter
      do f=1,2 ! 1 => left interface, 2 => right interface
         IS(f,s1) = 0._R8P
         do s2=0,S-1
            do s3=0,S-1
               IS(f,s1) = IS(f,s1) + weno_d(s3,s2,s1,S)*v(f,s1-s3)*v(f,s1-s2)
            enddo
         enddo
      enddo
   enddo
   ! computing alfa coefficients
   a_tot = 0._R8P
   do s1=0,S-1
      do f=1,2 ! 1 => left interface, 2 => right interface
         ! the power written out for the exponents weno_weights_exponent returns (1..4): the integer power of a variable
         ! exponent is a library call (__powidf2, 3% of a CPU run); the products are the ones its binary exponentiation
         ! forms, in the same order, so the result is bitwise the same
         d = weno_zeps + IS(f,s1)*weno_rmu*weno_rmu
         select case(wexp)
         case(1)
            dp = d
         case(2)
            dp = d*d
         case(3)
            dp = d*(d*d)
         case(4)
            dp = (d*d)*(d*d)
         case default
            dp = d**wexp
         endselect
         a(f,s1) = weno_a(f,s1,S)*(1._R8P/dp) ; a_tot(f) = a_tot(f) + a(f,s1)
      enddo
   enddo
   ! computing weights
   do s1=0,S-1
      do f=1,2
         w(f,s1) = a(f,s1)/a_tot(f)
      enddo
   enddo
   endsubroutine weno_compute_weights_upwind

   pure subroutine weno_reconstruct_upwind(S, weno_a, weno_p, weno_d, weno_zeps, weno_rmu, v, vr)
   !< Reconstruct by WENO upwind method of 2S-1 order, non TBP.
   integer(I4P), intent(in)  :: S                   !< Number of stencils used.
   real(R8P),    intent(in)  :: weno_a(1:,0:,1:)    !< Optimal weights.
   real(R8P),    intent(in)  :: weno_p(1:,0:,0:,1:) !< Polinomials coefficients.
   real(R8P),    intent(in)  :: weno_d(0:,0:,0:,1:) !< Smoothness indicators coefficients.
   real(R8P),    intent(in)  :: weno_zeps           !< Parameter for avoiding division by zero in computing IS.
   real(R8P),    intent(in)  :: weno_rmu            !< Inverse of the descaler (1: Jiang-Shu weights).
   real(R8P),    intent(in)  :: v (1:2,1-S:-1+S)    !< Variables to be reconstructed.
   real(R8P),    intent(out) :: vr(1:2         )    !< Left and right (1,2) interface value of reconstructed v.
   real(R8P)                 :: vp(1:2,0:S_MAX_M1)  !< Polynomial reconstructions (0:S-1 used).
   real(R8P)                 :: w (1:2,0:S_MAX_M1)  !< Weights of the stencils (0:S-1 used).

#ifdef DEV_OMP
   !$omp declare target(weno_reconstruct_upwind)
#endif
   call weno_compute_polynomials_upwind(S=S, weno_p=weno_p, v=v(1:2,1-S:-1+S), vp=vp(1:2,0:S-1))
   call weno_compute_weights_upwind(S=S, weno_a=weno_a, weno_d=weno_d, weno_zeps=weno_zeps, weno_rmu=weno_rmu,     &
                                    v=v(1:2,1-S:-1+S), w=w(1:2,0:S-1))
   call weno_compute_convolution(S=S, vp=vp(1:2,0:S-1), w=w(1:2,0:S-1), vr=vr(1:2))
   endsubroutine weno_reconstruct_upwind

   pure subroutine weno_reconstruct_upwind_wratio(S, weno_a, weno_p, weno_d, weno_zeps, weno_rmu, v, vr, wr)
   !< Reconstruct (or interpolate) by WENO upwind method of 2S-1 order, returning also the smoothness ratio of each
   !< interface, `wr(f) = min_k w(f,k) / weno_a(f,k,S)`: 1 where the nonlinear weights equal the linear ones (smooth
   !< data), toward 0 where a stencil is discarded. The same values `vr` as `weno_reconstruct_upwind`; non TBP.
   integer(I4P), intent(in)  :: S                   !< Number of stencils used.
   real(R8P),    intent(in)  :: weno_a(1:,0:,1:)    !< Optimal weights.
   real(R8P),    intent(in)  :: weno_p(1:,0:,0:,1:) !< Polinomials coefficients.
   real(R8P),    intent(in)  :: weno_d(0:,0:,0:,1:) !< Smoothness indicators coefficients.
   real(R8P),    intent(in)  :: weno_zeps           !< Parameter for avoiding division by zero in computing IS.
   real(R8P),    intent(in)  :: weno_rmu            !< Inverse of the descaler (1: Jiang-Shu weights).
   real(R8P),    intent(in)  :: v (1:2,1-S:-1+S)    !< Variables to be reconstructed.
   real(R8P),    intent(out) :: vr(1:2         )    !< Left and right (1,2) interface value of reconstructed v.
   real(R8P),    intent(out) :: wr(1:2         )    !< Left and right (1,2) interface smoothness ratio.
   real(R8P)                 :: vp(1:2,0:S_MAX_M1)  !< Polynomial reconstructions (0:S-1 used).
   real(R8P)                 :: w (1:2,0:S_MAX_M1)  !< Weights of the stencils (0:S-1 used).
   integer(I4P)              :: k, f                !< Counters.

#ifdef DEV_OMP
   !$omp declare target(weno_reconstruct_upwind_wratio)
#endif
   call weno_compute_polynomials_upwind(S=S, weno_p=weno_p, v=v(1:2,1-S:-1+S), vp=vp(1:2,0:S-1))
   call weno_compute_weights_upwind(S=S, weno_a=weno_a, weno_d=weno_d, weno_zeps=weno_zeps, weno_rmu=weno_rmu,     &
                                    v=v(1:2,1-S:-1+S), w=w(1:2,0:S-1))
   call weno_compute_convolution(S=S, vp=vp(1:2,0:S-1), w=w(1:2,0:S-1), vr=vr(1:2))
   do f=1, 2
      wr(f) = w(f,0) / weno_a(f,0,S)
      do k=1, S-1
         wr(f) = min(wr(f), w(f,k) / weno_a(f,k,S))
      enddo
   enddo
   endsubroutine weno_reconstruct_upwind_wratio

   pure function weno_descaling(S, weno_sigma, v) result(rmu)
   !< Return the inverse `1/mu` of the descaler of the smoothness indicators of one stencil pair, from its own values.
   !<
   !< Scale-invariant weights (Don, Li, Wang & Wang 2022, JCP 448, 110724; issue #49) read `(IS rmu) rmu` in place of
   !< `IS`, so that `zeps` acts relative to the magnitude of the data. Here `mu = mean |v|` over the 2(2S-1) values of
   !< both interfaces, floored at `tiny`: the generic descaler, for callers without a better magnitude. It is exact
   !< under `v -> 2**n v` (WENO(2**n v) = 2**n WENO(v) bit for bit while `mu` and `IS` stay normal), but it measures the
   !< values, so a field near zero by construction (a characteristic field projected with the local state) gets an
   !< effectively vanishing `zeps` and loses accuracy on smooth data: such callers pass the inverse of a magnitude of
   !< their own (FLUME: the absolute projection of the state, issue #49). `weno_sigma = 0` gives `rmu = 1` exactly, the
   !< Jiang-Shu weights bitwise; arithmetic, not a branch.
   integer(I4P), intent(in) :: S                !< Number of stencils used.
   real(R8P),    intent(in) :: weno_sigma       !< Descaler switch: 0 Jiang-Shu, 1 scale-invariant.
   real(R8P),    intent(in) :: v(1:2,1-S:-1+S)  !< Variable to be reconstructed.
   real(R8P)                :: rmu              !< Inverse of the descaler.
   real(R8P)                :: mu               !< Descaler.
   integer(I4P)             :: f, m             !< Counters.
   !$acc routine seq
   !$omp declare target

   mu = 0._R8P
   do m=1-S, S-1
      do f=1, 2
         mu = mu + abs(v(f,m))
      enddo
   enddo
   mu = max(mu/real(2*(2*S-1), R8P), tiny(1._R8P))
   rmu = weno_sigma*(1._R8P/mu) + (1._R8P - weno_sigma)
   endfunction weno_descaling

   pure function weno_weights_exponent(S) result(wexp)
   !< Return the exponent of the smoothness indicators in the WENO nonlinear weights, `a = d / (zeps + IS)**wexp`.
   !<
   !< The one definition used by the host weights, the device weights and `weno_object%wexp`: `wexp = S`, reduced to
   !< `S - 1` for `S > 4`. Before, the host used `S` and the device `2`, so the two backends ran different schemes for
   !< `S > 2`. Being a function of `S`, the exponent also follows the stencil of each ROR iteration.
   integer(I4P), intent(in) :: S    !< Number of stencils used.
   integer(I4P)             :: wexp !< Exponent of the smoothness indicators.
   !$acc routine seq
   !$omp declare target

   wexp = S
   if (S > 4) wexp = S - 1
   endfunction weno_weights_exponent

   pure subroutine compute_kep_flux(d, gamma, w, f)
   !< Compute the 6th-order KEP face flux of direction `d`.
   integer(I4P), intent(in)  :: d            !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: gamma        !< Specific heats ratio.
   real(R8P),    intent(in)  :: w(5,-2:3)    !< Cell states.
   real(R8P),    intent(out) :: f(5)         !< Face flux.
   real(R8P)                 :: a(3)         !< Coefficients.
   real(R8P)                 :: uu           !< u_l . u_r.
   real(R8P)                 :: mf, pf       !< Two-point mass flux and pressure.
   real(R8P)                 :: rg           !< 1 / (gamma - 1).
   integer(I4P)              :: l, m, n, c, r !< Counters.
   !$acc routine seq
   !$omp declare target

   a(1) = 0.75_R8P ; a(2) = -0.15_R8P ; a(3) = 1._R8P / 60._R8P
   rg = 1._R8P / (gamma - 1._R8P)
   do n=1, 5
      f(n) = 0._R8P
   enddo
   do l=1, 3
      do m=0, l-1
         c = -m ; r = -m + l
         mf = 0.25_R8P * (w(1,c) + w(1,r)) * (w(1+d,c) + w(1+d,r))
         pf = 0.5_R8P * (w(5,c) + w(5,r))
         f(1) = f(1) + 2._R8P * a(l) * mf
         do n=1, 3
            f(1+n) = f(1+n) + 2._R8P * a(l) * mf * 0.5_R8P * (w(1+n,c) + w(1+n,r))
         enddo
         f(1+d) = f(1+d) + 2._R8P * a(l) * pf
         uu = w(2,c) * w(2,r) + w(3,c) * w(3,r) + w(4,c) * w(4,r)
         f(5) = f(5) + 2._R8P * a(l) * (0.5_R8P * mf * uu + 0.5_R8P * rg * pf * (w(1+d,c) + w(1+d,r)) + &
                                        0.5_R8P * (w(1+d,c) * w(5,r) + w(1+d,r) * w(5,c)))
      enddo
   enddo
   endsubroutine compute_kep_flux

   pure subroutine compute_ducros(eps, jump, jump_max, g, psi)
   !< Compute the Ducros et al. (1999) sensor, `theta^2 / (theta^2 + |omega|^2 + eps)`, from the velocity gradient `g(a,b) = d u_a / d x_b`;
   !< 1 where the density jumps: the Ducros ratio cannot see a contact, where the velocity is uniform and a central flux rings.
   real(R8P), intent(in)  :: eps    !< Floor of the denominator.
   real(R8P), intent(in)  :: jump   !< Largest |rho_+ - rho_-| / (rho_+ + rho_-) of the cell along its active directions.
   real(R8P), intent(in)  :: jump_max !< Its bound.
   real(R8P), intent(in)  :: g(3,3) !< Velocity gradient.
   real(R8P), intent(out) :: psi    !< Sensor.
   real(R8P)              :: th, wx, wy, wz !< Dilatation and vorticity.
   !$acc routine seq
   !$omp declare target

   th = g(1,1) + g(2,2) + g(3,3)
   wx = g(3,2) - g(2,3) ; wy = g(1,3) - g(3,1) ; wz = g(2,1) - g(1,2)
   psi = th * th / (th * th + wx * wx + wy * wy + wz * wz + eps)
   if (jump > jump_max) psi = 1._R8P
   endsubroutine compute_ducros

   pure subroutine weno_reconstruct_centered(S, weno_a, weno_p, weno_d, weno_c, weno_zeps, weno_rmu, v, vr)
   !< Reconstruct by WENO centered method of 2S order, non TBP.
   integer(I4P), intent(in)  :: S                   !< Number of stencils used.
   real(R8P),    intent(in)  :: weno_a(1:,0:,1:)    !< Optimal weights.
   real(R8P),    intent(in)  :: weno_p(1:,0:,0:,1:) !< Polinomials coefficients.
   real(R8P),    intent(in)  :: weno_d(0:,0:,0:,1:) !< Smoothness indicators coefficients.
   real(R8P),    intent(in)  :: weno_c(1-S:,1:)     !< Centered polinomials coefficients.
   real(R8P),    intent(in)  :: weno_zeps           !< Parameter for avoiding division by zero in computing IS.
   real(R8P),    intent(in)  :: weno_rmu            !< Inverse of the descaler (1: Jiang-Shu weights).
   real(R8P),    intent(in)  :: v (1:2,1-S:S)       !< Variables to be reconstructed.
   real(R8P),    intent(out) :: vr(1:2)             !< Left and right (1,2) interface value of reconstructed v.
   real(R8P)                 :: vu(1:2,-S_MAX_M1:S_MAX_M1) !< Variables to be reconstructed, upwind stencils (1-S:S-1).
   real(R8P)                 :: vp(1:2,0:S_MAX_M1)  !< Polynomial reconstructions (0:S-1 used).
   real(R8P)                 :: w (1:2,0:S_MAX_M1)  !< Weights of the stencils (0:S-1 used).
   real(R8P)                 :: vc(1:2)             !< Centered reconstruction.
   real(R8P)                 :: IS(1:3),ISC(1:2)    !< Centered smoothness indicator.
   real(R8P)                 :: phi                 !< Centered-upwind blending coefficient.
   integer(I4P)              :: s1, f               !< Counter.

#ifdef DEV_OMP
   !$omp declare target(weno_reconstruct_centered)
#endif
   ! upwind reconstruction
   do s1=1-S, S
      if (s1>1-S) vu(1,s1-1) = v(1,s1)
      if (s1<S  ) vu(2,s1  ) = v(2,s1)
   enddo
   call weno_compute_polynomials_upwind(S=S, weno_p=weno_p, v=vu(1:2,1-S:-1+S), vp=vp(1:2,0:S-1))
   call weno_compute_weights_upwind(S=S, weno_a=weno_a, weno_d=weno_d, weno_zeps=weno_zeps, weno_rmu=weno_rmu,     &
                                    v=v(1:2,1-S:-1+S), w=w(1:2,0:S-1))
   call weno_compute_convolution(S=S, vp=vp(1:2,0:S-1), w=w(1:2,0:S-1), vr=vr(1:2))
   if (S>1) then
      ! centered reconstruction
      do f=1, 2
         vc(f) = 0._R8P
         do s1=1-S, S
            vc(f) = vc(f) + weno_c(s1,S) * v(f,s1)
         enddo
         ! blending centered reconstruction with upwind one for discontinuos solution
         IS(1) = (v(f,0) - v(f,-1))**2 + weno_zeps
         IS(2) = (v(f,1) - v(f, 0))**2 + weno_zeps
         IS(3) = (v(f,2) - v(f, 1))**2 + weno_zeps
         ISC(1) = abs(IS(2)-IS(1))/(IS(2)+IS(1))
         ISC(2) = abs(IS(3)-IS(2))/(IS(3)+IS(2))
         phi = max(ISC(1),ISC(2))
         vr(f) = vc(f) + phi*(vr(f)-vc(f))
      enddo
   endif
   endsubroutine weno_reconstruct_centered

   pure subroutine weno_reconstruct_centered_d1fix(S, weno_a, weno_p, weno_d, weno_c, weno_zeps, weno_rmu, v, vr)
   !< Reconstruct by WENO centered method of 2S order, non TBP.
   !< SETH: ADAM's routine with D1 fixed, the weights on the shifted cells vu (the one changed line).
   integer(I4P), intent(in)  :: S                   !< Number of stencils used.
   real(R8P),    intent(in)  :: weno_a(1:,0:,1:)    !< Optimal weights.
   real(R8P),    intent(in)  :: weno_p(1:,0:,0:,1:) !< Polinomials coefficients.
   real(R8P),    intent(in)  :: weno_d(0:,0:,0:,1:) !< Smoothness indicators coefficients.
   real(R8P),    intent(in)  :: weno_c(1-S:,1:)     !< Centered polinomials coefficients.
   real(R8P),    intent(in)  :: weno_zeps           !< Parameter for avoiding division by zero in computing IS.
   real(R8P),    intent(in)  :: weno_rmu            !< Inverse of the descaler (1: Jiang-Shu weights).
   real(R8P),    intent(in)  :: v (1:2,1-S:S)       !< Variables to be reconstructed.
   real(R8P),    intent(out) :: vr(1:2)             !< Left and right (1,2) interface value of reconstructed v.
   real(R8P)                 :: vu(1:2,-S_MAX_M1:S_MAX_M1) !< Variables to be reconstructed, upwind stencils (1-S:S-1).
   real(R8P)                 :: vp(1:2,0:S_MAX_M1)  !< Polynomial reconstructions (0:S-1 used).
   real(R8P)                 :: w (1:2,0:S_MAX_M1)  !< Weights of the stencils (0:S-1 used).
   real(R8P)                 :: vc(1:2)             !< Centered reconstruction.
   real(R8P)                 :: IS(1:3),ISC(1:2)    !< Centered smoothness indicator.
   real(R8P)                 :: phi                 !< Centered-upwind blending coefficient.
   integer(I4P)              :: s1, f               !< Counter.

#ifdef DEV_OMP
   !$omp declare target(weno_reconstruct_centered_d1fix)
#endif
   ! upwind reconstruction
   do s1=1-S, S
      if (s1>1-S) vu(1,s1-1) = v(1,s1)
      if (s1<S  ) vu(2,s1  ) = v(2,s1)
   enddo
   call weno_compute_polynomials_upwind(S=S, weno_p=weno_p, v=vu(1:2,1-S:-1+S), vp=vp(1:2,0:S-1))
   call weno_compute_weights_upwind(S=S, weno_a=weno_a, weno_d=weno_d, weno_zeps=weno_zeps, weno_rmu=weno_rmu,     &
                                    v=vu(1:2,1-S:-1+S), w=w(1:2,0:S-1))
   call weno_compute_convolution(S=S, vp=vp(1:2,0:S-1), w=w(1:2,0:S-1), vr=vr(1:2))
   if (S>1) then
      ! centered reconstruction
      do f=1, 2
         vc(f) = 0._R8P
         do s1=1-S, S
            vc(f) = vc(f) + weno_c(s1,S) * v(f,s1)
         enddo
         ! blending centered reconstruction with upwind one for discontinuos solution
         IS(1) = (v(f,0) - v(f,-1))**2 + weno_zeps
         IS(2) = (v(f,1) - v(f, 0))**2 + weno_zeps
         IS(3) = (v(f,2) - v(f, 1))**2 + weno_zeps
         ISC(1) = abs(IS(2)-IS(1))/(IS(2)+IS(1))
         ISC(2) = abs(IS(3)-IS(2))/(IS(3)+IS(2))
         phi = max(ISC(1),ISC(2))
         vr(f) = vc(f) + phi*(vr(f)-vc(f))
      enddo
   endif
   endsubroutine weno_reconstruct_centered_d1fix
endmodule seth_adam_reference
