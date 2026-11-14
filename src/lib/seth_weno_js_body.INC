   ! Body of a WENO-JS kernel of SETH, included by every routine of seth_weno_js after it declares the named constants
   ! of its scheme: S (candidates, of width S), WA(1:2,0:S-1), WP(1:2,0:S-1,0:S-1), WD(0:S-1,0:S-1,0:S-1), WEXP.
   ! With SETH_WR defined, the routine also returns the weights ratio of each state.
   !
   ! Stencil convention: cells by their offset from the face i+1/2 (cell 0 left of it, cell 1 right of it); vp holds
   ! the values reconstructed with a left bias, vm those reconstructed with a right bias, both unshifted. Internally the
   ! two are packed as ADAM's rows (row 2 = vp, row 1 = vm shifted one cell) and every operation follows ADAM's host
   ! routines in the same order, so the result is ADAM's.
   real(R8P), intent(in)  :: vp(1-S:S) !< Values reconstructed with a left bias (cells 1-S..S-1 used).
   real(R8P), intent(in)  :: vm(1-S:S) !< Values reconstructed with a right bias (cells 2-S..S used).
   real(R8P), intent(in)  :: eps       !< Regulariser of the smoothness indicators.
   real(R8P), intent(in)  :: rmu       !< Inverse descaler of the smoothness indicators (1: Jiang-Shu weights).
   real(R8P), intent(out) :: rp        !< Left state at the face, from vp.
   real(R8P), intent(out) :: rm        !< Right state at the face, from vm.
#ifdef SETH_WR
   real(R8P), intent(out) :: wrp       !< Weights ratio of the left state, min_k w_k/d_k: 1 smooth, toward 0 rough.
   real(R8P), intent(out) :: wrm       !< Weights ratio of the right state.
   real(R8P)              :: wr(1:2)   !< Weights ratios, ADAM's rows.
#endif
   real(R8P)              :: v(1:2,1-S:S-1) !< Stencils, ADAM's rows.
   real(R8P)              :: q(1:2,0:S-1)   !< Candidate reconstructions.
   real(R8P)              :: is(1:2,0:S-1)  !< Smoothness indicators.
   real(R8P)              :: al(1:2,0:S-1)  !< Non-normalised weights.
   real(R8P)              :: al_tot(1:2)    !< Sum of the non-normalised weights.
   real(R8P)              :: w(1:2,0:S-1)   !< Weights.
   real(R8P)              :: r(1:2)         !< Reconstructed states, ADAM's rows.
   real(R8P)              :: b, bp          !< Denominator base and its power.
   integer(I4P)           :: m, s1, s2, s3, f
   !$acc routine seq
   !$omp declare target

   do m=1-S, S-1
      v(1,m) = vm(m+1)
      v(2,m) = vp(m)
   enddo
   ! candidate polynomials
   q = 0._R8P
   do s1=0, S-1
      do s2=0, S-1
         do f=1, 2
            q(f,s1) = q(f,s1) + WP(f,s2,s1)*v(f,-s2+s1)
         enddo
      enddo
   enddo
   ! smoothness indicators
   do s1=0, S-1
      do f=1, 2
         is(f,s1) = 0._R8P
         do s2=0, S-1
            do s3=0, S-1
               is(f,s1) = is(f,s1) + WD(s3,s2,s1)*v(f,s1-s3)*v(f,s1-s2)
            enddo
         enddo
      enddo
   enddo
   ! weights: the power written out as ADAM's host does (WEXP is a constant: the branch folds)
   al_tot = 0._R8P
   do s1=0, S-1
      do f=1, 2
         b = eps + is(f,s1)*rmu*rmu
         select case(WEXP)
         case(1)
            bp = b
         case(2)
            bp = b*b
         case(3)
            bp = b*(b*b)
         case(4)
            bp = (b*b)*(b*b)
         case default
            bp = b**WEXP
         endselect
         al(f,s1) = WA(f,s1)*(1._R8P/bp) ; al_tot(f) = al_tot(f) + al(f,s1)
      enddo
   enddo
   do s1=0, S-1
      do f=1, 2
         w(f,s1) = al(f,s1)/al_tot(f)
      enddo
   enddo
   ! convolution
   r = 0._R8P
   do s1=0, S-1
      do f=1, 2
         r(f) = r(f) + w(f,s1)*q(f,s1)
      enddo
   enddo
   rm = r(1)
   rp = r(2)
#ifdef SETH_WR
   do f=1, 2
      wr(f) = w(f,0) / WA(f,0)
      do s1=1, S-1
         wr(f) = min(wr(f), w(f,s1) / WA(f,s1))
      enddo
   enddo
   wrm = wr(1)
   wrp = wr(2)
#endif
