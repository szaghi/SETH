!< SETH tests: the registry of every scheme and one dispatch, so that each test sweeps all of them.
module seth_test_schemes
!< SETH tests: the registry of every scheme and one dispatch, so that each test sweeps all of them.
!<
!< Every kernel reconstructs the two states at one face from two stencils with bounds 1-s:s (cells by offset from the
!< face). Here the stencils are always passed as `v(-4:5)` (the widest, WENO-9) and sliced to the scheme's bounds, so a
!< test writes its data once for all schemes. `reconstruct` calls the scheme's kernel with its default parameters, or
!< with the `eps`/`par` given (`par`: rmu for WENO-JS, p for WENO-Z, C_T for TENO); `wrp`/`wrm`, when present, come
!< from the `_wr` kernel (WENO-JS only). `interp` marks the interpolations (point values at the cell centres instead
!< of cell averages): WENO-JS (family FAM_JS_INTERP), WENO-Z and TENO (their own families, so that the family-driven
!< expectations apply to them too).
use seth
implicit none
private
public :: scheme_t, NSCHEMES, SCHEMES, FAM_JS, FAM_JS_INTERP, FAM_Z, FAM_TENO, VLO, VHI
public :: reconstruct

integer(I4P), parameter :: FAM_JS = 1, FAM_JS_INTERP = 2, FAM_Z = 3, FAM_TENO = 4
integer(I4P), parameter :: VLO = -4_I4P, VHI = 5_I4P !< Bounds of the test stencils.

type :: scheme_t
   character(16) :: name   = ''      !< Kernel name without the seth_ prefix.
   integer(I4P)  :: family = 0       !< FAM_*.
   integer(I4P)  :: s      = 0       !< Stencil bounds 1-s:s.
   real(R8P)     :: order  = 0._R8P  !< Order on smooth data (design order; the optimised TENO6 weights: 5 and 4).
   real(R8P)     :: eps    = 0._R8P  !< Default regulariser.
   real(R8P)     :: par    = 0._R8P  !< Default rmu (JS), p (Z), C_T (TENO).
   logical       :: interp = .false. !< Interpolation from point values (otherwise reconstruction from averages).
endtype scheme_t

integer(I4P), parameter :: NSCHEMES = 29_I4P
type(scheme_t), parameter :: SCHEMES(NSCHEMES) = [ &
   scheme_t('weno1_js',        FAM_JS,        1, 1._R8P, 1.e-6_R8P,  1._R8P), &
   scheme_t('weno3_js',        FAM_JS,        2, 3._R8P, 1.e-6_R8P,  1._R8P), &
   scheme_t('weno5_js',        FAM_JS,        3, 5._R8P, 1.e-6_R8P,  1._R8P), &
   scheme_t('weno7_js',        FAM_JS,        4, 7._R8P, 1.e-6_R8P,  1._R8P), &
   scheme_t('weno9_js',        FAM_JS,        5, 9._R8P, 1.e-6_R8P,  1._R8P), &
   scheme_t('weno1_js_interp', FAM_JS_INTERP, 1, 1._R8P, 1.e-6_R8P,  1._R8P, .true.), &
   scheme_t('weno3_js_interp', FAM_JS_INTERP, 2, 3._R8P, 1.e-6_R8P,  1._R8P, .true.), &
   scheme_t('weno5_js_interp', FAM_JS_INTERP, 3, 5._R8P, 1.e-6_R8P,  1._R8P, .true.), &
   scheme_t('weno7_js_interp', FAM_JS_INTERP, 4, 7._R8P, 1.e-6_R8P,  1._R8P, .true.), &
   scheme_t('weno9_js_interp', FAM_JS_INTERP, 5, 9._R8P, 1.e-6_R8P,  1._R8P, .true.), &
   scheme_t('wenoz3',          FAM_Z,         2, 3._R8P, 1.e-40_R8P, 1._R8P), &
   scheme_t('wenoz5',          FAM_Z,         3, 5._R8P, 1.e-40_R8P, 1._R8P), &
   scheme_t('wenoz7',          FAM_Z,         4, 7._R8P, 1.e-40_R8P, 1._R8P), &
   scheme_t('wenoz9',          FAM_Z,         5, 9._R8P, 1.e-40_R8P, 2._R8P), &
   scheme_t('teno5',           FAM_TENO,      3, 5._R8P, 1.e-40_R8P, SETH_TENO5_CT), &
   scheme_t('teno5_opt',       FAM_TENO,      3, 4._R8P, 1.e-40_R8P, SETH_TENO5_OPT_CT), &
   scheme_t('teno5_az',        FAM_TENO,      3, 5._R8P, 1.e-40_R8P, SETH_TENO5_AZ_CT), &
   scheme_t('teno6',           FAM_TENO,      3, 6._R8P, 1.e-40_R8P, SETH_TENO6_CT), &
   scheme_t('teno6_opt',       FAM_TENO,      3, 5._R8P, 1.e-40_R8P, SETH_TENO6_OPT_CT), &
   scheme_t('teno6_neno6',     FAM_TENO,      3, 4._R8P, 1.e-40_R8P, SETH_TENO6_NENO6_CT), &
   scheme_t('teno8',           FAM_TENO,      4, 6._R8P, 1.e-40_R8P, SETH_TENO8_CT), &
   scheme_t('teno5_a',         FAM_TENO,      3, 5._R8P, 1.e-40_R8P, SETH_TENO5_A_CR), &
   scheme_t('teno8_a',         FAM_TENO,      4, 6._R8P, 1.e-40_R8P, SETH_TENO8_A_CR), &
   scheme_t('wenoz3_interp',   FAM_Z,         2, 3._R8P, 1.e-40_R8P, 1._R8P, .true.), &
   scheme_t('wenoz5_interp',   FAM_Z,         3, 5._R8P, 1.e-40_R8P, 1._R8P, .true.), &
   scheme_t('wenoz7_interp',   FAM_Z,         4, 7._R8P, 1.e-40_R8P, 1._R8P, .true.), &
   scheme_t('wenoz9_interp',   FAM_Z,         5, 9._R8P, 1.e-40_R8P, 2._R8P, .true.), &
   scheme_t('teno5_interp',    FAM_TENO,      3, 5._R8P, 1.e-40_R8P, SETH_TENO5_CT, .true.), &
   scheme_t('teno6_interp',    FAM_TENO,      3, 6._R8P, 1.e-40_R8P, SETH_TENO6_CT, .true.)]

contains
   subroutine reconstruct(id, vp, vm, rp, rm, wrp, wrm, eps, par)
   !< The two states at the face from scheme `id`.
   integer(I4P), intent(in)            :: id               !< Scheme index in SCHEMES.
   real(R8P),    intent(in)            :: vp(VLO:VHI)      !< Left-biased values.
   real(R8P),    intent(in)            :: vm(VLO:VHI)      !< Right-biased values.
   real(R8P),    intent(out)           :: rp, rm           !< Left and right states.
   real(R8P),    intent(out), optional :: wrp, wrm         !< Weights ratios (WENO-JS: from the _wr kernel).
   real(R8P),    intent(in),  optional :: eps, par         !< Parameters (default: the registry's).
   real(R8P)                           :: e, q
   integer(I4P)                        :: p
   logical                             :: wr

   e = SCHEMES(id)%eps ; if (present(eps)) e = eps
   q = SCHEMES(id)%par ; if (present(par)) q = par
   p = 0 ; if (SCHEMES(id)%family == FAM_Z) p = nint(q) ! only WENO-Z reads an integer parameter
   wr = present(wrp) .and. present(wrm)
   select case(trim(SCHEMES(id)%name))
   case('weno1_js')
      if (wr) then ; call seth_weno1_js_wr(vp(0:1), vm(0:1), e, q, rp, rm, wrp, wrm)
      else ; call seth_weno1_js(vp(0:1), vm(0:1), e, q, rp, rm) ; endif
   case('weno3_js')
      if (wr) then ; call seth_weno3_js_wr(vp(-1:2), vm(-1:2), e, q, rp, rm, wrp, wrm)
      else ; call seth_weno3_js(vp(-1:2), vm(-1:2), e, q, rp, rm) ; endif
   case('weno5_js')
      if (wr) then ; call seth_weno5_js_wr(vp(-2:3), vm(-2:3), e, q, rp, rm, wrp, wrm)
      else ; call seth_weno5_js(vp(-2:3), vm(-2:3), e, q, rp, rm) ; endif
   case('weno7_js')
      if (wr) then ; call seth_weno7_js_wr(vp(-3:4), vm(-3:4), e, q, rp, rm, wrp, wrm)
      else ; call seth_weno7_js(vp(-3:4), vm(-3:4), e, q, rp, rm) ; endif
   case('weno9_js')
      if (wr) then ; call seth_weno9_js_wr(vp(-4:5), vm(-4:5), e, q, rp, rm, wrp, wrm)
      else ; call seth_weno9_js(vp(-4:5), vm(-4:5), e, q, rp, rm) ; endif
   case('weno1_js_interp')
      if (wr) then ; call seth_weno1_js_interp_wr(vp(0:1), vm(0:1), e, q, rp, rm, wrp, wrm)
      else ; call seth_weno1_js_interp(vp(0:1), vm(0:1), e, q, rp, rm) ; endif
   case('weno3_js_interp')
      if (wr) then ; call seth_weno3_js_interp_wr(vp(-1:2), vm(-1:2), e, q, rp, rm, wrp, wrm)
      else ; call seth_weno3_js_interp(vp(-1:2), vm(-1:2), e, q, rp, rm) ; endif
   case('weno5_js_interp')
      if (wr) then ; call seth_weno5_js_interp_wr(vp(-2:3), vm(-2:3), e, q, rp, rm, wrp, wrm)
      else ; call seth_weno5_js_interp(vp(-2:3), vm(-2:3), e, q, rp, rm) ; endif
   case('weno7_js_interp')
      if (wr) then ; call seth_weno7_js_interp_wr(vp(-3:4), vm(-3:4), e, q, rp, rm, wrp, wrm)
      else ; call seth_weno7_js_interp(vp(-3:4), vm(-3:4), e, q, rp, rm) ; endif
   case('weno9_js_interp')
      if (wr) then ; call seth_weno9_js_interp_wr(vp(-4:5), vm(-4:5), e, q, rp, rm, wrp, wrm)
      else ; call seth_weno9_js_interp(vp(-4:5), vm(-4:5), e, q, rp, rm) ; endif
   case('wenoz3')
      call seth_wenoz3(vp(-1:2), vm(-1:2), e, p, rp, rm)
   case('wenoz5')
      call seth_wenoz5(vp(-2:3), vm(-2:3), e, p, rp, rm)
   case('wenoz7')
      call seth_wenoz7(vp(-3:4), vm(-3:4), e, p, rp, rm)
   case('wenoz9')
      call seth_wenoz9(vp(-4:5), vm(-4:5), e, p, rp, rm)
   case('teno5')
      call seth_teno5(vp(-2:3), vm(-2:3), e, q, rp, rm)
   case('teno6')
      call seth_teno6(vp(-2:3), vm(-2:3), e, q, rp, rm)
   case('teno6_opt')
      call seth_teno6_opt(vp(-2:3), vm(-2:3), e, q, rp, rm)
   case('teno6_neno6')
      call seth_teno6_neno6(vp(-2:3), vm(-2:3), e, q, rp, rm)
   case('wenoz3_interp')
      call seth_wenoz3_interp(vp(-1:2), vm(-1:2), e, p, rp, rm)
   case('wenoz5_interp')
      call seth_wenoz5_interp(vp(-2:3), vm(-2:3), e, p, rp, rm)
   case('wenoz7_interp')
      call seth_wenoz7_interp(vp(-3:4), vm(-3:4), e, p, rp, rm)
   case('wenoz9_interp')
      call seth_wenoz9_interp(vp(-4:5), vm(-4:5), e, p, rp, rm)
   case('teno5_interp')
      call seth_teno5_interp(vp(-2:3), vm(-2:3), e, q, rp, rm)
   case('teno6_interp')
      call seth_teno6_interp(vp(-2:3), vm(-2:3), e, q, rp, rm)
   case('teno5_opt')
      call seth_teno5_opt(vp(-2:3), vm(-2:3), e, q, rp, rm)
   case('teno5_az')
      call seth_teno5_az(vp(-2:3), vm(-2:3), e, q, 1._R8P, rp, rm) ! Tang's form (rmu = 1): not scale invariant
   case('teno8')
      call seth_teno8(vp(-3:4), vm(-3:4), e, q, rp, rm)
   case('teno5_a') ! the adaptive cutoff, the papers' parameters (par: C_r)
      call seth_teno5_a(vp(-2:3), vm(-2:3), e, q, SETH_TENO5_A_XI, SETH_TENO5_A_A1, SETH_TENO5_A_A2, 1._R8P, rp, rm)
   case('teno8_a')
      call seth_teno8_a(vp(-3:4), vm(-3:4), e, q, SETH_TENO8_A_XI, SETH_TENO8_A_A1, SETH_TENO8_A_A2, 1._R8P, rp, rm)
   case default
      error stop 'reconstruct: unknown scheme'
   endselect
   if (.not.wr) then
      if (present(wrp)) wrp = -1._R8P
      if (present(wrm)) wrm = -1._R8P
   elseif (SCHEMES(id)%family /= FAM_JS .and. SCHEMES(id)%family /= FAM_JS_INTERP) then
      wrp = -1._R8P ; wrm = -1._R8P ! no weights ratio for this family
   endif
   endsubroutine reconstruct
endmodule seth_test_schemes
