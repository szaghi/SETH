!< Cookbook: the central branch of a hybrid for the Euler equations, the recipe chosen on the command line.
program cookbook_euler
!run cookbook-keep cookbook_euler keep
!run cookbook-ducros cookbook_euler ducros
use seth
use selah, only : stepper, box, box_interior, SCHEME_SSP_33, FORM_SHU_OSHER_2R   ! the time integrator: SELAH's
implicit none
real(R_P), parameter :: PI = acos(-1._R_P), GAMMA = 1.4_R_P
character(16)        :: recipe

call get_command_argument(1, recipe)
select case(trim(recipe))
case('keep')   ; call pressure_equilibrium
case('ducros') ; call ducros
case default   ; print '(A)', 'usage: cookbook_euler keep|ducros'
endselect

contains
   subroutine pressure_equilibrium
   !< A density wave in uniform velocity and pressure, advected one period by the sixth-order KEEP-PE flux.
   integer, parameter :: N = 64, NG = 3
   real(R_P)          :: q(5,N), h, t, dt, w(5), du, dp
   integer            :: i, k
   type(stepper)      :: st
   type(box)          :: bx
   real(R_P), pointer :: r(:,:)
   h = 1._R_P/N
   do i=1, N                                    ! rho = 1 + 0.5 sin(2 pi x), u = 1, p = 1
      q(:,i) = cons([1 + 0.5_R_P*sin(2*PI*(i - 0.5_R_P)*h), 1._R_P, 0._R_P, 0._R_P, 1._R_P])
   enddo
   bx = box_interior(lb=[1, 1], ub=[5, N], ngc=0, components=.true.)
   call st%initialize(SCHEME_SSP_33, FORM_SHU_OSHER_2R, lb=[1, 1], ub=[5, N], components=.true.)
   t = 0
   do while (t < 1._R_P)
      dt = min(0.4_R_P*h/(1 + sqrt(GAMMA*2)), 1._R_P - t)        ! |u| + c, c largest where rho is smallest
      call st%begin_step(bx, q)
      do k=1, st%stages()                          ! SELAH's SSP(3,3): the stage state is q itself
         call st%begin_stage(k, bx, dt, q)
         call st%rhs_register_view(r)
         call residual(q, r)
         call st%end_stage(k, bx, dt, q)
      enddo
      call st%finish_step(bx, dt, q)
      t = t + dt
   enddo
   call st%destroy
   du = 0 ; dp = 0
   do i=1, N
      w = prim(q(:,i)) ; du = max(du, abs(w(2) - 1)) ; dp = max(dp, abs(w(5) - 1))
   enddo
   print '(A)', 'a density wave in uniform u = 1, p = 1, one period with seth_keep_pe6, 64 cells'
   print '(A,L2)', '  velocity and pressure still uniform to 1e-12:', du < 1.e-12_R_P .and. dp < 1.e-12_R_P
   endsubroutine pressure_equilibrium

   subroutine residual(q, res)
   !< -(F_{i+1/2} - F_{i-1/2})/h with the KEEP-PE flux of the primitive state of cells i-2..i+3.
   real(R_P), intent(in)  :: q(:,:)
   real(R_P), intent(out) :: res(:,:)
   real(R_P)              :: wg(5,-2:size(q,2)+3), f(5,0:size(q,2)), h
   integer                :: i, n
   n = size(q, 2) ; h = 1._R_P/n
   do i=-2, n+3                                 ! periodic: the primitive state of every cell and ghost cell
      wg(:,i) = prim(q(:,modulo(i - 1, n) + 1))
   enddo
   !region keep
   do i=0, n                                    ! d = 1: the face's normal is x
      call seth_keep_pe6(1, GAMMA, wg(:,i-2:i+3), f(:,i))
   enddo
   !endregion keep
   res = -(f(:,1:n) - f(:,0:n-1))/h
   endsubroutine residual

   subroutine ducros
   !< The Ducros sensors on four velocity gradients.
   real(R_P) :: g(3,3), phi, psi
   integer   :: k
   character(22) :: label(4) = [character(22) :: 'compression', 'expansion', 'rotation', 'shear']
   print '(A)', 'velocity gradient            phi (Ducros)   psi (compression only)'
   do k=1, 4
      g = 0
      select case(k)
      case(1) ; g(1,1) = -1 ; g(2,2) = -1 ; g(3,3) = -1          ! div u = -3, curl u = 0
      case(2) ; g(1,1) =  1 ; g(2,2) =  1 ; g(3,3) =  1          ! div u = +3
      case(3) ; g(1,2) = -1 ; g(2,1) =  1                         ! a vortex: div u = 0, curl u = 2 e_z
      case(4) ; g(1,2) =  1                                       ! shear: div u = 0, curl u = -e_z
      endselect
      !region ducros
      call seth_ducros(1.e-12_R_P, g, phi)
      call seth_ducros_compression(1.e-12_R_P, g, psi)
      !endregion ducros
      print '(A22,F14.4,F18.4)', label(k), phi, psi
   enddo
   endsubroutine ducros

   pure function cons(w) result(q)
   real(R_P), intent(in) :: w(5)                ! rho, u, v, w, p
   real(R_P)             :: q(5)                ! rho, rho u, rho v, rho w, E
   q = [w(1), w(1)*w(2:4), w(5)/(GAMMA - 1) + 0.5_R_P*w(1)*sum(w(2:4)**2)]
   endfunction cons

   pure function prim(q) result(w)
   real(R_P), intent(in) :: q(5)
   real(R_P)             :: w(5)
   w(1) = q(1) ; w(2:4) = q(2:4)/q(1) ; w(5) = (GAMMA - 1)*(q(5) - 0.5_R_P*sum(q(2:4)**2)/q(1))
   endfunction prim
endprogram cookbook_euler
