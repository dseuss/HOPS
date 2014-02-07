module hierarchy_interface
use hierarchy
use iso_c_binding, only: c_double, c_double_complex, c_int, c_bool
implicit none

contains

subroutine c_init(tLength, tSteps, depth, modes, hs_dim, g, gamma, Omega, h, &
         Lmap, with_terminator, populated_modes) bind(c)
   real(c_double), intent(in)            :: tLength
   integer(c_int), intent(in)            :: tSteps
   integer(c_int), intent(in)            :: depth
   integer(c_int), intent(in)            :: modes
   integer(c_int), intent(in)            :: hs_dim
   complex(c_double_complex), intent(in) :: g(modes)
   real(c_double), intent(in)            :: gamma(modes)
   real(c_double), intent(in)            :: Omega(modes)
   complex(c_double_complex), intent(in) :: h(hs_dim, hs_dim)
   integer(c_int), intent(in)            :: Lmap(modes)
   logical(c_bool), intent(in)           :: with_terminator
   integer(c_int), intent(in)            :: populated_modes
   call init(tLength, tSteps, depth, g, gamma, omega, h, Lmap, &
         logical(with_terminator), populated_modes)
end subroutine c_init


subroutine c_free() bind(c)
   implicit none
   call free()
end subroutine c_free


subroutine c_run_trajectory_z0_rk4(hs_dim, tSteps, psi0, psi) bind(c)
   implicit none
   integer(c_int), intent(in) :: hs_dim
   integer(c_int), intent(in) :: tSteps
   complex(c_double_complex), intent(in) :: psi0(hs_dim)
   complex(c_double_complex), intent(out) :: psi(tSteps, hs_dim)
   psi = run_trajectory_z0_rk4(psi0)
end subroutine c_run_trajectory_z0_rk4


subroutine c_run_trajectory_z0_zvode(hs_dim, tSteps, psi0, psi) bind(c)
   implicit none
   integer(c_int), intent(in) :: hs_dim
   integer(c_int), intent(in) :: tSteps
   complex(c_double_complex), intent(in) :: psi0(hs_dim)
   complex(c_double_complex), intent(out) :: psi(tSteps, hs_dim)
   psi = run_trajectory_z0_zvode(psi0)
end subroutine c_run_trajectory_z0_zvode


end module hierarchy_interface
