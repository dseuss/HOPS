module noisegen_interface
use system
use noisegen, only: ExponentialNoiseGenerator, init_random_seed
use iso_c_binding, only: c_double, c_double_complex, c_int
implicit none


public c_test
contains

subroutine c_test(dt, tSteps, modes, g, gamma, Omega, realizations, EZ, EZZ, &
         EZccZ) bind(c)
   implicit none
   real(c_double), intent(in)             :: dt
   integer(c_int), intent(in)             :: tSteps
   integer(c_int), intent(in)             :: modes
   complex(c_double_complex), intent(in)  :: g(modes)
   real(c_double), intent(in)             :: gamma(modes)
   real(c_double), intent(in)             :: Omega(modes)
   integer(c_int), intent(in)             :: realizations
   complex(c_double_complex), intent(out) :: EZ(tSteps)
   complex(c_double_complex), intent(out) :: EZZ(tSteps)
   complex(c_double_complex), intent(out) :: EZccZ(tSteps)

   type(ExponentialNoiseGenerator) :: gen
   complex(c_double_complex) :: Z(tSteps)
   integer(c_int) :: i

   call gen%init(dt, tSteps, g, gamma, Omega)

   EZ = (0._dp, 0._dp)
   EZZ = (0._dp, 0._dp)
   EZccZ = (0._dp, 0._dp)

   do i = 1, realizations
      Z = gen%get_realization()
      EZ = EZ + Z
      EZZ = EZZ + Z*Z(1)
      EZccZ = EZccZ + Z*conjg(Z(1))
   end do

   EZ = EZ / realizations
   EZZ = EZZ / realizations
   EZccZ = EZccZ / realizations

   call gen%free()
end subroutine c_test


subroutine c_init_random_seed() bind(c)
   implicit none
   call init_random_seed()
end subroutine c_init_random_seed

end module noisegen_interface
