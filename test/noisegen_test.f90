module noisegen_test
use fruit
use system
use noisegen, only: ExponentialNoiseGenerator, init_random_seed
implicit none

real(dp), parameter :: tol = 0.05_dp

contains

subroutine test_noise_generator()
   implicit none
   integer, parameter :: &
         tSteps = 1000, &
         Realizations = 10000
   real(dp), parameter :: &
         dt = 0.025_dp, &
         gamma(1) = [1._dp], &
         Omega(1) = [4._dp]
   complex(dp), parameter :: g(1) = [1._dp]

   type(ExponentialNoiseGenerator) :: gen
   complex(dp) :: &
         Z(tSteps), &
         EZ(tSteps), &
         EZZ(tSteps), &
         EZccZ(tSteps)
   real(dp) :: diff(tSteps)
   integer :: N

   call init_random_seed()
   call gen%init(dt, tSteps, g, gamma, Omega)
   EZ = 0._dp
   EZZ = 0._dp
   EZccZ = 0._dp

   do N = 1, Realizations
      Z = gen%get_realization()
      EZ = EZ + Z
      EZZ = EZZ + Z * Z(1)
      EZccZ = EZccZ + Z * conjg(Z(1))
   end do

   EZ = EZ / Realizations
   EZZ = EZZ / Realizations
   EZccZ = EZccZ / Realizations

   do N = 1, tSteps
      diff(N) = abs(sum(g * exp(-(gamma + ii * Omega) * (N-1)*dt) - EZccZ(N)))
   end do

   call assert_true(all(abs(EZ) <= tol), "EZ")
   call assert_true(all(abs(EZZ) <= tol), "EZZ")
   call assert_true(all(diff <= tol), "EZccZ")
end subroutine test_noise_generator

end module noisegen_test
