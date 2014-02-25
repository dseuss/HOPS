program main
use system
use hierarchy
implicit none

real(dp)    , parameter :: tLength = 400
integer     , parameter :: tSteps = 10000
integer     , parameter :: depth = 2
complex(dp) , parameter :: g(2) = [1._dp, 1._dp]
real(dp)    , parameter :: gamma(2) = [1._dp, 1._dp]
real(dp)    , parameter :: Omega(2) = [4._dp, 4._dp]
complex(dp) , parameter :: h(2, 2) = reshape([0._dp, .5_dp, .5_dp, 0._dp], &
                                             [2, 2])
integer     , parameter :: Lmap(2) = [1, 2]
logical     , parameter :: with_terminator = .true.
integer     , parameter :: populated_modes = 2
complex(dp) , parameter :: psi0(2) = [1._dp, 0._dp]
complex(dp) :: psi(tSteps, 2)

call init(tLength, tSteps, depth, g, gamma, Omega, h, Lmap, with_terminator, &
      populated_modes)
psi = run_trajectory_rk4(psi0)
call free()
end program main
