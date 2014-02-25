module randomsparse_test
use system
use fruit
use randomsparse, only: RandomSparseMatrix
implicit none

contains

subroutine test_randomsparse_resetval()
   implicit none
   type(RandomSparseMatrix) :: A
   complex(dp) :: x(2), y(2), Z(2)

   call A%init(2, 2)
   call A%add(1, 1, 1, (0.5_dp, 0._dp))
   call A%add(1, 2, 1, (1._dp, 0._dp))
   call A%add(2, 1, 2, (1.5_dp, 0._dp))
   call A%finalize()
   x = [1._dp, 2._dp]
   y = 0._dp

   Z = [1._dp, 4._dp]
   call A%multiply(Z, x, y)
   call assert_equals((2.5_dp, 0._dp), y(1), "Randomsparse (1)")
   call assert_equals((6._dp, 0._dp), y(2), "Randomsparse (2)")

   Z = [5._dp, 0._dp]
   call A%multiply(Z, x, y, (1., 0._dp))
   call assert_equals((12.5_dp, 0._dp), y(1), "Randomsparse update (1)")
   call assert_equals((0._dp, 0._dp), y(2), "Randomsparse update (2)")
end subroutine test_randomsparse_resetval

end module randomsparse_test
