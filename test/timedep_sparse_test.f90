module timedep_sparse_test
use system
use fruit
use timedep_sparse, only: TimeDepSparseMatrix
implicit none

contains

subroutine test_timedep_sparse_resetval()
   implicit none
   type(TimeDepSparseMatrix) :: A
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
   call assert_equals((2.5_dp, 0._dp), y(1), "TimeDepSparse (1)")
   call assert_equals((6._dp, 0._dp), y(2), "TimeDepSparse (2)")

   Z = [5._dp, 0._dp]
   call A%multiply(Z, x, y, (1., 0._dp))
   call assert_equals((12.5_dp, 0._dp), y(1), "TimeDepSparse update (1)")
   call assert_equals((0._dp, 0._dp), y(2), "TimeDepSparse update (2)")
end subroutine test_timedep_sparse_resetval

end module timedep_sparse_test
