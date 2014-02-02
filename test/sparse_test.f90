module sparse_test
use system
use fruit
use sparse, only: SparseMatrix

implicit none

contains

subroutine test_sparse_matmul()
   implicit none
   type(SparseMatrix) :: A
   complex(dp) :: x(2), y(2)
   call A%init(2, 1)

   call A%add(1, 1, (1._dp, 0._dp))
   call A%add(1, 2, (2._dp, 0._dp))
   call A%add(2, 2, (4._dp, 1._dp))
   call A%add(2, 1, (3._dp, 0._dp))
   call A%add(2, 2, (4._dp, 1._dp))
   call A%finalize()
   x = [1._dp, 2._dp]
   call A%multiply(x, y, (1._dp, 0._dp))

   call assert_equals((5._dp, 0._dp), y(1))
   call assert_equals((19._dp, 4._dp), y(2))

   call A%free()
end subroutine test_sparse_matmul


subroutine test_sparse_setblock()
   implicit none
   type(SparseMatrix) :: A
   complex(dp) :: x(2), y(2)
   ! Recall we are using column major order!
   complex(dp), parameter :: h(2, 2) = reshape([1, 3, 2, 4], [2, 2])
   call A%init(2, 2)

   call A%add_block(1, 1, h)
   call A%finalize()
   x = [1._dp, 2._dp]
   call A%multiply(x, y, (1._dp, 0._dp))

   call assert_equals((5._dp, 0._dp), y(1))
   call assert_equals((11._dp, 0._dp), y(2))

   call A%free()
end subroutine test_sparse_setblock

end module sparse_test
