module dynarray_test
use fruit
use dynarray_int
implicit none

contains

subroutine test_dynarray_reshape()
   implicit none
   type(IntDynamicArray) :: A

   call A%init(4)
   call A%add(1)
   call assert_equals(1, A%size())
   call A%add(2)
   call A%add(3)
   call A%add(4)
   call A%add(5)
   call assert_equals(5, A%size())

   call assert_equals(1, A%get(1))
   call assert_equals(2, A%get(2))
   call assert_equals(3, A%get(3))
   call assert_equals(4, A%get(4))
   call assert_equals(5, A%get(5))

   call A%free()
end subroutine test_dynarray_reshape

end module dynarray_test
