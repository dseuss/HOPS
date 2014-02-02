module hstruct_test
use fruit
use hstruct, only: HStructure
use hstructtab, only: INVALID_INDEX

implicit none

contains

subroutine test_triangular_completeness()
   implicit none
   type(HStructure) :: h
   call h%init(2, 3)

   call assert_not_equals(INVALID_INDEX, h%find_index([0, 0]))

   call assert_not_equals(INVALID_INDEX, h%find_index([1, 0]))
   call assert_not_equals(INVALID_INDEX, h%find_index([0, 1]))

   call assert_not_equals(INVALID_INDEX, h%find_index([2, 0]))
   call assert_not_equals(INVALID_INDEX, h%find_index([1, 1]))
   call assert_not_equals(INVALID_INDEX, h%find_index([0, 2]))

   call assert_not_equals(INVALID_INDEX, h%find_index([3, 0]))
   call assert_not_equals(INVALID_INDEX, h%find_index([2, 1]))
   call assert_not_equals(INVALID_INDEX, h%find_index([1, 2]))
   call assert_not_equals(INVALID_INDEX, h%find_index([0, 3]))

   call assert_equals(INVALID_INDEX, h%find_index([1, 3]))
   call assert_equals(INVALID_INDEX, h%find_index([0, 8]))
   call assert_equals(INVALID_INDEX, h%find_index([4, 3]))

   call assert_equals(10, h%entries())

   call h%free()
end subroutine test_triangular_completeness


subroutine test_triangular_limited()
   implicit none
   type(HStructure) :: h
   call h%init(3, 3, 2)

   call assert_not_equals(INVALID_INDEX, h%find_index([0, 0, 0]))

   call assert_not_equals(INVALID_INDEX, h%find_index([0, 0, 1]))
   call assert_not_equals(INVALID_INDEX, h%find_index([0, 0, 2]))
   call assert_not_equals(INVALID_INDEX, h%find_index([0, 0, 3]))
   call assert_not_equals(INVALID_INDEX, h%find_index([0, 1, 0]))
   call assert_not_equals(INVALID_INDEX, h%find_index([0, 2, 0]))
   call assert_not_equals(INVALID_INDEX, h%find_index([0, 3, 0]))
   call assert_not_equals(INVALID_INDEX, h%find_index([1, 0, 0]))
   call assert_not_equals(INVALID_INDEX, h%find_index([2, 0, 0]))
   call assert_not_equals(INVALID_INDEX, h%find_index([3, 0, 0]))

   call assert_not_equals(INVALID_INDEX, h%find_index([0, 1, 1]))
   call assert_not_equals(INVALID_INDEX, h%find_index([0, 1, 2]))
   call assert_not_equals(INVALID_INDEX, h%find_index([0, 2, 1]))

   call assert_not_equals(INVALID_INDEX, h%find_index([1, 0, 1]))
   call assert_not_equals(INVALID_INDEX, h%find_index([1, 0, 2]))
   call assert_not_equals(INVALID_INDEX, h%find_index([2, 0, 1]))

   call assert_not_equals(INVALID_INDEX, h%find_index([1, 1, 0]))
   call assert_not_equals(INVALID_INDEX, h%find_index([1, 2, 0]))
   call assert_not_equals(INVALID_INDEX, h%find_index([2, 1, 0]))

   call assert_equals(INVALID_INDEX, h%find_index([1, 1, 1]))
   call assert_equals(INVALID_INDEX, h%find_index([2, 1, 1]))
   call assert_equals(INVALID_INDEX, h%find_index([0, 1, 4]))

   call assert_equals(19, h%entries())

   call h%free()
end subroutine test_triangular_limited


subroutine test_triangular_coupling()
   implicit none
   type(HStructure) :: h
   integer i
   call h%init(2, 3)

   i = h%find_index([1, 1]);
   call assert_true(all([1, 1] == h%vecind(i)), "1) Check")
   call assert_true(all([2, 1] == h%vecind(h%indab(i, 1))), "1) indab(1)")
   call assert_true(all([1, 2] == h%vecind(h%indab(i, 2))), "1) indab(2)")
   call assert_true(all([0, 1] == h%vecind(h%indbl(i, 1))), "1) indbl(1)")
   call assert_true(all([1, 0] == h%vecind(h%indbl(i, 2))), "1) indbl(2)")

   i = h%find_index([3, 0]);
   call assert_true(all([3, 0] == h%vecind(i)), "2) Check")
   call assert_equals(INVALID_INDEX, h%indab(i, 1), "2) indab(1)")
   call assert_equals(INVALID_INDEX, h%indab(i, 2), "2) indab(2)")
   call assert_true(all([2, 0] == h%vecind(h%indbl(i, 1))), "2) indbl(1)")
   call assert_equals(INVALID_INDEX, h%indbl(i, 2), "2) indbl(2)")

   call h%free()
end subroutine test_triangular_coupling

end module hstruct_test
