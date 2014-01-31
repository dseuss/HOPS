module hstructtab_test
use fruit
implicit none

contains

subroutine test_hashtab_add_and_lookup()
   use hstructtab, only: HStructureTable, INVALID_INDEX
   implicit none
   type(HStructureTable) :: ht

   call ht%init(3, 3)
   call ht%add([0, 0, 0], 1)
   call ht%add([1, 0, 0], 2)
   call ht%add([1, 3, 0], 3)
   call ht%add([3, 0, 8], 4)
   call ht%add([3, 0, 6], 5)
   call ht%add([4, 3, 6], 6)
   call ht%add([9, 3, 1], 7)
   call ht%add([2, 7, 0], 8)
   call ht%add([1, 1, 1], 9)
   print *, "Done with adding!"

   call assertEquals(1, ht%get([0, 0, 0]))
   call assertEquals(2, ht%get([1, 0, 0]))
   call assertEquals(3, ht%get([1, 3, 0]))
   call assertEquals(4, ht%get([3, 0, 8]))
   call assertEquals(5, ht%get([3, 0, 6]))
   call assertEquals(6, ht%get([4, 3, 6]))
   call assertEquals(7, ht%get([9, 3, 1]))
   call assertEquals(8, ht%get([2, 7, 0]))
   call assertEquals(9, ht%get([1, 1, 1]))
end subroutine test_hashtab_add_and_lookup

end module hstructtab_test
