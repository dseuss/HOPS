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

   call assert_not_equals(INVALID_INDEX, h%find_index([0, 0]), "T 00")

   call assert_not_equals(INVALID_INDEX, h%find_index([1, 0]), "T 10")
   call assert_not_equals(INVALID_INDEX, h%find_index([0, 1]), "T 01")

   call assert_not_equals(INVALID_INDEX, h%find_index([2, 0]), "T 20")
   call assert_not_equals(INVALID_INDEX, h%find_index([1, 1]), "T 11")
   call assert_not_equals(INVALID_INDEX, h%find_index([0, 2]), "T 02")

   call assert_not_equals(INVALID_INDEX, h%find_index([3, 0]), "T 30")
   call assert_not_equals(INVALID_INDEX, h%find_index([2, 1]), "T 21")
   call assert_not_equals(INVALID_INDEX, h%find_index([1, 2]), "T 12")
   call assert_not_equals(INVALID_INDEX, h%find_index([0, 3]), "T 03")

   call assert_equals(INVALID_INDEX, h%find_index([1, 3]), "T 13")
   call assert_equals(INVALID_INDEX, h%find_index([0, 8]), "T 08")
   call assert_equals(INVALID_INDEX, h%find_index([4, 3]), "T 43")

   call assert_equals(10, h%entries(), "T Length")

   call h%free()
end subroutine test_triangular_completeness


subroutine test_triangular_limited()
   implicit none
   type(HStructure) :: h
   call h%init(3, 3, 2)

   call assert_not_equals(INVALID_INDEX, h%find_index([0, 0, 0]), "TL 000")

   call assert_not_equals(INVALID_INDEX, h%find_index([0, 0, 1]), "TL 001")
   call assert_not_equals(INVALID_INDEX, h%find_index([0, 0, 2]), "TL 002")
   call assert_not_equals(INVALID_INDEX, h%find_index([0, 0, 3]), "TL 003")
   call assert_not_equals(INVALID_INDEX, h%find_index([0, 1, 0]), "TL 010")
   call assert_not_equals(INVALID_INDEX, h%find_index([0, 2, 0]), "TL 020")
   call assert_not_equals(INVALID_INDEX, h%find_index([0, 3, 0]), "TL 030")
   call assert_not_equals(INVALID_INDEX, h%find_index([1, 0, 0]), "TL 100")
   call assert_not_equals(INVALID_INDEX, h%find_index([2, 0, 0]), "TL 200")
   call assert_not_equals(INVALID_INDEX, h%find_index([3, 0, 0]), "TL 300")

   call assert_not_equals(INVALID_INDEX, h%find_index([0, 1, 1]), "TL 011")
   call assert_not_equals(INVALID_INDEX, h%find_index([0, 1, 2]), "TL 012")
   call assert_not_equals(INVALID_INDEX, h%find_index([0, 2, 1]), "TL 021")

   call assert_not_equals(INVALID_INDEX, h%find_index([1, 0, 1]), "TL 101")
   call assert_not_equals(INVALID_INDEX, h%find_index([1, 0, 2]), "TL 102")
   call assert_not_equals(INVALID_INDEX, h%find_index([2, 0, 1]), "TL 201")

   call assert_not_equals(INVALID_INDEX, h%find_index([1, 1, 0]), "TL 110")
   call assert_not_equals(INVALID_INDEX, h%find_index([1, 2, 0]), "TL 121")
   call assert_not_equals(INVALID_INDEX, h%find_index([2, 1, 0]), "TL 210")

   call assert_equals(INVALID_INDEX, h%find_index([1, 1, 1]), "TL 111")
   call assert_equals(INVALID_INDEX, h%find_index([2, 1, 1]), "TL 211")
   call assert_equals(INVALID_INDEX, h%find_index([0, 1, 4]), "TL 014")

   call assert_equals(19, h%entries(), "TL Length")

   call h%free()
end subroutine test_triangular_limited


subroutine test_triangular_cutoff()
   implicit none
   type(HStructure) :: h

   call h%init(3, 2, 0, [-1, 1, 0])

   call assert_not_equals(INVALID_INDEX, h%find_index([0, 0, 0]), "TC 000")

   call assert_not_equals(INVALID_INDEX, h%find_index([0, 1, 0]), "TC 010")

   call assert_not_equals(INVALID_INDEX, h%find_index([1, 0, 0]), "TC 100")
   call assert_not_equals(INVALID_INDEX, h%find_index([1, 1, 0]), "TC 110")
   call assert_not_equals(INVALID_INDEX, h%find_index([2, 0, 0]), "TC 200")

   call assert_equals(INVALID_INDEX, h%find_index([3, 0, 0]), "TC 300")
   call assert_equals(INVALID_INDEX, h%find_index([0, 0, 1]), "TC 001")
   call assert_equals(INVALID_INDEX, h%find_index([1, 1, 1]), "TC 111")
   call assert_equals(INVALID_INDEX, h%find_index([2, 1, 1]), "TC 211")
   call assert_equals(INVALID_INDEX, h%find_index([0, 1, 4]), "TC 014")
   call assert_equals(INVALID_INDEX, h%find_index([1, 0, 2]), "TC 102")

   call assert_equals(5, h%entries(), "TC Length")

   call h%free()
end subroutine test_triangular_cutoff


subroutine test_triangular_cutoff_limited()
   implicit none
   type(HStructure) :: h
   call h%init(3, 3, 2, [-1, -1, 1])

   call assert_not_equals(INVALID_INDEX, h%find_index([0, 0, 0]), "TCL 000")

   call assert_not_equals(INVALID_INDEX, h%find_index([0, 0, 1]), "TCL 001")
   call assert_not_equals(INVALID_INDEX, h%find_index([0, 1, 0]), "TCL 010")
   call assert_not_equals(INVALID_INDEX, h%find_index([0, 2, 0]), "TCL 020")
   call assert_not_equals(INVALID_INDEX, h%find_index([0, 3, 0]), "TCL 030")
   call assert_not_equals(INVALID_INDEX, h%find_index([1, 0, 0]), "TCL 100")
   call assert_not_equals(INVALID_INDEX, h%find_index([2, 0, 0]), "TCL 200")
   call assert_not_equals(INVALID_INDEX, h%find_index([3, 0, 0]), "TCL 300")

   call assert_not_equals(INVALID_INDEX, h%find_index([0, 1, 1]), "TCL 011")
   call assert_not_equals(INVALID_INDEX, h%find_index([0, 2, 1]), "TCL 021")

   call assert_not_equals(INVALID_INDEX, h%find_index([1, 0, 1]), "TCL 101")
   call assert_not_equals(INVALID_INDEX, h%find_index([2, 0, 1]), "TCL 201")

   call assert_not_equals(INVALID_INDEX, h%find_index([1, 1, 0]), "TCL 110")
   call assert_not_equals(INVALID_INDEX, h%find_index([1, 2, 0]), "TCL 120")
   call assert_not_equals(INVALID_INDEX, h%find_index([2, 1, 0]), "TCL 210")

   call assert_equals(INVALID_INDEX, h%find_index([1, 1, 1]), "TCL 111")
   call assert_equals(INVALID_INDEX, h%find_index([2, 1, 1]), "TCL 211")
   call assert_equals(INVALID_INDEX, h%find_index([0, 1, 4]), "TCL 014")
   call assert_equals(INVALID_INDEX, h%find_index([1, 0, 2]), "TCL 102")
   call assert_equals(INVALID_INDEX, h%find_index([1, 1, 2]), "TCL 112")

   call assert_equals(15, h%entries(), "TCL Length")

   call h%free()
end subroutine test_triangular_cutoff_limited



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
