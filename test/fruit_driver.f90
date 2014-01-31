program fruit_driver
use fruit
use hstructtab_test
use hstruct_test

implicit none

call init_fruit
! Hashtab test suite
call test_hashtab_add_and_lookup

! HStructure test suite
call test_triangular_completeness
call test_triangular_limited
call test_triangular_coupling

call fruit_summary

end program fruit_driver
