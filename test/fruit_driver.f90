program fruit_driver
use fruit
use hstructtab_test
use hstruct_test
use noisegen_test

implicit none

call init_fruit

! HStructureTable test suite
call test_hashtab_add_and_lookup

! HStructure test suite
call test_triangular_completeness
call test_triangular_limited
call test_triangular_coupling

! ExponentialNoiseGenerator test suite
call test_noise_generator

call fruit_summary

end program fruit_driver
