program fruit_driver
use fruit
use hstructtab_test
use hstruct_test
use noisegen_test
use dynarray_test
use sparse_test
use randomsparse_test

implicit none

call init_fruit

! HStructureTable test suite
call test_hashtab_add_and_lookup

! HStructure test suite
call test_triangular_completeness
call test_triangular_limited
call test_triangular_cutoff
call test_triangular_cutoff_limited
call test_triangular_coupling

! ExponentialNoiseGenerator test suite
call test_noise_generator

! DynArray test
call test_dynarray_reshape

! SparseMatrix test suite
call test_sparse_matmul
call test_sparse_setblock
call test_randomsparse_resetval

call fruit_summary

end program fruit_driver
