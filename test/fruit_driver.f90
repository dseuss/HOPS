program fruit_driver
use fruit
use hstructtab_test
implicit none

call init_fruit
call test_hashtab_add_and_lookup
call fruit_summary

end program fruit_driver
