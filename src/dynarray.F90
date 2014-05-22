! Module interface for dynamical arrays
!
! :TYPENAME: Name of the newly created type (dynamical array with entries of
!            type DTYPE)
! :DTYPE: Type of the entries of the array
! :PADVAL: Default value of the entries


module dynarray_int
#define TYPENAME IntDynamicArray
#define DTYPE integer
#define PADVAL 0
#include "dynarray_template.src"
#undef TYPENAME
#undef DTYPE
#undef PADVAL
end module dynarray_int

module dynarray_cmplx
use system, only: dp
#define TYPENAME CmplxDynamicArray
#define DTYPE complex(dp)
#define PADVAL cmplx(0., 0., dp)
#include "dynarray_template.src"
#undef TYPENAME
#undef DTYPE
#undef PADVAL
end module dynarray_cmplx
