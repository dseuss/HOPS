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
