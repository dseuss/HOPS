module system
implicit none

integer, parameter, public :: &
      sp = kind(1.0), &
      dp = kind(1.d0)
complex(dp), parameter, public :: &
      ii = (0._dp, 1._dp)

public array_free

interface array_free
   procedure :: array_free_int1D
   procedure :: array_free_int2D
   procedure :: array_free_complex1D
end interface array_free

contains

subroutine array_free_int1D(array)
   implicit none
   integer, allocatable :: array(:)
   if (allocated(array)) then
      deallocate(array)
   end if
end subroutine array_free_int1D

subroutine array_free_int2D(array)
   implicit none
   integer, allocatable :: array(:, :)
   if (allocated(array)) then
      deallocate(array)
   end if
end subroutine array_free_int2D

subroutine array_free_complex1D(array)
   implicit none
   complex(dp), allocatable :: array(:)
   if (allocated(array)) then
      deallocate(array)
   end if
end subroutine array_free_complex1D

end module system
