module hstructent
use system
implicit none

private
type, public :: HStructureEntry
   integer, public, allocatable :: vind(:)
   integer, public              :: iind

   contains
   procedure, public :: init
   procedure, public :: free
end type HStructureEntry

contains

subroutine init(self, vind, iind)
   implicit none
   class(HStructureEntry) :: self
   integer, intent(in)    :: vind(:), iind

   if (allocated(self%vind)) then
      deallocate(self%vind)
   end if

   allocate(self%vind(size(vind)))
   self%vind = vind
   self%iind = iind
end subroutine init

subroutine free(self)
   implicit none
   class(HStructureEntry) :: self
   call array_free(self%vind)
end subroutine free

end module hstructent
