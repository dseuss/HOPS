module hstructlist

use hstructent, only: HStructureEntry

private
type, public :: HStructureList
   type(HStructureList), public, pointer :: next
   type(HStructureEntry), public         :: val

   contains
   procedure, public :: add
   procedure, public :: free
end type HStructureList

! Needed for allocatable array
type, public :: HStructureListPointer
   type(HStructureList), pointer :: p
end type HStructureListPointer


contains

subroutine add(self, k, ind)
   implicit none
   class(HStructureList) :: self
   integer, intent(in)   :: k(:), ind

   type(HStructureList), pointer :: old

   old => self%next
   allocate(self%next)
   self%next%next => old
   call self%next%val%init(k, ind)
end subroutine add

recursive subroutine free(self)
   implicit none
   class(HStructureList) :: self

   call self%val%free()
   if (associated(self%next)) then
      call self%next%free()
      deallocate(self%next)
   end if

end subroutine free

end module hstructlist
