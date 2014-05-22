! Linked list of hierarchy-structure entries HStructureEntry (more precise: base
! node of a linked list). Only contains basic functions to add an entry and the
! free ALL of the space allocated by the entries in the list.  The retrieval of
! entries needs to be done manually by traversing the list:
!
! Usage:
!     type(HStructureListPointer) :: list
!     call p%add([0, 0], 0)
!     call p%add([0, 1], 1)
!     call p%add([1, 0], 2)
!
!     list%p%val                    ! == ([1, 0], 2)
!     list%p%next%val               ! == ([0, 1], 1)
!     list%p%next%next%val           ! == ([0, 0], 0)
!     list%p%next%next%next         ! <- unassociated
!
! Note: This type is intended to be used in the HStructureTable hash table and
!       does not provide enough functionality for a stand-alone usage.

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
