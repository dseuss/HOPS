! Entry of the hierarchy-structure linked list. Each is characterized in two
! equivalent ways:
!     - its vector index (k), which is used in defining the hierarchy. Its
!       entries k_i describe the order of the corresponding auxilliary state
!       with respect to the functional derivative d/dZ_i
!     - its integer index, which arises from some enumeration of all auxilliary
!       states involved

module hstructent
use system
implicit none

private
type, public :: HStructureEntry
   integer, public, allocatable :: vind(:)         ! Vector index
   integer, public              :: iind            ! Integer index

   contains
   procedure, public :: init
   procedure, public :: free
end type HStructureEntry

contains

subroutine init(self, vind, iind)
   ! Intializes the hierarchy-structure entry with given vector index and
   ! integer index
   !
   ! :vind[:]: Vector index
   ! :iind: Integer index

   implicit none
   class(HStructureEntry) :: self
   integer, intent(in)    :: vind(:), iind

   call array_free(self%vind)
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
