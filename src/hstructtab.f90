! NOTE: All final keywords are incompatible with gfortran 4.8
module hstructtab

use hstructlist, only: HStructureList, HStructureListPointer
implicit none

integer, parameter :: INVALID_INDEX = 0

type, public :: HStructureTable

   private
   integer :: buckets_
   integer :: modes_
   type(HStructureListPointer), allocatable :: data_(:)

   contains
   private
   procedure hash
   final :: free

   procedure, public :: init
   procedure, public :: add
   procedure, public :: get

end type HStructureTable

contains

subroutine init(self, buckets, modes)
   implicit none
   class(HStructureTable) :: self
   integer, intent(in)   :: buckets
   integer, intent(in)   :: modes

   integer :: i

   self%buckets_ = buckets
   self%modes_ = modes
   allocate(self%data_(buckets))
   do i=1, self%buckets_
      allocate(self%data_(i)%p)
   end do
end subroutine init


subroutine free(self)
   implicit none
   type(HStructureTable) :: self

   integer :: i

   if (allocated(self%data_)) then
      do i=1, self%buckets_
         call self%data_(i)%p%free()
      end do
      deallocate(self%data_)
   end if
end subroutine free


! Convert the vector-indices vind to valid hash
!  :vind:
function hash(self, k)
   implicit none
   class(HStructureTable) :: self
   integer, intent(in)    :: k(self%modes_)
   integer                :: hash

   integer :: i

   ! FIXME This will not work for many modes, find better hash function
   hash = 0
   do i=1, self%modes_
      hash = hash + k(i) * 7**i
   end do
   hash = modulo(hash, self%buckets_) + 1
end function hash

!
! Add an entry with vector-index k and index ind
subroutine add(self, k, ind)
   use hstructent
   implicit none
   class(HStructureTable) :: self
   integer, intent(in)    :: k(self%modes_)
   integer                :: ind

   integer :: hash

   hash = self%hash(k)
   call self%data_(hash)%p%add(k, ind)
end subroutine add


! Lookup vector-index k; returns INVALID_INDEX if it was not found
function get(self, k) result(ind)
   implicit none
   class(HStructureTable) :: self
   integer, intent(in)    :: k(self%modes_)
   integer                :: ind

   type(HStructureList), pointer :: current
   integer hashval

   hashval = self%hash(k)

   current => self%data_(hashval)%p%next
   do while(associated(current))
      if(all(k == current%val%vind)) then
         ind = current%val%iind
         return
      end if
      current => current%next
   end do
   ind = 0
end function get

end module hstructtab
