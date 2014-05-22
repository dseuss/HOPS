! Hash table for hierarchy-structure entries HStructureEntry. Stores the
! integer indices of a structure indexed by the vector index.
!
! Usage:
!     type(HStructureTable) :: tab
!     call tab%init(3, 2)
!     call tab%add([0, 0], 1)
!     call tab%add([1, 0], 2)
!     call tab%add([0, 1], 3)
!
!     call tab%get([1, 0])    ! == 2
!     call tab%get([3, 0])    ! == INVALID_INDEX
!
! Algorithm:
!     Create an array of linked lists where the length of the array is equal
!     to buckets_. The hash-function maps each vector index to a number between
!     1 and buckets_. Adding an vector index (k) with integer index i to the
!     hash table results in adding the entry ((k), i) to the linked list at
!     position hash(k).
!
!     Retrieving the integer index corresponding to an vector index (k)
!     is performed by first a linear search in the linked list at position
!     hash(k).
!
!     In order to be effective, the linked lists should be as short as possible
!     and evenly populated. This can be done by a "good" choice for the hash
!     function.

module hstructtab

use hstructlist, only: HStructureList, HStructureListPointer
implicit none

integer, parameter :: INVALID_INDEX = 0

type, public :: HStructureTable

   private
   integer :: buckets_                          ! Number of buckets
   integer :: modes_                            ! Length of the vector indices
   type(HStructureListPointer), allocatable :: data_(:)

   contains
   private
   procedure hash

   procedure, public :: init
   procedure, public :: free
   procedure, public :: add
   procedure, public :: get

end type HStructureTable

contains

subroutine init(self, buckets, modes)
   ! :buckets: Number of buckets (individual linked lists)
   ! :modes: Length of the vector indices to be added to the list

   implicit none
   class(HStructureTable) :: self
   integer, intent(in)    :: buckets
   integer, intent(in)    :: modes

   integer :: i

   self%buckets_ = buckets
   self%modes_ = modes
   allocate(self%data_(buckets))

   ! Initialize each linked list by setting the base node to point to NULL
   do i=1, self%buckets_
      allocate(self%data_(i)%p)
      self%data_(i)%p%next => null()
   end do
end subroutine init


subroutine free(self)
   implicit none
   class(HStructureTable) :: self

   integer :: i

   if (allocated(self%data_)) then
      do i=1, self%buckets_
         call self%data_(i)%p%free()
         deallocate(self%data_(i)%p)
      end do
      deallocate(self%data_)
   end if

end subroutine free


pure function hash(self, k)
   ! Convert the vector index (k) to a valid hash ( = 1,...,self%buckets)
   !
   !     hash = sum_i  k_i * 7^i
   !
   ! Taking the base 7 is random. Taking base self%modes_ would be better, but
   ! generates too large numbers for large problems.
   !
   ! :k[self%modes_]: Vector index to hashed
   !
   !
   ! FIXME Better choice of hash function? Check performance!

   implicit none
   class(HStructureTable), intent(in) :: self
   integer, intent(in)                :: k(self%modes_)
   integer                            :: hash

   integer :: i

   hash = 0
   do i=1, self%modes_
      hash = hash + k(i) * 7**i
   end do
   hash = modulo(hash, self%buckets_) + 1
end function hash


subroutine add(self, k, ind)
   ! Add hierarchy structure entry ((k), ind) to the table
   !
   ! :k[self%modes]: Vector index
   ! :ind: Integer index

   implicit none
   class(HStructureTable) :: self
   integer, intent(in)    :: k(self%modes_)
   integer                :: ind

   integer :: hash

   hash = self%hash(k)
   call self%data_(hash)%p%add(k, ind)
end subroutine add


function get(self, k) result(ind)
   ! Looks up the vector index k in the hash table
   !
   ! :k[self%modes_]: Vector index
   ! :result: Integer index of k or INVALID_INDEX if k is not present

   implicit none
   class(HStructureTable), intent(in) :: self
   integer, intent(in)                :: k(self%modes_)
   integer                            :: ind

   type(HStructureList), pointer :: current
   integer hash

   hash = self%hash(k)

   current => self%data_(hash)%p%next        ! Base node of list number hash
   ! Linear search through list
   do while(associated(current))
      if(all(k == current%val%vind)) then
         ind = current%val%iind
         return
      end if
      current => current%next
   end do
   ind = INVALID_INDEX
end function get

end module hstructtab
