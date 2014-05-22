! Module to manage the structure of HOPS. The hierarchy structure is determined
! by the set of all valid vector indices combined with an enumeration of the
! latter.
!
! The valid entries are determined by the truncation condition.
module hstruct
use system
use hstructtab, only: HStructureTable
implicit none
private

integer, parameter :: BUCKETFILL = 10

type, public :: HStructure
   private
   integer :: &
         entries_, &
         modes_, &
         depth_, &
         populated_modes_

   type(HStructureTable) :: intind_
   integer, allocatable :: &
         cutoff_(:), &
         vecind_(:, :), &
         indab_(:, :), &
         indbl_(:, :)

   contains
   private
   procedure :: recursive_indices_
   procedure :: add_entry_

   procedure, public :: init
   procedure, public :: free
   procedure, public :: find_index
   procedure, public :: indab
   procedure, public :: indbl
   procedure, public :: vecind
   procedure, public :: entries
   procedure, public :: print
end type HStructure

contains

recursive subroutine recursive_indices_(self, at_mode, indices, currpop)
   implicit none

   class(HStructure)      :: self
   integer, intent(in)    :: at_mode
   integer, intent(in)    :: indices(self%modes_)
   integer, intent(in)    :: currpop

   integer i, indices_inc(self%modes_), at_mode_inc, currpop_inc

   ! first we treat the at_mode - unpopulated case
   if (at_mode >= self%modes_) then
      call self%add_entry_(indices)
   else
      at_mode_inc = at_mode + 1
      call self%recursive_indices_(at_mode_inc, indices, currpop)
   end if

   if (currpop >= self%populated_modes_) then
      return
   end if

   indices_inc = indices
   currpop_inc = currpop + 1
   do i = 1, self%depth_ - sum(indices)
      indices_inc(at_mode) = i
      if ((at_mode >= self%modes_)) then
         call self%add_entry_(indices_inc)
      else
         call self%recursive_indices_(at_mode_inc, indices_inc, currpop_inc)
      end if
   end do
end subroutine recursive_indices_


subroutine add_entry_(self, k)
   implicit none
   class(HStructure) :: self
   integer, intent(in) :: k(self%modes_)

   if (all(k <= self%cutoff_)) then
      self%entries_ = self%entries_ + 1
      self%vecind_(self%entries_, :) = k
   end if
end subroutine add_entry_

subroutine init(self, modes, depth, populated_modes, cutoff)
   implicit none
   class(HStructure)             :: self
   integer, intent(in)           :: modes
   integer, intent(in)           :: depth
   integer, intent(in), optional :: populated_modes
   integer, intent(in), optional :: cutoff(:)

   !----------------------------------------------------------------------------
   integer i, j, max_entries, k(modes), currpop

   self%modes_ = modes
   self%depth_ = depth
   if (present(populated_modes)) then
      if (populated_modes < 1) then
         self%populated_modes_ = modes
      else
         self%populated_modes_ = populated_modes
      endif
   else
      self%populated_modes_ = modes
   end if

   allocate(self%cutoff_(modes))
   if (present(cutoff)) then
      do i = 1, modes
         if (cutoff(i) < 0) then
            self%cutoff_(i) = depth
         else
            self%cutoff_(i) = cutoff(i)
         end if
      end do
   else
      self%cutoff_ = depth
   end if

   ! Calculate the maximum number of entries
   max_entries = 0
   do i = 0, depth
      max_entries = max_entries + choose(i + modes - 1, modes - 1)
   end do

   call self%intind_%init(max_entries / BUCKETFILL + 1, modes)
   allocate(self%vecind_(max_entries, modes))

   k = 0
   self%entries_ = 0
   currpop = 0
   call self%recursive_indices_(1, k, currpop)

   ! FIXME Resize vecind_
   allocate(self%indab_(self%entries_, self%modes_))
   allocate(self%indbl_(self%entries_, self%modes_))
   do i = 1, self%entries_
      call self%intind_%add(self%vecind_(i, :), i)
   end do

   ! FIXME Parallelize this
   do i = 1, self%entries_
      do j = 1, self%modes_
         k = self%vecind_(i, :)
         k(j) = k(j) - 1
         self%indbl_(i, j) = self%find_index(k)
         k(j) = k(j) + 2
         self%indab_(i, j) = self%find_index(k)
      end do
   end do

end subroutine init


subroutine free(self)
   implicit none
   class(HStructure) :: self
   call array_free(self%cutoff_)
   call array_free(self%indab_)
   call array_free(self%indbl_)
   call array_free(self%vecind_)
   call self%intind_%free()
end subroutine free


function find_index(self, k) result(i)
   implicit none
   class(HStructure), intent(in) :: self
   integer                       :: i
   integer, intent(in)           :: k(self%modes_)
   i = self%intind_%get(k)
end function find_index


elemental function indab(self, n, mode)
   implicit none
   class(HStructure), intent(in) :: self
   integer, intent(in)           :: n
   integer, intent(in)           :: mode
   integer                       :: indab
   ! FIXME Add some checks
   indab = self%indab_(n, mode)
end function indab


elemental function indbl(self, n, mode)
   implicit none
   class(HStructure), intent(in) :: self
   integer, intent(in)           :: n
   integer, intent(in)           :: mode
   integer                       :: indbl
   ! FIXME Add some checks
   indbl = self%indbl_(n, mode)
end function indbl

pure function vecind(self, n)
   implicit none
   class(HStructure), intent(in) :: self
   integer, intent(in)           :: n
   integer                       :: vecind(self%modes_)
   ! FIXME Add some checks
   vecind = self%vecind_(n, :)
end function vecind


elemental function entries(self)
   implicit none
   class(HStructure), intent(in) :: self
   integer                       :: entries
   entries = self%entries_
end function entries


subroutine print(self)
   implicit none
   class(HStructure) :: self
   integer :: ind
   print *, ''
   print *, '------------------------------------------------------------------'
   print *, 'Depth=', self%depth_, '  Populated Modes=', self%populated_modes_
   print*, 'Cutoff=', self%cutoff_
   do ind = 1, self%entries_
      print *, ind, "=>", self%vecind_(ind, :)
   end do
   print *, '------------------------------------------------------------------'
end subroutine print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                               HELPER FUNCTIONS                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure function choose(n, k) result (res)
   ! Returns the binomial coefficient (n over k) by calculating Pascals
   ! triangle.
   implicit none
   integer, intent (in) :: n, k
   integer              :: res

   integer b(0:n+1), i, j

   if (k > n) then
      res = 0
      return
   end if

   b(0) = 1
   do i = 1, n
      b(i) = 1
      do j = i-1, 1, -1
         b(j) = b(j) + b(j-1)
      end do !j
   end do !i

   res = b(k)
end function choose

end module hstruct
