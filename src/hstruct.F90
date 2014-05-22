! Module to manage the structure of HOPS. The hierarchy structure is determined
! by the set of all valid vector indices combined with an enumeration of the
! latter.
!
! A hierarchy structure is characterized by the number of modes N (the dimension
! of the vector index (k)) and the truncation condition. Currently, there are
! three different criteria, which determine valid (k):
!
!     - depth D of the hierarchy (triangular truncation): (k) is a valid vector
!     index of the hierarchy, iff
!                             k_1 + ... + k_N <= D
!
!     - populated modes P: (k) is a valid vector index of the hierarchy, iff
!     at most P components of (k) are non-zero.
!     Example: P=2, N=3, D=4: valid [0,0,0], [4,0,0], [3,1,0]
!                         not valid: [1,1,1] (3 non-zero components)
!
!     - manual cutoff (c): (k) is a valid vector index of the hierarchy, iff
!               for all i=1,...,N:   k_i <= c_i
!
! FIXME vecind_ should be column-major order!

module hstruct
use system
use hstructtab, only: HStructureTable
implicit none
private

! Average number of entries in each linked list in the hash table
!     => number of buckets for hash table is given by int(entries / BUCKETFILL)
integer, parameter :: BUCKETFILL = 10

type, public :: HStructure
   private
   integer :: &
         entries_, &                   ! Number of total entries in hierarchy
         modes_, &                     ! Number of modes = length(k)
         depth_, &                     ! Depth of the hierarchy
         populated_modes_              ! Max. number of populated modes

   type(HStructureTable) :: intind_    ! Hash-table to lookup integer-index
   integer, allocatable :: &
         cutoff_(:), &                 ! Manual cutoff vector
         vecind_(:, :), &              ! vecind_(i, :) = i-th vector index
         indab_(:, :), &               ! indab_(i, j) = integer-index of
                                       ! (k) + (e_j), where (k) = vecind_(i, :)
         indbl_(:, :)                  ! Same as indab_, but with "- (e_j)"

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
   ! Recursively calculate all valid vector indices. To start call as
   !
   !     call self%recursive_indices_(1, [0,0,...,0], 0)
   !
   !  Fills in the member variables of self.
   !
   ! :at_mode: Fill in all valid combinations for "at_mode"-th component of (k)
   ! :indices[self%modes_]: Components indices_i with i < at_mode are filled in
   ! :currpop: Number of currently populated states

   implicit none

   class(HStructure)      :: self
   integer, intent(in)    :: at_mode
   integer, intent(in)    :: indices(self%modes_)
   integer, intent(in)    :: currpop

   integer i, indices_inc(self%modes_), at_mode_inc, currpop_inc

   ! First we treat the at_mode-th component with no population
   !     -> this contributes irrespective of currpop
   ! Have we reached the last mode? -> Fill in!
   if (at_mode >= self%modes_) then
      call self%add_entry_(indices)
   ! no? go to next mode and call with same indices (as these are 0 by default)
   else
      at_mode_inc = at_mode + 1
      call self%recursive_indices_(at_mode_inc, indices, currpop)
   end if

   if (currpop >= self%populated_modes_) then
      return
   end if

   indices_inc = indices
   currpop_inc = currpop + 1
   ! Fill in at_mode-th mode with all valid components
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
   ! Adds the vector index k unless the manual cutoff prevents it
   !
   ! :k[self%modes_]: vector index to add

   implicit none
   class(HStructure) :: self
   integer, intent(in) :: k(self%modes_)

   if (all(k <= self%cutoff_)) then
      self%entries_ = self%entries_ + 1
      self%vecind_(self%entries_, :) = k
   end if
end subroutine add_entry_


subroutine init(self, modes, depth, populated_modes, cutoff)
   ! :modes: Number of modes in the hierarchy == length of the vector indices
   ! :depth: Depth of the hierarchy for the triangular truncation condition
   ! :populated_modes (modes): Max. number of modes with excitation
   ! :cutoff[modes] ([depth,...,depth]): Manual cutoff condition

   implicit none
   class(HStructure)             :: self
   integer, intent(in)           :: modes
   integer, intent(in)           :: depth
   integer, intent(in), optional :: populated_modes
   integer, intent(in), optional :: cutoff(modes)

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
   ! Look up the vector index k and return the corresponding integer index
   !
   ! :k[self%modes_]: Vector index to look up
   ! :result: Integer index corresponding to k

   implicit none
   class(HStructure), intent(in) :: self
   integer                       :: i
   integer, intent(in)           :: k(self%modes_)
   i = self%intind_%get(k)
end function find_index

#ifdef NDEBUG
function indab(self, n, mode)
#else
elemental function indab(self, n, mode)
#endif
   ! Look up the integer index belonging to (k) + (e_mode), where
   ! (k) = vecind_(n, :).
   !
   ! :n: Integer index of the "current" entry
   ! :mode: Mode number where are coupling to (above)
   ! :result: Integer index belonging to (k) + (e_mode)

   implicit none
   class(HStructure), intent(in) :: self
   integer, intent(in)           :: n
   integer, intent(in)           :: mode
   integer                       :: indab
#ifdef NDEBUG
   if ((mode < 0) .or. (mode > self%modes_)) then
      print *, 'ERROR in indab: mode out of bounds. Is', mode, &
            ', number of modes', self%modes_
      stop -1
   end if

   if ((n < 0) .or. (n > self%entries_)) then
      print *, 'ERROR in indab: entry out of bounds. Is', n, &
            ', number of entries', self%entries_
      stop -1
   end if
#endif

   indab = self%indab_(n, mode)
end function indab


#ifdef NDEBUG
function indbl(self, n, mode)
#else
elemental function indbl(self, n, mode)
#endif
   ! Look up the integer index belonging to (k) - (e_mode), where
   ! (k) = vecind_(n, :).
   !
   ! :n: Integer index of the "current" entry
   ! :mode: Mode number where are coupling to (below)
   ! :result: Integer index belonging to (k) - (e_mode)

   implicit none
   class(HStructure), intent(in) :: self
   integer, intent(in)           :: n
   integer, intent(in)           :: mode
   integer                       :: indbl
#ifdef NDEBUG
   if ((mode < 0) .or. (mode > self%modes_)) then
      print *, 'ERROR in indbl: mode out of bounds. Is', mode, &
            ', number of modes', self%modes_
      stop -1
   end if

   if ((n < 0) .or. (n > self%entries_)) then
      print *, 'ERROR in indbl: entry out of bounds. Is', n, &
            ', number of entries', self%entries_
      stop -1
   end if
#endif

   indbl = self%indbl_(n, mode)
end function indbl


#ifdef NDEBUG
function vecind(self, n)
#else
pure function vecind(self, n)
#endif
   ! Returns the vector index corresponding to the integer index n
   !
   ! :n: Integer index
   ! :result[self%modes_]: Vector index corresponding to n
   implicit none
   class(HStructure), intent(in) :: self
   integer, intent(in)           :: n
   integer                       :: vecind(self%modes_)
#ifdef NDEBUG
   if ((n < 0) .or. (n > self%entries_)) then
      print *, 'ERROR in vecind: entry out of bounds. Is', n, &
            ', number of entries', self%entries_
      stop -1
   end if
#endif

   vecind = self%vecind_(n, :)
end function vecind


elemental function entries(self)
   ! Getter function for the member variable self%entries_
   !
   ! :returns: The number of entries in the hierarchy
   implicit none
   class(HStructure), intent(in) :: self
   integer                       :: entries
   entries = self%entries_
end function entries


subroutine print(self)
   ! Prints the hierarchy structure on screen, just for debugging...

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
   !
   ! :n:
   ! :k:
   ! :result: (n over k)

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
