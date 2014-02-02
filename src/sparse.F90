module sparse
use system
use dynarray_int
use dynarray_cmplx
implicit none
include 'mkl.fi'
private

type, public :: SparseMatrix
   integer :: size_
   integer :: block_
   integer :: nnz_

   ! integer, allocatable :: cooI_(:)
   type(IntDynamicArray)   :: cooI_
   type(IntDynamicArray)   :: cooJ_
   type(CmplxDynamicArray) :: cooA_

   integer, allocatable     :: csrI_(:)
   integer, allocatable     :: csrJ_(:)
   complex(dp), allocatable :: csrA_(:)

   contains
   private

   procedure, public :: init
   procedure, public :: free
   procedure, public :: finalize
   procedure, public :: add
   procedure, public :: add_block
   procedure, public :: multiply
end type SparseMatrix

contains

subroutine init(self, size, block, chunk)
   implicit none
   class(SparseMatrix)           :: self
   integer, intent(in)           :: size
   integer, intent(in)           :: block
   integer, intent(in), optional :: chunk

   self%size_ = size
   self%block_ = block
   self%nnz_ = 0
   if (present(chunk)) then
      call self%cooI_%init(chunk)
      call self%cooJ_%init(chunk)
      call self%cooA_%init(chunk)
   else
      call self%cooI_%init(1)
      call self%cooJ_%init(1)
      call self%cooA_%init(1)
   end if
end subroutine init


subroutine free(self)
   implicit none
   class(SparseMatrix) :: self

   call self%cooI_%free()
   call self%cooJ_%free()
   call self%cooA_%free()
   call array_free(self%csrI_)
   call array_free(self%csrJ_)
   call array_free(self%csrA_)
end subroutine free


subroutine add(self, i, j, val)
   implicit none
   class(SparseMatrix)     :: self
   integer, intent(in)     :: i
   integer, intent(in)     :: j
   complex(dp), intent(in) :: val

   if ((i < 1) .or. (i > self%size_)) then
      print *, "OUT OF BOUND ERROR: Adding i to SparseMatrix failed"
      print *, "0 < ", i, " < ", self%size_
   end if
   if ((j < 1) .or. (j > self%size_)) then
      print *, "OUT OF BOUND ERROR: Adding j to SparseMatrix failed"
      print *, "0 < ", j, " < ", self%size_
   end if

   self%nnz_ = self%nnz_ + 1
   call self%cooI_%add(i)
   call self%cooJ_%add(j)
   call self%cooA_%add(val)
end subroutine add


subroutine add_block(self, i, j, vals)
   implicit none
   class(SparseMatrix) :: self
   integer, intent(in) :: i
   integer, intent(in) :: j
   complex(dp), intent(in) :: vals(self%block_, self%block_)

   integer :: m, n

   do m = 1, self%block_
      do n = 1, self%block_
         call self%add(i + m - 1, j + n - 1, vals(m, n))
      end do
   end do
end subroutine add_block


! subroutine add_block_flat(self, i, j, vals)
!    implicit none
!    class(SparseMatrix) :: self
!    integer, intent(in) :: i
!    integer, intent(in) :: j
!    complex(dp), intent(in) :: vals(self%block_*self%block_)

!    call self%add_block(i, j, reshape(vals, [self%block_, self%block_]))
! end subroutine add_block_flat


subroutine finalize(self)
   implicit none
   class(SparseMatrix) :: self

   integer :: job(8), info
   ! NOTE Sorting (setting job(1)=2) does not work properly
   ! see http://software.intel.com/en-us/forums/topic/375484
   job = [1, 0, 1, 0, self%nnz_, 0, 0, 0]

   allocate(self%csrI_(self%size_+1), self%csrJ_(self%nnz_), &
         self%csrA_(self%nnz_))
   call mkl_zcsrcoo(job, self%size_, self%csrA_, self%csrJ_, self%csrI_, &
         self%nnz_, self%cooA_%data, self%cooI_%data, self%cooJ_%data, info)

   call self%cooI_%free()
   call self%cooJ_%free()
   call self%cooA_%free()

   ! call csr_sort_indices(self%csrI_, self%csrJ_, self%csrA_)
   ! call csr_sum_duplicates(self%csrI_, self%csrJ_, self%csrA_)

   print *, ""
   print *, "csrA_= ", self%csrA_
   print *, "csrI_= ", self%csrI_
   print *, "csrJ_= ", self%csrJ_

   ! TODO Remove duplicates by adding them together
end subroutine finalize


! Calculate y = alpha * dot(SELF,x)
subroutine multiply(self, x, y, alpha)
   implicit none
   class(SparseMatrix), intent(in)   :: self
   complex(dp), intent(in)           :: x(self%size_)
   complex(dp), intent(out)          :: y(self%size_)
   complex(dp), intent(in), optional :: alpha

   character, parameter :: trans = 'n'

   call mkl_cspblas_zcsrgemv(trans, self%size_, self%csrA_, self%csrI_, &
         self%csrJ_, x, y)

   if (present(alpha)) then
      y = alpha * y
   end if
end subroutine multiply


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                               HELPER FUNCTIONS                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define DTYPE complex(dp)
#include "sparse_template.src"
#undef DTYPE

end module sparse
