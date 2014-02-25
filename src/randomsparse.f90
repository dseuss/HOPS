! TODO Make another one not based on MKL with direct multiplication (instead
! of resetting data.

module randomsparse
use system
use dynarray_int
use dynarray_cmplx
! include 'mkl_spblas.fi'
implicit none
private

type, public :: RandomSparseMatrix
   integer :: size_
   integer :: num_proc_
   integer :: nnz_

   type(IntDynamicArray)   :: cooI_
   type(IntDynamicArray)   :: cooJ_
   type(IntDynamicArray)   :: cooID_
   type(CmplxDynamicArray) :: cooC_

   integer, allocatable     :: csrI_(:)
   integer, allocatable     :: csrJ_(:)
   integer, allocatable     :: csrID_(:)
   complex(dp), allocatable :: csrC_(:)

   contains
   private

   procedure, public :: init
   procedure, public :: free
   procedure, public :: finalize
   procedure, public :: multiply
   procedure, public :: add
end type RandomSparseMatrix

contains

subroutine init(self, size, num_proc, chunk)
   implicit none
   class(RandomSparseMatrix)           :: self
   integer, intent(in)           :: size
   integer, intent(in)           :: num_proc
   integer, intent(in), optional :: chunk

   integer :: N

   self%size_ = size
   self%nnz_ = 0
   self%num_proc_ = num_proc
   if (present(chunk)) then
      N = chunk
   else
      N = size
   end if

   call self%cooI_%init(N)
   call self%cooJ_%init(N)
   call self%cooID_%init(N)
   call self%cooC_%init(N)
end subroutine init


subroutine free(self)
   implicit none
   class(RandomSparseMatrix) :: self

   call self%cooI_%free()
   call self%cooJ_%free()
   call self%cooID_%free()
   call self%cooC_%free()
   call array_free(self%csrI_)
   call array_free(self%csrJ_)
   call array_free(self%csrID_)
   call array_free(self%csrC_)
end subroutine free


subroutine add(self, i, j, id, coeff)
   implicit none
   class(RandomSparseMatrix)               :: self
   integer, intent(in)               :: i
   integer, intent(in)               :: j
   integer, intent(in)               :: id
   complex(dp), intent(in), optional :: coeff

   if ((i < 1) .or. (i > self%size_)) then
      print *, "OUT OF BOUND ERROR: Adding i to RandomSparseMatrix failed"
      print *, "0 < ", i, " < ", self%size_
   end if
   if ((j < 1) .or. (j > self%size_)) then
      print *, "OUT OF BOUND ERROR: Adding j to RandomSparseMatrix failed"
      print *, "0 < ", j, " < ", self%size_
   end if
   if ((id < 1) .or. (id > self%num_proc_)) then
      print *, "OUT OF BOUND ERROR: Wrong process id in RandomSparseMatrix"
      print *, "0 < ", id, " < ", self%num_proc_
   end if

   self%nnz_ = self%nnz_ + 1
   call self%cooI_%add(i)
   call self%cooJ_%add(j)
   call self%cooID_%add(id)
   if (present(coeff)) then
      call self%cooC_%add(coeff)
   else
      call self%cooC_%add((1._dp, 0._dp))
   end if
end subroutine add


subroutine finalize(self)
   implicit none
   class(RandomSparseMatrix) :: self

   integer :: job(8), info, i, j, N
   ! TODO Remove duplicates by adding them together
   ! NOTE Sorting (setting job(1)=2) does not work properly with pre 11.1 MKL
   ! see http://software.intel.com/en-us/forums/topic/375484
   job = [2, 0, 1, 0, self%nnz_, 0, 0, 0]

   allocate(self%csrI_(self%size_ + 1))
   allocate(self%csrJ_(self%nnz_))
   allocate(self%csrID_(self%nnz_))
   allocate(self%csrC_(self%nnz_))

   call mkl_scsrcoo(job, self%size_, self%csrID_, self%csrJ_, self%csrI_, &
         self%nnz_, self%cooID_%data, self%cooI_%data, self%cooJ_%data, info)
   call mkl_zcsrcoo(job, self%size_, self%csrC_, self%csrJ_, self%csrI_, &
         self%nnz_, self%cooC_%data, self%cooI_%data, self%cooJ_%data, info)

   call self%cooI_%free()
   call self%cooJ_%free()
   call self%cooID_%free()
   call self%cooC_%free()
end subroutine finalize


! Calculate y = alpha * dot(SELF,x)
! TODO Check if its not better to pre-allocate multiplier
subroutine multiply(self, Z, x, y, alpha)
   implicit none
   class(RandomSparseMatrix), intent(in)   :: self
   complex(dp), intent(in)           :: Z(self%num_proc_)
   complex(dp), intent(in)           :: x(self%size_)
   complex(dp), intent(out)          :: y(self%size_)
   complex(dp), intent(in), optional :: alpha

   complex(dp) :: multiplier(self%nnz_)
   character, parameter :: trans = 'n'

   multiplier = self%csrC_ * Z(self%csrID_)
   call mkl_cspblas_zcsrgemv(trans, self%size_, multiplier, self%csrI_, &
         self%csrJ_, x, y)

   if (present(alpha)) then
      y = alpha * y
   end if
end subroutine multiply

end module randomsparse
