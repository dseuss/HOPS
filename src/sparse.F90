module sparse
use system
use dynarray_int
use dynarray_cmplx
implicit none
private

type, public :: SparseMatrix
   integer :: size_
   integer :: block_
   integer :: nnz_

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
   procedure, public :: print
end type SparseMatrix

contains

subroutine init(self, size, block, chunk)
   implicit none
   class(SparseMatrix)           :: self
   integer, intent(in)           :: size
   integer, intent(in)           :: block
   integer, intent(in), optional :: chunk

   integer :: N

   self%size_ = size
   self%block_ = block
   self%nnz_ = 0
   if (present(chunk)) then
      N = chunk
   else
      N = size
   end if

   call self%cooI_%init(N)
   call self%cooJ_%init(N)
   call self%cooA_%init(N)
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
         call self%add((i-1) * self%block_ + m, (j-1) * self%block_ + n, &
               vals(m, n))
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

   integer :: job(8), info, i, j, N
   ! TODO Remove duplicates by adding them together
   ! NOTE Sorting (setting job(1)=2) does not work properly with pre 11.1 MKL
   ! see http://software.intel.com/en-us/forums/topic/375484
   job = [2, 0, 1, 0, self%nnz_, 0, 0, 0]

   allocate(self%csrI_(self%size_ + 1))
   allocate(self%csrJ_(self%nnz_))
   allocate(self%csrA_(self%nnz_))

   call mkl_zcsrcoo(job, self%size_, self%csrA_, self%csrJ_, self%csrI_, &
         self%nnz_, self%cooA_%data, self%cooI_%data, self%cooJ_%data, info)

   call self%cooI_%free()
   call self%cooJ_%free()
   call self%cooA_%free()
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


! This only works if all duplicates have been summed and the array is sorted!
subroutine print(self)
   implicit none
   class(SparseMatrix) :: self
   integer :: i, j, counter

   counter = 1
   do i = 1, self%size_
      do j = 1, self%size_
         ! Remember: csrI_/csrJ_ are zero based!
         if ((counter <= self%nnz_) .and. (self%csrJ_(counter) == j - 1) &
               .and. (counter <= self%csrI_(i+1))) then
         write(*, fmt='(a,f5.2,a,f5.2,a)', advance='no') &
               '(', real(self%csrA_(counter)), ',', &
               aimag(self%csrA_(counter)), ')'
         counter = counter + 1
         else
            write(*, fmt='(a)', advance='no') '      X      '
         end if

         if (mod(j, self%block_) == 0) then
            write(*, fmt='(a)', advance='no') ' '
         end if
      end do
      write(*, fmt='(a)', advance='yes') ''
      if (mod(i, self%block_) == 0) then
         write(*, fmt='(a)', advance='yes') ''
      end if
   end do

   print *, 'csrA_'
   do i = 1, size(self%csrA_)
      write(*, fmt='(a, f6.2, a, f6.2, a)', advance='no') &
            '(', real(self%csrA_(i)), ', ', aimag(self%csrA_(i)), ') '
   end do
   print *, ''
   print *, "csrI:", self%csrI_
   print *, "csrJ:", self%csrJ_
end subroutine print

end module sparse
