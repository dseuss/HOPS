! Unified interface for sparse matrix libraries. Supports only basic actions
! such as building the sparse matrix and matrix-vector multiplication
! (sparse-dense) for complex quadratic sparse matrices.
!
! Usage:
!     type(SparseMatrix) :: A
!     complex(dp) :: x(2), y(2)
!
!     call A%init(2, 1)
!     call A%add(1, 1, 0.5)
!     call A%finalize()
!
!     x = [1, 0]
!     call A%multiply(x, y)   ! y = Ax
!
! Currently supported backends:
!     - Intel MKL
!
! WARNING: Module may not function properly with pre 11.1 MKL
!          (see http://software.intel.com/en-us/forums/topic/375484).

module sparse
use system
use dynarray_int
use dynarray_cmplx
implicit none
include 'mkl_spblas.fi'

private

type, public :: SparseMatrix
   integer :: size_                    ! total size of matrix
   integer :: block_                   ! block size of matrix
   integer :: nnz_                     ! number non-zero elements

   ! Representation of sparse matrix in coordinate form
   type(IntDynamicArray)   :: cooI_    ! Row-indices
   type(IntDynamicArray)   :: cooJ_    ! Column-indices
   type(CmplxDynamicArray) :: cooA_    ! Entries

   ! Representation of sparse matrix in compressed sparse row format
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
   ! :size: total dimension of sparse matrix
   ! :block(1): block size of matrix, used for add_block
   ! :chunk(size): chunk size used for dynamical arrays, number of additional
   !               entries allocated in each resizing. Size should be fine
   !               except for very large or dense matrices.

   implicit none
   class(SparseMatrix)           :: self
   integer, intent(in)           :: size
   integer, intent(in), optional :: block
   integer, intent(in), optional :: chunk

   integer :: N

   self%size_ = size
   self%nnz_ = 0

   if (present(block)) then
      self%block_ = block
   else
      self%block_ = 1
   end if

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
   ! Adds the value val at the (i,j) position in the matrix. If A[i,j] is not
   ! zero (default value), val is added to the current value.
   ! :i: row index
   ! :j: column index
   ! :val: value

   implicit none
   class(SparseMatrix)     :: self
   integer, intent(in)     :: i
   integer, intent(in)     :: j
   complex(dp), intent(in) :: val

#ifdef NDEBUG
   if ((i < 1) .or. (i > self%size_)) then
      print *, "OUT OF BOUND ERROR: Adding i to SparseMatrix failed"
      print *, "0 < ", i, " < ", self%size_
      stop -1
   end if
   if ((j < 1) .or. (j > self%size_)) then
      print *, "OUT OF BOUND ERROR: Adding j to SparseMatrix failed"
      print *, "0 < ", j, " < ", self%size_
      stop -1
   end if
#endif

   self%nnz_ = self%nnz_ + 1
   call self%cooI_%add(i)
   call self%cooJ_%add(j)
   call self%cooA_%add(val)
end subroutine add


subroutine add_block(self, i, j, vals)
   ! Adds the submatrix vals of size block*block to the matrix at the (i,j)-
   ! submatrix position, i.e. starting at ((i-1)*block, (j-1)*block) component.
   !
   ! :i: submatrix row index
   ! :j: submatrix column index
   ! :vals[self%block_, self%block_]: submatrix values to add

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


subroutine finalize(self)
   ! Converts the matrix from coo to csr format. Afterwards, no more components
   ! can be added!

   implicit none
   class(SparseMatrix) :: self

   integer :: job(8), info, i, j, N
   ! TODO Remove duplicates by adding them together
   ! NOTE Sorting (setting job(1)=2) does not work properly with pre 11.1 MKL
   ! see http://software.intel.com/en-us/forums/topic/375484
   job = [2, 0, 1, 0, self%nnz_, 0, 0, 0]
   ! job(1)=2: the matrix in the coordinate format is converted to the CSR
   !           format, and the column indices in CSR representation are sorted
   !           in the increasing order within each row.
   ! job(2)=0: zero-based indexing for the matrix in CSR format is used;
   ! job(3)=1: one-based indexing for the matrix in coordinate format is used.
   ! job(5)=nnz: sets number of the non-zero elements of the matrix if job(1)=2
   ! job(6)=0: all arrays acsr, ja, ia are filled in for the output storage.

   allocate(self%csrI_(self%size_ + 1))
   allocate(self%csrJ_(self%nnz_))
   allocate(self%csrA_(self%nnz_))

   call mkl_zcsrcoo(job, self%size_, self%csrA_, self%csrJ_, self%csrI_, &
         self%nnz_, self%cooA_%data, self%cooI_%data, self%cooJ_%data, info)

   call self%cooI_%free()
   call self%cooJ_%free()
   call self%cooA_%free()

   ! self%nnz_ = csr_sum_duplicates(self%csrI_, self%csrJ_, self%csrA_)
end subroutine finalize


subroutine multiply(self, x, y, alpha)
   ! Calculate matrix multiplication of instance with vector x:  y = alpha * Ax
   !
   ! :x[self%size_]: Vector to multiply matrix with
   ! :y[self%size_]: Result
   ! :alpha(1): Scalar multiplier

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


subroutine print(self)
   ! Prints the current matrix to screen.
   !
   ! Note: Only works after finalizing and only if duplicates have been summed,
   !       i.e. there is at most one elmenent for each row-column-index
   !       combination.

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                               Helper Functions                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! function csr_sum_duplicates(Ai, Aj, Ax) result(nnz)
!    integer, intent(inout)     :: Ai(:)
!    integer, intent(inout)     :: Aj(:)
!    complex(dp), intent(inout) :: Ax(:)
!    integer                    :: nnz

!    integer :: r1, r2, i, j, jj
!    complex(dp) :: x

!    nnz = 1
!    r2 = 1
!    do i = 1, size(Ai) - 1
!       r1 = r2
!       r2 = Ai(i+1)
!       jj = r1
!       do while (jj < r2)
!          j = Aj(jj)
!          x = Ax(jj)
!          jj = jj + 1
!          do while (jj < r2)
!             if (Aj(jj) == j) then
!                x = x + Ax(jj)
!                jj = jj + 1
!             else
!                exit
!             end if
!          end do
!          Aj(nnz) = j
!          Ax(nnz) = x
!          nnz = nnz + 1
!       end do
!       Ai(i+1) = nnz
!    end do
! end function

end module sparse
