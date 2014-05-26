! Module for time-dependent sparse matrices. Instead of knowing the values of
! the matrix during its preparation, we only know "indices", which are used
! as a placeholder. Later, we can determine (and update) the value of the
! (i,j)-th component of the matrix by passing in an array of values.
!
! Usage:
!     type(TimeDepSparseMatrix) :: A
!     complex(dp) :: x(2), y(2), Z(2)
!
!     call A%init(2, 2)
!     call A%add(1, 1, 1)
!     call A%add(1, 2, 1)
!     call A%add(2, 1, 2)
!     call A%finalize()

!     ! Calculate y = Ax, where A = [[1., 1.], [4., 0.]]
!     x = [1._dp, 2._dp]
!     Z = [1._dp, 4._dp]
!     call A%multiply(Z, x, y)
!
! WARNING: Module may not function properly with pre 11.1 MKL
!          (see http://software.intel.com/en-us/forums/topic/375484).
!
! TODO More effective to do multiplication by hand (than setting data for MKL)
! TODO add_multiply function to calculate y += Ax

module timedep_sparse
use system
use dynarray_int
use dynarray_cmplx
implicit none

! FIXME Can't use header file since we abuse ordering-function from float.
! include 'mkl_spblas.fi'

private

type, public :: TimeDepSparseMatrix
   integer :: size_                          ! total size of matrix
   integer :: num_elem_                      ! Number of distinct values
   integer :: nnz_                           ! Number of non-zero elements

   ! Representation in coordinate form
   type(IntDynamicArray)   :: cooI_          ! Row indices
   type(IntDynamicArray)   :: cooJ_          ! Column indices
   type(IntDynamicArray)   :: cooID_         ! Index variables for values
   type(CmplxDynamicArray) :: cooC_          ! Coefficients

   ! Representation in compressed sparse row format
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
end type TimeDepSparseMatrix

contains

subroutine init(self, size, num_elem, chunk)
   ! :size: total dimension of sparse matrix
   ! :num_elem: Max. number of distinct time-dependent components
   ! :chunk(size): chunk size used for dynamical arrays, number of additional
   !               entries allocated in each resizing. Size should be fine
   !               except for very large or dense matrices.

   implicit none
   class(TimeDepSparseMatrix)    :: self
   integer, intent(in)           :: size
   integer, intent(in)           :: num_elem
   integer, intent(in), optional :: chunk

   integer :: N

   self%size_ = size
   self%nnz_ = 0
   self%num_elem_ = num_elem
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
   class(TimeDepSparseMatrix) :: self

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
   ! Adds the index id with coefficient coeff at the (i,j) position in the
   ! matrix. This means that to the final matrix for multiplication (with data-
   ! array Z) an entry coeff * Z[id] is added at the (i.j)-position.
   ! :i: row index
   ! :j: column index
   ! :id: Placeholder index
   ! :coeff(1): Time-independent Coefficient

   implicit none
   class(TimeDepSparseMatrix)        :: self
   integer, intent(in)               :: i
   integer, intent(in)               :: j
   integer, intent(in)               :: id
   complex(dp), intent(in), optional :: coeff

#ifdef NDEBUG
   if ((i < 1) .or. (i > self%size_)) then
      print *, "OUT OF BOUND ERROR: Adding i to TimeDepSparseMatrix failed"
      print *, "0 < ", i, " < ", self%size_
   end if
   if ((j < 1) .or. (j > self%size_)) then
      print *, "OUT OF BOUND ERROR: Adding j to TimeDepSparseMatrix failed"
      print *, "0 < ", j, " < ", self%size_
   end if
   if ((id < 1) .or. (id > self%num_elem_)) then
      print *, "OUT OF BOUND ERROR: Wrong process id in TimeDepSparseMatrix"
      print *, "0 < ", id, " < ", self%num_elem_
   end if
#endif

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
   ! Converts the matrix from coo to csr format. Afterwards, no more components
   ! can be added.

   implicit none
   class(TimeDepSparseMatrix) :: self

   integer :: job(8), info, i, j, N
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
   allocate(self%csrID_(self%nnz_))
   allocate(self%csrC_(self%nnz_))

   ! Note the abuse of scsrcoo: Supposed to convert single precision float
   ! matrices, but here we convert integer data (same number of bits for 32bit
   ! integers).
   call mkl_scsrcoo(job, self%size_, self%csrID_, self%csrJ_, self%csrI_, &
         self%nnz_, self%cooID_%data, self%cooI_%data, self%cooJ_%data, info)
   call mkl_zcsrcoo(job, self%size_, self%csrC_, self%csrJ_, self%csrI_, &
         self%nnz_, self%cooC_%data, self%cooI_%data, self%cooJ_%data, info)

   call self%cooI_%free()
   call self%cooJ_%free()
   call self%cooID_%free()
   call self%cooC_%free()
end subroutine finalize


subroutine multiply(self, Z, x, y, alpha)
   ! Calculates the matrix multiplication of the instance sparse matrix with
   ! vector x:  y = alpha * Ax. Here, A is the matrix defined in its coordinate
   ! form
   !
   !                   A[i,j] = cooC[n] * Z[cooID[n]]
   !
   ! where i is the n-th element of cooI and j is the n-th element of cooJ.
   ! If there are more than one entries at the (i,j) position, they are summed.
   !
   ! :Z[self%num_elem_]: The values of the entries
   ! :x[self%size_]: Vector to multiply matrix with
   ! :y[self%size_]: Result
   ! :alpha(1): Scalar multiplier

   implicit none
   class(TimeDepSparseMatrix), intent(in)   :: self
   complex(dp), intent(in)           :: Z(self%num_elem_)
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

end module timedep_sparse
