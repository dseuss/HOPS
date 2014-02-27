module hierarchy
use system
use hstruct, only: HStructure
use hstructtab, only: INVALID_INDEX
use sparse, only: SparseMatrix
use randomsparse, only: RandomSparseMatrix
use noisegen, only: ExponentialNoiseGenerator, init_random_seed
implicit none
private
public init, free, run_trajectory_z0_rk4, run_trajectory_z0_zvode, &
      trajectory_step_z0, run_trajectory_rk4


integer, parameter :: &
      RK_COPIES = 7, &
      RK_OLD = 1, &
      RK_NEW = 2, &
      RK_1 = 3, &
      RK_2 = 4, &
      RK_3 = 5, &
      RK_4 = 6, &
      RK_BUF = 7

integer :: &
      depth_, &
      modes_, &
      dim_, &
      tSteps_
integer, public :: &
      size_
integer, allocatable :: &
      Lmap_(:)
complex(dp), allocatable :: &
      g_(:), &
      h_(:, :)
real(dp), allocatable :: &
      gamma_(:), &
      Omega_(:)
real(dp) :: &
      tLength_, &
      dt_

type(SparseMatrix) :: linProp_
type(RandomSparseMatrix) :: noiseProp_
type(RandomSparseMatrix) :: nonlinProp_

type(ExponentialNoiseGenerator), allocatable :: noisegen_(:)

contains


!!!!!!!!!!!!!!!!!!!!!!
!  Setup & Managing  !
!!!!!!!!!!!!!!!!!!!!!!

subroutine init(tLength, tSteps, depth, g, gamma, Omega, h, Lmap, &
         with_terminator, populated_modes)
   implicit none
   real(dp), intent(in)    :: tLength
   integer, intent(in)     :: tSteps
   integer, intent(in)     :: depth
   complex(dp), intent(in) :: g(:)
   real(dp), intent(in)    :: gamma(:)
   real(dp), intent(in)    :: Omega(:)
   complex(dp), intent(in) :: h(:, :)
   integer, intent(in)     :: Lmap(:)
   logical, intent(in)     :: with_terminator
   integer, intent(in)     :: populated_modes
   !----------------------------------------------------------------------------

   type(HStructure) :: struct

   ! TODO Error checks
   tLength_ = tLength
   tSteps_ = tSteps
   dt_ = tLength / (tSteps_ - 1)
   dim_ = size(h(1, :))
   modes_ = size(g)

   allocate(h_(dim_, dim_), g_(modes_), gamma_(modes_), Omega_(modes_), &
         Lmap_(modes_))
   h_ = h
   ! TODO Why is this conj(g)
   g_ = conjg(g)
   gamma_ = gamma
   Omega_ = Omega
   Lmap_ = Lmap

   call struct%init(modes_, depth, populated_modes)
   size_ = struct%entries() * dim_

   call linProp_%init(size_, dim_)
   call setup_propagator_simple(struct, with_terminator)
   call setup_noiseprop_simple(struct, with_terminator)
   call setup_nonlinprop_simple(struct, with_terminator)

   call init_random_seed()

   !----------------------------------------------------------------------------
   !TODO Make this optional
   ! call struct%print()
   ! call linProp_%print()

   ! print *, "Hsys= ", h
   ! print *, "g=", g
   ! print *, "gamma=", gamma
   ! print *, "Omega=", Omega
   ! print *, "Lmap=", Lmap
   !----------------------------------------------------------------------------

   call struct%free()
end subroutine init


subroutine free()
   implicit none
   integer :: i

   call linProp_%free()
   call noiseProp_%free()
   call nonlinProp_%free()
   do i = 1, size(noisegen_)
      call noisegen_(i)%free()
   end do

   deallocate(h_, g_, gamma_, Omega_, Lmap_, noisegen_)
end subroutine free


subroutine setup_propagator_simple(struct, with_terminator)
   implicit none
   type(HStructure), intent(in) :: struct
   logical, intent(in)          :: with_terminator
   !----------------------------------------------------------------------------
   integer :: entry, mode, mode_term, k(modes_), k_term(modes_), i
   integer :: indbl, indab, indbl_term, ind_term
   complex(dp) :: identity(dim_, dim_), k_dot_w
   ! TODO Error checks

   identity = 0._dp
   do i = 1, dim_
      identity(i, i) = (1._dp, 0._dp)
   end do

   do entry = 1, struct%entries()
      k = struct%vecind(entry)
      k_dot_w = dot_product(k, gamma_) + ii*dot_product(k, Omega_)
      call linProp_%add_block(entry, entry, -ii*h_ - k_dot_w*identity)
      do mode = 1, modes_
         indbl = struct%indbl(entry, mode)
         indab = struct%indab(entry, mode)

         if (indbl /= INVALID_INDEX) then
            call linProp_%add(dim_ * (entry - 1) + Lmap_(mode), &
                  dim_ * (indbl - 1) + Lmap_(mode), &
                  k(mode) * g_(mode))
         end if

         if (indab /= INVALID_INDEX) then
            call linProp_%add(dim_ * (entry - 1) + Lmap_(mode), &
                  dim_ * (indab - 1) + Lmap_(mode), &
                  (-1._dp, 0._dp))
         else if (with_terminator) then
            k_term = k
            k_term(mode) = k_term(mode) + 1
            k_dot_w = dot_product(k_term, gamma_) &
                  + ii*dot_product(k_term, Omega_)

            do mode_term = 1, modes_
               if (mode_term == mode) then
                  call linProp_%add(dim_ * (entry - 1) + Lmap_(mode), &
                        dim_ * (entry - 1) + Lmap_(mode), &
                        -k_term(mode_term) * g_(mode_term) / k_dot_w)
                  cycle
               end if

               if (struct%indbl(entry ,mode_term) == INVALID_INDEX) then
                  cycle
               end if

               ind_term = struct%indab(struct%indbl(entry, mode_term), mode)
               if (ind_term == INVALID_INDEX) then
                  cycle
               end if

               call linProp_%add(dim_ * (entry - 1) + Lmap_(mode_term), &
                     dim_ * (ind_term - 1) + Lmap_(mode_term), &
                     -k_term(mode_term) * g_(mode_term) / k_dot_w)
            end do
         end if
      end do
   end do

   call linProp_%finalize()
end subroutine setup_propagator_simple


subroutine setup_noiseprop_simple(struct, with_terminator)
   implicit none
   type(HStructure), intent(in) :: struct
   logical, intent(in)          :: with_terminator
   !----------------------------------------------------------------------------

   integer :: i, j, counter, entry
   ! Copies of bcf-parameters for noise generator
   complex(dp) :: gc(modes_)
   real(dp) :: gammac(modes_), Omegac(modes_)

   call noiseProp_%init(size_, dim_)
   do entry = 1, struct%entries()
      do i = 1, dim_
         call noiseProp_%add((entry - 1)*dim_ + i, (entry - 1)*dim_ + i, i)
      end do
   end do
   call noiseProp_%finalize()

   allocate(noisegen_(dim_))
   do i = 1, dim_
      counter = 1
      gc = (0._dp, 0._dp)
      gammac = 0._dp
      Omegac = 0._dp
      do j = 1, modes_
         if (Lmap_(j) == i) then
            gc(counter) = g_(j)
            gammac(counter) = gamma_(j)
            Omegac(counter) = Omega_(j)
            counter = counter + 1
         end if
      end do
      call noisegen_(i)%init(dt_/2., 2*tSteps_, gc, gammac, Omegac)
   end do
end subroutine setup_noiseprop_simple


subroutine setup_nonlinprop_simple(struct, with_terminator)
   implicit none
   type(HStructure), intent(in) :: struct
   logical, intent(in)          :: with_terminator
   !----------------------------------------------------------------------------

   integer :: entry, mode, j, mode_term, ind_term, indab, k(modes_), &
         k_term(modes_)
   complex(dp) :: k_dot_w

   call nonlinProp_%init(size_, dim_)

   do entry = 1, struct%entries()
      k = struct%vecind(entry)
      do mode = 1, modes_
         indab = struct%indab(entry, mode)

         if (indab /= INVALID_INDEX) then
            do j = 1, dim_
               call nonlinProp_%add((entry - 1)*dim_ + j, &
                     (indab - 1)*dim_ + j, &
                     Lmap_(mode))
            end do

         else if (with_terminator) then
            k_term = k
            k_term(mode) = k_term(mode) + 1
            k_dot_w = dot_product(k_term, gamma_) &
                  + ii*dot_product(k_term, Omega_)

            do mode_term = 1, modes_
               if (mode_term == mode) then
                  call nonlinProp_%add(dim_ * (entry - 1) + Lmap_(mode_term), &
                        dim_ * (entry - 1) + Lmap_(mode_term), &
                        Lmap_(mode), &
                        k_term(mode_term) * g_(mode_term) / k_dot_w)
                  cycle
               end if

               if (struct%indbl(entry ,mode_term) == INVALID_INDEX) then
                  cycle
               end if

               ind_term = struct%indab(struct%indbl(entry, mode_term), mode)
               if (ind_term == INVALID_INDEX) then
                  cycle
               end if

               call nonlinProp_%add(dim_ * (entry - 1) + Lmap_(mode_term), &
                     dim_ * (entry - 1) + Lmap_(mode_term), &
                     Lmap_(mode), &
                     k_term(mode_term) * g_(mode_term) / k_dot_w)
            end do
         end if
      end do
   end do

   call nonlinProp_%finalize()
end subroutine setup_nonlinprop_simple

!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Transfer Calculation  !
!!!!!!!!!!!!!!!!!!!!!!!!!!
! TODO Make this more compliant with usual integrator-rules
subroutine trajectory_step_rk4(dt, hierarchy, memTerms, Zt, rk_step)
   implicit none
   real(dp), intent(in)       :: dt
   complex(dp), intent(inout) :: hierarchy(size_, RK_COPIES)
   complex(dp), intent(inout) :: memTerms(modes_, RK_COPIES)
   complex(dp), intent(in)    :: Zt(dim_)
   integer, intent(in)        :: rk_step

   complex(dp) :: Z_cp(dim_)
   complex(dp) :: ExpL(dim_)
   integer :: i

   ExpL = Abs(hierarchy(1:dim_, RK_NEW))**2 / &
         dot_product(hierarchy(1:dim_, RK_NEW), hierarchy(1:dim_, RK_NEW))

   Z_cp = conjg(Zt)
   ! Z_cp = (Zt)
   do i = 1, modes_
      Z_cp(Lmap_(i)) = Z_cp(Lmap_(i)) + memTerms(i, RK_NEW)
   end do

   memTerms(:, rk_step) = dt * (conjg(g_) * ExpL(Lmap_) &
         + (-gamma_ + ii*Omega_) * memTerms(:, RK_NEW))

   call linProp_%multiply(hierarchy(:, RK_NEW), hierarchy(:, rk_step))
   ! TODO Do this without extra buffer?
   call noiseProp_%multiply(Z_cp, hierarchy(:, RK_NEW), &
         hierarchy(:, RK_BUF))
   hierarchy(:, rk_step) = hierarchy(:, rk_step) +  hierarchy(:, RK_BUF)
   call nonlinprop_%multiply(expl, hierarchy(:, rk_new), hierarchy(:, RK_BUF))
   hierarchy(:, rk_step) = dt_ * (hierarchy(:, rk_step) + hierarchy(:, RK_BUF))
end subroutine trajectory_step_rk4


function run_trajectory_rk4(psi0) result(psi)
   implicit none
   complex(dp), intent(in) :: psi0(dim_)
   complex(dp) :: psi(tSteps_, dim_)

   complex(dp) :: hierarchy(size_, RK_COPIES), Z(dim_, 2*tSteps_)
   complex(dp) :: memTerms(modes_, RK_COPIES)
   integer :: t, i
   !TODO Error Check

   psi(1, :) = psi0
   hierarchy = 0._dp
   hierarchy(1:dim_, RK_OLD) = psi0
   memTerms = (0._dp, 0._dp)

   do i = 1, dim_
      Z(i, :) = noisegen_(i)%get_realization()
   end do

   do t = 2, tSteps_
      hierarchy(:, RK_NEW) = hierarchy(:, RK_OLD)
      memTerms(:, RK_NEW) = memTerms(:, RK_OLD)
      call trajectory_step_rk4(dt_, hierarchy, memTerms, Z(:, 2*t-3), RK_1)

      hierarchy(:, RK_NEW) = hierarchy(:, RK_OLD) + .5*hierarchy(:, RK_1)
      memTerms(:, RK_NEW) = memTerms(:, RK_OLD) + .5*memTerms(:, RK_1)
      call trajectory_step_rk4(dt_, hierarchy, memTerms, Z(:, 2*t-2), RK_2)

      hierarchy(:, RK_NEW) = hierarchy(:, RK_OLD) + .5*hierarchy(:, RK_2)
      memTerms(:, RK_NEW) = memTerms(:, RK_OLD) + .5*memTerms(:, RK_2)
      call trajectory_step_rk4(dt_, hierarchy, memTerms, Z(:, 2*t-2), RK_3)

      hierarchy(:, RK_NEW) = hierarchy(:, RK_OLD) + hierarchy(:, RK_3)
      memTerms(:, RK_NEW) = memTerms(:, RK_OLD) + memTerms(:, RK_3)
      call trajectory_step_rk4(dt_, hierarchy, memTerms, Z(:, 2*t-1), RK_4)

      hierarchy(:, RK_OLD) = hierarchy(:, RK_OLD) &
            + 1./6. * hierarchy(:, RK_1) &
            + 2./6. * hierarchy(:, RK_2) &
            + 2./6. * hierarchy(:, RK_3) &
            + 1./6. * hierarchy(:, RK_4)
      memTerms(:, RK_OLD) = memTerms(:, RK_OLD) &
            + 1./6. * memTerms(:, RK_1) &
            + 2./6. * memTerms(:, RK_2) &
            + 2./6. * memTerms(:, RK_3) &
            + 1./6. * memTerms(:, RK_4)
      psi(t, :) = hierarchy(1:dim_, RK_OLD)
   end do
end function run_trajectory_rk4


!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Spectrum Calculation  !
!!!!!!!!!!!!!!!!!!!!!!!!!!
function run_trajectory_z0_rk4(psi0) result(psi)
   implicit none
   complex(dp), intent(in) :: psi0(dim_)
   complex(dp)             :: psi(tSteps_, dim_)

   complex(dp) :: hierarchy(size_, RK_COPIES), dt
   integer :: t
   !TODO Error Check

   dt = tLength_ / (tSteps_ - 1)

   psi(1, :) = psi0
   hierarchy = 0._dp
   hierarchy(1:dim_, RK_OLD) = psi0

   do t = 2, tSteps_
      call linProp_%multiply(hierarchy(:, RK_OLD), hierarchy(:, RK_1), dt)

      hierarchy(:, RK_NEW) = hierarchy(:, RK_OLD) + .5*hierarchy(:, RK_1)
      call linProp_%multiply(hierarchy(:, RK_NEW), hierarchy(:, RK_2), dt)

      hierarchy(:, RK_NEW) = hierarchy(:, RK_OLD) + .5*hierarchy(:, RK_2)
      call linProp_%multiply(hierarchy(:, RK_NEW), hierarchy(:, RK_3), dt)

      hierarchy(:, RK_NEW) = hierarchy(:, RK_OLD) + hierarchy(:, RK_3)
      call linProp_%multiply(hierarchy(:, RK_NEW), hierarchy(:, RK_4), dt)

      hierarchy(:, RK_OLD) = hierarchy(:, RK_OLD) &
            + 1./6. * hierarchy(:, RK_1) &
            + 2./6. * hierarchy(:, RK_2) &
            + 2./6. * hierarchy(:, RK_3) &
            + 1./6. * hierarchy(:, RK_4)
      psi(t, :) = hierarchy(1:dim_, RK_OLD)
   end do
end function run_trajectory_z0_rk4


subroutine trajectory_step_z0(NEQ, T, Y, YDOT)
   implicit none
   integer, intent(in) :: NEQ
   real(dp), intent(in)     :: T
   complex(dp), intent(in)  :: Y(NEQ)
   complex(dp), intent(out) :: YDOT(NEQ)

   ! FIXME Disable this outside debug mode
   if (NEQ /= size_) then
      print *, 'Error in trajectory_step_z0: input/output has wrong shape.'
   end if

   call linProp_%multiply(Y, YDOT)
end subroutine trajectory_step_z0


! Dummy routine
subroutine trajectory_jac_z0()
end subroutine trajectory_jac_z0


function run_trajectory_z0_zvode(psi0, method, rtol, atol) result(psi)
   implicit none
   complex(dp), intent(in) :: psi0(dim_)
   complex(dp)             :: psi(tSteps_, dim_)
   integer, intent(in)     :: method
   real(dp), intent(in)    :: rtol
   real(dp), intent(in)    :: atol

   integer :: t

   complex(dp), allocatable :: psi_hierarchy(:), ZWORK(:)
   real(dp), allocatable :: RWORK(:)
   integer, allocatable :: IWORK(:)
   integer :: LZW, LRW, LIW, flag
   real(dp), parameter :: RPAR(1) = [0._dp]
   integer, parameter :: IPAR(1) = [0]
   real(dp) :: t_current

   select case (method)
      case(10)
         LZW = 15*size_*2
         LIW = 30
      case(22)
         LZW = 8*size_ * 2*size_**2
         LIW = 30 + size_
      end select
   LRW = 20 + size_

   allocate(ZWORK(LZW), RWORK(LRW), IWORK(LIW))
   allocate(psi_hierarchy(size_))
   psi_hierarchy = 0._dp
   psi_hierarchy(1:dim_) = psi0

   t_current = 0._dp
   flag = 1
   psi(1, :) = psi0
   do t = 1, tSteps_ - 1
      call ZVODE (&
            trajectory_step_z0, & ! F
            size_,              & ! NEQ
            psi_hierarchy,      & ! Y
            t_current,          & ! T
            t_current+dt_,      & ! TOUT
            1,                  & ! ITOL (ATOL is scalar)
            rtol,               & ! RTOL
            atol,               & ! ATOL
            1,                  & ! ITASK
            flag,               & ! ISTATE
            0,                  & ! IOPT
            ZWORK,              & ! ZWORK
            LZW,                & ! LZW
            RWORK,              & ! RWORK
            LRW,                & ! LRW
            IWORK,              & ! IWORK
            LIW,                & ! LIW
            trajectory_jac_z0,  & ! JAC (this is just dummy variable)
            method              & ! MF
            )
      ! FIXME Check the error flag
      psi(t+1, :) = psi_hierarchy(1:dim_)
   end do
   deallocate(psi_hierarchy, ZWORK, RWORK, IWORK)
end function run_trajectory_z0_zvode

end module hierarchy
