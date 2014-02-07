module hierarchy
use system
use hstruct, only: HStructure
use hstructtab, only: INVALID_INDEX
use sparse, only: SparseMatrix
implicit none
private
public init, free, run_trajectory_z0_rk4, run_trajectory_z0_zvode, &
      trajectory_step_z0


integer, parameter :: &
      RK_COPIES = 6, &
      RK_OLD = 1, &
      RK_NEW = 2, &
      RK_1 = 3, &
      RK_2 = 4, &
      RK_3 = 5, &
      RK_4 = 6

integer :: &
      depth_, &
      modes_, &
      dim_, &
      size_, &
      tSteps_
real(dp) :: &
      tLength_, &
      dt_

type(SparseMatrix) :: linProp_

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

   call struct%init(modes_, depth, populated_modes)
   size_ = struct%entries() * dim_
   call linProp_%init(size_, dim_)

   call setup_propagator_simple(struct, modes_, dim_, g, gamma, Omega, h, &
         Lmap, with_terminator)

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
   call linProp_%free()
end subroutine free


subroutine setup_propagator_simple(struct, modes, hs_dim, g, gamma, Omega, h, &
         Lmap, with_terminator)
   implicit none
   type(HStructure), intent(in) :: struct
   integer, intent(in)          :: modes
   integer, intent(in)          :: hs_dim
   complex(dp), intent(in)      :: g(modes)
   real(dp), intent(in)         :: gamma(modes)
   real(dp), intent(in)         :: Omega(modes)
   complex(dp), intent(in)      :: h(hs_dim, hs_dim)
   integer, intent(in)          :: Lmap(modes)
   logical, intent(in)          :: with_terminator
   !----------------------------------------------------------------------------
   integer :: entry, mode, mode_term, k(modes), k_term(modes), i
   integer :: indbl, indab, indbl_term, ind_term
   complex(dp) :: identity(hs_dim, hs_dim), k_dot_w
   ! TODO Error checks

   identity = 0._dp
   do i = 1, hs_dim
      identity(i, i) = (1._dp, 0._dp)
   end do

   do entry = 1, struct%entries()
      k = struct%vecind(entry)
      k_dot_w = dot_product(k, gamma) + ii*dot_product(k, Omega)
      call linProp_%add_block(entry, entry, -ii*h - k_dot_w*identity)
      do mode = 1, modes
         indbl = struct%indbl(entry, mode)
         indab = struct%indab(entry, mode)

         if (indbl /= INVALID_INDEX) then
            call linProp_%add(hs_dim * (entry - 1) + Lmap(mode), &
                  hs_dim * (indbl - 1) + Lmap(mode), &
                  k(mode) * g(mode))
         end if

         if (indab /= INVALID_INDEX) then
            call linProp_%add(hs_dim * (entry - 1) + Lmap(mode), &
                  hs_dim * (indab - 1) + Lmap(mode), &
                  (-1._dp, 0._dp))
         ! else if (with_terminator) then
         !    k_term = k
         !    k_term(mode) = k_term(mode) + 1
         !    k_dot_w = dot_product(k_term, gamma) + ii*dot_product(k_term, Omega)

         !    do mode_term = 1, modes
         !       if (mode_term == mode) then
         !          call linProp_%add(hs_dim * (entry - 1) + Lmap(mode), &
         !                hs_dim * (entry - 1) + Lmap(mode), &
         !                -k_term(mode_term) * g(mode_term) / k_dot_w)
         !          cycle
         !       end if

         !       if (struct%indbl(entry ,mode_term) == INVALID_INDEX) then
         !          cycle
         !       end if

         !       ind_term = struct%indab(struct%indbl(entry, mode_term), mode)
         !       if (ind_term == INVALID_INDEX) then
         !          cycle
         !       end if

         !       call linProp_%add(hs_dim * (entry - 1) + Lmap(mode_term), &
         !             hs_dim * (ind_term - 1) + Lmap(mode_term), &
         !             -k_term(mode_term) * g(mode_term) / k_dot_w)
         !    end do
         end if
      end do
   end do

   call linProp_%finalize()
end subroutine setup_propagator_simple


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

subroutine trajectory_jac_z0(NEQ, T, Y, ML, MU, PD, NROWPD, RPAR, IPAR)
   integer :: NEQ
   real(dp) :: T
   complex(dp) :: Y(NEQ)
   integer :: ML
   integer :: MU
   complex(dp) :: PD(NROWPD, NEQ)
   integer :: NROWPD
   real(dp) :: RPAR(*)
   integer :: IPAR(*)
end subroutine trajectory_jac_z0

function run_trajectory_z0_zvode(psi0) result(psi)
   implicit none
   complex(dp), intent(in) :: psi0(dim_)
   complex(dp)             :: psi(tSteps_, dim_)

   integer :: t

   integer, parameter :: MF = 10 ! (or 22 for stiff)
   complex(dp), allocatable :: psi_hierarchy(:), ZWORK(:)
   real(dp), allocatable :: RWORK(:)
   integer, allocatable :: IWORK(:)
   integer :: LZW, LRW, LIW, flag
   real(dp), parameter :: RPAR(1) = [0._dp]
   integer, parameter :: IPAR(1) = [0]
   real(dp) :: t_current

   print *, 'Using ZVODE'
   allocate(psi_hierarchy(size_))
   psi_hierarchy = 0._dp
   psi_hierarchy(1:dim_) = psi0

   select case (MF)
      case(10)
         LZW = 15*size_*2
         LIW = 30
      case(22)
         LZW = 8*size_ * 2*size_**2
         LIW = 30 + size_
      end select
   LRW = 20 + size_

   allocate(ZWORK(LZW), RWORK(LRW), IWORK(LIW))

   t_current = 0._dp
   flag = 1
   do t = 1, tSteps_ - 1
      call ZVODE (&
            trajectory_step_z0, & ! F
            size_,              & ! NEQ
            psi_hierarchy,      & ! Y
            t_current,          & ! T
            t*dt_,              & ! TOUT
            1,                  & ! ITOL (ATOL is scalar)
            1.D-7,              & ! RTOL
            0._dp,             & ! ATOL
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
            MF                  & ! MF
            )
      print *, t_current
      psi(t+1, :) = psi_hierarchy(1:dim_)
   end do
   deallocate(psi_hierarchy, ZWORK, RWORK, IWORK)
end function run_trajectory_z0_zvode

end module hierarchy
