module hierarchy
use system
use hstruct, only: HStructure
use hstructtab, only: INVALID_INDEX
use sparse, only: SparseMatrix
use iso_c_binding, only: c_double, c_double_complex, c_int, c_bool
implicit none
private
public init, run_trajectory_z0_rk4, free, &
      c_init, c_run_trajectory_z0_rk4, c_free

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

subroutine c_init(tLength, tSteps, depth, modes, hs_dim, g, gamma, Omega, h, &
         Lmap, with_terminator, populated_modes) bind(c)
   real(c_double), intent(in)            :: tLength
   integer(c_int), intent(in)            :: tSteps
   integer(c_int), intent(in)            :: depth
   integer(c_int), intent(in)            :: modes
   integer(c_int), intent(in)            :: hs_dim
   complex(c_double_complex), intent(in) :: g(modes)
   real(c_double), intent(in)            :: gamma(modes)
   real(c_double), intent(in)            :: Omega(modes)
   complex(c_double_complex), intent(in) :: h(hs_dim, hs_dim)
   integer(c_int), intent(in)            :: Lmap(modes)
   logical(c_bool), intent(in)           :: with_terminator
   integer(c_int), intent(in)            :: populated_modes
   call init(tLength, tSteps, depth, g, gamma, omega, h, Lmap, &
         logical(with_terminator), populated_modes)
end subroutine c_init

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


subroutine c_free() bind(c)
   implicit none
   call free()
end subroutine c_free

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


subroutine c_run_trajectory_z0_rk4(hs_dim, tSteps, psi0, psi) bind(c)
   implicit none
   integer(c_int), intent(in) :: hs_dim
   integer(c_int), intent(in) :: tSteps
   complex(c_double_complex), intent(in) :: psi0(hs_dim)
   complex(c_double_complex), intent(out) :: psi(tSteps, hs_dim)

   if ((hs_dim /= dim_) .or. (tSteps /= tSteps_)) then
      print *, 'ERROR: Wrong parameters in C-wrapper.'
      call Exit(1)
   end if
   psi = run_trajectory_z0_rk4(psi0)
end subroutine c_run_trajectory_z0_rk4

function run_trajectory_z0_rk4(psi0) result(psi)
   implicit none
   complex(dp), intent(in) :: psi0(dim_)
   complex(dp)             :: psi(tSteps_, dim_)
   !----------------------------------------------------------------------------
   !F2PY INTENT(IN) psi0, hs_dim, tSteps
   !F2PY INTENT(OUT) psi
   !----------------------------------------------------------------------------

   complex(dp) :: hierarchy(size_, RK_COPIES), dt
   integer :: t
   !TODO Error Check
   print *, "psi0 =", psi0

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

end module hierarchy
