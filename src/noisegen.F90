! Functions and classes related to generating colored noise
!
! ExponentialNoiseGenerator contains a generator for noise with correlation
! function
!
!          alpha(t) = sum_j  g_j * exp(-gamma_j * |t| - ii * Omega_j * t).
!
! In order to obtain good results, two conditions have to be met:
!     - alpha(0) should be real, otherwise the resulting correlation function
!     is distorted around t=0 with larger real part and zero imaginary part
!     - the propagation time should be enough for alpha to be decayed to
!     approximately zero.
!
! Usage:
!     type(ExponentialNoiseGenerator) :: ng
!     complex                         :: Z(tSteps)
!
!     call ng%init(dt, tSteps, g, gamma, Omega)
!     Z = ng%get_realization()
!     call ng%free()
!
! The Noise Generator is based on the Fourier-Filter algorithm described in
! details in Garcia-Ojalvo, Sancho: Noise in Spacially Extended Systems.

module noisegen
use system
use, intrinsic :: iso_c_binding
include 'fftw3.f03'
private

real(dp), parameter :: twopi = 8._dp * atan(1._dp)

type, public :: ExponentialNoiseGenerator
   integer :: N_
   ! Square root of the spectral density
   complex(C_DOUBLE_COMPLEX), allocatable :: sqrtJ_(:)
   ! Spaces allocated for the Fourier Transform routines
   complex(C_DOUBLE_COMPLEX), allocatable :: Z_(:)

   ! FFTW3 plans
   type(C_PTR) :: Zt_to_w_
   type(C_PTR) :: Zw_to_t_

   contains
   private
   procedure, public :: init
   procedure, public :: free
   procedure, public :: get_realization
end type ExponentialNoiseGenerator

public init_random_seed

contains

subroutine init(self, dt, tSteps, g, gamma, Omega)
   ! :dt: Time step size
   ! :tSteps: Number of time steps
   ! :g[*]: Coupling strenghts in bcf
   ! :gamma[*]: Dampings in bcf
   ! :Omega[*]: Center frequencies in bcf
   !
   ! Note: g, gamma, Omega should be of same size.

   implicit none
   class(ExponentialNoiseGenerator) :: self
   real(dp), intent(in)             :: dt
   integer, intent(in)              :: tSteps
   complex(dp), intent(in)          :: g(:)
   real(dp), intent(in)             :: gamma(:)
   real(dp), intent(in)             :: Omega(:)

   complex(C_DOUBLE_COMPLEX) :: alpha(2*tSteps), a
   type(C_PTR) :: plan
   integer :: t

   if (size(g) /= size(gamma) .or. (size(g) /= size(Omega))) then
      print *, "ERROR in ExponentialNoiseGenerator.init: Arrays do not match"
   endif

   self%N_ = tSteps
   ! Setup the Fourier-transformation plans used in noise creation
   ! FIXME Is this simd-valid?
   allocate(self%Z_(2*tSteps), self%sqrtJ_(2*tSteps))
   self%Zt_to_w_ = fftw_plan_dft_1d(2*tSteps, self%Z_, self%Z_, FFTW_FORWARD, &
         FFTW_ESTIMATE)
   self%Zw_to_t_ = fftw_plan_dft_1d(2*tSteps, self%Z_, self%Z_, FFTW_BACKWARD, &
         FFTW_ESTIMATE)

   ! Calculate the spectral density by Inverse Fourier transforming bcf
   ! FIXME Can we do that directly without alpha?
   do t = -tSteps, tSteps - 1
      a = sum(g * exp((-gamma - ii * Omega) * abs(t * dt)))
      if (t < 0) then
         alpha(2 * tSteps + t + 1) = conjg(a)
      else
         alpha(t + 1) = a
      end if
   end do

   plan = fftw_plan_dft_1d(2*tSteps, alpha, self%sqrtJ_, FFTW_FORWARD, &
         FFTW_ESTIMATE)
   call fftw_execute_dft(plan, alpha, self%sqrtJ_)
   call fftw_destroy_plan(plan)

   ! We actually need square root of spectral density (properly normalized)
   self%sqrtJ_ = sqrt(self%sqrtJ_) / (2*tSteps)
end subroutine init


subroutine free(self)
   implicit none
   class(ExponentialNoiseGenerator) :: self
   call array_free(self%Z_)
   call array_free(self%sqrtJ_)
   call fftw_destroy_plan(self%Zt_to_w_)
   call fftw_destroy_plan(self%Zw_to_t_)
end subroutine free


function get_realization(self) result(Z)
   ! Create a realization of the random process, call after initialization
   !
   ! :result[tSteps]: Realization of the noise process

   implicit none
   class(ExponentialNoiseGenerator) :: self
   complex(dp)                      :: Z(self%N_)

   call white_noise_fill(self%Z_)
   ! Calculate power spectrum of white noise
   call fftw_execute_dft(self%Zt_to_w_, self%Z_, self%Z_)
   ! Multiply in freq. domain = folding of white noise with sqrtJ in time domain
   self%Z_ = self%Z_ * self%sqrtJ_
   call fftw_execute_dft(self%Zw_to_t_, self%Z_, self%Z_)

   Z = self%Z_(1:self%N_)
end function get_realization

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                               HELPER FUNCTIONS                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine white_noise_fill(Z)
   ! Fill the array Z with "discrete white noise", i.e. indepent, complex
   ! Gaussian random variables with unit covariance in each component using the
   ! Box-Mueller algorithm.
   !
   ! :Z[*]: Array to be filled
   !
   ! Note: We don't write this as function, since space for Z is already
   ! allocated during initialization.

   implicit none
   complex(dp), intent(out) :: Z(:)

   integer :: t
   real(dp) :: x, y

   do t = 1, size(Z)
      call random_number(x)
      call random_number(y)
      Z(t) = cmplx( sqrt(-log(x)) * cos(twopi*y), sqrt(-log(x)) * sin(twopi*y),&
            dp)
   end do
end subroutine white_noise_fill


subroutine init_random_seed()
   ! Sets the random seed of the machine in a very "random" fashion in order
   ! to yield indepent results even in multicore calculations on the same
   ! machine.
   !
   ! Stolen from  http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html

#ifdef __INTEL_COMPILER
   use ifport
#endif
   implicit none
   integer, allocatable :: myseed(:)
   integer :: i, n, un, istat, dt(8), pid, t(2), s
   integer(8) :: count, tms

   call random_seed(size = n)
   allocate(myseed(n))
   ! First try if the OS provides a random number generator
   open(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
   if (istat == 0) then
      read(un) myseed
      close(un)
   else
      ! Fallback to XOR:ing the current time and pid. The PID is
      ! useful in case one launches multiple instances of the same
      ! program in parallel.
      call system_clock(count)
      if (count /= 0) then
         t = transfer(count, t)
      else
         call date_and_time(values=dt)
         tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24 * 60 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
         t = transfer(tms, t)
      end if
      s = ieor(t(1), t(2))
      pid = getpid() + 1099279 ! Add a prime
      s = ieor(s, pid)
      if (n >= 3) then
         myseed(1) = t(1) + 36269
         myseed(2) = t(2) + 72551
         myseed(3) = pid
         if (n > 3) then
            myseed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
         end if
      else
         myseed = s + 37 * (/ (i, i = 0, n - 1 ) /)
      end if
   end if
   call random_seed(put=myseed)
end subroutine init_random_seed

end module noisegen
