module rho_interface
  !! Define a scalar function over 2-space.
  implicit none

  private
  public :: rho, dp

  integer, parameter :: precision = 15, range = 307
  integer, parameter :: dp = selected_real_kind(precision, range)

  interface 
    pure real(dp) module function rho(x,y)
      !! Poisson equation inhomogeneous term.
      implicit none
      real(dp), intent(in) :: x,y
    end function
  end interface

end module 

submodule(rho_interface) rho_definition
  implicit none
contains
  module procedure rho
      !! To Do: give meaningful names to the magic numbers.
      !! See https://en.wikipedia.org/wiki/Magic_number_(programming)#Unnamed_numerical_constants.
      if (all([x,y]>0.6_dp .and. [x,y]<0.8_dp)) then
        rho = 1._dp
      else
        rho = merge(-1., 0., all([x,y]>0.2_dp .and. [x,y]<0.4_dp))
      end if
  end procedure
end submodule

program poisson
  !! Solve the 2D Poisson equation using a smoothing operation to iterate.
  use rho_interface, only: rho, dp
  implicit none

  integer i,j
  integer, parameter :: M=150
  real(dp), parameter :: dx=0.01
  real(dp) :: t_start, rho_sampled(M,M)
  
  call cpu_time(t_start)

  associate( dy => (dx) ) ! Associating with an expression provides immutable state so dy cannot be inadvertently redefined.

    do concurrent(i=1:M, j=1:M)
      rho_sampled(i,j) = rho(i*dx,j*dy)
    end do
 
    block ! Tighten the scoping to declutter the code above.
      real(dp) :: delta_phi, t_end
      real(dp), parameter :: epsilon0=8.85E-12_dp, tolerance=1E-6_dp
      real(dp), dimension(M,M) :: phi_prime, phi
      integer iteration

      phi = 0.

      phi_prime([1,M], 2:M-1) = phi([1,M], 2:M-1) ! Initialize only boundary values except corners. (Internal values will
      phi_prime(2:M-1, [1,M]) = phi(2:M-1, [1,M]) ! be overwritten in the first iteration. Corners will never be used.)

      delta_phi = tolerance + epsilon(tolerance) ! Ensure at least 1 iteration.
      iteration = 0
      do while (delta_phi > tolerance )
        iteration = iteration + 1
        do concurrent(i=2:M-1, j=2:M-1) ! Compute updated solution estimate at internal points.
          phi_prime(i,j) =   (phi(i+1,j) + phi(i-1,j) + phi(i,j+1) + phi(i,j-1))/4._dp &
                           + (dx/2._dp)*(dy/2._dp)/epsilon0*rho_sampled(i,j)
        end do
        delta_phi = maxval(abs(phi_prime - phi))
        phi(2:M-1, 2:M-1) = phi_prime(2:M-1, 2:M-1) ! Update internal values.
      end do

      call cpu_time(t_end)
      print *, t_end-t_start, iteration

    end block

  end associate

end program
