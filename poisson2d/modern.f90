module modern_func
  !! Code by Damian Rouson, @rouson.
  !! Modified by ARC, @arunningcroc

  !! Define a scalar function over 2-space.
  implicit none

  private
  public :: dp, run_modern

  integer, parameter :: precision = 15, range = 50
  integer, parameter :: dp = selected_real_kind(precision, range)
contains
  pure real(dp) function rho(x,y)
    real(dp), intent(in) :: x,y
    rho = -(6.0_dp*x*y*(1.0_dp-y)-2.0_dp*x**3.0_dp)
  end function

  pure elemental real(dp) function boundary(y)
    real(dp),intent(in)   :: y
    boundary = y*(1-y)
  end function

  pure real(dp) function analytical_sol(x,y)
    real(dp), intent(in) :: x,y
    analytical_sol = y*(1-y)*x**3
  end function

  integer function run_modern(M)
    real(dp)    :: dx, rho_sampled(M,M)
    integer     :: M,i,j

    dx = 1.0/(M-1)
    associate( dy => (dx) ) ! Associating with an expression provides immutable state so dy cannot be inadvertently redefined.

      do concurrent(i=1:M, j=1:M)
        rho_sampled(i,j) = rho(i*dx,j*dy)
      end do
   
      block ! Tighten the scoping to declutter the code above.
        real(dp)                         :: delta_phi, avgerror
        real(dp), parameter              ::  tolerance=1E-8_dp, analytical_tol = 5.0e-3_dp
        real(dp), dimension(0:M-1,0:M-1) :: phi_prime, phi, analytical
        integer                          :: iteration
  
        phi = 0.
        phi(M-1,:) = boundary([(i*dx,i=0,M-1)])

        do concurrent(i=0:M-1, j=0:M-1)
          analytical(i,j) = analytical_sol(i*dx,j*dx)
        end do

        phi_prime([0,M-1], 1:M-2) = phi([0,M-1], 1:M-2) ! Initialize only boundary values except corners. (Internal values will
        phi_prime(1:M-2, [0,M-1]) = phi(1:M-2, [0,M-1]) ! be overwritten in the first iteration. Corners will never be used.)
        
        delta_phi = 1.0_dp !tolerance + epsilon(tolerance) ! Ensure at least 1 iteration.
        iteration = 0
        do while (delta_phi > tolerance )
          iteration = iteration + 1
          do concurrent(i=1:M-2, j=1:M-2) ! Compute updated solution estimate at internal points.
            phi_prime(i,j) =   (phi(i+1,j) + phi(i-1,j) + phi(i,j+1) + phi(i,j-1))/4._dp &
                             + (dx/2._dp)*(dy/2._dp)*rho_sampled(i,j)
          end do
          delta_phi = maxval(abs(phi_prime - phi))
          phi(1:M-2, 1:M-2) = phi_prime(1:M-2, 1:M-2) ! Update internal values.
          !print *, iteration
        end do
        phi = phi/maxval(abs(phi))
        analytical = analytical/maxval(abs(analytical))
        avgerror = sum(abs(phi-analytical))/(M*M)
        if(avgerror > analytical_tol) then
            print *, "Failed to match analytical solution", avgerror
        end if
        run_modern = iteration
      end block
      

    end associate
  end function
end module 


program poisson
  !! Solve the 2D Poisson equation using a smoothing operation to iterate.
  use modern_func, only: run_modern, dp
  implicit none

  integer             :: iter
  integer, parameter  :: M=170
  real(dp)            :: t_start, t_end

  
  call cpu_time(t_start)
  iter = run_modern(M)
  call cpu_time(t_end)
  print *, t_end-t_start, iter
end program