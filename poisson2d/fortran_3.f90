module poisson_m
  use common_m
  implicit none

contains

  pure real(dp) function rho(x,y)
    real(dp), intent(in) :: x,y
    if (x > 0.6_dp .and. x < 0.8_dp .and. y > 0.6_dp .and. y<0.8) then
      rho = 1.0_dp
    else if (x> 0.2_dp .and. x<0.4_dp .and. y>0.2_dp .and. y<0.4_dp) then
      rho = -1.0_dp
    else
      rho = 0.0_dp
    end if
  end function rho

  subroutine run(M, N_ITER)
    integer, intent(in) :: M, N_ITER

    integer :: i,j, iter
    real(dp) :: tic, toc
    real(dp) :: delta, a2
    real(dp) :: phiprime(M,M), phi(M,M), rhoa(M,M), temp(M,M)

    ! Since arrays are allocated on the stack
    call cpu_time(tic)
    
    phiprime(:,:) = 0.0_dp
    phi(:,:) = 0.0_dp
    do j=2, M-1
      do i=2, M-1
        rhoa(i,j) = rho((i-1)*a,(j-1)*a)
      end do
    end do

    a2 = a**2

    iter = 0
    do while (iter < N_ITER )
      
      iter = iter + 1
      do j=2, M-1
        do i=2, M-1
          phiprime(i,j) = (phi(i+1,j) + phi(i-1,j) + phi(i,j+1) + phi(i,j-1) + &
              a2 / epsilon0 * rhoa(i,j))*0.25_dp
        end do
      end do

      delta = maxval(abs(phiprime - phi))
      temp = phi
      phi = phiprime 
      phiprime = temp
      
    end do
    call cpu_time(toc)

    call write_timing(M, iter, delta, toc - tic)

  end subroutine run

end module poisson_m

program poisson
  use common_m, only : parse_args
  use poisson_m, only: run
  implicit none
  integer :: M, N_ITER

  call parse_args(M, N_ITER)

  call run(M, N_iter)

end program poisson
