program poisson
  use common_m

  implicit none

  integer :: M, N_ITER
  integer :: i,j, iter
  real(dp) :: tic, toc
  real(dp) :: delta, a2

  real(dp), allocatable :: rhoa(:,:)
  real(dp), allocatable, target :: phiprime(:,:), phi(:,:)
  real(dp), pointer :: p(:,:), pnew(:,:), tmp(:,:)

  call parse_args(M, N_ITER)
  
  allocate(phi(M,M))
  allocate(phiprime(M,M))
  allocate(rhoa(M,M))

  ! we don't time allocation/deallocation
  call cpu_time(tic)

  ! Fortran doesn't care too much about pow
  ! since it figures out if the exponent is an integer or not
  a2 = a*a / epsilon0
  do j=2, M-1
    do i=2, M-1
      rhoa(i,j) = rho((i-1)*a,(j-1)*a) * a2
    end do
  end do

  phi(:,:) = 0.0_dp
  phiprime(:,:) = 0.0_dp
  p => phi
  pnew => phiprime
  iter = 0
  do while ( iter < N_ITER )
    iter = iter + 1

    call iterate(rhoa, p, pnew, delta)

    ! Easer to swap pointer than copying stuff around
    ! In fortran we could even use pointers
    ! But we can just call
    tmp => p
    p => pnew
    pnew => tmp

  end do
  
  call cpu_time(toc)

  deallocate(phi, phiprime, rhoa)

  call write_timing(M, iter, delta, toc - tic)

contains

  subroutine iterate(rhoa, p, pnew, delta)
    real(dp), intent(in) :: p(:,:), rhoa(:,:)
    real(dp), intent(out) :: pnew(:,:)
    real(dp), intent(out) :: delta

    integer :: i, j

    delta = 0._dp
    do j=2, M-1
      do i=2, M-1
        pnew(i,j) = (p(i+1,j) + p(i-1,j) + p(i,j+1) + p(i,j-1) + rhoa(i,j)) * 0.25_dp
        delta = max(delta, abs(pnew(i,j) - p(i,j)))
      end do
    end do

  end subroutine iterate

  pure real(dp) function rho(x,y)
    real(dp), intent(in) :: x,y
    if ( 0.6_dp < x .and. x < 0.8_dp .and. 0.6_dp < y .and. y < 0.8_dp ) then
      rho = 1.0_dp
    else if ( 0.2_dp < x .and. x < 0.4_dp .and. 0.2_dp < y .and. y < 0.4_dp ) then
      rho = -1.0_dp
    else
      rho = 0.0_dp
    end if
  end function rho

end program poisson
