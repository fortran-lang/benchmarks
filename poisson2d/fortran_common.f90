!< Module for handling common things such as input arguments and timings
module common_m

  implicit none

  public
  integer, parameter :: sp = kind(0.)
  integer, parameter :: dp = kind(0.d0)

  real(dp), parameter :: a = 0.01_dp
  real(dp), parameter :: epsilon0 = 8.85E-12_dp

contains

  subroutine parse_args(N, N_ITER)
    integer, intent(out) :: N, N_ITER

    character(len=128) :: input
    integer :: length

    if ( command_argument_count() < 2 ) then
      stop 'Not enough arguments [N N_ITER]'
    end if

    call get_command_argument(1,input,length)
    read(input,*) N
    call get_command_argument(2,input,length)
    read(input,*) N_ITER

  end subroutine parse_args

  subroutine write_timing(N, iter, delta, timing)
    integer, intent(in) :: N, iter
    real(dp), intent(in) :: delta, timing

    write(*,'(2(i16,tr1),2(e22.15,tr1))') N, iter, delta, timing/real(iter,dp)

  end subroutine write_timing

end module common_m
  

  
