module rhofunc
implicit none
public
integer, parameter :: dp=kind(0.d0)
contains
    pure real(dp) function rho(x,y)
        real(dp), intent(in) :: x,y
        if (x > 0.6_dp .and. x < 0.8_dp .and. y > 0.6_dp .and. y<0.8_dp) then
            rho = 1.0_dp
        else if (x> 0.2_dp .and. x<0.4_dp .and. y>0.2_dp .and. y<0.4_dp) then
            rho = -1.0_dp
        else
            rho = 0.0_dp
        end if
    end function

end module

program poisson
use rhofunc, only: rho
implicit none
integer, parameter :: dp=kind(0.d0), M=300
integer            :: i,j, iter
real(dp),parameter :: epsilon0=8.85E-12_dp, target=1E-6_dp, a=0.01_dp
real(dp)           :: delta, b, e, phiprime(M,M), phi(M,M), a2, rhoarr(M,M), temp(M,M)


delta = 1.0_dp
iter = 0
call cpu_time(b)
phiprime(:,:) = 0.0_dp
phi(:,:) = 0.0_dp

do while (delta > target )
    iter = iter + 1
    a2 = a**2.0_dp
    do i=2, M-1
        do j=2, M-1
            phiprime(i,j) = (phi(i+1,j) + phi(i-1,j) + phi(i,j+1) + phi(i,j-1))/4.0_dp &
            + a2/4.0_dp/epsilon0*rho(i*a,j*a)
        end do
    end do
    delta = maxval(abs(phiprime - phi))
    temp = phi
    phi = phiprime 
    phiprime = temp

end do
call cpu_time(e)
print *, e-b, iter
end program
