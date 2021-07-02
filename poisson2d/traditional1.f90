module naiveimpl
implicit none
public
integer, parameter :: dp=kind(0.d0)
contains

pure real(dp) function rho(x,y)
    real(dp), intent(in) :: x,y
    rho = -(6.0_dp*x*y*(1.0_dp-y)-2.0_dp*x**3.0_dp)
end function

pure elemental real(dp) function boundary(y)
    real(dp), intent(in) :: y

    boundary = y*(1-y)
end function

pure real(dp) function analytical_sol(x,y)
        real(dp), intent(in) :: x,y
        analytical_sol = y*(1-y)*x**3
end function

integer function run_naive(M)
    integer, intent(in):: M
    integer, parameter :: dp=kind(0.d0)
    real(dp)           :: phiprime(0:M-1,0:M-1), phi(0:M-1,0:M-1), a2, analytical(0:M-1,0:M-1), temp(0:M-1,0:M-1),delta,a, avgerror
    real(dp),parameter :: target=1E-8_dp, analytical_tol = 5.0e-3
    integer            :: i,j, iter
    a = 1.0/(M-1)
    delta = 1.0_dp
    iter = 0
    phi = 0
    phiprime = 0
    phiprime(M-1,:) = boundary([(i*a,i=0,M-1)])
    phi(M-1,:) = boundary([(i*a,i=0,M-1)])
    do i=0,M-1
        do j=0, M-1
            analytical(i,j) = analytical_sol(i*a,j*a)
        end do
    end do
    do while (delta > target )
        iter = iter + 1
        a2 = a**2.0_dp
        do i=1, M-2
            do j=1, M-2
                phiprime(i,j) = (phi(i+1,j) + phi(i-1,j) + phi(i,j+1) + phi(i,j-1))/4.0_dp &
                + a2*rho(i*a,j*a)/4.0_dp
            end do
        end do
        delta = maxval(abs(phiprime - phi))
        temp = phi
        phi = phiprime 
        phiprime = temp
    
    end do
    phi = phi/maxval(abs(phi))
    analytical = analytical/maxval(abs(analytical))
    avgerror = sum(abs(phi-analytical))/(M*M)
    print *, avgerror
    if(avgerror > analytical_tol) then
        print *, "Failed to agree with analytical solution", avgerror
    end if
    run_naive = iter
end function
end module naiveimpl


program poisson
use naiveimpl, only: run_naive
implicit none
integer, parameter :: dp=kind(0.d0), M=100 
integer            :: iter
real(dp)           :: b,e


call cpu_time(b)
iter = run_naive(M)
call cpu_time(e)
print *, e-b, iter
end program