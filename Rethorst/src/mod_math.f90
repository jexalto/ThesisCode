module MOD_MATH
    use iso_fortran_env
    
    implicit none

    public :: factorial, gamma

    contains

    subroutine factorial(fac, facout)

        integer, parameter :: LargeInt_K = selected_int_kind (32)
        integer(kind=LargeInt_K), intent(inout) :: facout, fac
        ! integer, parameter :: MyLongIntType = selected_int_kind (12) (kind=MyLongIntType) 
        integer(kind=LargeInt_K) :: ret
        integer :: m
        
        ret = 1
        
        if(fac==0)then
            facout = 1
        else
            do m = 1, fac
                ret = ret * m
            end do
            facout = ret
        endif

    end subroutine factorial

    subroutine gamma(z, gam)

        integer :: discretisation, i
        real(8) :: dx, x
        real(8), intent(inout) :: z
        real(8), intent(inout) :: gam

        x = 0.0001
        discretisation = 100
        dx = 9./discretisation
        gam = 0.
        do i = 1, discretisation
            gam = gam + x**(z-1) * exp(-x) * dx
            x = x + dx
        end do
        
    end subroutine gamma

end module MOD_MATH