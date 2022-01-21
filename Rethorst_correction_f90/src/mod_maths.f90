module MDO_MATHS
    use iso_fortran_env
    
    implicit none

    public :: factorial

    contains

    subroutine factorial(fac, facout)

        integer, parameter :: LargeInt_K = selected_int_kind (32)
        integer(kind=LargeInt_K), intent(inout) :: facout
        ! integer, parameter :: MyLongIntType = selected_int_kind (12) (kind=MyLongIntType) 
        integer(kind=LargeInt_K) :: ret
        integer :: m, fac
        
        ret = 1

        if(fac==0)then
            facout = 1
        else
            do m = 1, fac
                ret = ret * m
                ! if (ret<0)then
                !     print*,'negative ret: ', m
                !     exit
                ! endif
            end do
            facout = ret
        endif

        ! deallocate(ret)
        ! deallocate(m)
        ! deallocate(fac)

    end subroutine factorial

end module MDO_MATHS