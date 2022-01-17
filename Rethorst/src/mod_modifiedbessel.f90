module MOD_MODIFIEDBESSEL
    use MDO_MATHS, only: &
        factorial
    
    implicit none

    public :: firstkind_bessel, secondkind_bessel

    contains

    subroutine firstkind_bessel(order, x, fk_bessel)

        real, intent(inout) :: x, fk_bessel
        real :: bessel, addition
        integer, intent(in) :: order
        integer :: m, m_fac, m_fac2
        integer, parameter :: LargeInt_K = selected_int_kind (32)
        integer(kind=LargeInt_K) :: facout, facout2

        fk_bessel = 0.
        facout = 1
        facout2 = 1

        do m = 0, 10
            if(m==0)then
                facout = 1
            else
                facout = facout * m
            endif
            facout2 = facout2* (m + order)
            addition = (x/2)**(2*m + order)/(facout * facout2)
            fk_bessel = fk_bessel + addition
        end do

        ! deallocate(facout)
        ! deallocate(facout2)
        ! deallocate(bessel)
        ! deallocate(addition)
        ! deallocate(m)
        ! deallocate(m_fac2)

    end subroutine firstkind_bessel

    subroutine secondkind_bessel(order, x, sk_bessel)
        
        real, intent(inout) :: x, Iv, Iv_min
        real :: bessel, addition, pi
        integer, intent(in) :: order
        integer :: m, m_fac, m_fac2
        integer, parameter :: LargeInt_K = selected_int_kind (32)
        integer(kind=LargeInt_K) :: facout, facout2

        pi = 3.141592653589793

        do m = 0, 10
            if(m==0)then
                facout = 1
            else
                facout = facout * m
            endif
            facout2 = facout2* (m + order)
            addition = (x/2)**(2*m + order)/(facout * facout2)
            Iv = Iv + addition
        end do

        do m = 0, 10
            if(m==0)then
                facout = 1
            else
                facout = facout * m
            endif
            facout2 = facout2* (m - order)
            addition = (x/2)**(2*m - order)/(facout * facout2)
            Iv_min = Iv_min + addition
        end do

        sk_bessel = pi/2 * (Iv_min - Iv)/sin(order*pi)

    end subroutine secondkind_bessel

end module MOD_MODIFIEDBESSEL