module MOD_MODIFIEDBESSEL_LST
    use MOD_MATH, only: &
        factorial
    
    implicit none

    public :: firstkind_bessel, secondkind_bessel

    contains

    subroutine firstkind_bessel(order, x, fk_bessel)

        real, intent(inout), DIMENSION(1:50) :: x, fk_bessel
        real :: addition
        integer, intent(in) :: order
        integer :: m
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
        
        real, intent(inout) :: x, sk_bessel
        real :: pi, t, dt
        integer, intent(in) :: order
        integer :: m

        pi = 3.14159265
        sk_bessel = 0.

        t = 0.05
        dt = 0.1

        do m = 0, 50
            sk_bessel = sk_bessel + (EXP(-x * cosh(t)) * cosh(order*t))*dt
            t = t + dt
        end do   

    end subroutine secondkind_bessel

end module MOD_MODIFIEDBESSEL_LST