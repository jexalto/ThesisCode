module MOD_MODIFIEDBESSEL
    use MOD_MATH, only: &
        gamma, factorial
    
    implicit none

    public :: firstkind_bessel, secondkind_bessel, firstkind_bessel_der, secondkind_bessel_der

    contains

    subroutine firstkind_bessel(order, x, fk_bessel)

        real(8), intent(inout) :: fk_bessel, x
        real(8) :: addition, start, end
        integer, intent(inout) :: order
        integer :: m
        integer, parameter :: LargeInt_K = selected_int_kind (32)
        integer(kind=LargeInt_K) :: facout, facout2, gammainput

        fk_bessel = 0.
    
        facout = 1
        facout2 = 1

        do m = 0, 7
            if(m==0)then
                facout = 1
            else
                facout = facout * m
            endif

            gammainput = order + m
            call factorial(gammainput, facout2)
            addition = (x/2)**(2*m + order) / (facout * facout2)

            fk_bessel = fk_bessel + addition

        end do

    end subroutine firstkind_bessel

    subroutine secondkind_bessel(order, x, sk_bessel)
        
        real(8), intent(inout) :: x, sk_bessel
        real(8) :: pi, t, dt
        integer, intent(inout) :: order
        integer :: m, discretisation

        pi = 3.141592653_8
        sk_bessel = 0._8
        discretisation = 1000_8
        t = 0.001_8
        dt = 10./discretisation

        do m = 0, discretisation
            sk_bessel = sk_bessel + exp(-x*cosh(t)) * cosh(order*t) * dt

            t = t + dt
        end do

    end subroutine secondkind_bessel

    subroutine firstkind_bessel_der(order, x, fk_bessel_der)

        ! Derivatives found on: https://math.stackexchange.com/questions/2846402/derivative-of-modified-bessel-function-of-second-kind

        real(8) :: x, fk_bessel_der, fk_bessel_1, fk_bessel_2, input
        integer :: order, order_low, order_high

        order_low = order - 1
        ordeR_high = order + 1

        input = x

        call firstkind_bessel(order, input, fk_bessel_1)
        call firstkind_bessel(order, input, fk_bessel_2)
        
        fk_bessel_der = 0.5 * (fk_bessel_1 + fk_bessel_2)

    end subroutine firstkind_bessel_der

    subroutine secondkind_bessel_der(order, x, sk_bessel_der)

        ! Derivatives found on: https://math.stackexchange.com/questions/2846402/derivative-of-modified-bessel-function-of-second-kind

        real(8) :: x, sk_bessel_der, sk_bessel_1, sk_bessel_2, input
        integer :: order, order_low, order_high

        order_low = order - 1
        ordeR_high = order + 1

        input = x

        call secondkind_bessel(order_low, input, sk_bessel_1)
        call secondkind_bessel(order_high, input, sk_bessel_2)
        
        sk_bessel_der = -0.5 * (sk_bessel_1 + sk_bessel_2)

    end subroutine secondkind_bessel_der

end module MOD_MODIFIEDBESSEL