module MOD_INTEGRALS
    ! call cpu_time(start)
    ! call cpu_time(finish)
    ! print*, "Time = ", finish-start, "seconds"
    ! real(8) :: start, finish
    
    use MOD_MODIFIEDBESSEL, only: &
        firstkind_bessel, secondkind_bessel, firstkind_bessel_der, secondkind_bessel_der

    implicit none

    public :: Iv_integral, Kv_integral, Wjj_int, Woj_int, Wjo_int, Woo_int

    contains

    subroutine Iv_integral(order, lambda, c, d, int_Iv)

        real(8) :: dlambda, maxlambda, startlambda, lambda_run, fk_bessel     ! Running variables
        real(8), intent(out) :: int_Iv                                             ! Integral quantities
        real(8), intent(inout) :: lambda, c, d
        integer :: order, m, integration

        integration = 200
        maxlambda = lambda*d
        startlambda = lambda*c
        dlambda = (maxlambda - startlambda)/integration

        ! --- Iv_lambdaB integral
        lambda_run = startlambda
        fk_bessel = 0.
        int_Iv = 0.

        do m = 1, integration
            call firstkind_bessel(order, lambda_run, fk_bessel)
            int_Iv = int_Iv + fk_bessel/lambda_run * dlambda
            lambda_run = lambda_run + dlambda
        end do
    
    end subroutine Iv_integral


    subroutine Kv_integral(order, lambda, c, d, int_Kv)

        real(8) :: dlambda, lambda, maxlambda, startlambda, lambda_run    ! Running variables
        real(8) :: sk_bessel, int_Kv                                      ! Integral quantities
        real(8) :: c, d
        integer :: order, m, integration

        integration = 25
        maxlambda = lambda*d
        startlambda = lambda*c
        dlambda = (maxlambda - startlambda)/integration

        ! --- Kv_lambdaB integral
        lambda_run = startlambda
        sk_bessel = 0.
        int_Kv = 0.

        do m = 1, integration
            call secondkind_bessel(order, lambda_run, sk_bessel)
            int_Kv = int_Kv + sk_bessel/lambda_run * dlambda
            lambda_run = lambda_run + dlambda
        end do
    
    end subroutine Kv_integral

    subroutine Wjj_int(order, eta, ksi, mu, c, d, BIG_int)
        ! Wjj, vortex INSIDE jet
        real(8), intent(inout) :: eta, ksi, mu, c, d, BIG_int
        real(8) :: Iv, Kv, Kv_p, Iv_etalambda, int_Iv, etalambda
        real(8) :: lambda, dlambda, multiplier
        integer :: integration, m, order

        integration = 30

        lambda = 0.001
        BIG_int = 0.
        dlambda = (5.-lambda)/integration
        etalambda = lambda * eta

        do m = 1, integration-1
            ! --- Get modified bessel functions ---
            call secondkind_bessel(order, lambda, Kv)            ! Kv
            call firstkind_bessel(order, lambda, Iv)             ! Iv
            call secondkind_bessel_der(order, lambda, Kv_p)      ! Kv_prime
            call firstkind_bessel(order, etalambda, Iv_etalambda)  ! Iv(mu*lambda)
            ! --- Get nested integral ---
            call Iv_integral(order, lambda, c, d, int_Iv)

            multiplier = (Kv * Kv_p * Iv_etalambda)/(1/(lambda * (1/mu**2 - 1)) - Iv * Kv_p)
            
            BIG_int = BIG_int + ( multiplier * sin(ksi*lambda)/lambda * int_Iv ) * dlambda

            lambda = lambda + dlambda
            etalambda = lambda * eta
        end do

    end subroutine Wjj_int

    subroutine Woj_int(order, eta, ksi, mu, c, d, BIG_int)
        ! Woj, vortex INSIDE jet
        real(8), intent(inout) :: eta, ksi, mu, c, d
        real(8) :: Iv, Kv, Kv_p, Kv_mulambda, int_Iv, mulambda
        real(8) :: lambda, dlambda, BIG_int, multiplier
        integer :: integration, m, order

        integration = 50

        lambda = 0.001
        BIG_int = 0.
        dlambda = 5./integration
        mulambda = lambda * mu

        do m = 1, integration
            ! --- Get modified bessel functions ---
            call secondkind_bessel(order, lambda, Kv)            ! Kv
            call firstkind_bessel(order, lambda, Iv)             ! Iv
            call secondkind_bessel_der(order, lambda, Kv_p)      ! Kv_prime
            call secondkind_bessel(order, lambda, Kv_mulambda)   !  Iv(mu*lambda)
            ! --- Get nested integral ---
            call Iv_integral(order, lambda, c, d, int_Iv)

            multiplier = (1 / (mu - lambda*(1/mu - mu) * Iv * Kv_p) - 1) 

            BIG_int = BIG_int + (multiplier * (Kv_mulambda * sin(ksi*lambda))/lambda * int_Iv ) * dlambda

            lambda = lambda + dlambda
            mulambda = lambda * mu
        end do

    end subroutine Woj_int

    subroutine Wjo_int(order, eta, ksi, mu, c, d, BIG_int)
        ! Wjo, vortex OUTSIDE jet
        real(8), intent(inout) :: eta, ksi, mu, c, d
        real(8) :: Iv, Kv, Kv_p, Iv_mulambda, int_Kv, mulambda
        real(8) :: lambda, dlambda, BIG_int, multiplier
        integer :: integration, m, order

        integration = 50

        lambda = 0.001
        BIG_int = 0.
        dlambda = 5./integration
        mulambda = lambda * mu

        do m = 1, integration
            ! --- Get modified bessel functions ---
            call secondkind_bessel(order, lambda, Kv)            ! Kv
            call firstkind_bessel(order, lambda, Iv)             ! Iv
            call firstkind_bessel_der(order, lambda, Kv_p)       ! Iv_prime
            call firstkind_bessel(order, lambda, Iv_mulambda)    ! Iv(mu*lambda)
            ! --- Get nested integral ---
            call Kv_integral(order, lambda, c, d, int_Kv)

            multiplier = 1 / (mu - lambda*(1/mu - mu) * Iv * Kv_p) - 1

            BIG_int = BIG_int + multiplier * (Iv_mulambda * sin(ksi*lambda))/lambda * int_Kv * dlambda

            lambda = lambda + dlambda
            mulambda = lambda * mu
        end do

    end subroutine Wjo_int

    subroutine Woo_int(order, eta, ksi, mu, c, d, BIG_int)
        ! Woo, vortex OUTSIDE jet
        real(8), intent(inout) :: eta, ksi, mu, c, d
        real(8) :: Iv, Kv, Kv_p, Kv_mulambda, int_Kv, Iv_p, mulambda
        real(8) :: lambda, dlambda, BIG_int, multiplier
        integer :: integration, m, order

        integration = 50

        lambda = 0.001
        BIG_int = 0.
        dlambda = 5./integration
        mulambda = lambda * mu

        do m = 1, integration
            ! --- Get modified bessel functions ---
            call secondkind_bessel(order, lambda, Kv)            ! Kv
            call firstkind_bessel(order, lambda, Iv)             ! Iv
            call secondkind_bessel_der(order, lambda, Kv_p)      ! Kv_prime
            call firstkind_bessel_der(order, lambda, Iv_p)       ! Kv_prime
            call secondkind_bessel(order, lambda, Kv_mulambda)   ! Kv(mu*lambda)
            ! --- Get nested integral ---
            call Kv_integral(order, lambda, c, d, int_Kv)

            multiplier = (Iv * Iv_p * Kv_mulambda)/(1/(lambda * (1/mu**2 - 1)) - Iv * Kv_p)

            BIG_int = BIG_int + multiplier * sin(ksi*lambda)/lambda * int_Kv * dlambda

            lambda = lambda + dlambda
            mulambda = lambda * mu
        end do

    end subroutine Woo_int

end module MOD_INTEGRALS