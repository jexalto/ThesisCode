module MOD_GVALUES_ODD
    use MOD_INTEGRALS, only: &
    Wjj_int, Woj_int, Wjo_int, Woo_int

    implicit none

    public :: Gjj_odd, Goj_odd, Gjo_odd, Goo_odd

    contains

    subroutine Gjj_odd(mu, ksi, eta, c, d, r, Gjj_out)
        ! Input: mu (tip speed ratio), ksi (location x/r), eta (location y/r), c-d (location vortex), r (radius)
        ! Output; Gjj (correction factor)
        real(8) :: mu, ksi, eta, c, d, pi, r
        real(8) :: Wjj_int_out, sum_jj, Gjj_out, start, end
        integer :: order, summation, m
            
        order = 1
        summation = 4
        pi = 3.141592653589793_8

        ! --- Wjj value ---
        call cpu_time(start)
        do m = 1, summation
            call Wjj_int(order, eta, ksi, mu, c, d, Wjj_int_out)
            sum_jj = sum_jj + (order)**2 * Wjj_int_out
            order = 2 * m + 1
        end do

        Gjj_out = 8/(r * eta * pi) * sum_jj
        call cpu_time(end)
        print*, 'Gjj value: ', Gjj_out
        print*, 'Time need to converge: ', (end-start), ' seconds'
    
    end subroutine Gjj_odd

    subroutine Goj_odd(mu, ksi, eta, c, d, r, Goj_out)
        ! Input: mu (tip speed ratio), ksi (location x/r), eta (location y/r), c-d (location vortex), r (radius)
        ! Output; Gjj (correction factor)
        real(8) :: mu, ksi, eta, c, d, pi, r
        real(8), intent(inout) :: Woj_int_out, sum_oj, Goj_out, start, end
        integer :: order, summation, m

        order = 1
        summation = 4
        pi = 3.141592653589793_8

        ! --- Woj value ---
        call cpu_time(start)
        do m = 1, summation
            call Woj_int(order, eta, ksi, mu, c, d, Woj_int_out)
            ! print*, Woj_int_out
            sum_oj= sum_oj + (order)**2 * Woj_int_out
            order = 2 * m + 1
        end do

        Goj_out = 8/(r * eta * pi) * sum_oj
        call cpu_time(end)
        print*, 'Goj value: ', Goj_out
        print*, 'Time need to converge: ', (end-start), ' seconds'
    
    end subroutine Goj_odd

    subroutine Gjo_odd(mu, ksi, eta, c, d, r, Gjo_out)
        ! Input: mu (tip speed ratio), ksi (location x/r), eta (location y/r), c-d (location vortex), r (radius)
        ! Output; Gjj (correction factor)
        real(8) :: mu, ksi, eta, c, d, pi, r
        real(8) :: Wjo_int_out, sum_jo, Gjo_out, start, end
        integer :: order, summation, m

        order = 1
        summation = 4
        pi = 3.141592653589793_8

        ! --- Woj value
        call cpu_time(start)
        do m = 1, summation
            call Wjo_int(order, eta, ksi, mu, c, d, Wjo_int_out)
            sum_jo = sum_jo + (order)**2 * Wjo_int_out
            order = 2 * m + 1
        end do

        Gjo_out = 8/(r * eta * pi) * sum_jo
        call cpu_time(end)
        print*, 'Gjo value: ', Gjo_out
        print*, 'Time need to converge: ', (end-start), ' seconds'

    end subroutine Gjo_odd

    subroutine Goo_odd(mu, ksi, eta, c, d, r, Goo_out)
        ! Input: mu (tip speed ratio), ksi (location x/r), eta (location y/r), c-d (location vortex), r (radius)
        ! Output; Gjj (correction factor)
        real(8) :: mu, ksi, eta, c, d, pi, r
        real(8) :: Woo_int_out, sum_oo, Goo_out, start, end
        integer :: order, summation, m

        order = 1
        summation = 4
        pi = 3.141592653589793_8

        ! --- Woo value ---
        call cpu_time(start)
        do m = 1, summation
            call Woo_int(order, eta, ksi, mu, c, d, Woo_int_out)
            sum_oo = sum_oo + (order)**2 * Woo_int_out
            order = 2 * m + 1
        end do

        Goo_out = 8/(r * eta * pi) * sum_oo
        call cpu_time(end)
        print*, 'Goo value: ', Goo_out
        print*, 'Time need to converge: ', (end-start), ' seconds'

    end subroutine Goo_odd
    
end module MOD_GVALUES_ODD