program Gmatrix
    use MOD_INTEGRALS, only: &
    Iv_integral, Kv_integral, Wjj_int, Woj_int, Wjo_int, Woo_int

    use MOD_MODIFIEDBESSEL, only:&
    firstkind_bessel, secondkind_bessel, secondkind_bessel_der

    use MOD_MATH, only:&
    factorial

    implicit none

    real(8) :: ksi, eta, mu, c, d, pi, r, lambda, sum_jj, Wjj_int_out, Gjj
    ! real(8) :: Wjj_int_out, Woj_int_out, Wjo_int_out, Woo_int_out
    ! real(8) :: sum_jj, sum_jo, sum_oj, sum_oo
    ! real(8) :: Gjj, Gjo, Goj, Goo
    ! real(8) :: int_Iv, fk_bessel
    integer :: order, summation, m
    integer, parameter :: LargeInt_K = selected_int_kind (32)
    integer(kind=LargeInt_K) :: output, gammainput
    real(8) :: end, start, Kv, Iv, Kv_p, Iv_etalambda

    call cpu_time(start)
    c = 1.2
    d = 1.3
    mu = 0.9
    ksi = 0.5
    eta = 0.1
    r = 1.0
    order = 3
    lambda = 4.8

    pi = 3.1415926
    sum_jj = 0.
    summation = 5
    gammainput = 5

    ! call cpu_time(start)
    ! call secondkind_bessel(order, lambda, Kv)            ! Kv
    ! call cpu_time(end)
    ! ! print*, 'Second kind bessel value: ', Kv
    ! print*, 'Second kind bessel: ', (end-start), ' seconds'

    ! call cpu_time(start)
    ! call factorial(gammainput, output)
    ! call cpu_time(end)
    ! print*, 'Factorial function: ', (end-start), ' seconds'

    ! ! --- Get nested integral ---
    ! call cpu_time(start)
    ! call Iv_integral(order, lambda, c, d, int_Iv)
    ! call cpu_time(end)
    ! print*, 'Iv integral: ', (end-start), ' seconds'

    call cpu_time(start)
    call firstkind_bessel(order, lambda, Iv)             ! Iv
    call cpu_time(end)
    ! print*, 'First kind bessel value: ', Iv
    print*, 'First kind bessel: ', (end-start), ' seconds'

    ! call cpu_time(start)
    ! call secondkind_bessel_der(order, lambda, Kv_p)      ! Kv_prime
    ! call cpu_time(end)
    ! print*, 'Second kind bessel derivative: ', (end-start), ' seconds'

    call cpu_time(start)
    do m = 1, summation
        call Wjj_int(order, eta, ksi, mu, c, d, Wjj_int_out)
        sum_jj = sum_jj + (order)**2 * Wjj_int_out
        ! print*, 'Order, Wjj ', order, Wjj_int_out * (order)**2
        order = 2 * m + 1
    end do

    Gjj = 8/(r * eta * pi) * sum_jj
    call cpu_time(end)
    print*, '--- FORTRAN ---'
    print*, 'Gjj value: ', Gjj
    print*, 'Time need to converge: ', (end-start), ' seconds'
    
end program Gmatrix