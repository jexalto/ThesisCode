program Gmatrix
    use MOD_INTEGRALS, only: &
    Iv_integral, Kv_integral, Wjj_int, Woj_int, Wjo_int, Woo_int

    use MOD_MODIFIEDBESSEL, only:&
    firstkind_bessel, secondkind_bessel

    implicit none

    real(8) :: ksi, eta, mu, c, d, pi, r, lambda, output
    real(8) :: Wjj_int_out, Woj_int_out, Wjo_int_out, Woo_int_out
    real(8) :: sum_jj, sum_jo, sum_oj, sum_oo
    real(8) :: Gjj, Gjo, Goj, Goo
    integer :: order, summation, m
    real(8) :: end, start

    call cpu_time(start)
    c = 1.2
    d = 1.3
    mu = 0.9
    ksi = 0.5
    eta = 0.1
    r = 1.0
    order = 1
    lambda = 0.5

    pi = 3.1415926
    sum_jj = 0.
    summation = 10

    do m = 1, summation
        call Wjj_int(order, eta, ksi, mu, c, d, Wjj_int_out)
        sum_jj = sum_jj + (order)**2 * Wjj_int_out

        order = 2 * m + 1
    end do

    Gjj = 8/(r * eta * pi) * sum_jj
    call cpu_time(end)
    print*, '--- FORTRAN ---'
    print*, 'Gjj value: ', Gjj
    print*, 'Time need to converge: ', (end-start), ' seconds'
    
end program Gmatrix