program Gmatrix_lst
    ! call cpu_time(start)
    ! call cpu_time(finish)
    ! print*, "Time = ", finish-start, "seconds"
    ! real :: start, finish
    
    use MOD_MODIFIEDBESSEL_LST, only: &
        firstkind_bessel, secondkind_bessel

    implicit none

    real :: dlambda, lambda, lambdab                        ! Running variables
    real :: fk_bessel, sk_bessel, int_Iv, int_Kv            ! Integral quantities
    integer :: order, m
    real :: start, finish

    order = 1
    dlambda = 0.1

    print*, '--- FORTRAN ---'
    call cpu_time(start)
    ! Iv_lambdaB integral
    fk_bessel = 0.
    int_Iv = 0.
    lambda = dlambda
    do m = 1, 50-1
        lambdab = lambda
        call firstkind_bessel(order, lambdab, fk_bessel)
        lambda = lambda + dlambda
        int_Iv = int_Iv + fk_bessel*(dlambda)/(lambdab)
    end do
    call cpu_time(finish)
    print*, "Iv integral: time = ", finish-start, "seconds"

end program Gmatrix_lst