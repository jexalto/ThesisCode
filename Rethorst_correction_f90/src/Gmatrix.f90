program Gmatrix
    ! call cpu_time(start)
    ! call cpu_time(finish)
    ! print*, "Time = ", finish-start, "seconds"
    ! real :: start, finish
    
    use MOD_MODIFIEDBESSEL, only: &
        firstkind_bessel

    implicit none

    character(len=*), parameter :: name = 'hello world'
    real :: dlambda, lambda, lambdab, fk_bessel, int_Iv, beta, int_Kv
    integer :: order, m

    order = 1
    lambda = 0.
    fk_bessel = 0.
    int_Iv = 0.
    int_Kv = 0.
    beta = 1.1

    ! Iv_lambdaB integral
    do m = 1, 50+1
        lambdab = lambda * beta
        call firstkind_bessel(order, lambdab, fk_bessel)
        lambda = lambda + 0.1
        int_Iv = int_Iv + fk_bessel*dlambda/(lambdab)
    end do

    ! Kv_lambdaB integral
    lambda = 0.
    do m = 1, 50+1
        lambdab = lambda * beta
        call secondkind_bessel(order, lambdab, sk_bessel)
        lambda = lambda + 0.1
        int_Iv = int_Kv + sk_bessel*dlambda/(lambdab)
    end do

end program Gmatrix