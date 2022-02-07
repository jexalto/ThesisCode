program main
    use MOD_GVALUES_ODD, only:&
    Gjj_odd, Goj_odd, Gjo_odd, Goo_odd

    implicit none

    real(8) :: mu, ksi, eta, c, d, r
    real(8) :: Gjj_out, Goj_out, Gjo_out, Goo_out

    mu = 0.95
    ksi = 1.01
    eta = 1.
    r = 1.0
    c = 0.4
    d = 0.5

    ! Subscript: first letter is control point, second letter is vortex
    ! o = outside, j = jet - inside jet

    call Gjj_odd(mu, ksi, eta, c, d, r, Gjj_out)

    call Goj_odd(mu, ksi, eta, c, d, r, Goj_out)

    call Gjo_odd(mu, ksi, eta, c, d, r, Gjo_out)

    call Goo_odd(mu, ksi, eta, c, d, r, Goo_out)

end program main