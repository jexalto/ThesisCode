module MOD_GVALUES_EVEN

    implicit none

    public :: Gjj_even, Goj_even, Gjo_even, Goo_even

    contains

    subroutine Gjj_even(mu, eta, c, d, r, Gjj_out)
        ! Input: mu (tip speed ratio), ksi (location x/r), eta (location y/r), c-d (location vortex), r (radius)
        ! Output; Gjj (correction factor), even
        real(8), intent(inout) :: mu, eta, c, d, r, Gjj_out

        Gjj_out = 1/r * (1 - mu*mu) / (1 + mu*mu) * (1/(1/d - eta) - 1/(1/c - eta) + 1/(1/d + eta) + 1/(1/d + eta))
    
    end subroutine Gjj_even

    subroutine Goj_even(mu, eta, c, d, r, Goj_out)
        ! Input: mu (tip speed ratio), ksi (location x/r), eta (location y/r), c-d (location vortex), r (radius)
        ! Output; Goj (correction factor), even
        real(8), intent(inout) :: mu, eta, c, d, r, Goj_out

        Goj_out = 1/r * (1 - mu)*(1 - mu) / (1 + mu*mu) * (1/(eta - c) - 1/(eta - d) + 1/(eta + d) + 1/(eta + c))

    end subroutine Goj_even

    subroutine Gjo_even(mu, eta, c, d, r, Gjo_out)
        ! Input: mu (tip speed ratio), ksi (location x/r), eta (location y/r), c-d (location vortex), r (radius)
        ! Output; Gjj (correction factor), even
        real(8), intent(inout) :: mu, eta, c, d, r, Gjo_out

        Gjo_out = 1/r * (1 - mu)*(1 - mu) / (1 + mu*mu) * (1/(eta - c) - 1/(eta - d) + 1/(eta + d) + 1/(eta + c))

    end subroutine Gjo_even

    subroutine Goo_even(mu, eta, c, d, r, Goo_out)
        ! Input: mu (tip speed ratio), ksi (location x/r), eta (location y/r), c-d (location vortex), r (radius)
        ! Output; Goo (correction factor), even
        real(8), intent(inout) :: mu, eta, c, d, r, Goo_out

        Goo_out = -1/r * (1 - mu*mu) / (1 + mu*mu) * (1/(1/d - eta) - 1/(1/c - eta) + 1/(1/d + eta) + 1/(1/d + eta))

    end subroutine Goo_even
    
end module MOD_GVALUES_EVEN