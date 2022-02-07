module MATRIX
    use MOD_GVALUES_EVEN, only:&
    Gjj_even, Goj_even, Gjo_even, Goo_even

    use MOD_GVALUES_odd, only:&
    Gjj_odd, Goj_odd, Gjo_odd, Goo_odd

    public :: corr_matrix, Gmatrix

    contains

    subroutine Gmatrix(nx, ny, nr_prop, loc_prop, r_prop, Vj, V0, span, chord)
        integer, intent(inout) :: nx, ny, nr_prop
        real(8), intent(inout) :: loc_prop, r_prop, Vj
        real(8), intent(inout) :: V0, span, chord
        real(8) :: coll_pts
        real(8), dimension(1:ny+1) :: coll_pt

        coll_pts = span/(ny+1)


    end subroutine Gmatrix
    
    subroutine corr_matrix(nx, ny, nr_prop, loc_prop, r_prop, Vj, V0, span, chord)
        integer, intent(inout) :: nx, ny, nr_prop
        real(8), dimension(1:nr_prop) intent(inout) :: loc_prop, r_prop, Vj
        real(8), intent(inout) :: V0


    end subroutine corr_matrix