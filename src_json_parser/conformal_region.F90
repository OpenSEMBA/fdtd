module conformal_region_mod

    use geometry_mod, only: triangle_t
    use cells_mod, only: cell_interval_t

    implicit none
    
    integer, parameter :: REGION_TYPE_VOLUME = 3
    integer, parameter :: REGION_TYPE_SURFACE = 2


    type, public :: conformal_region_t
        type(triangle_t), dimension(:), allocatable :: triangles
        integer :: type
    end type

end module