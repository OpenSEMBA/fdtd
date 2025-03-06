module conformal_types_mod

    use geometry_mod, only: triangle_t
    use cells_mod, only: cell_interval_t

    implicit none
    
    type, public :: conformal_region_t
        type(triangle_t), dimension(:), allocatable :: triangles
        type(cell_interval_t), dimension(:), allocatable :: intervals
        integer :: offgrid_points
    end type
    


end module