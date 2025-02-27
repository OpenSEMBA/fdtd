module conformal_types_mod

    use cells_mod

    implicit none
    
    integer, parameter :: FACE_X = 1
    integer, parameter :: FACE_Y = 2
    integer, parameter :: FACE_Z = 3
    integer, parameter :: NOT_ON_FACE = -1
    
    integer, parameter :: EDGE_X = 1
    integer, parameter :: EDGE_Y = 2
    integer, parameter :: EDGE_Z = 3
    integer, parameter :: NOT_ON_EDGE = -1
 
    type, public :: coord_t
        real, dimension(3) :: position
        integer :: id
    end type

    type, public :: side_t
        real, dimension(3) :: direction
    contains 
        procedure :: getEdge
        procedure :: isSideOnEdge
    end type

    type, public :: triangle_t
        type(coord_t), dimension(3) :: vertices
    contains
        procedure :: getNormal
        procedure :: getFace
        procedure :: getSides
        procedure :: getCell
        procedure :: isOnFace

        procedure, private :: centroid
    end type

        
    type, public :: conformal_region_t
        type(triangle_t), dimension(:), allocatable :: triangles
        type(cell_interval_t), dimension(:), allocatable :: intervals
    end type
        

contains


    function getCell(this) result(res)
        class(triangle_t) :: this
        integer(kind=4), dimension(3) :: res
        res = floor(this%centroid())
    end function

    logical function isSideOnEdge(this, edge)
        class(side_t) :: this
        integer :: edge
        isSideOnEdge = this%direction(edge) /= 0 .and. &
                       this%direction(mod(edge,3) + 1) == 0 .and. &
                       this%direction(mod(edge+1,3) + 1) == 0
    end function

    function getEdge(this) result(res)
        class(side_t) :: this
        integer :: res, edge
        res = NOT_ON_EDGE
        do edge = EDGE_X, EDGE_Z
            if (this%isSideOnEdge(edge))  then 
                res = edge
            end if
        end do
    end function

    function getSides(this) result(res)
        class(triangle_t) :: this
        type(side_t), dimension(3) :: res
        res(1)%direction = this%vertices(2)%position - this%vertices(1)%position
        res(2)%direction = this%vertices(3)%position - this%vertices(2)%position
        res(3)%direction = this%vertices(1)%position - this%vertices(3)%position
    end function

    function getFace(this) result(res)
        class(triangle_t) :: this
        integer :: face, res
        res = NOT_ON_FACE
        do face = FACE_X, FACE_Z
            if (this%isOnFace(face)) res = face
        end do
    end function

    function centroid(this) result(res)
        class(triangle_t) :: this
        real, dimension(3) :: res
        integer :: i
        res(:) = 0.0
        do i = 1,3
            res = res + this%vertices(i)%position
        end do
        res = res/3.0
    end function
    
    function getNormal(this) result(res)
        class(triangle_t) :: this
        real, dimension(3) :: v1,v2, res
        v1 = this%vertices(2)%position - this%vertices(1)%position
        v2 = this%vertices(3)%position - this%vertices(2)%position
        res = [ v1(2)*v2(3)-v1(3)*v2(2), &
                -(v1(1)*v2(3)-v1(3)*v2(1)), &
                v1(1)*v2(2)-v1(2)*v2(1)]
        res = res/norm2(res)
    end function
    
    logical function isOnFace(this, dir)
        class(triangle_t) :: this
        integer :: dir
        real, dimension(3) :: n
        n = this%getNormal()
        isOnFace = abs(n(dir)) == 1 .and. &
                   abs(n(mod(dir,3) + 1)) == 0 .and. &
                   abs(n(mod(dir+1,3) + 1)) == 0
    end function

end module