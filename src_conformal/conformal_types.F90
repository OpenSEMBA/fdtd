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
        logical :: isOnEdge = .false.
        integer :: edge = NOT_ON_EDGE
    contains 
        procedure :: getEdge
        procedure, private :: isSideOnEdge
    end type

    type, public :: triangle_t
        ! real, dimension(3) :: vertices(3)
        type(coord_t), dimension(3) :: vertices
        ! private
        real, dimension(3) :: normal
        integer :: face = NOT_ON_FACE
        type(side_t), dimension(3) :: sides
        integer(kind=4), dimension(3) :: cell
    contains
        procedure :: buildTriangle
        procedure :: getNormal
        procedure :: getFace
        procedure :: getSides
        procedure :: getCell

        procedure, private :: centroid
        procedure, private :: setCell
        procedure, private :: computeNormal
        procedure, private :: setFace
        procedure, private :: buildSides
        procedure, private :: setEdges
        procedure, private :: isOnFace
    end type

        
    type, public :: conformal_region_t
        type(triangle_t), dimension(:), allocatable :: triangles
    end type
        

contains

    subroutine buildTriangle(this)
        class(triangle_t) :: this
        call this%computeNormal()
        call this%setFace()
        call this%setCell()
        call this%buildSides()
    end subroutine

    function getNormal(this) result(res)
        class(triangle_t) :: this
        real, dimension(3) :: res
        res = this%normal
    end function

    integer function getFace(this)
        class(triangle_t) :: this
        getFace = this%face
    end function

    function getCell(this) result(res)
        class(triangle_t) :: this
        integer(kind=4), dimension(3) :: res
        res = this%cell
    end function

    function getSides(this) result(res)
        class(triangle_t) :: this
        type(side_t), dimension(3) :: res
        res = this%sides
    end function

    integer function getEdge(this)
        class(side_t) :: this
        getEdge = this%edge
    end function

    logical function isSideOnEdge(this, edge)
        class(side_t) :: this
        integer :: edge
        isSideOnEdge = this%direction(edge) /= 0 .and. &
                       this%direction(mod(edge,3) + 1) == 0 .and. &
                       this%direction(mod(edge+1,3) + 1) == 0
    end function

    subroutine buildSides(this)
        class(triangle_t) :: this
        this%sides(1)%direction = this%vertices(2)%position - this%vertices(1)%position
        this%sides(2)%direction = this%vertices(3)%position - this%vertices(2)%position
        this%sides(3)%direction = this%vertices(1)%position - this%vertices(3)%position
        call this%setEdges()
    end subroutine

    subroutine setEdges(this)
        class(triangle_t) :: this
        integer :: edge, i
        do i = 1,3
            this%sides(i)%edge = NOT_ON_EDGE
            this%sides(i)%isOnEdge = .false.
            do edge = EDGE_X, EDGE_Z
                if (this%sides(i)%isSideOnEdge(edge))  then 
                    this%sides(i)%edge = edge
                    this%sides(i)%isOnEdge = .true.
                end if
            end do
        end do

    end subroutine

    subroutine setFace(this)
        class(triangle_t) :: this
        integer :: face
        this%face = NOT_ON_FACE
        do face = FACE_X, FACE_Z
            if (this%isOnFace(face)) this%face = face
        end do
    end subroutine

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
    
    subroutine setcell(this)
        class(triangle_t) :: this
        this%cell = floor(this%centroid())
    end subroutine   

    subroutine computeNormal(this)
        class(triangle_t) :: this
        real, dimension(3) :: v1,v2, cross
        v1 = this%vertices(2)%position - this%vertices(1)%position
        v2 = this%vertices(3)%position - this%vertices(2)%position
        cross = [ v1(2)*v2(3)-v1(3)*v2(2), &
                 -(v1(1)*v2(3)-v1(3)*v2(1)), &
                   v1(1)*v2(2)-v1(2)*v2(1)]
        this%normal = cross/norm2(cross)
    end subroutine
    
    logical function isOnFace(this, dir)
        class(triangle_t) :: this
        integer :: dir
        isOnFace = abs(this%normal(dir)) == 1 .and. &
                   abs(this%normal(mod(dir,3) + 1)) == 0 .and. &
                   abs(this%normal(mod(dir+1,3) + 1)) == 0
    end function

end module