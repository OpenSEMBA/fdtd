module conformal_types_mod

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
        integer :: id = -1
    contains
        procedure :: isOnVertex => coord_isOnVertex
        procedure :: isOnEdge => coord_isOnEdge
        procedure :: isOnAnyEdge => coord_isOnAnyEdge
        procedure :: isOnFace => coord_isOnFace
        procedure :: isOnAnyFace => coord_isOnAnyFace
        procedure :: getEdge => coord_getEdge
    end type

    type, public :: side_t
        type(coord_t) :: init, end
        real, dimension(3) :: normal = [0.0,0.0,0.0]
    contains 
        procedure :: getEdge => side_getEdge
        procedure :: getCell => side_getCell
        procedure :: getFace => side_getFace
        procedure :: isOnEdge => side_isOnEdge
        procedure :: isInCell => side_isInCell
        procedure :: isOnAnyEdge => side_isOnAnyEdge
        procedure :: isOnFace => side_isOnFace
        procedure :: isOnAnyFace => side_isOnAnyFace
        procedure :: length
        procedure :: isEquiv
    end type

    type, public :: side_list_t
        type(side_t), dimension(:), allocatable :: sides
    end type

    type, public :: triangle_t
        type(coord_t), dimension(3) :: vertices
    contains
        procedure :: getNormal => triangle_getNormal
        procedure :: getFace => triangle_getFace
        procedure :: getSides
        procedure :: getCell => triangle_getCell
        procedure :: isOnFace => triangle_isOnFace
        procedure :: isOnAnyFace => triangle_isOnAnyFace
        procedure, private :: centroid
    end type

    type :: cell_t
        integer, dimension(3) :: cell
    end type

    type :: interval_t
        type(cell_t) :: ini, end
    end type


contains

    logical function coord_isOnEdge(this, edge) 
        class(coord_t) :: this
        integer :: edge
        real, dimension(3) :: delta
        delta = (this%position - floor(this%position))
        coord_isOnEdge = delta(edge) /= 0 .and. &
                         delta(mod(edge,3) + 1) == 0 .and. &
                         delta(mod(edge+1,3) + 1) == 0
    end function

    logical function coord_isOnFace(this, face)
        class(coord_t) :: this
        integer :: face
        real, dimension(3) :: delta
        delta = (this%position - floor(this%position))
        coord_isOnFace = delta(face) == 0 .and. &
                         delta(mod(face,3) + 1) /= 0 .and. &
                         delta(mod(face+1,3) + 1) /= 0
    end function


    function coord_getEdge(this) result(res)
        class(coord_t) :: this
        integer :: res, edge
        res = NOT_ON_EDGE
        if (.not. this%isOnVertex()) then
            do edge = EDGE_X, EDGE_Z
                if (this%isOnEdge(edge))  then 
                    res = edge
                end if
            end do
        end if
    end function

    logical function coord_isOnVertex(this) 
        class(coord_t) :: this
        real, dimension(3) :: delta
        delta = (this%position - floor(this%position))
        coord_isOnVertex = all(delta .eq. [0,0,0])
    end function

    logical function coord_isOnAnyFace(this)
        class(coord_t) :: this
        coord_isOnAnyFace = (this%isOnFace(FACE_X) .or. this%isOnFace(FACE_Y) .or. this%isOnFace(FACE_Z))
    end function

    logical function coord_isOnAnyEdge(this)
        class(coord_t) :: this
        coord_isOnAnyEdge = (this%isOnEdge(EDGE_X) .or. this%isOnEdge(EDGE_Y) .or. this%isOnEdge(EDGE_Z))
    end function

    function side_getCell(this) result(res)
        class(side_t) :: this
        real, dimension(3) :: c
        integer(kind=4), dimension(3) :: res
        c = 0.5*(this%init%position + this%end%position)
        res = floor(c)
    end function


    logical function side_isInCell(this, cell)
        class(side_t) :: this
        integer(kind=4), dimension(3) :: cell
        side_isInCell = all(this%getCell() .eq. cell)
    end function

    logical function side_isOnEdge(this, edge)
        class(side_t) :: this
        integer :: edge
        type(coord_t) :: c
        c = coord_t(position=0.5*(this%end%position + this%init%position))
        side_isOnEdge = c%isOnEdge(edge)
    end function

    function side_getEdge(this) result(res)
        class(side_t) :: this
        integer :: res, edge
        res = NOT_ON_EDGE
        do edge = EDGE_X, EDGE_Z
            if (this%isOnEdge(edge))  then 
                res = edge
            end if
        end do
    end function

    function side_getFace(this) result(res)
        class(side_t) :: this
        integer :: face, res
        res = NOT_ON_FACE
        do face = FACE_X, FACE_Z
            if (this%isOnFace(face)) res = face
        end do
    end function

    logical function side_isOnFace(this, face)
        class(side_t) :: this
        integer :: face
        type(coord_t) :: mean
        mean%position = 0.5*(this%init%position+this%end%position)
        side_isOnFace = mean%getEdge() == NOT_ON_EDGE .and. &
                        mean%isOnFace(face)
    end function

    logical function side_isOnAnyFace(this)
        class(side_t) :: this
        side_isOnAnyFace = (this%isOnFace(FACE_X) .or. this%isOnFace(FACE_Y) .or. this%isOnFace(FACE_Z))
    end function


    logical function side_isOnAnyEdge(this)
        class(side_t) :: this
        side_isOnAnyEdge = (this%isOnEdge(EDGE_X) .or. this%isOnEdge(EDGE_Y) .or. this%isOnEdge(EDGE_Z))
    end function

    function length(this) result(res)
        class(side_t) :: this
        real :: res
        res = norm2(this%init%position - this%end%position)
    end function


    logical function isEquiv(this, side)
        class(side_t) :: this
        type(side_t), intent(in) :: side
        logical :: eq, eq_inv
        integer :: i
        eq = .true.
        eq_inv = .true.
        do i = 1, 3
            eq = eq .and. & 
                (abs(this%init%position(i) - side%init%position(i)) < 0.01) .and. &
                (abs(this%end%position(i) - side%end%position(i)) < 0.01)
            eq_inv = eq_inv .and. & 
                (abs(this%init%position(i) - side%end%position(i)) < 0.01) .and. &
                (abs(this%end%position(i) - side%init%position(i)) < 0.01)
        end do
        isEquiv = eq .or. eq_inv
    end function

    function triangle_getCell(this) result(res)
        class(triangle_t) :: this
        integer(kind=4), dimension(3) :: res
        res = floor(this%centroid())
    end function



    function getSides(this) result(res)
        class(triangle_t) :: this
        type(side_t), dimension(3) :: res
        integer :: i
        do i = 1, 3
            res(i)%init%position = this%vertices(i)%position
            res(i)%end%position = this%vertices(mod(i,3)+1)%position
            res(i)%init%id = this%vertices(i)%id
            res(i)%end%id = this%vertices(mod(i,3)+1)%id
            res(i)%normal = this%getNormal()
        end do
    end function

    function triangle_getFace(this) result(res)
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
    
    function getNormal(contour) result(res)
        type(side_t), dimension(:), allocatable :: contour
        real, dimension(3) :: v1, v2, res
        v1 = contour(1)%end%position - contour(1)%init%position
        v2 = contour(2)%end%position - contour(2)%init%position
        res = [ v1(2)*v2(3)-v1(3)*v2(2), &
                -(v1(1)*v2(3)-v1(3)*v2(1)), &
                v1(1)*v2(2)-v1(2)*v2(1)]
        res = res/norm2(res)
        
    end function

    function triangle_getNormal(this) result(res)
        class(triangle_t) :: this
        real, dimension(3) :: v1,v2, res
        v1 = this%vertices(2)%position - this%vertices(1)%position
        v2 = this%vertices(3)%position - this%vertices(2)%position
        res = [ v1(2)*v2(3)-v1(3)*v2(2), &
                -(v1(1)*v2(3)-v1(3)*v2(1)), &
                v1(1)*v2(2)-v1(2)*v2(1)]
        res = res/norm2(res)
    end function

    logical function triangle_isOnFace(this, face)
        class(triangle_t) :: this
        integer :: face
        type(coord_t) :: c
        c = coord_t(position=this%centroid())
        triangle_isOnFace = c%isOnFace(face)
    end function

    logical function triangle_isOnAnyFace(this)
        class(triangle_t) :: this
        triangle_isOnAnyFace = (this%isOnFace(FACE_X) .or. this%isOnFace(FACE_Y) .or. this%isOnFace(FACE_Z))
    end function


end module