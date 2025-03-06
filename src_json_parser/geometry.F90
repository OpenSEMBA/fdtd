module geometry_mod

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
    contains
        procedure :: isOnVertex => coord_isOnVertex
        procedure :: isOnEdge => coord_isOnEdge
        procedure :: isOnFace => coord_isOnFace
        procedure :: getEdge => coord_getEdge
    end type

    type, public :: side_t
        type(coord_t) :: init, end
        real, dimension(3) :: normal
    contains 
        procedure :: getEdge => side_getEdge
        procedure :: getCell => side_getCell
        procedure :: isOnEdge => side_isOnEdge
        procedure :: isOnFace => side_isOnFace
    end type

    type, public :: triangle_t
        type(coord_t), dimension(3) :: vertices
    contains
        procedure :: getNormal
        procedure :: getFace
        procedure :: getSides
        procedure :: getCell => triangle_getCell
        procedure :: isOnFace => triangle_isOnFace
        procedure :: isOnAnyFace
        procedure, private :: centroid
    end type
    
contains


    function triangle_getCell(this) result(res)
        class(triangle_t) :: this
        integer(kind=4), dimension(3) :: res
        res = floor(this%centroid())
    end function

    function side_getCell(this) result(res)
        class(side_t) :: this
        real, dimension(3) :: c
        integer(kind=4), dimension(3) :: res
        c = 0.5*(this%init%position + this%end%position)
        res = floor(c)
    end function

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

    logical function side_isOnEdge(this, edge)
        class(side_t) :: this
        integer :: edge
        real, dimension(3) :: direction
        direction = this%end%position - this%init%position
        side_isOnEdge = direction(edge) /= 0 .and. &
                        direction(mod(edge,3) + 1) == 0 .and. &
                        direction(mod(edge+1,3) + 1) == 0
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

    function getSides(this) result(res)
        class(triangle_t) :: this
        type(side_t), dimension(3) :: res
        res(1)%init%position = this%vertices(1)%position
        res(1)%end%position = this%vertices(2)%position
        res(1)%normal = this%getNormal()
        res(2)%init%position = this%vertices(2)%position
        res(2)%end%position = this%vertices(3)%position 
        res(2)%normal = this%getNormal()
        res(3)%init%position = this%vertices(3)%position
        res(3)%end%position = this%vertices(1)%position 
        res(3)%normal = this%getNormal()
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
    
    logical function side_isOnFace(this, face)
        class(side_t) :: this
        integer :: face
        type(coord_t) :: mean
        mean%position = 0.5*(this%init%position+this%end%position)
        side_isOnFace = mean%getEdge() == NOT_ON_EDGE .and. &
                        mean%isOnFace(face)
    end function

    logical function triangle_isOnFace(this, face)
        class(triangle_t) :: this
        integer :: face
        real, dimension(3) :: n
        n = this%getNormal()
        triangle_isOnFace = abs(n(face)) == 1 .and. &
                            abs(n(mod(face,3) + 1)) == 0 .and. &
                            abs(n(mod(face+1,3) + 1)) == 0
    end function

    logical function isOnAnyFace(this)
        class(triangle_t) :: this
        isOnAnyFace = .not. (this%isOnFace(FACE_X) .or. this%isOnFace(FACE_Y) .or. this%isOnFace(FACE_Z))
    end function


    function buildSidesContour(sides, face) result (res)
        type(side_t), dimension(:), allocatable, intent(in) :: sides
        integer, intent(in) :: face
        type(side_t), dimension(:), allocatable :: inner_path
        type(side_t), dimension(:), allocatable :: res
        type(coord_t) :: init, end
        integer :: i, n

        allocate(inner_path(size(sides)))
        n = 0
        do while (n < size(inner_path))
            do i = i, size(sides)
                if (n == 0 .and. (sides(i)%init%isOnEdge(1) .or. sides(i)%init%isOnEdge(2) .or. sides(i)%init%isOnEdge(3))) then 
                    n = n + 1
                    inner_path(n) = sides(i)
                else if (n /= 0 .and. all(sides(i)%init%position .eq. inner_path(n)%end%position)) then 
                    n = n + 1
                    inner_path(n) = sides(i)
                end if
            end do
        end do

        init = inner_path(1)%init
        end = inner_path(size(inner_path))%end
        if (init%isOnVertex() .and. end%isOnVertex()) then 
            res = buildVertexToVertexContour(init, end, inner_path, face)
        end if

        if (init%isOnVertex() .and. .not. end%isOnVertex()) then 
            res = buildVertexToSideContour(init, end, inner_path, face)
        else if (.not. init%isOnVertex() .and. end%isOnVertex()) then 
            res = buildVertexToSideContour(end, init, inner_path, face)
        end if
        
        if (.not. init%isOnVertex() .and. .not. end%isOnVertex()) then 
            res = buildSideToSideContour(init, end, inner_path, face)
        end if

    end function

    function buildVertexToVertexContour(init, end, inner_path, face) result(res)
        type(coord_t), intent(in) :: init, end
        type(side_t), dimension(:), allocatable, intent(in) :: inner_path
        integer, intent(in) :: face
        type(side_t), dimension(:), allocatable :: res
        integer :: idx
        real, dimension(3,4) :: corners
        
        corners = buildCorners(floor(init%position), face)
        
        allocate(res(size(inner_path) + 2))
        res(1:size(inner_path)) = inner_path
        idx = mod(cornerIndex(corners, end%position),4) + 1
        res(size(inner_path) + 1) = buildSide(end%position, corners(:,idx))
        res(size(inner_path) + 2) = buildSide(corners(:,idx), init%position)
    end function
    
    function buildVertexToSideContour(init, end, inner_path, face) result(res)
        type(coord_t), intent(in) :: init, end
        type(side_t), dimension(:), allocatable, intent(in) :: inner_path
        integer, intent(in) :: face
        type(side_t), dimension(:), allocatable :: res
        integer :: i, idx
        real, dimension(3,4) :: corners
        type(side_t) :: cell_side

        corners = buildCorners(floor(init%position), face)
        allocate(res(size(inner_path)))
        res = inner_path
        do i = 1, 4
            cell_side%init%position = corners(:,i)
            cell_side%end%position  = corners(:,mod(i,4) + 1)
            if (all(cell_side%getCell() .eq. floor(end%position)) .and. &
                (cell_side%getEdge() == end%getEdge()) ) then 
                idx = i
                exit
            end if
        end do
        call addSide(res, buildSide(end%position, corners(:,mod(idx,4) + 1)))
        do while (.not. all(corners(:,mod(idx,4) + 1) .eq. init%position))
            call addSide(res, buildSide(corners(:,mod(idx,4) + 1), corners(:,mod(idx + 1,4) + 1)))
            idx = idx + 1
        end do
    end function
    
    function buildSideToSideContour(init, end, inner_path, face) result(res)
        type(coord_t), intent(in) :: init, end
        type(side_t), dimension(:), allocatable, intent(in) :: inner_path
        integer, intent(in) :: face
        type(side_t) :: cell_side
        integer :: i, idx_i, idx_e
        real, dimension(3,4) :: corners

        type(side_t), dimension(:), allocatable :: res
        corners = buildCorners(floor(init%position), face)
        allocate(res(size(inner_path)))
        res = inner_path
        do i = 1, 4
            cell_side%init%position = corners(:,i)
            cell_side%end%position  = corners(:,mod(i,4) + 1)
            if (all(cell_side%getCell() .eq. floor(init%position)) .and. &
                (cell_side%getEdge() == init%getEdge()) ) then 
                idx_i = i
            end if
            if (all(cell_side%getCell() .eq. floor(end%position)) .and. &
                (cell_side%getEdge() == end%getEdge()) ) then 
                idx_e = i
            end if
        end do
        call addSide(res, buildSide(end%position, corners(:,mod(idx_e,4) + 1)))
        i = idx_e
        do while (mod(i,4) + 1 < idx_i)
            call addSide(res, buildSide(corners(:,mod(i,4) + 1), corners(:,mod(i + 1,4) + 1)))
            i = i + 1
        end do
        call addSide(res, buildSide(corners(:,mod(i,4) + 1), init%position))

    end function

    subroutine addSide(sides, side) 
        type(side_t), dimension(:), allocatable, intent(inout) :: sides
        type(side_t), dimension(:), allocatable :: aux
        type(side_t) :: side
        allocate(aux(size(sides) + 1))
        aux(1:size(sides)) = sides
        aux(size(sides) + 1) = side
        deallocate(sides)
        allocate(sides(size(aux)))
        sides = aux        
    end subroutine

    function buildSide(c1, c2) result(res)
        real, dimension(3), intent(in) :: c1, c2
        type(side_t) :: res
        res%init%position = c1
        res%end%position= c2
    end function

    integer function cornerIndex(corners, vertex)
        real, dimension(3,4),intent(in) :: corners
        real, dimension(3), intent(in) :: vertex
        integer :: i
        do i = 1,4
            if (all(vertex .eq. corners(:,i))) then 
                cornerIndex = i
                exit
            end if
        end do

    end function

    function buildCorners(cell, face) result(res)
        integer, dimension(3), intent(in) :: cell
        integer, intent(in) :: face
        integer, dimension(3,4) :: res
        if (face == FACE_X) then 
            res(:,1) = cell(:)
            res(:,2) = cell(:) + [0,1,0]
            res(:,3) = cell(:) + [0,1,1]
            res(:,4) = cell(:) + [0,0,1]
        else if (face == FACE_Y) then 
            res(:,1) = cell(:)
            res(:,2) = cell(:) + [0,0,1]
            res(:,3) = cell(:) + [1,0,1]
            res(:,4) = cell(:) + [1,0,0]
        else if (face == FACE_Z) then 
            res(:,1) = cell(:)
            res(:,2) = cell(:) + [1,0,0]
            res(:,4) = cell(:) + [1,1,0]
            res(:,3) = cell(:) + [0,1,0]
        end if
    end function

    function getContourOnFace(sides) result (res)
        type(side_t), dimension(:), allocatable, intent(in) :: sides
        type(coord_t), dimension(:), allocatable :: res
        type(coord_t) :: init, end
        integer :: i, n
        allocate(res(size(sides) + 1))
        n = 0
        do while (n < size(res))
           do i = 1, size(sides) 
              if (n == 0 .and. (sides(i)%init%isOnEdge(1) .or. sides(i)%init%isOnEdge(2) .or. sides(i)%init%isOnEdge(3))) then 
                 n = n + 1
                 res(n) = sides(i)%init
              else if (n /= 0 .and. all(sides(i)%init%position .eq. res(n)%position)) then 
                 n = n + 1
                 res(n) = sides(i)%end
              end if
           end do
        end do
    end function

    function getTrianglesOffFaces(triangles) result (res)
        type(triangle_t), dimension(:), allocatable, intent(in) :: triangles
        type(triangle_t), dimension(:), allocatable :: res
        integer :: i, j, n
        n = 0
        do i = 1, size(triangles)
           if (.not. (triangles(i)%isOnFace(1) .or. triangles(i)%isOnFace(2) .or. triangles(i)%isOnFace(3))) n = n + 1
        end do
        allocate(res(n))
        n = 0
        do i = 1, size(triangles)
           if (.not. (triangles(i)%isOnFace(1) .or. triangles(i)%isOnFace(2) .or. triangles(i)%isOnFace(3))) then 
              n = n + 1
              res(n) = triangles(j)
           end if
        end do
 
    end function   
 
    function getSidesOnFace(sides, face) result (res)
        type(side_t), dimension(:), allocatable, intent(in) :: sides
        integer, intent(in) :: face
        type(side_t), dimension(:), allocatable :: res
        integer :: i, n
        n = 0
        do i = 1, size(sides)
            if (sides(i)%isOnFace(face)) then
                n = n + 1
            end if
        end do
        do i = 1, size(sides)
            if (sides(i)%isOnFace(face)) then
                n = n + 1
                res(n) = sides(i)
            end if
        end do
    end function
 
    function getSidesOnEdges(sides, face) result (res)
        type(side_t), dimension(:), allocatable, intent(in) :: sides
        integer, intent(in) :: face
        type(side_t), dimension(:), allocatable :: res
        integer :: i, n
        n = 0
        do i = 1, size(sides)
            if (sides(i)%isOnEdge(mod(face,3) + 1) .or. & 
                sides(i)%isOnEdge(mod(face+1,3) + 1)) then
                n = n + 1
            end if
        end do
        do i = 1, size(sides)
            if (sides(i)%isOnEdge(mod(face,3) + 1) .or. & 
                sides(i)%isOnEdge(mod(face+1,3) + 1)) then
                n = n + 1
                res(n) = sides(i)
            end if
        end do
    end function

end module