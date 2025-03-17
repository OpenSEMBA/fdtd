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
        procedure :: isOnAnyEdge => side_isOnAnyEdge
        procedure :: isOnFace => side_isOnFace
        procedure :: isOnAnyFace => side_isOnAnyFace
        procedure :: length
        procedure :: isEquiv
    end type

    type, public :: triangle_t
        type(coord_t), dimension(3) :: vertices
    contains
        procedure :: getNormal
        procedure :: getFace => triangle_getFace
        procedure :: getSides
        procedure :: getCell => triangle_getCell
        procedure :: isOnFace => triangle_isOnFace
        procedure :: isOnAnyFace => triangle_isOnAnyFace
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
        type(coord_t) :: c
        c = coord_t(position=0.5*(this%end%position + this%init%position))
        side_isOnEdge = c%isOnEdge(edge)
    end function

    function length(this) result(res)
        class(side_t) :: this
        real :: res
        res = norm2(this%init%position - this%end%position)
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

    function side_getFace(this) result(res)
        class(side_t) :: this
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
        res = cross(v1,v2)
        ! res = [ v1(2)*v2(3)-v1(3)*v2(2), &
        !         -(v1(1)*v2(3)-v1(3)*v2(1)), &
        !         v1(1)*v2(2)-v1(2)*v2(1)]
        res = res/norm2(res)
    end function
    
    function cross(v1,v2) result(res)
        real, dimension(3) :: v1,v2, res
        res = [ v1(2)*v2(3)-v1(3)*v2(2), &
                -(v1(1)*v2(3)-v1(3)*v2(1)), &
                v1(1)*v2(2)-v1(2)*v2(1)]
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
        type(coord_t) :: c
        c = coord_t(position=this%centroid())
        triangle_isOnFace = c%isOnFace(face)
    end function

    logical function triangle_isOnAnyFace(this)
        class(triangle_t) :: this
        triangle_isOnAnyFace = (this%isOnFace(FACE_X) .or. this%isOnFace(FACE_Y) .or. this%isOnFace(FACE_Z))
    end function

    logical function side_isOnAnyFace(this)
        class(side_t) :: this
        side_isOnAnyFace = (this%isOnFace(FACE_X) .or. this%isOnFace(FACE_Y) .or. this%isOnFace(FACE_Z))
    end function

    logical function coord_isOnAnyFace(this)
        class(coord_t) :: this
        coord_isOnAnyFace = (this%isOnFace(FACE_X) .or. this%isOnFace(FACE_Y) .or. this%isOnFace(FACE_Z))
    end function

    logical function coord_isOnAnyEdge(this)
        class(coord_t) :: this
        coord_isOnAnyEdge = (this%isOnEdge(EDGE_X) .or. this%isOnEdge(EDGE_Y) .or. this%isOnEdge(EDGE_Z))
    end function

    logical function side_isOnAnyEdge(this)
        class(side_t) :: this
        side_isOnAnyEdge = (this%isOnEdge(EDGE_X) .or. this%isOnEdge(EDGE_Y) .or. this%isOnEdge(EDGE_Z))
    end function


    function buildCellSideSet(sides) result(res)
        type(side_t), dimension(:), allocatable, intent(in) :: sides
        type(side_t), dimension(:), allocatable :: contour, res
        integer :: face, edge
        allocate(res(0))
        do face = FACE_X, FACE_Z
            contour = buildSidesContour(getSidesOnFace(sides, face))
            call addNewSides(res, contour)
        end do
        do edge = EDGE_X, EDGE_Z
            contour = getSidesOnEdge(sides, edge)
            call addNewSides(res, contour)
        end do
    end function 

    subroutine addNewSides(sides, new_sides)
        type(side_t), dimension(:), allocatable, intent(inout) :: sides
        type(side_t), dimension(:), allocatable, intent(in) :: new_sides
        integer :: i
        do i = 1, size(new_sides)
            if (isNewSide(sides, new_sides(i))) call addNewSide(sides, new_sides(i))
        end do
    end subroutine

    logical function isNewSide(sides, side)
        type(side_t), dimension(:), allocatable, intent(in) :: sides
        type(side_t), intent(in) :: side
        integer :: i
        isNewSide = .true.
        do i = 1, size(sides)
            if (sides(i)%isEquiv(side)) isNewSide = .false.
        end do
    end function    

    subroutine addNewSide(sides, side)
        type(side_t), dimension(:), allocatable, intent(inout) :: sides
        type(side_t), intent(in) :: side
        type(side_t), dimension(:), allocatable :: aux
        if (size(sides) == 0) then 
            deallocate(sides)
            allocate(sides(1))
            sides(1) = side
        else 
            allocate(aux(size(sides) + 1))
            aux(1:size(sides)) = sides
            aux(size(sides) + 1) = side
            deallocate(sides)
            allocate(sides(size(aux)))
            sides = aux
        end if
    end subroutine

    function buildSidesContour(sides) result (res)
        type(side_t), dimension(:), allocatable, intent(in) :: sides
        type(side_t), dimension(:), allocatable :: inner_path
        type(side_t), dimension(:), allocatable :: res
        type(coord_t) :: init, end

        if (size(sides) == 0) then 
            allocate(res(0))
        else
            inner_path = getPathOnFace(sides)
            init = inner_path(1)%init
            end = inner_path(size(inner_path))%end
            if (init%isOnVertex() .and. end%isOnVertex()) then 
                res = buildVertexToVertexContour(inner_path)
            else if (init%isOnVertex() .and. .not. end%isOnVertex()) then 
                res = buildVertexToSideContour(inner_path)
            else if (.not. init%isOnVertex() .and. end%isOnVertex()) then 
                res = buildSideToVertexContour(inner_path)
            elseif (.not. init%isOnVertex() .and. .not. end%isOnVertex()) then 
                res = buildSideToSideContour(inner_path)
            end if
        end if
    end function

    function buildVertexToVertexContour(inner_path) result(res)
        type(side_t), dimension(:), allocatable, intent(in) :: inner_path
        type(side_t), dimension(:), allocatable :: res
        integer :: mid_corner_idx
        real, dimension(3,4) :: corners

        corners = buildCorners(inner_path(1), inner_path(1)%getFace())
        
        allocate(res(size(inner_path) + 2))
        res(1:size(inner_path)) = inner_path
        mid_corner_idx = mod(cornerIndex(corners, inner_path(size(inner_path))%end%position),4) + 1
        res(size(inner_path) + 1) = buildSide(inner_path(size(inner_path))%end%position, corners(:,mid_corner_idx))
        res(size(inner_path) + 2) = buildSide(corners(:,mid_corner_idx), inner_path(1)%init%position)
    end function
    
    function buildVertexToSideContour(inner_path) result(res)
        type(side_t), dimension(:), allocatable, intent(in) :: inner_path
        type(side_t), dimension(:), allocatable :: res
        type(coord_t) :: init, end
        integer :: i, idx
        real, dimension(3,4) :: corners
        type(side_t) :: cell_side

        init = inner_path(1)%init
        end = inner_path(size(inner_path))%end
        corners = buildCorners(inner_path(1), inner_path(1)%getFace())

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

    function buildSideToVertexContour(inner_path) result(res)
        type(side_t), dimension(:), allocatable, intent(in) :: inner_path
        type(side_t), dimension(:), allocatable :: res
        type(coord_t) :: init, end
        integer :: idx
        real, dimension(3,4) :: corners
        type(side_t) :: cell_side

        init = inner_path(1)%init
        end = inner_path(size(inner_path))%end
        corners = buildCorners(inner_path(1), inner_path(1)%getFace())

        allocate(res(size(inner_path)))
        res = inner_path

        idx = cornerIndex(corners, end%position)

        cell_side%init%position = corners(:,idx)
        cell_side%end%position  = corners(:,mod(idx,4) + 1)
        do while (.not. (all(cell_side%getCell() .eq. floor(init%position)) .and. &
                 (cell_side%getEdge() == init%getEdge())))
                 call addSide(res, buildSide(cell_side%init%position, cell_side%end%position))
                cell_side%init%position = corners(:,mod(idx,4) + 1)
                cell_side%end%position  = corners(:,mod(idx + 1,4) + 1)
                idx = idx + 1
        end do
        call addSide(res, buildSide(cell_side%init%position, init%position))

    end function
    
    function buildSideToSideContour(inner_path) result(res)
        type(side_t), dimension(:), allocatable, intent(in) :: inner_path
        type(side_t), dimension(:), allocatable :: res
        type(side_t) :: cell_side
        integer :: i, idx_i, idx_e
        real, dimension(3,4) :: corners
        type(coord_t) :: init, end


        init = inner_path(1)%init
        end = inner_path(size(inner_path))%end

        corners = buildCorners(inner_path(1), inner_path(1)%getFace())
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
        idx_e = mod(idx_e,4) + 1
        do while (mod(idx_e,4) + 1 <= idx_i)
            call addSide(res, buildSide(corners(:,idx_e), corners(:,mod(idx_e,4) + 1)))
            idx_e = mod(idx_e, 4) + 1
        end do
        call addSide(res, buildSide(corners(:,idx_e), init%position))

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

    function buildCorners(side, face) result(res)
        type(side_t), intent(in) :: side
        integer, intent(in) :: face
        integer, dimension(3,4) :: res
        integer, dimension(3) :: cell, aux
        cell = side%getCell()

        if (face == FACE_X) then 
            res(:,1) = cell(:)
            res(:,2) = cell(:) + [0,1,0]
            res(:,3) = cell(:) + [0,1,1]
            res(:,4) = cell(:) + [0,0,1]
        else if (face == FACE_Y) then 
            res(:,1) = cell(:)
            res(:,4) = cell(:) + [1,0,0]
            res(:,3) = cell(:) + [1,0,1]
            res(:,2) = cell(:) + [0,0,1]
        else if (face == FACE_Z) then 
            res(:,1) = cell(:)
            res(:,2) = cell(:) + [1,0,0]
            res(:,3) = cell(:) + [1,1,0]
            res(:,4) = cell(:) + [0,1,0]
        end if

        if (isClockwise(side, face)) then 
            aux = res(:,2)
            res(:,2) = res(:,4)
            res(:,4) = aux
        end if
    end function

    logical function isClockwise(side, face)
        type(side_t), intent(in) :: side
        integer, intent(in) :: face
        real, dimension(3) :: x_prod
        isClockwise = .true.
        x_prod = cross(side%end%position - side%init%position, side%normal)
        if (x_prod(face) < 0) isClockwise = .false.
    end function

    function contourArea(contour) result(res)
        type(side_t), dimension(:), allocatable, intent(in) :: contour
        type(side_t), dimension(:), allocatable :: aux_contour
        ! type(side_t) :: aux_side
        real :: res
        integer :: face, i, dir1,dir2
        face = contour(1)%getFace()
        allocate(aux_contour(size(contour)))
        aux_contour = contour
        if (isClockwise(contour(1), face)) then 
            do i = 1, size(contour)
                aux_contour(size(contour) + 1 - i)%init = contour(i)%end
                aux_contour(size(contour) + 1 - i)%end = contour(i)%init
            end do
        end if

        dir1 = mod(face,3)+1
        dir2 = mod(face+1,3)+1
        res = 0
        do i = 1, size(aux_contour)
           res = res + aux_contour(i)%init%position(dir1)*aux_contour(i)%end%position(dir2) - & 
                       aux_contour(i)%end%position(dir1)*aux_contour(i)%init%position(dir2)
        end do
        res = 0.5*res
    end function
  

    function getPathOnFace(sides) result (res)
        type(side_t), dimension(:), allocatable, intent(in) :: sides
        type(side_t), dimension(:), allocatable :: res
        integer :: i, n
        allocate(res(size(sides)))
        n = 0
        do while (n < size(res))
           do i = 1, size(sides) 
              if (n == 0) then 
                if(.not. sides(i)%init%isOnAnyFace()) then 
                 n = n + 1
                 res(n) = sides(i)
                end if
              else if (n /= 0 .and. all(sides(i)%init%position .eq. res(n)%end%position)) then 
                 n = n + 1
                 res(n) = sides(i)
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
        allocate(res(n))
        n = 0
        do i = 1, size(sides)
            if (sides(i)%isOnFace(face)) then
                n = n + 1
                res(n) = sides(i)
            end if
        end do
    end function
 
    function getSidesOnEdge(sides, edge) result (res)
        type(side_t), dimension(:), allocatable, intent(in) :: sides
        integer, intent(in) :: edge
        type(side_t), dimension(:), allocatable :: res
        integer :: i, n
        n = 0
        do i = 1, size(sides)
            if (sides(i)%isOnEdge(edge)) then
                n = n + 1
            end if
        end do
        allocate(res(n))
        n = 0
        do i = 1, size(sides)
            if (sides(i)%isOnEdge(edge)) then
                n = n + 1
                res(n) = sides(i)
            end if
        end do
    end function
 
    function getSidesOnAdjacentEdges(sides, face) result (res)
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
        allocate(res(n))
        n = 0
        do i = 1, size(sides)
            if (sides(i)%isOnEdge(mod(face,3) + 1) .or. & 
                sides(i)%isOnEdge(mod(face+1,3) + 1)) then
                n = n + 1
                res(n) = sides(i)
            end if
        end do
    end function

end module