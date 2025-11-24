module geometry_mod

    use conformal_types_mod

contains

    function getArea(triangle) result(res)
        type(triangle_t) :: triangle
        real :: res
        type(side_t), dimension(:), allocatable :: sides
        sides = triangle%getSides()
        res = contourArea(sides)
    end function

    function cross(v1,v2) result(res)
        real, dimension(3) :: v1,v2, res
        res = [ v1(2)*v2(3)-v1(3)*v2(2), &
                -(v1(1)*v2(3)-v1(3)*v2(1)), &
                v1(1)*v2(2)-v1(2)*v2(1)]
    end function

    function mergeSides(sides, edge) result(res)
        type(side_t), dimension(:), allocatable, intent(in) :: sides
        type(side_t), dimension(:), allocatable :: sides_copy
        integer(kind=4), intent(in) :: edge
        real :: c
        type(side_t) :: res
        integer :: i
        sides_copy = sides
        do i = 1, size(sides_copy)
            if (sides_copy(i)%init%position(edge) > sides_copy(i)%end%position(edge)) then 
                c = sides_copy(i)%init%position(edge)
                sides_copy(i)%init%position(edge) = sides_copy(i)%end%position(edge)
                sides_copy(i)%end%position(edge) = c
            end if
        end do
        res = sides_copy(1)
        do i = 2, size(sides_copy)
            if (sides_copy(i)%init%position(edge) < res%init%position(edge)) then 
                res%init%position(edge) = sides_copy(i)%init%position(edge)
            end if
            if (sides_copy(i)%end%position(edge) > res%end%position(edge)) then 
                res%end%position(edge) = sides_copy(i)%end%position(edge)
            end if
        end do
    end function

    function findContourCell(contour) result(res)
        type(side_t), dimension(:), allocatable, intent(in) :: contour
        integer, dimension(3) :: res
        integer :: i 
        do i = 1, size(contour)
            if (contour(i)%isOnAnyFace()) then 
                res = contour(i)%getCell()
            end if
        end do
    end function

    function findContourFace(contour) result(res)
        type(side_t), dimension(:), allocatable, intent(in) :: contour
        integer :: res
        integer :: i 
        res = NOT_ON_FACE
        do i = 1, size(contour)
            if (contour(i)%isOnAnyFace()) then 
                res = contour(i)%getFace()
            end if
        end do
        if (res == NOT_ON_FACE) error stop 'Contour face could not be identified'
    
    end function

    function buildCellSideSet(sides, on_sides) result(res)
        type(side_t), dimension(:), allocatable, intent(in) :: sides, on_sides
        type(side_t), dimension(:), allocatable :: sides_on_face
        type(side_t), dimension(:), allocatable :: contour, aux, res
        integer :: face, edge
        allocate(res(0))
        do face = FACE_X, FACE_Z
            sides_on_face = getSidesOnFace(sides, face)
            contour = buildSidesContour(sides_on_face)
            call addNewSides(res, contour)
        end do
        do edge = EDGE_X, EDGE_Z
            aux = getSidesOnEdge(sides, edge)
            call addNewSides(res, aux)
            aux = getSidesOnEdge(on_sides, edge)
            call addNewSides(res, aux)
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
                 (cell_side%getEdge() == init%getEdge()) ))
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
        integer :: i, idx_i, idx_e, idx
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
        idx = mod(idx_e,4) + 1
        call addSide(res, buildSide(end%position, corners(:,idx)))

        
        cell_side%init%position = corners(:,idx)
        cell_side%end%position  = corners(:,mod(idx,4) + 1)
        do while (.not. (all(cell_side%getCell() .eq. floor(init%position)) .and. &
                 (cell_side%getEdge() == init%getEdge()) ))
                 call addSide(res, buildSide(cell_side%init%position, cell_side%end%position))
                cell_side%init%position = corners(:,mod(idx,4) + 1)
                cell_side%end%position  = corners(:,mod(idx + 1,4) + 1)
                idx = idx + 1
        end do
        call addSide(res, buildSide(cell_side%init%position, init%position))

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

    function contourArea(contour, orientation) result(res)
        type(side_t), dimension(:), allocatable, intent(in) :: contour
        type(side_t), dimension(:), allocatable :: aux_contour
        integer, optional :: orientation
        real :: res
        integer :: face, i, dir1,dir2
        if (present(orientation)) then 
            face = orientation
        else
            do i = 1, size(contour)
                face = contour(i)%getFace()
                if (face /= NOT_ON_FACE) exit
            end do
        end if
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

    function getCellDistance(ref_cell, cell,edge) result(res)
        integer(kind=4), dimension(3), intent(in) :: ref_cell, cell
        integer, intent(in) :: edge
        integer :: res
        integer :: i 
        res = 0
        do i = 1, 3
            res = res + i*abs(ref_cell(i)-cell(i))
        end do
        if (edge == EDGE_X) then 
            if (res == 2) res = 1
            if (res == 3) res = 2
            if (res == 5) res = 3
        else if (edge == EDGE_Y) then 
            if (res == 3) res = 2
            if (res == 4) res = 3
        end if
        res = res + 1
    end function

    function getTrianglesOffFaces(triangles) result (res) !delete
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
 
    function getTrianglesOnFace(tris, face) result(res)
        type(triangle_t), dimension(:), allocatable, intent(in) :: tris
        integer, intent(in) :: face
        type(triangle_t), dimension(:), allocatable :: res
        integer :: i, n
        n = 0
        do i = 1, size(tris)
            if (tris(i)%isOnFace(face)) n = n + 1
        end do
        allocate(res(n))
        n = 0
        do i = 1, size(tris)
            if (tris(i)%isOnFace(face)) then 
                n = n + 1
                res(n) = tris(i)
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