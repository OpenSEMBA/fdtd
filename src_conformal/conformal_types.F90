module conformal_types_mod

    use cells_mod, only: cell_interval_t
    use fhash, only: fhash_tbl_t, key=>fhash_key

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

        procedure, private :: centroid
    end type
    
    type, public :: conformal_region_t
        type(triangle_t), dimension(:), allocatable :: triangles
        type(cell_interval_t), dimension(:), allocatable :: intervals
        integer :: offgrid_points
    end type
    
    type :: tri_list_t
        type(triangle_t), dimension(:), allocatable :: triangles
    end type

    type :: cell_t
        integer, dimension(3) :: cell
     end type

    type, extends(fhash_tbl_t) :: cell_map_t
        type(cell_t), dimension(:), allocatable :: keys
    contains
        procedure :: hasKey
        procedure :: addTriangle
        procedure :: getTrianglesInCell
    end type

contains

    logical function hasKey(this, k)
        class(cell_map_t) :: this
        integer(kind=4), dimension(3), intent(in) :: k
        integer :: stat
        hasKey = .false.
        call this%check_key(key(k), stat)
        if (stat == 0) hasKey = .true.
    end function

    subroutine addTriangle(this, triangle)
        class(cell_map_t) :: this
        type(triangle_t) :: triangle
        class(*), allocatable :: alloc_list
        type(tri_list_t) :: aux_list
        integer (kind=4), dimension(3) :: cell
        type(cell_t), dimension(:), allocatable :: aux_keys
        cell = triangle%getCell()
        if (this%hasKey(cell)) then 

            call this%get_raw(key(cell), alloc_list)
            select type(alloc_list)
            type is(tri_list_t)
                allocate(aux_list%triangles(size(alloc_list%triangles) + 1))
                aux_list%triangles(1:size(alloc_list%triangles)) = alloc_list%triangles
                aux_list%triangles(size(alloc_list%triangles) + 1) = triangle
                deallocate(alloc_list%triangles)
                allocate(alloc_list%triangles(size(aux_list%triangles)))
                alloc_list%triangles = aux_list%triangles
                call this%set(key(cell), value = alloc_list)
            end select

        else 
            allocate(aux_list%triangles(1))
            aux_list%triangles(1) = triangle
            call this%set(key(cell), value = aux_list)

            if (.not. allocated(this%keys)) allocate(this%keys(0))

            allocate(aux_keys(size(this%keys) + 1))
            aux_keys(1:size(this%keys)) = this%keys
            aux_keys(size(this%keys) + 1)%cell = triangle%getCell()
            deallocate(this%keys)
            allocate(this%keys(size(aux_keys)))
            this%keys = aux_keys
        end if

    end subroutine

    function getTrianglesInCell(this, k) result(res)
        class(cell_map_t) :: this
        integer(kind=4), dimension(3) :: k
        class(*), allocatable :: alloc_list
        type(triangle_t), dimension(:), allocatable :: res

        if (this%hasKey(k)) then 
            call this%get_raw(key(k), alloc_list)
            select type(alloc_list)
            type is(tri_list_t)
                res = alloc_list%triangles
            end select
        else
            allocate(res(0))
        end if
    end function

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

end module