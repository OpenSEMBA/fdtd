module cell_map_mod

    use geometry_mod, only: triangle_t, side_t, FACE_X, FACE_Y, FACE_Z
    use fhash, only: fhash_tbl_t, key=>fhash_key

    implicit none


    type :: triangle_set_t
        type(triangle_t), dimension(:), allocatable :: triangles ! triangles on faces
    end type

    type :: side_set_t
        type(side_t), dimension(:), allocatable :: sides ! sides from triangles off faces
    end type

    type :: cell_t
        integer, dimension(3) :: cell
     end type


    type, extends(fhash_tbl_t) :: cell_map_t
        type(cell_t), dimension(:), allocatable :: keys
    contains
        procedure :: hasKey
    end type

    type, extends(cell_map_t) :: triangle_map_t
    contains
        procedure :: addTriangle
        procedure :: getTrianglesInCell
    end type

    type, extends(cell_map_t) :: side_map_t
    contains
        procedure :: addSide
        procedure :: getSidesInCell
    end type

contains

    function buildTriangleMap(triangles) result(res)
        type(triangle_t), dimension(:), allocatable :: triangles
        type(side_t), dimension(3) :: sides
        type(triangle_map_t) :: res
        integer :: i, j
        integer (kind=4), dimension(3) :: cell
        do i = 1, size(triangles)
            if (triangles(i)%isOnAnyFace()) then 
                call res%addTriangle(triangles(i))
            end if
        end do
    end function

    function buildSideMap(triangles) result(res)
        type(triangle_t), dimension(:), allocatable :: triangles
        type(side_t), dimension(3) :: sides
        type(side_map_t) :: res
        integer :: i, j
        integer (kind=4), dimension(3) :: cell
        do i = 1, size(triangles)
            if (.not. triangles(i)%isOnAnyFace()) then 
                sides = triangles(i)%getSides()
                do j = 1, 3
                    if (sides(j)%isOnAnyFace() .or. & 
                        sides(j)%isOnAnyEdge()) then 
                            call res%addSide(sides(j))
                    end if
                end do
            end if
        end do
    end function

    logical function hasKey(this, k)
        class(cell_map_t) :: this
        integer(kind=4), dimension(3), intent(in) :: k
        integer :: stat
        hasKey = .false.
        call this%check_key(key(k), stat)
        if (stat == 0) hasKey = .true.
    end function

    subroutine addTriangle(this, triangle)
        class(triangle_map_t) :: this
        type(triangle_t) :: triangle
        class(*), allocatable :: alloc_list
        type(triangle_set_t) :: aux_list
        integer (kind=4), dimension(3) :: cell
        type(cell_t), dimension(:), allocatable :: aux_keys
        cell = triangle%getCell()
        if (this%hasKey(cell)) then 

            call this%get_raw(key(cell), alloc_list)
            select type(alloc_list)
            type is(triangle_set_t)
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

    subroutine addSide(this, side)
        class(side_map_t) :: this
        type(side_t) :: side
        class(*), allocatable :: alloc_list
        type(side_set_t) :: aux_list
        integer (kind=4), dimension(3) :: cell
        type(cell_t), dimension(:), allocatable :: aux_keys
        cell = side%getCell()
        if (this%hasKey(cell)) then 

            call this%get_raw(key(cell), alloc_list)
            select type(alloc_list)
            type is(side_set_t)
                if (allocated(alloc_list%sides)) then 
                    allocate(aux_list%sides(size(alloc_list%sides) + 1))
                    aux_list%sides(1:size(alloc_list%sides)) = alloc_list%sides
                    aux_list%sides(size(alloc_list%sides) + 1) = side
                    deallocate(alloc_list%sides)
                    allocate(alloc_list%sides(size(aux_list%sides)))
                    alloc_list%sides = aux_list%sides
                    call this%set(key(cell), value = alloc_list)
                else 
                    allocate(aux_list%sides(1))
                    aux_list%sides(1) = side
                    call this%set(key(cell), value = aux_list)
                end if
            end select

        else 
            allocate(aux_list%sides(1))
            aux_list%sides(1) = side
            call this%set(key(cell), value = aux_list)

            if (.not. allocated(this%keys)) allocate(this%keys(0))

            allocate(aux_keys(size(this%keys) + 1))
            aux_keys(1:size(this%keys)) = this%keys
            aux_keys(size(this%keys) + 1)%cell = side%getCell()
            deallocate(this%keys)
            allocate(this%keys(size(aux_keys)))
            this%keys = aux_keys
        end if

    end subroutine

    function getTrianglesInCell(this, k) result(res)
        class(triangle_map_t) :: this
        integer(kind=4), dimension(3) :: k
        class(*), allocatable :: alloc_list
        type(triangle_t), dimension(:), allocatable :: res

        if (this%hasKey(k)) then 
            call this%get_raw(key(k), alloc_list)
            select type(alloc_list)
            type is(triangle_set_t)
                allocate(res(size(alloc_list%triangles)))
                res = alloc_list%triangles
            end select
        else
            allocate(res(0))
        end if
    end function

    function getSidesInCell(this, k) result(res)
        class(side_map_t) :: this
        integer(kind=4), dimension(3) :: k
        class(*), allocatable :: alloc_list
        type(side_t), dimension(:), allocatable :: res

        if (this%hasKey(k)) then 
            call this%get_raw(key(k), alloc_list)
            select type(alloc_list)
            type is(side_set_t)
                res = alloc_list%sides
            end select
        else
            allocate(res(0))
        end if
    end function



end module