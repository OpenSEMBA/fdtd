module cell_map_mod

    use geometry_mod, only: triangle_t, side_t, FACE_X, FACE_Y, FACE_Z, isNewSide
    use fhash, only: fhash_tbl_t, key=>fhash_key

    implicit none


    type :: element_set_t
        type(triangle_t), dimension(:), allocatable :: triangles ! triangles on faces
        type(side_t), dimension(:), allocatable :: sides ! sides from triangles off faces
        type(side_t), dimension(:), allocatable :: sides_on  ! sides from triangles on faces
    end type

    type :: cell_t
        integer, dimension(3) :: cell
    end type

    type, public :: cfl_info_t
        real, dimension(3) :: area = [1,1,1], length = [1,1,1]
    end type

    type, extends(fhash_tbl_t) :: cfl_map_t
        type(cell_t), dimension(:), allocatable :: keys
    contains
        procedure :: hasKey => cfl_hasKey
        procedure :: addFaceRatio
        procedure :: addEdgeRatio
        procedure :: getCFLInCell
    end type

    type, extends(fhash_tbl_t) :: cell_map_t
        type(cell_t), dimension(:), allocatable :: keys
    contains
        procedure :: hasKey => cell_hasKey
        procedure :: getTrianglesInCell
        procedure :: getSidesInCell
        procedure :: getOnSidesInCell
    end type

    type, extends(cell_map_t) :: triangle_map_t
    contains
        procedure :: addTriangle
    end type

    type, extends(cell_map_t) :: side_map_t
    contains
        procedure :: addSide
        procedure :: addSideOn
    end type

contains

    subroutine buildCellMap(res, triangles)
        type(cell_map_t), intent(inout) :: res
        type(triangle_t), dimension(:), allocatable :: triangles
        type(triangle_map_t) :: tri_map
        type(side_map_t) :: side_map, side_map_on
        type(cell_t), dimension(:), allocatable :: keys
        type(element_set_t) :: elems
        integer :: i
        call buildMapOfTrisOnFaces(tri_map, triangles)
        call buildMapOfSidesFromTrisNotOnFaces(side_map, triangles)
        call buildMapOfSidesOnEdgeFromTrisOnFaces(side_map_on, triangles)
        keys = mergeKeys(tri_map%keys, side_map%keys)
        keys = mergeKeys(keys, side_map_on%keys)
        do i = 1, size(keys)
            elems%triangles = tri_map%getTrianglesInCell(keys(i)%cell)
            elems%sides = side_map%getSidesInCell(keys(i)%cell)
            elems%sides_on = side_map_on%getOnSidesInCell(keys(i)%cell)
            call res%set(key(keys(i)%cell), value=elems)
        end do
        res%keys = keys
    end subroutine

    function mergeKeys(tri_keys, side_keys) result(res)
        type(cell_t), dimension(:), allocatable, intent(in) :: tri_keys, side_keys
        type(cell_t), dimension(:), allocatable :: res, aux
        integer :: i
        if (size(tri_keys) == 0) then 
            allocate(res(0))
        else
            allocate(res(size(tri_keys)))
        end if
        res = tri_keys
        if (size(side_keys) /= 0) then 
            do i = 1, size(side_keys)
                if (isNewKey(tri_keys, side_keys(i))) call addKey(res, side_keys(i))
            end do
        end if
    end function

    subroutine addKey(keys, new_key)
        type(cell_t), dimension(:), allocatable, intent(inout) :: keys
        type(cell_t), intent(in) :: new_key
        type(cell_t), dimension(:), allocatable :: aux
        allocate(aux(size(keys) + 1))
        aux(1:size(keys)) = keys
        aux(size(keys) + 1) = new_key
        deallocate(keys)
        allocate(keys(size(aux)))
        keys = aux
    end subroutine

    logical function isNewKey(keys, k) 
        type(cell_t), dimension(:), allocatable, intent(in) :: keys
        type(cell_t), intent(in) :: k
        integer :: i
        isNewKey = .true.
        if (size(keys) /= 0) then 
            do i = 1, size(keys)
                if (all(keys(i)%cell .eq. k%cell)) isNewKey = .false.
            end do
        end if
    end function    

    subroutine buildMapOfTrisOnFaces(res, triangles)
        type(triangle_map_t), intent(inout) :: res
        type(triangle_t), dimension(:), allocatable :: triangles
        type(side_t), dimension(3) :: sides
        integer :: i, j
        integer (kind=4), dimension(3) :: cell
        if (.not. allocated(res%keys)) allocate(res%keys(0))
        do i = 1, size(triangles)
            if (triangles(i)%isOnAnyFace()) then 
                call res%addTriangle(triangles(i))
            end if
        end do
    end subroutine
    
    subroutine buildMapOfSidesFromTrisNotOnFaces(res, triangles)
        type(side_map_t), intent(inout) :: res
        type(triangle_t), dimension(:), allocatable :: triangles
        type(side_t), dimension(3) :: sides
        integer :: i, j
        integer (kind=4), dimension(3) :: cell
        if (.not. allocated(res%keys)) allocate(res%keys(0))
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
    end subroutine

    subroutine buildMapOfSidesOnEdgeFromTrisOnFaces(res, triangles)
        type(side_map_t), intent(inout) :: res
        type(triangle_t), dimension(:), allocatable :: triangles
        type(side_t), dimension(3) :: sides
        integer :: i, j
        integer (kind=4), dimension(3) :: cell
        if (.not. allocated(res%keys)) allocate(res%keys(0))
        do i = 1, size(triangles)
            if (triangles(i)%isOnAnyFace()) then 
                sides = triangles(i)%getSides()
                do j = 1, 3
                    if (sides(j)%isOnAnyEdge()) then 
                        call res%addSideOn(sides(j))
                    end if
                end do
            end if
        end do
    end subroutine

    logical function cell_hasKey(this, k)
        class(cell_map_t) :: this
        integer(kind=4), dimension(3), intent(in) :: k
        integer :: stat
        cell_hasKey = .false.
        call this%check_key(key(k), stat)
        if (stat == 0) cell_hasKey = .true.
    end function

    logical function cfl_hasKey(this, k)
        class(cfl_map_t) :: this
        integer(kind=4), dimension(3), intent(in) :: k
        integer :: stat
        cfl_hasKey = .false.
        call this%check_key(key(k), stat)
        if (stat == 0) cfl_hasKey = .true.
    end function

    subroutine addTriangle(this, triangle)
        class(triangle_map_t) :: this
        type(triangle_t) :: triangle
        class(*), allocatable :: alloc_list
        type(element_set_t) :: aux_list
        integer (kind=4), dimension(3) :: cell
        type(cell_t), dimension(:), allocatable :: aux_keys
        cell = triangle%getCell()
        if (this%hasKey(cell)) then 

            call this%get_raw(key(cell), alloc_list)
            select type(alloc_list)
            type is(element_set_t)
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
        type(element_set_t) :: aux_list
        integer (kind=4), dimension(3) :: cell
        type(cell_t), dimension(:), allocatable :: aux_keys
        cell = side%getCell()
        if (this%hasKey(cell)) then 

            call this%get_raw(key(cell), alloc_list)
            select type(alloc_list)
            type is(element_set_t)
            
                if (isNewSide(alloc_list%sides, side)) then 
                    allocate(aux_list%sides(size(alloc_list%sides) + 1))
                    aux_list%sides(1:size(alloc_list%sides)) = alloc_list%sides
                    aux_list%sides(size(alloc_list%sides) + 1) = side
                    deallocate(alloc_list%sides)
                    allocate(alloc_list%sides(size(aux_list%sides)))
                    alloc_list%sides = aux_list%sides
                    call this%set(key(cell), value = alloc_list)
                end if
            end select

        else 
            allocate(aux_list%sides(1))
            aux_list%sides(1) = side
            call this%set(key(cell), value = aux_list)

            allocate(aux_keys(size(this%keys) + 1))
            aux_keys(1:size(this%keys)) = this%keys
            aux_keys(size(this%keys) + 1)%cell = side%getCell()
            deallocate(this%keys)
            allocate(this%keys(size(aux_keys)))
            this%keys = aux_keys
        end if

    end subroutine

    subroutine addSideOn(this, side)
        class(side_map_t) :: this
        type(side_t) :: side
        class(*), allocatable :: alloc_list
        type(element_set_t) :: aux_list
        integer (kind=4), dimension(3) :: cell
        type(cell_t), dimension(:), allocatable :: aux_keys
        cell = side%getCell()
        if (this%hasKey(cell)) then 

            call this%get_raw(key(cell), alloc_list)
            select type(alloc_list)
            type is(element_set_t)
            
                allocate(aux_list%sides_on(size(alloc_list%sides_on) + 1))
                aux_list%sides_on(1:size(alloc_list%sides_on)) = alloc_list%sides_on
                aux_list%sides_on(size(alloc_list%sides_on) + 1) = side
                deallocate(alloc_list%sides_on)
                allocate(alloc_list%sides_on(size(aux_list%sides_on)))
                alloc_list%sides_on = aux_list%sides_on
                call this%set(key(cell), value = alloc_list)
            end select

        else 
            allocate(aux_list%sides_on(1))
            aux_list%sides_on(1) = side
            call this%set(key(cell), value = aux_list)

            allocate(aux_keys(size(this%keys) + 1))
            aux_keys(1:size(this%keys)) = this%keys
            aux_keys(size(this%keys) + 1)%cell = side%getCell()
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
            type is(element_set_t)
                allocate(res(size(alloc_list%triangles)))
                res = alloc_list%triangles
            end select
        else
            allocate(res(0))
        end if
    end function

    function getSidesInCell(this, k) result(res)
        class(cell_map_t) :: this
        integer(kind=4), dimension(3) :: k
        class(*), allocatable :: alloc_list
        type(side_t), dimension(:), allocatable :: res

        if (this%hasKey(k)) then 
            call this%get_raw(key(k), alloc_list)
            select type(alloc_list)
            type is(element_set_t)
                res = alloc_list%sides
            end select
        else
            allocate(res(0))
        end if
    end function


    function getOnSidesInCell(this, k) result(res)
        class(cell_map_t) :: this
        integer(kind=4), dimension(3) :: k
        class(*), allocatable :: alloc_list
        type(side_t), dimension(:), allocatable :: res

        if (this%hasKey(k)) then 
            call this%get_raw(key(k), alloc_list)
            select type(alloc_list)
            type is(element_set_t)
                res = alloc_list%sides_on
            end select
        else
            allocate(res(0))
        end if
    end function


    subroutine addFaceRatio(this, cell, direction, ratio)
        class(cfl_map_t) :: this
        class(*), allocatable :: alloc_list
        type(cfl_info_t) :: aux_cfl
        integer (kind=4), dimension(3) :: cell
        integer :: direction
        real :: ratio
        type(cell_t), dimension(:), allocatable :: aux_keys

        if (this%hasKey(cell)) then 

            call this%get_raw(key(cell), alloc_list)
            select type(alloc_list)
            type is(cfl_info_t)
                alloc_list%area(direction) = ratio
                call this%set(key(cell), value = alloc_list)
            end select

        else 
            aux_cfl%area(direction) = ratio
            call this%set(key(cell), value = aux_cfl)

            allocate(aux_keys(size(this%keys) + 1))
            aux_keys(1:size(this%keys)) = this%keys
            aux_keys(size(this%keys) + 1)%cell = cell
            deallocate(this%keys)
            allocate(this%keys(size(aux_keys)))
            this%keys = aux_keys
        end if

    end subroutine
    subroutine addEdgeRatio(this, cell, direction, ratio)
        class(cfl_map_t) :: this
        
        class(*), allocatable :: alloc_list
        type(cfl_info_t) :: aux_cfl
        integer (kind=4), dimension(3) :: cell
        integer :: direction
        real :: ratio

        type(cell_t), dimension(:), allocatable :: aux_keys
        if (this%hasKey(cell)) then 

            call this%get_raw(key(cell), alloc_list)
            select type(alloc_list)
            type is(cfl_info_t)
                alloc_list%length(direction) = ratio
                call this%set(key(cell), value = alloc_list)
            end select

        else 
            aux_cfl%length(direction) = ratio
            call this%set(key(cell), value = aux_cfl)

            allocate(aux_keys(size(this%keys) + 1))
            aux_keys(1:size(this%keys)) = this%keys
            aux_keys(size(this%keys) + 1)%cell = cell
            deallocate(this%keys)
            allocate(this%keys(size(aux_keys)))
            this%keys = aux_keys
        end if

    end subroutine

    function getCFLInCell(this, k) result(res)
        class(cfl_map_t) :: this
        integer(kind=4), dimension(3) :: k
        class(*), allocatable :: alloc_list
        type(cfl_info_t) :: res

        if (this%hasKey(k)) then 
            call this%get_raw(key(k), alloc_list)
            select type(alloc_list)
            type is(cfl_info_t)
                res = alloc_list
            end select
        end if
    end function


end module