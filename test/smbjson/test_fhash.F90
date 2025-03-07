integer function test_fhash_coords() bind(C) result(err)
    use cell_map_mod
    use fhash, only: fhash_tbl_t, key=>fhash_key
    use iso_fortran_env, only: int32, int64

    implicit none

    type(fhash_tbl_t) :: tbl
    type(triangle_t) :: t_set, t_get 
    integer (int32), dimension(3) :: cell_set
    class(*), allocatable :: t_alloc
    integer :: stat, n_cond
    err = 0

    t_set%vertices(1)%position = [0,0,0]
    t_set%vertices(2)%position = [0,1,0]
    t_set%vertices(3)%position = [1,0,1]
    cell_set = t_set%getCell()
    call tbl%set(key(cell_set), value = t_set)
    
    call tbl%get_raw(key(cell_set), t_alloc, stat)
    select type(t_alloc)
    type is (triangle_t)
        if (.not. all(t_alloc%vertices(1)%position .eq. [0,0,0])) err = err + 1
        if (.not. all(t_alloc%vertices(2)%position .eq. [0,1,0])) err = err + 1
        if (.not. all(t_alloc%vertices(3)%position .eq. [1,0,1])) err = err + 1
    end select

end function

integer function test_fhash_array() bind(C) result(err)
    use cell_map_mod
    use fhash, only: fhash_tbl_t, key=>fhash_key
    use iso_fortran_env, only: int32, int64

    implicit none
    
    type(triangle_map_t) :: tbl
    type(triangle_t) :: t1, t2
    type(triangle_t), dimension(:), allocatable :: t_array, t_array_aux
    type(triangle_set_t) :: tri_list, tri_list_aux
    integer (kind=4), dimension(3) :: cell_set
    class(*), allocatable :: t_alloc
    integer :: stat, n_cond
    err = 0

    

    t1%vertices(1)%position = [0,0,0]
    t1%vertices(2)%position = [0,1,0]
    t1%vertices(3)%position = [1,0,1]
    cell_set = t1%getCell()
    allocate(tri_list%triangles(1))
    tri_list%triangles(1) = t1
    call tbl%set(key(cell_set), value = tri_list)

    t1%vertices(1)%position = [0,0,0]
    t1%vertices(2)%position = [0,1,0]
    t1%vertices(3)%position = [1,0,0]
    cell_set = t1%getCell()

    if (.not. tbl%hasKey(cell_set)) err = err +1

    if (tbl%hasKey(cell_set)) then 
        call tbl%get_raw(key(cell_set), t_alloc, stat)
        select type(t_alloc)
        type is (triangle_set_t)
            if (size(t_alloc%triangles) /= 1) err = err +1

            allocate(tri_list_aux%triangles(size(t_alloc%triangles) + 1))
            
            tri_list_aux%triangles(1:size(t_alloc%triangles)) = t_alloc%triangles
            tri_list_aux%triangles(size(t_alloc%triangles) + 1) = t1
            
            deallocate(t_alloc%triangles)
            allocate(t_alloc%triangles(size(tri_list_aux%triangles)))
            t_alloc%triangles = tri_list_aux%triangles
            call tbl%set(key(cell_set), value = t_alloc)

        end select
    end if

    if (.not. tbl%hasKey(cell_set)) err = err +1
    if (tbl%hasKey(cell_set)) then 
        call tbl%get_raw(key(cell_set), t_alloc, stat)
        select type(t_alloc)
        type is (triangle_set_t)
            if (size(t_alloc%triangles) /= 2) err = err + 1
            write(*,*) t_alloc%triangles(1)%vertices
            if (.not. all(t_alloc%triangles(1)%vertices(1)%position .eq. [0,0,0])) err = err + 1
            if (.not. all(t_alloc%triangles(1)%vertices(2)%position .eq. [0,1,0])) err = err + 1
            if (.not. all(t_alloc%triangles(1)%vertices(3)%position .eq. [1,0,1])) err = err + 1

            if (.not. all(t_alloc%triangles(2)%vertices(1)%position .eq. [0,0,0])) err = err + 1
            if (.not. all(t_alloc%triangles(2)%vertices(2)%position .eq. [0,1,0])) err = err + 1
            if (.not. all(t_alloc%triangles(2)%vertices(3)%position .eq. [1,0,0])) err = err + 1


        end select
    end if

end function

integer function test_fhash_add_triangle() bind(C) result(err)
    use cell_map_mod
    use fhash, only: fhash_tbl_t, key=>fhash_key
    use iso_fortran_env, only: int32, int64

    implicit none
    
    type(triangle_map_t) :: tbl
    type(triangle_t) :: t1, t2
    type(triangle_t), dimension(:), allocatable :: t_array, t_array_aux
    type(triangle_set_t) :: tri_list, tri_list_aux
    integer (kind=4), dimension(3) :: cell_set
    class(*), allocatable :: t_alloc
    integer :: stat, n_cond
    type(cell_t), dimension(:), allocatable :: keys
    err = 0

    

    t1%vertices(1)%position = [0,0,0]
    t1%vertices(2)%position = [0,1,0]
    t1%vertices(3)%position = [1,0,1]
    cell_set = t1%getCell()
    call tbl%addTriangle(t1)
    
    if (.not. tbl%hasKey(cell_set)) err = err +1
    call tbl%get_raw(key(cell_set), t_alloc, stat)
    select type(t_alloc)
    type is (triangle_set_t)
        if (size(t_alloc%triangles) /= 1) err = err + 1
        write(*,*) t_alloc%triangles(1)%vertices
        if (.not. all(t_alloc%triangles(1)%vertices(1)%position .eq. [0,0,0])) err = err + 1
        if (.not. all(t_alloc%triangles(1)%vertices(2)%position .eq. [0,1,0])) err = err + 1
        if (.not. all(t_alloc%triangles(1)%vertices(3)%position .eq. [1,0,1])) err = err + 1
    end select
    
    t1%vertices(1)%position = [0,0,0]
    t1%vertices(2)%position = [0,1,0]
    t1%vertices(3)%position = [1,0,0]
    cell_set = t1%getCell()
    call tbl%addTriangle(t1)

    call tbl%get_raw(key(cell_set), t_alloc, stat)
    select type(t_alloc)
    type is (triangle_set_t)
        if (size(t_alloc%triangles) /= 2) err = err + 1
        write(*,*) t_alloc%triangles(1)%vertices
        if (.not. all(t_alloc%triangles(1)%vertices(1)%position .eq. [0,0,0])) err = err + 1
        if (.not. all(t_alloc%triangles(1)%vertices(2)%position .eq. [0,1,0])) err = err + 1
        if (.not. all(t_alloc%triangles(1)%vertices(3)%position .eq. [1,0,1])) err = err + 1

        if (.not. all(t_alloc%triangles(2)%vertices(1)%position .eq. [0,0,0])) err = err + 1
        if (.not. all(t_alloc%triangles(2)%vertices(2)%position .eq. [0,1,0])) err = err + 1
        if (.not. all(t_alloc%triangles(2)%vertices(3)%position .eq. [1,0,0])) err = err + 1
    end select

    if (size(tbl%keys) /= 1) err = err + 1
    if (.not. all(tbl%keys(1)%cell == [0,0,0])) err = err + 1
end function

integer function test_fhash_cellmap_set_get() bind(C) result(err)
    use cell_map_mod
    use geometry_mod
    implicit none 

    type(triangle_map_t) :: tri_map
    type(side_t), dimension(:), allocatable :: sides
    type(triangle_t), dimension(:), allocatable :: tri_set, tri_get
    type(triangle_t) :: t1, t2, t3
    type(coord_t) :: c1, c2, c3, c4, c5
    err = 0

    c1 = coord_t(position = [0,0,0], id = 1)
    c2 = coord_t(position = [1,0,0], id=  2)
    c3 = coord_t(position = [0,0,1], id=  3)
    c4 = coord_t(position = [0,1,1], id=  4)
    c5 = coord_t(position = [0,1,0], id=  5)
    t1 = triangle_t(vertices = [c1,c5,c3])
    t2 = triangle_t(vertices = [c1,c2,c5])
    t3 = triangle_t(vertices = [c1,c2,c4])

    allocate(tri_set(1))
    tri_set(1) = t1
    
    tri_map = buildTriangleMap(tri_set)
    tri_get = tri_map%getTrianglesInCell(floor(c1%position))
    if (size(tri_get) /= 1) err = err + 1
    if (tri_get(1)%vertices(1)%id /= tri_set(1)%vertices(1)%id) err = err + 1
    if (tri_get(1)%vertices(2)%id /= tri_set(1)%vertices(2)%id) err = err + 1
    if (tri_get(1)%vertices(3)%id /= tri_set(1)%vertices(3)%id) err = err + 1
    
    deallocate(tri_set)
    allocate(tri_set(2))
    tri_set(1) = t1
    tri_set(2) = t2

    tri_map = buildTriangleMap(tri_set)
    tri_get = tri_map%getTrianglesInCell(floor(c1%position))
    if (size(tri_get) /= 2) err = err + 1
    if (tri_get(1)%vertices(1)%id /= tri_set(1)%vertices(1)%id) err = err + 1
    if (tri_get(1)%vertices(2)%id /= tri_set(1)%vertices(2)%id) err = err + 1
    if (tri_get(1)%vertices(3)%id /= tri_set(1)%vertices(3)%id) err = err + 1
    if (tri_get(2)%vertices(1)%id /= tri_set(2)%vertices(1)%id) err = err + 1
    if (tri_get(2)%vertices(2)%id /= tri_set(2)%vertices(2)%id) err = err + 1
    if (tri_get(2)%vertices(3)%id /= tri_set(2)%vertices(3)%id) err = err + 1

end function
