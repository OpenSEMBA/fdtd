integer function test_triangle_normal() bind(C) result(err)

    use conformal_types_mod
    implicit none
    type(triangle_t) :: t
    err = 0

    
    t%vertices(1)%position = [0,0,0]; t%vertices(2)%position = [0,0,1]; t%vertices(3)%position = [0,1,0]
    if (t%getFace() /= FACE_X) err = err + 1
    
    t%vertices(1)%position = [0,0,0]; t%vertices(2)%position = [0,0,1]; t%vertices(3)%position = [1,0,0]
    if (t%getFace() /= FACE_Y) err = err + 1
    
    t%vertices(1)%position = [0,0,0]; t%vertices(2)%position = [0,1,0]; t%vertices(3)%position = [1,0,0]
    if (t%getFace() /= FACE_Z) err = err + 1

    t%vertices(1)%position = [0,0,0]; t%vertices(2)%position = [0,1,0]; t%vertices(3)%position = [1,0,1]
    if (t%getFace() /= NOT_ON_FACE) err = err + 1


end function

integer function test_triangle_edges() bind(C) result(err)

    use conformal_types_mod
    implicit none
    type(triangle_t) :: t
    type(side_t), dimension(3) :: sides
    err = 0
    
    t%vertices(1)%position = [0,0,0]; t%vertices(2)%position = [0,0,1]; t%vertices(3)%position = [0,1,0]
    sides = t%getSides()
    if (sides(1)%getEdge() /= EDGE_Z) err = err + 1
    if (sides(2)%getEdge() /= NOT_ON_EDGE) err = err + 1
    if (sides(3)%getEdge() /= EDGE_Y) err = err + 1
    
    t%vertices(1)%position = [0,0,0]; t%vertices(2)%position = [0,0,1]; t%vertices(3)%position = [1,0,0]
    sides = t%getSides()
    if (sides(1)%getEdge() /= EDGE_Z) err = err + 1
    if (sides(2)%getEdge() /= NOT_ON_EDGE) err = err + 1
    if (sides(3)%getEdge() /= EDGE_X) err = err + 1
    
    t%vertices(1)%position = [0,0,0]; t%vertices(2)%position = [0,1,0]; t%vertices(3)%position = [1,0,0]
    sides = t%getSides()
    if (sides(1)%getEdge() /= EDGE_Y) err = err + 1
    if (sides(2)%getEdge() /= NOT_ON_EDGE) err = err + 1
    if (sides(3)%getEdge() /= EDGE_X) err = err + 1

    t%vertices(1)%position = [0,0,0]; t%vertices(2)%position = [0,1,0]; t%vertices(3)%position = [1,0,1]
    sides = t%getSides()
    if (sides(1)%getEdge() /= EDGE_Y) err = err + 1
    if (sides(2)%getEdge() /= NOT_ON_EDGE) err = err + 1
    if (sides(3)%getEdge() /= NOT_ON_EDGE) err = err + 1


end function

integer function test_triangle_cell() bind(C) result(err)

    use conformal_types_mod
    implicit none
    type(triangle_t) :: t
    type(side_t), dimension(3) :: sides
    integer, dimension(3) :: cell
    err = 0
    
    t%vertices(1)%position = [0,0,0]; t%vertices(2)%position = [0,0,1]; t%vertices(3)%position = [0,1,0]
    cell = [0,0,0]
    if (all(t%getCell() .eq. cell) .eqv. .false.) err = err + 1
    
    t%vertices(1)%position = [1,0,0]; t%vertices(2)%position = [1,0,1]; t%vertices(3)%position = [1,1,0]
    cell = [1,0,0]
    if (all(t%getCell() .eq. cell) .eqv. .false.) err = err + 1
    
    t%vertices(1)%position = [1,0,1]; t%vertices(2)%position = [1,0,2]; t%vertices(3)%position = [1,1,1]
    cell = [1,0,1]
    if (all(t%getCell() .eq. cell) .eqv. .false.) err = err + 1

    t%vertices(1)%position = [0,0,1]; t%vertices(2)%position = [1,0,1]; t%vertices(3)%position = [1,1,1]
    cell = [0,0,1]
    if (all(t%getCell() .eq. cell) .eqv. .false.) err = err + 1


end function

integer function test_fhash_coords() bind(C) result(err)
    use conformal_types_mod
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
    use conformal_types_mod
    use fhash, only: fhash_tbl_t, key=>fhash_key
    use iso_fortran_env, only: int32, int64

    implicit none
    
    type(cell_map_t) :: tbl
    type(triangle_t) :: t1, t2
    type(triangle_t), dimension(:), allocatable :: t_array, t_array_aux
    type(tri_list_t) :: tri_list, tri_list_aux
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
        type is (tri_list_t)
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
        type is (tri_list_t)
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
    use conformal_types_mod
    use fhash, only: fhash_tbl_t, key=>fhash_key
    use iso_fortran_env, only: int32, int64

    implicit none
    
    type(cell_map_t) :: tbl
    type(triangle_t) :: t1, t2
    type(triangle_t), dimension(:), allocatable :: t_array, t_array_aux
    type(tri_list_t) :: tri_list, tri_list_aux
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
    type is (tri_list_t)
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
    type is (tri_list_t)
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