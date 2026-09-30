module conformal_timestepping_m
    use conformal_types_m, only: map_key_t
    use geometry_m, only: buildEdgesOnFace, buildFacesOnEdge
    use FDETypes_m, only: face_t, edge_t, FACE_X, FACE_Y, FACE_Z, EDGE_X, EDGE_Y, EDGE_Z
    use NFDETypes_m, only: rkind
    use fhash, only: fhash_tbl_t, key=>fhash_key


    type, extends(fhash_tbl_t) :: face_map_t
        type(map_key_t), dimension(:), allocatable :: keys
    contains
        procedure :: hasKey  => face_hasKey
        procedure :: getFace => face_getFace
        procedure :: addFace => face_addFace
    end type

    type, extends(fhash_tbl_t) :: edge_map_t
        type(map_key_t), dimension(:), allocatable :: keys
    contains
        procedure :: hasKey  => edge_hasKey
        procedure :: getEdge => edge_getEdge
        procedure :: addEdge => edge_addEdge
    end type

    type, public :: conformal_maps_t
        type(face_map_t) :: face_map
        type(edge_map_t) :: edge_map
    end type

    type, public :: conformal_fields_t
        type(face_t), dimension(:), allocatable :: faces
        type(edge_t), dimension(:), allocatable :: edges
    end type

contains

    subroutine addNewEdge(edge, edge_map, conformal_edges)
        type(edge_t), intent(in) :: edge
        type(edge_map_t), intent(inout) :: edge_map
        type(edge_t), dimension(:), allocatable, intent(inout) :: conformal_edges
        type(edge_t), dimension(:), allocatable :: aux_edges
        allocate(aux_edges(size(conformal_edges) + 1))
        aux_edges(1:size(conformal_edges)) = conformal_edges
        aux_edges(size(conformal_edges) + 1) = edge
        deallocate(conformal_edges)
        allocate(conformal_edges(size(aux_edges)))
        conformal_edges = aux_edges
        call edge_map%addEdge(conformal_edges(size(conformal_edges)))
    end subroutine

    subroutine addAdditionalConformalFeatures(face_map, edge_map, conformal_edges)
        type(edge_map_t), intent(inout) :: edge_map
        type(face_map_t), intent(inout) :: face_map
        type(edge_t), dimension(:), allocatable, intent(inout) :: conformal_edges
        type(face_t), pointer :: face
        type(edge_t), pointer :: edge
        type(edge_t), target :: new_edge
        type(map_key_t), dimension(4) :: edges_on_face
        integer(kind=4) :: i, j, k

        do i = 1, size(face_map%keys)
            face => face_map%getFace(face_map%keys(i)%key)
            if (.not. face%is_two_sided) cycle
            edges_on_face = buildEdgesOnFace(face%cell, face%direction)
            do j = 1, size(edges_on_face) 
                if (.not. edge_map%hasKey(edges_on_face(j)%key)) then 
                    new_edge = edge_t(cell=edges_on_face(j)%key(1:3), & 
                                      ratio=1.0, & 
                                      direction=edges_on_face(j)%key(4), &
                                      material_coords = [0.0,1.0])
                    call addNewEdge(new_edge, edge_map, conformal_edges)
                    edge => new_edge
                else
                    edge=> edge_map%getEdge(edges_on_face(j)%key)
                end if
            end do
        end do
    end subroutine


    subroutine assignEdgeFieldsOnFace(Ex, Ey, Ez, face, j, cell)
        type(face_t), pointer :: face
        integer(kind=4), intent(in) :: j
        integer(kind=4), dimension(3), intent(in) :: cell
        integer(kind=4), dimension(3) :: c
        real(kind=rkind), pointer, dimension(:,:,:) :: Ex, Ey, Ez, E
        integer :: dir
        c = cell
        select case (face%direction)
        case (FACE_X)
        if (mod(j,2)==0) then 
            E => Ey
            dir = EDGE_Z
        else if (mod(j,2)/=0) then 
            E => Ez
            dir = EDGE_Y
        end if
        case (FACE_Y)
        if (mod(j,2)==0) then 
            E => Ez
            dir = EDGE_X
        else if (mod(j,2)/=0) then 
            E => Ex
            dir = EDGE_Z
        end if
        case (FACE_Z)
        if (mod(j,2)==0) then 
            E => Ex
            dir = EDGE_Y
        else if (mod(j,2)/=0) then 
            E => Ey
            dir = EDGE_X
        end if
        end select

        if (j==1) then 
        if (face%normal(dir) > 0) then 
            face%region_I_fields%E1 = 0.0
            face%region_II_fields%E1 => E(c(1),c(2),c(3))
        else if (face%normal(dir) < 0) then 
            face%region_I_fields%E1 => E(c(1),c(2),c(3))
            face%region_II_fields%E1 = 0.0
        end if
        else if (j == 2) then 
        c(dir) = cell(dir) + 1
        if (face%normal(dir) > 0) then 
            face%region_I_fields%E2=> E(c(1),c(2),c(3))
            face%region_II_fields%E2 = 0.0
        else if (face%normal(dir) < 0) then 
            face%region_I_fields%E2 = 0.0
            face%region_II_fields%E2 => E(c(1),c(2),c(3))
        end if
        else if (j == 3) then 
        c(dir) = cell(dir) + 1
        if (face%normal(dir) > 0) then 
            face%region_I_fields%E3 => E(c(1),c(2),c(3))
            face%region_II_fields%E3 = 0.0
        else if (face%normal(dir) < 0) then 
            face%region_I_fields%E3 = 0.0
            face%region_II_fields%E3 => E(c(1),c(2),c(3))
        end if
        else if (j == 4) then 
        if (face%normal(dir) > 0) then 
            face%region_I_fields%E4 = 0.0
            face%region_II_fields%E4 => E(c(1),c(2),c(3))
        else if (face%normal(dir) < 0) then 
            face%region_I_fields%E4 => E(c(1),c(2),c(3))
            face%region_II_fields%E4 = 0.0
        end if
        end if

    end subroutine

    subroutine assignSplitEdgeFieldsOnFace(face, edge, j)
        type(face_t), pointer :: face
        type(edge_t), pointer :: edge
        integer(kind=4), intent(in) :: j
        if (j==1) then 
            face%region_I_fields%E1 => edge%region_I_fields%E
            face%region_II_fields%E1 => edge%region_II_fields%E
        else if (j==2) then 
            face%region_I_fields%E2 => edge%region_I_fields%E
            face%region_II_fields%E2 => edge%region_II_fields%E
        else if (j==3) then 
            face%region_I_fields%E3 => edge%region_I_fields%E
            face%region_II_fields%E3 => edge%region_II_fields%E
        else if (j==4) then 
            face%region_I_fields%E4 => edge%region_I_fields%E
            face%region_II_fields%E4 => edge%region_II_fields%E
        end if
    end subroutine

    subroutine assignFaceFieldsOnEdge(Hx, Hy, Hz, edge, j, cell)
        type(edge_t), pointer :: edge
        integer(kind=4), intent(in) :: j
        integer(kind=4), dimension(3), intent(in) :: cell
        integer(kind=4), dimension(3) :: c
        real(kind=rkind), pointer, dimension(:,:,:) :: Hx, Hy, Hz, H
        integer :: dir
        c = cell
        select case (edge%direction)
        case (EDGE_X)
            if (mod(j,2)==0) then 
                H => Hy
                dir = FACE_Z
            else if (mod(j,2)/=0) then 
                H => Hz
                dir = FACE_Y
            end if
        case (EDGE_Y)
            if (mod(j,2)==0) then 
                H => Hz
                dir = FACE_X
            else if (mod(j,2)/=0) then 
                H => Hx
                dir = FACE_Z
            end if
        case (EDGE_Z)
            if (mod(j,2)==0) then 
                H => Hx
                dir = FACE_Y
            else if (mod(j,2)/=0) then 
                H => Hy
                dir = FACE_X
            end if
        end select
        if (k==1) then 
            if (edge%ratio == 0.0) then 
                edge%region_II_fields%H1 => H(c(1),c(2),c(3))
                edge%region_I_fields%H1 = 0.0
            else if (edge%ratio == 1.0) then 
                edge%region_I_fields%H1  => H(c(1),c(2),c(3))
                edge%region_II_fields%H1 = 0.0
            else
                edge%region_I_fields%H1  => H(c(1),c(2),c(3))
                edge%region_II_fields%H1 => H(c(1),c(2),c(3))
            end if
        else if (k==2) then 
            c(dir) = cell(dir) - 1
            if (edge%ratio == 0.0) then 
                edge%region_II_fields%H2 => H(c(1),c(2),c(3)-1)
                edge%region_I_fields%H2 = 0.0
            else if (edge%ratio == 1.0) then 
                edge%region_I_fields%H2  => H(c(1),c(2),c(3)-1)
                edge%region_II_fields%H2 = 0.0
            else
                edge%region_I_fields%H2  => H(c(1),c(2),c(3)-1)
                edge%region_II_fields%H2 => H(c(1),c(2),c(3)-1)
            end if
        else if (k==3) then 
            c(dir) = cell(dir) - 1
            if (edge%ratio == 0.0) then 
                edge%region_II_fields%H3 => H(c(1),c(2)-1,c(3))
                edge%region_I_fields%H3 = 0.0
            else if (edge%ratio == 1.0) then 
                edge%region_I_fields%H3 => H(c(1),c(2)-1,c(3))
                edge%region_II_fields%H3 = 0.0
            else
                edge%region_I_fields%H3  => H(c(1),c(2)-1,c(3))
                edge%region_II_fields%H3 => H(c(1),c(2)-1,c(3))
            end if
        else if (k==4) then 
            if (edge%ratio == 0.0) then 
                edge%region_II_fields%H4 => H(c(1),c(2),c(3))
                edge%region_I_fields%H4 = 0.0
            else if (edge%ratio == 1.0) then 
                edge%region_I_fields%H4  => H(c(1),c(2),c(3))
                edge%region_II_fields%H4 = 0.0
            else
                edge%region_I_fields%H4   => H(c(1),c(2),c(3))
                edge%region_II_fields%H4  => H(c(1),c(2),c(3))
            end if
        end if


    end subroutine

    subroutine assignSplitFaceFieldsOnEdge(face, edge, j)
        type(face_t), pointer :: face
        type(edge_t), pointer :: edge
        integer(kind=4), intent(in) :: j
        if (j==1) then 
            if (edge%ratio == 0.0) then 
                edge%region_I_fields%H1  = 0.0
                edge%region_II_fields%H1 => face%region_II_fields%H
            else if (edge%ratio == 1.0) then 
                edge%region_I_fields%H1  => face%region_I_fields%H
                edge%region_II_fields%H1 = 0.0
            else 
                edge%region_I_fields%H1  => face%region_I_fields%H
                edge%region_II_fields%H1 => face%region_II_fields%H
            end if
        else if (j==2) then 
            if (edge%ratio == 0.0) then 
                edge%region_I_fields%H2  = 0.0
                edge%region_II_fields%H2 => face%region_II_fields%H
            else if (edge%ratio == 1.0) then 
                edge%region_I_fields%H2  => face%region_I_fields%H
                edge%region_II_fields%H2 = 0.0
            else
                edge%region_I_fields%H2  => face%region_I_fields%H
                edge%region_II_fields%H2 => face%region_II_fields%H
            end if
        else if (j==3) then 
            if (edge%ratio == 0.0) then 
                edge%region_I_fields%H3 = 0.0
                edge%region_II_fields%H3 => face%region_II_fields%H
            else if (edge%ratio == 1.0) then 
                edge%region_I_fields%H3  => face%region_I_fields%H
                edge%region_II_fields%H3 = 0.0
            else 
                edge%region_I_fields%H3  => face%region_I_fields%H
                edge%region_II_fields%H3 => face%region_II_fields%H
            end if
        else if (j==4) then 
            if (edge%ratio == 0.0) then 
                edge%region_I_fields%H4 = 0.0
                edge%region_II_fields%H4 => face%region_II_fields%H
            else if (edge%ratio == 1.0) then 
                edge%region_I_fields%H4  => face%region_I_fields%H
                edge%region_II_fields%H4 = 0.0
            else
                edge%region_I_fields%H4  => face%region_I_fields%H
                edge%region_II_fields%H4 => face%region_II_fields%H
            end if
        end if
    end subroutine


    subroutine buildConformalMaps(face_map, edge_map, faces, edges)
        type(face_map_t), intent(inout) :: face_map
        type(edge_map_t), intent(inout) :: edge_map
        ! type(edges_on_face_map_t), intent(inout) :: edges_on_face_map
        type(face_t), dimension(:), allocatable :: faces
        type(edge_t), dimension(:), allocatable :: edges
        ! type(map_key_t), dimension(4) :: edges_on_face
        integer :: i,j
        if (.not. allocated(face_map%keys)) allocate(face_map%keys(0))
        if (.not. allocated(edge_map%keys)) allocate(edge_map%keys(0))
        ! if (.not. allocated(edges_on_face_map%keys)) allocate(edges_on_face_map%keys(0))
        do i = 1, size(faces)
            call face_map%addFace(faces(i))
            ! edges_on_face = buildEdgesOnFace(faces(i))
            ! call edges_on_face_map%addEdges(edges_on_face)
        end do
        do i = 1, size(edges)
            call edge_map%addEdge(edges(i))
            ! edges_on_face = buildEdgesOnFace(faces(i))
            ! call edges_on_face_map%addEdges(edges_on_face)
        end do
    end subroutine


    subroutine face_addFace(this, face)
        class(face_map_t) :: this
        type(face_t), intent(in), target :: face
        integer(kind=4), dimension(4) :: face_key
        ! type(map_key_t), dimension(:), allocatable :: aux_keys
        face_key(1:3) = face%cell
        face_key(4) = face%direction
        if (.not. this%hasKey(face_key)) then 
            call this%set_ptr(key(face_key), value = face)
            ! allocate(aux_keys(size(this%keys) + 1))
            ! aux_keys(1:size(this%keys)) = this%keys
            ! aux_keys(size(this%keys) + 1)%key = face_key
            ! deallocate(this%keys)
            ! allocate(this%keys(size(aux_keys)))
            ! this%keys = aux_keys
        end if
    end subroutine

    function face_getFace(this, k, found) result(res)
        class(face_map_t) :: this
        integer(kind=4), dimension(4) :: k
        logical, intent(inout), optional :: found
        class(*), pointer :: alloc_val
        type(face_t), pointer :: res
        integer :: stat
        if (present(found)) found = .false.
        if (this%hasKey(k)) then 
            call this%get_raw_ptr(key = key(k), value = alloc_val, stat = stat)
            select type(alloc_val)
            type is(face_t)
                if (present(found)) found = .true.
                res = alloc_val
            end select
        end if
    end function

    subroutine edge_addEdge(this, edge)
        class(edge_map_t) :: this
        type(edge_t), target, intent(in) :: edge
        integer(kind=4), dimension(4) :: edge_key
        ! type(map_key_t), dimension(:), allocatable :: aux_keys
        edge_key(1:3) = edge%cell
        edge_key(4) = edge%direction
        if (.not. this%hasKey(edge_key)) then 
            call this%set_ptr(key(edge_key), value = edge)
            ! allocate(aux_keys(size(this%keys) + 1))
            ! aux_keys(1:size(this%keys)) = this%keys
            ! aux_keys(size(this%keys) + 1)%key = edge_key
            ! deallocate(this%keys)
            ! allocate(this%keys(size(aux_keys)))
            ! this%keys = aux_keys
        end if
    end subroutine

    function edge_getedge(this, k, found) result(res)
        class(edge_map_t) :: this
        integer(kind=4), dimension(4) :: k
        logical, intent(inout), optional :: found
        class(*), pointer :: alloc_val
        type(edge_t), pointer :: res
        integer :: stat
    
        if (present(found)) found = .false.
        if (this%hasKey(k)) then 
            call this%get_raw_ptr(key = key(k), value = alloc_val, stat = stat)
            select type(alloc_val)
            type is(edge_t)
                if (present(found)) found = .true.
                res = alloc_val
            end select
        end if
    end function

    logical function face_hasKey(this, k)
        class(face_map_t) :: this
        integer(kind=4), dimension(4), intent(in) :: k
        integer :: stat
        face_hasKey = .false.
        call this%check_key(key(k), stat)
        if (stat == 0) face_hasKey = .true.
    end function

    logical function edge_hasKey(this, k)
        class(edge_map_t) :: this
        integer(kind=4), dimension(4), intent(in) :: k
        integer :: stat
        edge_hasKey = .false.
        call this%check_key(key(k), stat)
        if (stat == 0) edge_hasKey = .true.
    end function

end module