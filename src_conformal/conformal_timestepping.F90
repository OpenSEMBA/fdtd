module conformal_timestepping_m
    use conformal_types_m, only: map_key_t
    use geometry_m, only: buildEdgesOnFace, buildFacesOnEdge
    use FDETypes_m, only: face_t, edge_t, FACE_X, FACE_Y, FACE_Z, EDGE_X, EDGE_Y, EDGE_Z, EDGE_REGION_CONFORMAL, EDGE_REGION_I, EDGE_REGION_II
    use NFDETypes_m, only: rkind
    use fhash, only: fhash_tbl_t, key=>fhash_key
    use Report_m

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
        procedure :: allocateEdges => edge_allocateEdges
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

    subroutine addAdditionalConformalFeatures(conf_maps, conformal_edges)
        type(conformal_maps_t), intent(inout) :: conf_maps
        type(edge_t), dimension(:), allocatable, intent(inout) :: conformal_edges
        type(edge_map_t) :: empty_map
        type(fhash_tbl_t) :: seen
        type(edge_t), dimension(:), allocatable :: aux_conformal_edges
        type(face_t), pointer :: face
        type(map_key_t), dimension(4) :: edges_on_face
        integer(kind=4), dimension(:,:), allocatable :: missing_keys
        logical, dimension(:), allocatable :: missing_in_region_II
        integer(kind=4) :: i, j, n_missing, prev_size, stat

        allocate(missing_keys(4, 4*size(conf_maps%face_map%keys)))
        allocate(missing_in_region_II(4*size(conf_maps%face_map%keys)))
        n_missing = 0
        do i = 1, size(conf_maps%face_map%keys)
            face => conf_maps%face_map%getFace(conf_maps%face_map%keys(i)%key)
            if (.not. face%is_two_sided) cycle
            edges_on_face = buildEdgesOnFace(face%cell, face%direction)
            do j = 1, size(edges_on_face)
                if (conf_maps%edge_map%hasKey(edges_on_face(j)%key)) cycle
                call seen%check_key(key(edges_on_face(j)%key), stat)
                if (stat == 0) cycle
                call seen%set(key(edges_on_face(j)%key), 1)
                n_missing = n_missing + 1
                missing_keys(:, n_missing) = edges_on_face(j)%key
                missing_in_region_II(n_missing) = isInRegionII(face, edges_on_face(j))
            end do
        end do

        prev_size = size(conformal_edges)
        allocate(aux_conformal_edges(prev_size + n_missing))
        aux_conformal_edges(1:prev_size) = conformal_edges
        do i = 1, n_missing
            call initAdditionalEdge(aux_conformal_edges(prev_size + i), &
                                    missing_keys(:, i), missing_in_region_II(i))
        end do
        call move_alloc(aux_conformal_edges, conformal_edges)

        ! Replaces the table, the finalizer frees the previous nodes. 
        ! The map must be rebuilt because conformal_edges has been reallocated. 
        conf_maps%edge_map = empty_map
        call conf_maps%edge_map%allocateEdges(4*size(conformal_edges) + 16)
        call buildConformalEdgeMap(conf_maps%edge_map, conformal_edges)
    end subroutine

    subroutine initAdditionalEdge(edge, edge_key, in_region_II)
        type(edge_t), intent(out) :: edge
        integer(kind=4), dimension(4), intent(in) :: edge_key
        logical, intent(in) :: in_region_II
        edge%cell = edge_key(1:3)
        edge%direction = edge_key(4)
        edge%material_coords = [0.0, 1.0]
        allocate(edge%region_I_fields%E)
        allocate(edge%region_II_fields%E)
        if (in_region_II) then
            edge%ratio = 0.0
            edge%edge_region = EDGE_REGION_II
            edge%region_I_fields%E => null()
            edge%region_II_fields%E = 0.0
        else
            edge%ratio = 1.0
            edge%edge_region = EDGE_REGION_I
            edge%region_I_fields%E = 0.0
            edge%region_II_fields%E => null()
        end if
    end subroutine

    logical function isInRegionII(face, edge_on_face)
        type(face_t), intent(in) :: face
        type(map_key_t) :: edge_on_face
        integer(kind=4) :: l
        isInRegionII = .false.
        do l = 1, size(face%contour_inside_edges)
            isInRegionII = isInRegionII .or. ((all(edge_on_face%key(1:3) == face%contour_inside_edges(l)%cell)) .and. &
                                              edge_on_face%key(4) == face%contour_inside_edges(l)%direction)
        end do
    end function

    subroutine assignEdgeFieldsOnFace(Ex, Ey, Ez, face, edge_on_face, j)
        real(kind=rkind), pointer, dimension(:,:,:) :: Ex, Ey, Ez
        real(kind=rkind), pointer :: E
        type(face_t), pointer :: face
        integer(kind=4), intent(in) :: j
        type(map_key_t), intent(in) ::edge_on_face
        integer(kind=4), dimension(3) :: c
        integer :: dir
        c = edge_on_face%key(1:3)
        select case (face%direction)
        case (FACE_X)
        if (mod(j,2)==0) then 
            E => Ey(c(1),c(2),c(3))
            dir = EDGE_Z
        else if (mod(j,2)/=0) then 
            E => Ez(c(1),c(2),c(3))
            dir = EDGE_Y
        end if
        case (FACE_Y)
        if (mod(j,2)==0) then 
            E => Ez(c(1),c(2),c(3))
            dir = EDGE_X
        else if (mod(j,2)/=0) then 
            E => Ex(c(1),c(2),c(3))
            dir = EDGE_Z
        end if
        case (FACE_Z)
        if (mod(j,2)==0) then 
            E => Ex(c(1),c(2),c(3))
            dir = EDGE_Y
        else if (mod(j,2)/=0) then 
            E => Ey(c(1),c(2),c(3))
            dir = EDGE_X
        end if
        end select

        if (isInRegionII(face,edge_on_face)) then 
            if (j == 1) then 
                allocate(face%region_I_fields%E1)
                face%region_I_fields%E1 = 0.0
                face%region_II_fields%E1 => E
            else if (j==2) then 
                allocate(face%region_I_fields%E2)
                face%region_I_fields%E2 = 0.0
                face%region_II_fields%E2 => E
            else if (j==3) then 
                allocate(face%region_I_fields%E3)
                face%region_I_fields%E3 = 0.0
                face%region_II_fields%E3 => E
            else if (j==4) then 
                allocate(face%region_I_fields%E4)
                face%region_I_fields%E4 = 0.0
                face%region_II_fields%E4 => E
            end if
        ! if (isInRegionI(face,edge_on_face)) then 
        else
            if (j == 1) then 
                face%region_I_fields%E1 => E
                allocate(face%region_II_fields%E1)
                face%region_II_fields%E1 = 0.0
            else if (j == 2) then 
                face%region_I_fields%E2 => E
                allocate(face%region_II_fields%E2)
                face%region_II_fields%E2 = 0.0
            else if (j == 3) then 
                face%region_I_fields%E3 => E
                allocate(face%region_II_fields%E3)
                face%region_II_fields%E3 = 0.0
            else if (j == 4) then 
                face%region_I_fields%E4 => E
                allocate(face%region_II_fields%E4)
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

    subroutine assignFaceFieldsOnEdge(Hx, Hy, Hz, edge, j, c)
        type(edge_t), pointer :: edge
        integer(kind=4), intent(in) :: j
        integer(kind=4), intent(in), dimension(3) :: c
        real(kind=rkind), pointer, dimension(:,:,:) :: Hx, Hy, Hz
        real(kind=rkind), pointer :: H
        integer :: dir
        select case (edge%direction)
        case (EDGE_X)
            if (mod(j,2)==0) then 
                H => Hy(c(1),c(2),c(3))
                dir = FACE_Z
            else if (mod(j,2)/=0) then 
                H => Hz(c(1),c(2),c(3))
                dir = FACE_Y
            end if
        case (EDGE_Y)
            if (mod(j,2)==0) then 
                H => Hz(c(1),c(2),c(3))
                dir = FACE_X
            else if (mod(j,2)/=0) then 
                H => Hx(c(1),c(2),c(3))
                dir = FACE_Z
            end if
        case (EDGE_Z)
            if (mod(j,2)==0) then 
                H => Hx(c(1),c(2),c(3))
                dir = FACE_Y
            else if (mod(j,2)/=0) then 
                H => Hy(c(1),c(2),c(3))
                dir = FACE_X
            end if
        end select

        if (edge%edge_region == EDGE_REGION_I) then 
            if (j == 1) then 
                edge%region_I_fields%H1  => H
                allocate(edge%region_II_fields%H1)
                edge%region_II_fields%H1 = 0.0
            else if (j == 2) then 
                edge%region_I_fields%H2  => H
                allocate(edge%region_II_fields%H2)
                edge%region_II_fields%H2 = 0.0
            else if (j == 3) then 
                edge%region_I_fields%H3  => H
                allocate(edge%region_II_fields%H3)
                edge%region_II_fields%H3 = 0.0
            else if (j == 4) then
                edge%region_I_fields%H4  => H
                allocate(edge%region_II_fields%H4)
                edge%region_II_fields%H4 = 0.0
            end if 
        else if (edge%edge_region == EDGE_REGION_II) then 
            if (j == 1) then 
                allocate(edge%region_I_fields%H1)
                edge%region_I_fields%H1 = 0.0
                edge%region_II_fields%H1 => H
            else if (j == 2) then 
                allocate(edge%region_I_fields%H2)
                edge%region_I_fields%H2 = 0.0
                edge%region_II_fields%H2 => H
            else if (j == 3) then 
                allocate(edge%region_I_fields%H3)
                edge%region_I_fields%H3 = 0.0
                edge%region_II_fields%H3 => H
            else if (j == 4) then
                allocate(edge%region_I_fields%H4)
                edge%region_I_fields%H4 = 0.0
                edge%region_II_fields%H4 => H
            end if 

        else if (edge%edge_region == EDGE_REGION_CONFORMAL) then 
            if (j == 1) then 
                edge%region_I_fields%H1  => H
                edge%region_II_fields%H1 => H
            else if (j == 2) then 
                edge%region_I_fields%H2  => H
                edge%region_II_fields%H2 => H
            else if (j == 3) then 
                edge%region_I_fields%H3  => H
                edge%region_II_fields%H3 => H
            else if (j == 4) then
                edge%region_I_fields%H4  => H
                edge%region_II_fields%H4 => H
            end if 

        end if
    end subroutine

    subroutine assignSplitFaceFieldsOnEdge(face, edge, j)
        type(face_t), pointer :: face
        type(edge_t), pointer :: edge
        integer(kind=4), intent(in) :: j

        if (edge%edge_region == EDGE_REGION_I) then 
            if (j == 1) then 
                edge%region_I_fields%H1  => face%region_I_fields%H
                allocate(edge%region_II_fields%H1)
                edge%region_II_fields%H1 = 0.0
            else if (j == 2) then 
                edge%region_I_fields%H2  => face%region_I_fields%H
                allocate(edge%region_II_fields%H2)
                edge%region_II_fields%H2 = 0.0
            else if (j == 3) then 
                edge%region_I_fields%H3  => face%region_I_fields%H
                allocate(edge%region_II_fields%H3)
                edge%region_II_fields%H3 = 0.0
            else if (j == 4) then 
                edge%region_I_fields%H4  => face%region_I_fields%H
                allocate(edge%region_II_fields%H4)
                edge%region_II_fields%H4 = 0.0
            end if
        else if (edge%edge_region == EDGE_REGION_II) then 
            if (j == 1) then 
                allocate(edge%region_I_fields%H1)
                edge%region_I_fields%H1  = 0.0
                edge%region_II_fields%H1 => face%region_II_fields%H
            else if (j == 2) then 
                allocate(edge%region_I_fields%H2)
                edge%region_I_fields%H2  = 0.0
                edge%region_II_fields%H2 => face%region_II_fields%H
            else if (j == 3) then 
                allocate(edge%region_I_fields%H3)
                edge%region_I_fields%H3  = 0.0
                edge%region_II_fields%H3 => face%region_II_fields%H
            else if (j == 4) then 
                allocate(edge%region_I_fields%H4)
                edge%region_I_fields%H4  = 0.0
                edge%region_II_fields%H4 => face%region_II_fields%H
            end if
        else if (edge%edge_region == EDGE_REGION_CONFORMAL) then
            if (j == 1) then 
                edge%region_I_fields%H1  => face%region_I_fields%H
                edge%region_II_fields%H1 => face%region_II_fields%H
            else if (j == 2) then 
                edge%region_I_fields%H2  => face%region_I_fields%H
                edge%region_II_fields%H2 => face%region_II_fields%H
            else if (j == 3) then 
                edge%region_I_fields%H3  => face%region_I_fields%H
                edge%region_II_fields%H3 => face%region_II_fields%H
            else if (j == 4) then 
                edge%region_I_fields%H4  => face%region_I_fields%H
                edge%region_II_fields%H4 => face%region_II_fields%H
            end if
        end if 
    end subroutine


    subroutine buildConformalMaps(face_map, edge_map, faces, edges)
        type(face_map_t), intent(inout) :: face_map
        type(edge_map_t), intent(inout) :: edge_map
        type(face_t), dimension(:), allocatable :: faces
        type(edge_t), dimension(:), allocatable :: edges
        integer :: i,j
        call buildConformalFaceMap(face_map, faces)
        call allocateConformalEdgeFields(edges)
        call buildConformalEdgeMap(edge_map, edges)
    end subroutine

    subroutine buildConformalFaceMap(face_map, faces)
        type(face_map_t), intent(inout) :: face_map
        type(face_t), dimension(:), allocatable :: faces
        integer :: i,j
        if (.not. allocated(face_map%keys)) allocate(face_map%keys(0))
        do i = 1, size(faces)
            allocate(faces(i)%region_I_fields%H)
            faces(i)%region_I_fields%H = 0.0
            allocate(faces(i)%region_II_fields%H)
            faces(i)%region_II_fields%H = 0.0
            call face_map%addFace(faces(i))
        end do
    end subroutine

    subroutine buildConformalEdgeMap(edge_map, edges)
        type(edge_map_t), intent(inout) :: edge_map
        type(edge_t), dimension(:), allocatable :: edges
        integer :: i
        if (.not. allocated(edge_map%keys)) allocate(edge_map%keys(0))
        do i = 1, size(edges)
            if (edges(i)%ratio /= -1) call edge_map%addEdge(edges(i))
        end do
    end subroutine

    subroutine allocateConformalEdgeFields(edges)
        type(edge_t), dimension(:), allocatable :: edges
        integer :: i
        do i = 1, size(edges)
            if (edges(i)%ratio /= -1) then
                allocate(edges(i)%region_I_fields%E)
                edges(i)%region_I_fields%E = 0.0
                allocate(edges(i)%region_II_fields%E)
                edges(i)%region_II_fields%E = 0.0
            end if
        end do
    end subroutine

    subroutine face_addFace(this, face)
        class(face_map_t) :: this
        type(face_t), intent(in), target :: face
        integer(kind=4), dimension(4) :: face_key
        type(map_key_t), dimension(:), allocatable :: aux_keys
        face_key(1:3) = face%cell
        face_key(4) = face%direction
        if (.not. allocated(this%keys)) allocate(this%keys(0))
        if (.not. this%hasKey(face_key)) then 
            call this%set_ptr(key(face_key), value = face)
            allocate(aux_keys(size(this%keys) + 1))
            aux_keys(1:size(this%keys)) = this%keys
            aux_keys(size(this%keys) + 1)%key = face_key
            deallocate(this%keys)
            allocate(this%keys(size(aux_keys)))
            this%keys = aux_keys
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
                res => alloc_val
            end select
        end if
    end function

    subroutine edge_addEdge(this, edge)
        class(edge_map_t) :: this
        type(edge_t), target, intent(in) :: edge
        integer(kind=4), dimension(4) :: edge_key
        type(map_key_t), dimension(:), allocatable :: aux_keys
        edge_key(1:3) = edge%cell
        edge_key(4) = edge%direction
        if (.not. allocated(this%keys)) allocate(this%keys(0))
        if (.not. this%hasKey(edge_key)) then 
            call this%set_ptr(key(edge_key), value = edge)
            allocate(aux_keys(size(this%keys) + 1))
            aux_keys(1:size(this%keys)) = this%keys
            aux_keys(size(this%keys) + 1)%key = edge_key
            deallocate(this%keys)
            allocate(this%keys(size(aux_keys)))
            this%keys = aux_keys
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
                res => alloc_val
            end select
        else
            write(*,*) 'key: ', k
            call WarnErrReport('Key not found in map ', .true.)
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

    subroutine addConformalFaces(faces, media_faces)
        type(face_t), dimension(:), allocatable, intent(inout) :: faces
        type(face_t), dimension(:), allocatable, intent(in) :: media_faces
        type(face_t), dimension(:), allocatable :: aux_faces
        integer :: i, n
        if (.not. allocated(faces)) allocate(faces(0))
        if (.not. allocated(media_faces)) return
        if (size(media_faces) == 0) return
        n = 0
        do i = 1, size(media_faces)
            if (media_faces(i)%is_two_sided) n = n + 1
        end do

        allocate(aux_faces(size(faces) + n))
        aux_faces(1:size(faces)) = faces
        n = 0
        do i = 1, size(media_faces)
            if (media_faces(i)%is_two_sided) then 
                n = n + 1
                aux_faces(size(faces) + n) = media_faces(i)
            end if
        end do
        deallocate(faces)
        allocate(faces(size(aux_faces)))
        faces = aux_faces
    end subroutine

    subroutine addConformalEdges(edges, media_edges)
        type(edge_t), dimension(:), allocatable, intent(inout) :: edges
        type(edge_t), dimension(:), allocatable, intent(in) :: media_edges
        type(edge_t), dimension(:), allocatable :: aux_edges
        integer :: i, n
        if (.not. allocated(edges)) allocate(edges(0))
        if (.not. allocated(media_edges)) return
        if (size(media_edges) == 0) return
        n = 0
        do i = 1, size(media_edges) 
            if (media_edges(i)%ratio /= 0 .and. media_edges(i)%is_two_sided) n = n + 1
        end do
        allocate(aux_edges(size(edges) + n))
        aux_edges(1:size(edges)) = edges
        n = 0
        do i = 1, size(media_edges)
            if (media_edges(i)%ratio /= 0 .and. media_edges(i)%is_two_sided) then 
                n = n + 1
                aux_edges(size(edges) + n) = media_edges(i)
            end if
        end do
        deallocate(edges)
        allocate(edges(size(aux_edges)))
        edges = aux_edges
    end subroutine

    subroutine edge_allocateEdges(this, buck)
        class(edge_map_t) :: this
        integer :: buck
        call this%allocate(buck)
    end subroutine

    function computeEdgeRegionIICorrection(edge) result(res)
        type(edge_t), intent(in) :: edge
        real(kind=rkind) :: res
        res = 1.0
        if (edge%edge_region == EDGE_REGION_CONFORMAL) then 
            res = (1.0-edge%ratio)/edge%ratio
        end if
    end function

    function computeFaceRegionIICorrection(face) result(res)
        type(face_t), intent(in) :: face
        real(kind=rkind) :: res
        res = face%ratio/(1.0-face%ratio)
    end function


end module