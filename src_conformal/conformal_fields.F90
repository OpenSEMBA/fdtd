module conformal_fields_m
    use conformal_types_m, only: map_key_t
    use geometry_m, only: buildEdgesOnFace, buildFacesOnEdge
    use FDETypes_m, only: face_t, edge_t, FACE_X, FACE_Y, FACE_Z, EDGE_X, EDGE_Y, EDGE_Z, EDGE_REGION_CONFORMAL, EDGE_REGION_I, EDGE_REGION_II, MediaData_t, SGGFDTDINFO_t, sim_control_t, SGBCMaterialProfile_t, RKIND_TIME
    use conformal_sgbc_m, only: sgbc_t
    use NFDETypes_m, only: rkind
    use fhash, only: fhash_tbl_t, key=>fhash_key
    use Report_m

    implicit none

    type, extends(fhash_tbl_t) :: sgbc_map_t
        type(map_key_t), dimension(:), allocatable :: keys
    contains
        procedure :: hasKey  => sgbc_hasKey
        procedure :: getSGBC => sgbc_getSGBC
        procedure :: addSGBC => sgbc_addSGBC
    end type

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
        type(sgbc_map_t) :: sgbc_map
    end type

    type, public :: conformal_fields_t
        type(face_t), dimension(:), allocatable :: faces
        type(edge_t), dimension(:), allocatable :: edges
        type(sgbc_t), dimension(:), allocatable :: sgbcs
        integer, dimension(:,:), allocatable :: edges_id_list
        integer, dimension(:,:), allocatable :: faces_id_list
        integer, dimension(:,:), allocatable :: sgbcs_id_list
    end type



contains

    subroutine initializeConformalStates(sgg, conformal_fields, Ex, Ey, Ez, Hx, Hy, Hz, idx, idy, idz, control)
        type(SGGFDTDINFO_t), intent(in) :: sgg
        type(conformal_fields_t), intent(inout) :: conformal_fields
        real(kind=rkind), pointer, dimension(:,:,:) :: Ex, Ey, Ez, Hx, Hy, Hz
        real(kind=rkind), dimension(:), pointer :: idx, idy, idz
        type(sim_control_t), intent(in) :: control
        type(face_t), pointer :: face, face_of_edge
        type(edge_t), pointer :: edge
        type(sgbc_t), pointer :: sgbc
        type(conformal_maps_t) :: conformal_maps
        type(map_key_t), dimension(4) :: edges_on_face, faces_on_edge
        integer(kind=4), dimension(3) :: cell
        integer :: i,j,k,n, n_sgbc

        n = 0
        n_sgbc = 0
        do i = 1, sgg%NumMedia
            if (isConformalPECSurface(sgg%Med(i))) n = n + 1
            if (isConformalSGBCSurface(sgg%Med(i))) n_sgbc = n_sgbc + 1
        end do
        if (n == 0) then 
            allocate(conformal_fields%faces(0))
            allocate(conformal_fields%edges(0))
            allocate(conformal_fields%edges_id_list(3,0), conformal_fields%faces_id_list(3,0))            
        end if
        if (n_sgbc == 0) then 
            allocate(conformal_fields%sgbcs(0))
            allocate(conformal_fields%sgbcs_id_list(3,0)) 
        end if
        if (n == 0 .and. n_sgbc == 0) return 
!!!
        do i = 1, sgg%NumMedia
            if (.not. (isConformalSurface(sgg%Med(i)))) cycle
            call addConformalFaces(conformal_fields%faces, sgg%Med(i)%ConformalFace)
            call addConformalEdges(conformal_fields%edges, sgg%Med(i)%ConformalEdge)
            if (sgg%Med(i)%is%ConformalSGBC) call addConformalSGBCs(conformal_fields%sgbcs, sgg%Med(i)%ConformalFace)
        end do

        call buildConformalMaps(conformal_maps%face_map, & 
                                conformal_maps%edge_map, & 
                                conformal_maps%sgbc_map, & 
                                conformal_fields%faces, &
                                conformal_fields%edges, &
                                conformal_fields%sgbcs)

        call addAdditionalConformalFeatures(conformal_maps, &
                                            conformal_fields%edges)


        call initConformalSGBCStates(conformal_fields%sgbcs, conformal_maps%edge_map, sgg%ConformalSGBCProfiles, sgg%dt, control, idx, idy, idz)

                                            
        call buildIndexLists(conformal_fields)

        do i = 1, size(conformal_maps%face_map%keys)
            face => conformal_maps%face_map%getFace(conformal_maps%face_map%keys(i)%key)
            allocate(face%region_I_fields%H)
            face%region_I_fields%H = 0.0
            allocate(face%region_II_fields%H)
            face%region_II_fields%H = 0.0
            ! if (face%is_sgbc) call assignTangentialFields(face, idx, idy, idz)

            if (.not. face%is_two_sided) cycle
            edges_on_face = buildEdgesOnFace(face%cell, face%direction)
            do j = 1, size(edges_on_face)

                edge => conformal_maps%edge_map%getEdge(edges_on_face(j)%key)
                if (.not. edge%is_two_sided) then 
                    call assignEdgeFieldsOnFace(Ex, Ey, Ez, face, edges_on_face(j), j)
                else
                    call assignSplitEdgeFieldsOnFace(face, edge, j)
                end if

                faces_on_edge = buildFacesOnEdge(edge%cell, edge%direction)
                do k = 1, size(faces_on_edge)
                    if (.not. conformal_maps%face_map%hasKey(faces_on_edge(k)%key)) then 
                        cell = faces_on_edge(k)%key(1:3)
                        call assignFaceFieldsOnEdge(Hx,Hy,Hz, edge, k, cell)
                    else
                        face_of_edge => conformal_maps%face_map%getFace(faces_on_edge(k)%key)
                        call assignSplitFaceFieldsOnEdge(face_of_edge, edge, k)
                    end if
                end do

            end do

            ! link sgbc fields to face fields
            if (conformal_maps%sgbc_map%hasKey(conformal_maps%face_map%keys(i)%key)) then 
                sgbc => conformal_maps%sgbc_map%getSGBC(conformal_maps%face_map%keys(i)%key)
                call assignTangentialFields(face, sgbc)
            else 
                call voidTangentialFields(face)
            end if

        end do
    
    end subroutine

    subroutine assignTangentialFields(face, sgbc)
        type(face_t), pointer :: face
        type(sgbc_t), pointer :: sgbc
        face%region_I_fields%Etan => sgbc%state%e_1
        face%region_II_fields%Etan => sgbc%state%e_n
        sgbc%state%h_0 => face%region_I_fields%H
        sgbc%state%h_n => face%region_II_fields%H
    end subroutine

    subroutine voidTangentialFields(face)
        type(face_t), pointer :: face
        face%region_I_fields%Etan  = 0.0
        face%region_II_fields%Etan = 0.0
    end subroutine

    logical function isConformalSurface(medium)
        type(MediaData_t), intent(in) :: medium
        isConformalSurface = (medium%Is%ConformalPEC .or. medium%Is%ConformalSGBC) .and. & 
                            medium%Is%Surface .and. &
                            .not. medium%is%pml
    end function

    logical function isConformalPECSurface(medium)
        type(MediaData_t), intent(in) :: medium
        isConformalPECSurface = (medium%Is%ConformalPEC) .and. & 
                            medium%Is%Surface .and. &
                            .not. medium%is%pml
    end function

    logical function isConformalSGBCSurface(medium)
        type(MediaData_t), intent(in) :: medium
        isConformalSGBCSurface = (medium%Is%ConformalSGBC) .and. & 
                            medium%Is%Surface .and. &
                            .not. medium%is%pml
    end function

    subroutine initConformalSGBCStates(sgbcs, edge_map, profiles, dt, control, idx, idy, idz)
        type(sgbc_t), dimension(:), allocatable, intent(inout) :: sgbcs
        type(edge_map_t), intent(in) :: edge_map
        type(SGBCMaterialProfile_t), dimension(:), allocatable, intent(in) :: profiles
        real(kind=RKIND_TIME), intent(in) :: dt
        type(sim_control_t), intent(in) :: control
        real(kind=RKIND), dimension(:), pointer :: idx, idy, idz
        integer :: i, j, key(4)
        type(edge_t), pointer :: edge
        real(kind=rkind) :: lower, upper, step
        logical :: reverse = .true.
        lower = 0.0
        upper = 0.0
        do i = 1, size(sgbcs)
            do j = 1, size(sgbcs(i)%contour_inside_edges)
                step = 0.0
                key(1:3) = sgbcs(i)%contour_inside_edges(j)%cell(1:3)
                key(4) = sgbcs(i)%contour_inside_edges(j)%direction
                edge => edge_map%getEdge(key)
                if (.not. edge%is_two_sided) cycle
                if (edge%edge_region == EDGE_REGION_CONFORMAL) then 
                select case(edge%direction)
                    case(1); step = 1.0_RKIND/idx(edge%cell(1))
                    case(2); step = 1.0_RKIND/idy(edge%cell(2))
                    case(3); step = 1.0_RKIND/idz(edge%cell(3))
                end select
                lower = lower + 0.5*edge%ratio*step
                upper = upper + 0.5*(1.0-edge%ratio)*step
                end if
            end do
            call sgbcs(i)%state%init(profiles(sgbcs(i)%sgbc_profile_index), dt, lower, upper, control%sgbcFreq, control%sgbcresol, control%sgbcdepth, reverse)
            sgbcs(i)%state%length_ratio = sgbcs(i)%intersect%compute_length_ratio([1.0/idx(sgbcs(i)%cell(1)),1.0/idy(sgbcs(i)%cell(2)),1.0/idz(sgbcs(i)%cell(3))])
        end do
    end subroutine


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

    subroutine buildIndexLists(conformal_fields)
        type(conformal_fields_t), intent(inout) :: conformal_fields
        integer :: i, d
        integer, dimension(3) :: n_edges, n_faces, n_sgbcs

        n_edges = 0
        do i = 1, size(conformal_fields%edges)
            d = conformal_fields%edges(i)%direction
            n_edges(d) = n_edges(d) + 1
        end do
        n_faces = 0
        do i = 1, size(conformal_fields%faces)
            d = conformal_fields%faces(i)%direction
            n_faces(d) = n_faces(d) + 1
        end do
        n_sgbcs = 0
        do i = 1, size(conformal_fields%sgbcs)
            d = conformal_fields%sgbcs(i)%direction
            n_sgbcs(d) = n_sgbcs(d) + 1
        end do

        if (allocated(conformal_fields%edges_id_list)) deallocate(conformal_fields%edges_id_list)
        if (allocated(conformal_fields%faces_id_list)) deallocate(conformal_fields%faces_id_list)
        if (allocated(conformal_fields%sgbcs_id_list)) deallocate(conformal_fields%sgbcs_id_list)
        allocate(conformal_fields%edges_id_list(3, maxval(n_edges)))
        allocate(conformal_fields%faces_id_list(3, maxval(n_faces)))
        allocate(conformal_fields%sgbcs_id_list(3, maxval(n_sgbcs)))        
        conformal_fields%edges_id_list = 0
        conformal_fields%faces_id_list = 0
        conformal_fields%sgbcs_id_list = 0

        n_edges = 0
        do i = 1, size(conformal_fields%edges)
            d = conformal_fields%edges(i)%direction
            n_edges(d) = n_edges(d) + 1
            conformal_fields%edges_id_list(d, n_edges(d)) = i
        end do
        n_faces = 0
        do i = 1, size(conformal_fields%faces)
            d = conformal_fields%faces(i)%direction
            n_faces(d) = n_faces(d) + 1
            conformal_fields%faces_id_list(d, n_faces(d)) = i
        end do
        n_sgbcs = 0
        do i = 1, size(conformal_fields%sgbcs)
            d = conformal_fields%sgbcs(i)%direction
            n_sgbcs(d) = n_sgbcs(d) + 1
            conformal_fields%sgbcs_id_list(d, n_sgbcs(d)) = i
        end do    
    end subroutine

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


    subroutine buildConformalMaps(face_map, edge_map,sgbc_map, faces, edges, sgbcs)
        type(face_map_t), intent(inout) :: face_map
        type(edge_map_t), intent(inout) :: edge_map
        type(sgbc_map_t), intent(inout) :: sgbc_map
        type(face_t), dimension(:), allocatable :: faces
        type(edge_t), dimension(:), allocatable :: edges
        type(sgbc_t), dimension(:), allocatable :: sgbcs
        integer :: i,j
        call buildConformalFaceMap(face_map, faces)
        call allocateConformalEdgeFields(edges)
        call buildConformalEdgeMap(edge_map, edges)
        call buildConformalSGBCMap(sgbc_map, sgbcs)
    end subroutine

    subroutine buildConformalSGBCMap(sgbc_map, sgbcs)
        type(sgbc_map_t), intent(inout) :: sgbc_map
        type(sgbc_t), dimension(:), allocatable :: sgbcs
        integer :: i,j
        if (.not. allocated(sgbc_map%keys)) allocate(sgbc_map%keys(0))
        do i = 1, size(sgbcs)
            allocate(sgbcs(i)%state%e_1)
            sgbcs(i)%state%e_1 = 0.0
            allocate(sgbcs(i)%state%e_n)
            sgbcs(i)%state%e_n = 0.0
            allocate(sgbcs(i)%state%h_0)
            sgbcs(i)%state%h_0 = 0.0
            allocate(sgbcs(i)%state%h_n)
            sgbcs(i)%state%h_n = 0.0
            call sgbc_map%addSGBC(sgbcs(i))
        end do
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
            allocate(faces(i)%region_I_fields%Etan)
            allocate(faces(i)%region_II_fields%Etan)
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

    subroutine sgbc_addSGBC(this, sgbc)
        class(sgbc_map_t) :: this
        type(sgbc_t), intent(in), target :: sgbc
        integer(kind=4), dimension(4) :: sgbc_key
        type(map_key_t), dimension(:), allocatable :: aux_keys
        sgbc_key(1:3) = sgbc%cell
        sgbc_key(4) = sgbc%direction
        if (.not. allocated(this%keys)) allocate(this%keys(0))
        if (.not. this%hasKey(sgbc_key)) then 
            call this%set_ptr(key(sgbc_key), value = sgbc)
            allocate(aux_keys(size(this%keys) + 1))
            aux_keys(1:size(this%keys)) = this%keys
            aux_keys(size(this%keys) + 1)%key = sgbc_key
            deallocate(this%keys)
            allocate(this%keys(size(aux_keys)))
            this%keys = aux_keys
        end if
    end subroutine

    function sgbc_getSGBC(this, k, found) result(res)
        class(sgbc_map_t) :: this
        integer(kind=4), dimension(4) :: k
        logical, intent(inout), optional :: found
        class(*), pointer :: alloc_val
        type(sgbc_t), pointer :: res
        integer :: stat
        if (present(found)) found = .false.
        if (this%hasKey(k)) then 
            call this%get_raw_ptr(key = key(k), value = alloc_val, stat = stat)
            select type(alloc_val)
            type is(sgbc_t)
                if (present(found)) found = .true.
                res => alloc_val
            end select
        end if
    end function

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

    logical function sgbc_hasKey(this, k)
        class(sgbc_map_t) :: this
        integer(kind=4), dimension(4), intent(in) :: k
        integer :: stat
        sgbc_hasKey = .false.
        call this%check_key(key(k), stat)
        if (stat == 0) sgbc_hasKey = .true.
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

    subroutine addConformalSGBCs(sgbcs, media_faces)
        type(sgbc_t), dimension(:), allocatable, intent(inout) :: sgbcs
        type(face_t), dimension(:), allocatable, intent(in) :: media_faces
        type(sgbc_t), dimension(:), allocatable :: aux_sgbcs
        integer :: i, n
        if (.not. allocated(sgbcs)) allocate(sgbcs(0))
        if (.not. allocated(media_faces)) return
        if (size(media_faces) == 0) return
        n = 0
        do i = 1, size(media_faces)
            if (media_faces(i)%is_two_sided) n = n + 1
        end do

        allocate(aux_sgbcs(size(sgbcs) + n))
        aux_sgbcs(1:size(sgbcs)) = sgbcs
        n = 0
        do i = 1, size(media_faces)
            if (media_faces(i)%is_two_sided) then 
                n = n + 1
                
                aux_sgbcs(size(sgbcs) + n) = SGBCFace(media_faces(i))
            end if
        end do
        deallocate(sgbcs)
        allocate(sgbcs(size(aux_sgbcs)))
        sgbcs = aux_sgbcs
    end subroutine

    function SGBCFace(face) result(res)
        type(face_t), intent(in) :: face
        type(sgbc_t) :: res
        res%cell = face%cell
        res%direction = face%direction
        res%ratio = face%ratio
        res%is_two_sided = face%is_two_sided
        res%normal = face%normal
        res%intersect = face%intersect
        res%contour_inside_edges = face%contour_inside_edges
        res%sgbc_profile_index = face%sgbc_profile_index
        res%is_sgbc = .true.
    end function

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