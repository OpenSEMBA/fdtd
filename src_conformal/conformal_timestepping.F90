module conformal_mod

    ! use conformal_mod
    use fdetypes
    use fhash, only: fhash_tbl_t, key=>fhash_key

    integer, parameter :: CONFORMAL_FACE  =  0
    integer, parameter :: CONFORMAL_EDGE  =  1

    type, public :: media_maps_t
        type(fhash_tbl_t) :: Ex, Ey, Ez, Hz, Hy, Hz
    contains
        procedure :: getMedium => media_maps_getMedium
    end type

    type, public :: feature_map_t
        type(fhash_tbl_t) :: edge_map, face_map
    contains
        procedure :: getEdge
        procedure :: getFace
    end type


    ! advanceE avanza conformal+PEC+volumen
    ! conformal_timestepping comprueba PEC+surface y SIBC
    ! Si es PEC y superficie:
    !     crea 6 mapas de cell:medio (analogas a sggMi)


    !     itera sobre caras:
    !         en cada cara, conozco ratio
    !         ratio nuevo = 1 - ratio
    !         gm1 y gm2 dependen del ratio nuevo
    !         los E => a los edges o campos de la region I

    ! t_t alguno de los tipos conformal
    type :: conformal_field_t
        real(kind=rkind), pointer :: p => null()
        real(kind=rkind), allocatable :: owned
    end type

    type :: conformal_edge_fields_t
        real(kind=rkind), pointer :: E => null()
        real(kind=rkind), pointer :: H1 => null()
        real(kind=rkind), pointer :: H2 => null()
        real(kind=rkind), pointer :: H3 => null()
        real(kind=rkind), pointer :: H4 => null()
    end type

    type :: conformal_face_fields_t
        type(conformal_field_t) :: H
        type(conformal_field_t) :: E1
    ! real(kind=rkind), pointer :: H => null()
        ! real(kind=rkind), pointer :: E1 => null()
        real(kind=rkind), pointer :: E2 => null()
        real(kind=rkind), pointer :: E3 => null()
        real(kind=rkind), pointer :: E4 => null()
    end type

    type :: conformal_edge_t
        type(conformal_edge_fields_t) :: r1fields, r2fields
        ! type(constants_t) :: gI, gII
    end type
    type :: conformal_face_t
        type(conformal_face_fields_t) :: r1fields, r2fields
        ! type(constants_t) :: gI, gII
    end type
    
    type t_t
        type(conformal_edge_t), dimension(:), pointer :: edges
        type(conformal_face_t), dimension(:), pointer :: faces
    end type 

    type Conformal_t
        integer (kind=simple) :: nConformalMedia
        type(t_t), pointer, dimension(:) :: medium
        type(constants_t) :: g
        type(kind=rkind), pointer, dimension(:) :: ratio
        type(kind=rkind), pointer, dimension(:) :: id1, id2

    contains
        procedure :: setFields
    end type
    type (Conformal_t), save, target :: conformal

contains

    function media_maps_getMedium(this, cell, direction, found) result(res)
        class(media_maps_t) :: this
        integer (kind = 4), intent(in) :: cell
        integer, intent(in) :: direction
        logical, intent(inout), optional :: found
        integer :: stat
        integer :: medium

        if (present(found)) found = .false.
        select case (direction)
        case(EDGE_X)
            call this%Ex%get(key(cell), medium, stat)
        case(EDGE_Y)
            call this%Ey%get(key(cell), medium, stat)
        case(EDGE_Z)
            call this%Ez%get(key(cell), medium, stat)
        case(FACE_X)
            call this%Hx%get(key(cell), medium, stat)
        case(FACE_Y)
            call this%Hy%get(key(cell), medium, stat)
        case(FACE_Z)
            call this%Hz%get(key(cell), medium, stat)
        end select
        if (stat /= 0) return
        res = medium
        if (present(found)) found = .true.
    end function

    subroutine calc_conf_constants()

    end subroutine


    function buildFeatureMap(sgg) result(res)
        type(sggfdtdinfo), intent(in) :: sgg
        type(feature_map_t), intent(out) :: res
        integer :: i,j
        integer (kind=simple) :: k(4)
        do i = 1, sgg%NumMedia
            if (sgg%med(i)%Is%ConformalPEC .and. sgg%med(i)%Is%surface) then 
                if (allocated(sgg%med(i)%ConformalEdge)) then 
                    do j = 1, size(sgg%med(i)%ConformalEdge%edges) 
                        k(1:3) = sgg%med(i)%ConformalEdge%edges(j)%cell
                        k(4) = sgg%med(i)%ConformalEdge%edges(j)%direction
                        call res%edge_map%set(key(k), value = sgg%med(i)%ConformalEdge%edges(j))
                    end do
                else if (allocated(sgg%med(i)%ConformalFace)) then 
                    do j = 1, size(sgg%med(i)%ConformalFace%faces) 
                        k(1:3) = sgg%med(i)%ConformalFace%faces(j)%cell
                        k(4) = sgg%med(i)%ConformalFace%faces(j)%direction
                        call res%face_map%set(key(k), value = sgg%med(i)%ConformalFace%faces(j))
                    end do
                end if
            end if
        end do
    end function

    function getEdge(this, id, found) result(res)
        class(feature_map_t) :: this
        logical, intent(inout), optional :: found
        type(edge_t), intent(out) :: res
        integer(kind=simple) :: id(4)
        integer :: stat
        class(*), allocatable :: d

        if (present(found)) found = .false.
        call this%edge_map%get(key(k), id, stat)
        if (stat /= 0) return
        select type(d)
            type is (edge_t)
            res = d
            if (present(found)) found = .true.
        end select
    end function

    function getFace(this, id, found) result(res)
        class(feature_map_t) :: this
        logical, intent(inout), optional :: found
        type(face_t), intent(out) :: res
        integer(kind=simple) :: id(4)
        integer :: stat
        class(*), allocatable :: d

        if (present(found)) found = .false.
        call this%face_map%get(key(k), id, stat)
        if (stat /= 0) return
        select type(d)
            type is (face_t)
            res = d
            if (present(found)) found = .true.
        end select
    end function


    function buildMediaMaps(conformal_surface, edge_ratios, face_ratios, num_media) result(res)
        type(ConformalMedia_t), intent(in) :: conformal_surface
        real(kind=rkind), dimension(:), allocatable, intent(in) :: edge_ratios, face_ratios
        type(media_maps_t) :: res
        integer (kind=4) :: med_number
        integer (kind=4) :: cell(3)
        integer :: i, j

        ! el ratio que se busca deberia ser 1-ratio?
        ! y si no existe, crea medio nuevo?
        do i = 1, conformal_surface%n_edges_media
            med_number = findloc(edge_ratios, conformal_volumes%edge_media(i)%ratio, 1)
            do j = 1, conformal_volumes%edge_media(i)%size
                cell(:) = conformal_volumes%edge_media(i)%edges(j)%cell(:)
                select case(conformal_volumes%edge_media(i)%edges(j)%direction)
                case(EDGE_X)
                    call res%Ex%set(key(cell), value = med_number)
                case(EDGE_Y)
                    call res%Ey%set(key(cell), value = med_number)
                case(EDGE_Z)
                    call res%Ez%set(key(cell), value = med_number)
                end select
            end do
        end do

        do i = 1, conformal_surface%n_faces_media
            med_number = size(edge_ratios) + findloc(face_ratios, conformal_volumes%face_media(i)%ratio, 1)
            do j = 1, conformal_volumes%face_media(i)%size
                cell(:) = conformal_volumes%face_media(i)%faces(j)%cell(:)
                select case(conformal_volumes%face_media(i)%faces(j)%direction)
                case(FACE_X)
                    call res%Hx%set(key(cell), value = med_number)
                case(FACE_Y)
                    call res%Hy%set(key(cell), value = med_number)
                case(FACE_Z)
                    call res%Hz%set(key(cell), value = med_number)
                end select
            end do
        end do
    end function

    logical function featureIs(feature, featureType)
        class(*), intent(in) :: feature
        integer, intent(in) :: featureType
        integer :: foundType
        featureIs = .false.
        select type (feature)
        type is (conformal_face_media_t)
            foundType = CONFORMAL_FACE
        type is (conformal_edge_media_t)
            foundType = CONFORMAL_EDGE
        end select
        featureIs = (foundType == featureType)
    end function

    subroutine setFaceFields(this, n, i) 
        class(Conformal_t) :: this
        integer (kind=SINGLE) :: n, i
        allocate(this%medium(n)%faces(i)%r2fields%H%owned, 0.0)
        this%medium(n)%faces(i)%r2fields%H%p => this%medium(nmedia)%faces(j)%r2fields%H%owned

        allocate(this%medium(n)%faces(i)%r1fields%E1%owned, 0.0)
        allocate(this%medium(n)%faces(i)%r1fields%E2%owned, 0.0)
        allocate(this%medium(n)%faces(i)%r1fields%E3%owned, 0.0)
        allocate(this%medium(n)%faces(i)%r1fields%E4%owned, 0.0)
        this%medium(n)%faces(i)%r1fields%E1%p => this%medium(nmedia)%faces(j)%r1fields%E1%owned
        this%medium(n)%faces(i)%r1fields%E2%p => this%medium(nmedia)%faces(j)%r1fields%E2%owned
        this%medium(n)%faces(i)%r1fields%E3%p => this%medium(nmedia)%faces(j)%r1fields%E3%owned
        this%medium(n)%faces(i)%r1fields%E4%p => this%medium(nmedia)%faces(j)%r1fields%E4%owned

        allocate(this%medium(n)%faces(i)%r2fields%E1%owned, 0.0)
        allocate(this%medium(n)%faces(i)%r2fields%E2%owned, 0.0)
        allocate(this%medium(n)%faces(i)%r2fields%E3%owned, 0.0)
        allocate(this%medium(n)%faces(i)%r2fields%E4%owned, 0.0)
        this%medium(n)%faces(i)%r2fields%E1%p => this%medium(nmedia)%faces(j)%r2fields%E1%owned
        this%medium(n)%faces(i)%r2fields%E2%p => this%medium(nmedia)%faces(j)%r2fields%E2%owned
        this%medium(n)%faces(i)%r2fields%E3%p => this%medium(nmedia)%faces(j)%r2fields%E3%owned
        this%medium(n)%faces(i)%r2fields%E4%p => this%medium(nmedia)%faces(j)%r2fields%E4%owned
    end subroutine

    subroutine setEdgeFields(this, n, i) 
        class(Conformal_t) :: this
        integer (kind=SINGLE) :: n, i
        allocate(this%medium(n)%edges(i)%r2fields%H%owned, 0.0)
        this%medium(n)%edges(i)%r2fields%E%p => this%medium(nmedia)%edges(j)%r2fields%E%owned

        allocate(this%medium(n)%edges(i)%r1fields%H1%owned, 0.0)
        allocate(this%medium(n)%edges(i)%r1fields%H2%owned, 0.0)
        allocate(this%medium(n)%edges(i)%r1fields%H3%owned, 0.0)
        allocate(this%medium(n)%edges(i)%r1fields%H4%owned, 0.0)
        this%medium(n)%edges(i)%r1fields%H1%p => this%medium(nmedia)%edges(j)%r1fields%H1%owned
        this%medium(n)%edges(i)%r1fields%H2%p => this%medium(nmedia)%edges(j)%r1fields%H2%owned
        this%medium(n)%edges(i)%r1fields%H3%p => this%medium(nmedia)%edges(j)%r1fields%H3%owned
        this%medium(n)%edges(i)%r1fields%H4%p => this%medium(nmedia)%edges(j)%r1fields%H4%owned

        allocate(this%medium(n)%edges(i)%r2fields%H1%owned, 0.0)
        allocate(this%medium(n)%edges(i)%r2fields%H2%owned, 0.0)
        allocate(this%medium(n)%edges(i)%r2fields%H3%owned, 0.0)
        allocate(this%medium(n)%edges(i)%r2fields%H4%owned, 0.0)
        this%medium(n)%edges(i)%r2fields%H1%p => this%medium(nmedia)%edges(j)%r2fields%H1%owned
        this%medium(n)%edges(i)%r2fields%H2%p => this%medium(nmedia)%edges(j)%r2fields%H2%owned
        this%medium(n)%edges(i)%r2fields%H3%p => this%medium(nmedia)%edges(j)%r2fields%H3%owned
        this%medium(n)%edges(i)%r2fields%H4%p => this%medium(nmedia)%edges(j)%r2fields%H4%owned
    end subroutine

    subroutine initConformal(sgg, Ex, Ey, Ez, Hx, Hy, Hz, idxe, idye, idze, idxh, idyh, idzh)
        type (sggfdtdinfo), intent(inout) :: sgg
        real (kind=rkind), intent(in) , target :: &
        Ex(sgg%alloc(iEx)%XI : sgg%alloc(iEx)%XE,sgg%alloc(iEx)%YI : sgg%alloc(iEx)%YE,sgg%alloc(iEx)%ZI : sgg%alloc(iEx)%ZE),&
        Ey(sgg%alloc(iEy)%XI : sgg%alloc(iEy)%XE,sgg%alloc(iEy)%YI : sgg%alloc(iEy)%YE,sgg%alloc(iEy)%ZI : sgg%alloc(iEy)%ZE),&
        Ez(sgg%alloc(iEz)%XI : sgg%alloc(iEz)%XE,sgg%alloc(iEz)%YI : sgg%alloc(iEz)%YE,sgg%alloc(iEz)%ZI : sgg%alloc(iEz)%ZE),&
        Hx(sgg%alloc(iHx)%XI : sgg%alloc(iHx)%XE,sgg%alloc(iHx)%YI : sgg%alloc(iHx)%YE,sgg%alloc(iHx)%ZI : sgg%alloc(iHx)%ZE),&
        Hy(sgg%alloc(iHy)%XI : sgg%alloc(iHy)%XE,sgg%alloc(iHy)%YI : sgg%alloc(iHy)%YE,sgg%alloc(iHy)%ZI : sgg%alloc(iHy)%ZE),&
        Hz(sgg%alloc(iHz)%XI : sgg%alloc(iHz)%XE,sgg%alloc(iHz)%YI : sgg%alloc(iHz)%YE,sgg%alloc(iHz)%ZI : sgg%alloc(iHz)%ZE)

        real (kind=rkind), intent(in) :: 
        Idxe(sgg%alloc(iHx)%XI : sgg%ALLOC(iHx)%XE), &
        Idye(sgg%alloc(iHy)%YI : sgg%ALLOC(iHy)%YE), &
        Idze(sgg%alloc(iHz)%ZI : sgg%ALLOC(iHz)%ZE), &
        Idxh(sgg%alloc(iEx)%XI : sgg%ALLOC(iEx)%XE), &
        Idyh(sgg%alloc(iEy)%YI : sgg%ALLOC(iEy)%YE), &
        Idzh(sgg%alloc(iEz)%ZI : sgg%ALLOC(iEz)%ZE)


        integer (kind=simple) :: nmedia, cell(3), direction, k(4)
        integer :: i, j
        type(feature_map_t) :: feature_map
        type(edge_t) :: edge
        logical :: found
        nmedia = 0
        do i = 1, sgg%NumMedia
            if (sgg%med(i)%Is%ConformalPEC .and. sgg%med(i)%Is%surface) then 
                nmedia = nmedia + 1
            end if
        end do
        conformal%nConformalMedia = nmedia

        feature_map = buildFeatureMap(sgg)

        nmedia = 0
        allocate(conformal%medium(nmedia))
        allocate(conformal%g(nmedia))
        allocate(conformal%ratio(nmedia))
        allocate(conformal%id1(nmedia))
        allocate(conformal%id2(nmedia))

        do i = 1, sgg%NumMedia
            if (sgg%med(i)%Is%ConformalPEC .and. sgg%med(i)%Is%surface) then 
                if (allocated(sgg%med(i)%ConformalFace)) then 
                    nmedia = nmedia + 1
                    allocate(conformal%medium(nmedia)%faces,  size(sgg%med(i)%ConformalFace(1)%faces))
                    do j = 1, size(sgg%med(i)%ConformalFace(1)%faces)

                        call conformal%setFaceFields(nmedia, j)

                        cell(:) = sgg%med(i)%conformalFace(1)%faces(j)%cell(:)
                        direction = sgg%med(i)%conformalSurface(1)%faces(j)%direction
                        
                        conformal%ratio(nmedia) = sgg%med(i)%conformalSurface(1)%faces(j)%ratio
                        conformal%g(nmedia)%g1  = 1.0
                        conformal%g(nmedia)%g2  = sgg%dt / (Eps0*sgg%Med(nmedia)%Epr)
                        conformal%g(nmedia)%gm1 = 1.0
                        conformal%g(nmedia)%gm2 = sgg%dt / (Mu0*sgg%Med(nmedia)%Mur)
                        
                        select case(direction)
                        case(FACE_X)
                            conformal%medium(nmedia)%faces(j)%r1fields%H%p => Hx(cell(1), cell(2), cell(3))
                        case(FACE_Y)
                            conformal%medium(nmedia)%faces(j)%r1fields%H%p => Hy(cell(1), cell(2), cell(3))
                        case(FACE_Z)
                            
                            conformal%id1 = Idye(cell(2))
                            conformal%id2 = Idxe(cell(1))
                            conformal%medium(nmedia)%faces(j)%r1fields%H%p => Hz(cell(1), cell(2), cell(3))

                            k(1:3) = cell; k(4) = EDGE_Y
                            !del edge solo hace falta el ratio
                            ! map (n, j, dir) -> edge
                            ! de donde sale el ratio? hace falta?
                            ! quiza siempre es puntero a campo ppal, y a campo de region 2
                            edge = feature_map%getEdge(k, found)
                            if (found) then 
                                if (edge%ratio /= 0) then 
                                    ! deallocate(conformal%medium(nmedia)%faces(j)%r1fields%E1%p)
                                    conformal%medium(nmedia)%faces(j)%r1fields%E1%p => Ey(k(1), k(2), k(3))
                                    conformal%medium(nmedia)%faces(j)%r2fields%E1%p => 
                                        conformal%medium(xxx)%edges(yyy)%r2fields%E%p
                                end if
                            end if

                            k(1:3) = cell; k(2) = k(2) + 1; k(4) = EDGE_X
                            edge = feature_map%getEdge(k, found)
                            if (found) then 
                                if (edge%ratio /= 0) then 
                                    ! deallocate(conformal%medium(nmedia)%faces(j)%r1fields%E1%p)
                                    conformal%medium(nmedia)%faces(j)%r1fields%E2%p => Ex(k(1), k(2), k(3))
                                end if
                            end if

                            k(1:3) = cell; k(1) = k(1) + 1; k(4) = EDGE_Y
                            edge = feature_map%getEdge(k, found)
                            if (found) then 
                                if (edge%ratio /= 0) then 
                                    ! deallocate(conformal%medium(nmedia)%faces(j)%r1fields%E1%p)
                                    conformal%medium(nmedia)%faces(j)%r1fields%E3%p => Ey(k(1), k(2), k(3))
                                end if
                            end if

                            k(1:3) = cell; k(4) = EDGE_X
                            edge = feature_map%getEdge(k, found)
                            if (found) then 
                                if (edge%ratio /= 0) then 
                                    ! deallocate(conformal%medium(nmedia)%faces(j)%r1fields%E1%p)
                                    conformal%medium(nmedia)%faces(j)%r1fields%E4%p => Ex(k(1), k(2), k(3))
                                end if
                            end if


                        end select




                    end do
                else if (allocated(sgg%med(i)%ConformalEdge)) then 
                !
                end if

            end if
        end do


    end subroutine


    ! subroutine initConformal(sgg, Ex, Ey, Ez, Hx, Hy, Hz, conformal_media)
    !     type (sggfdtdinfo), intent(inout) :: sgg
    !     real (kind=rkind), intent(in) , target :: &
    !     Ex(sgg%alloc(iEx)%XI : sgg%alloc(iEx)%XE,sgg%alloc(iEx)%YI : sgg%alloc(iEx)%YE,sgg%alloc(iEx)%ZI : sgg%alloc(iEx)%ZE),&
    !     Ey(sgg%alloc(iEy)%XI : sgg%alloc(iEy)%XE,sgg%alloc(iEy)%YI : sgg%alloc(iEy)%YE,sgg%alloc(iEy)%ZI : sgg%alloc(iEy)%ZE),&
    !     Ez(sgg%alloc(iEz)%XI : sgg%alloc(iEz)%XE,sgg%alloc(iEz)%YI : sgg%alloc(iEz)%YE,sgg%alloc(iEz)%ZI : sgg%alloc(iEz)%ZE),&
    !     Hx(sgg%alloc(iHx)%XI : sgg%alloc(iHx)%XE,sgg%alloc(iHx)%YI : sgg%alloc(iHx)%YE,sgg%alloc(iHx)%ZI : sgg%alloc(iHx)%ZE),&
    !     Hy(sgg%alloc(iHy)%XI : sgg%alloc(iHy)%XE,sgg%alloc(iHy)%YI : sgg%alloc(iHy)%YE,sgg%alloc(iHy)%ZI : sgg%alloc(iHy)%ZE),&
    !     Hz(sgg%alloc(iHz)%XI : sgg%alloc(iHz)%XE,sgg%alloc(iHz)%YI : sgg%alloc(iHz)%YE,sgg%alloc(iHz)%ZI : sgg%alloc(iHz)%ZE)
    !     type(conformal_media_t), intent(in) :: conformal_media
    !     type (edge_t) :: edge
    !     integer (kind=4) :: cell(3), k(4), ratioII
    !     ! make fields in faces and edges point to the corresponding region I or II fields
    !     integer :: i, j

    !     do i = 1, sgg%NumMedia    
    !         if (sgg%med(i)%Is%conformal .and. sgg%med(i)%Is%surface) then 
    !             if (featureIs(sgg%med(i)%conformalSurface(1)), CONFORMAL_FACE) then 
    !                 do j = 1, sgg%med(i)%conformalSurface(1)%size 
    !                     cell(:) = sgg%med(i)%conformalSurface(1)%faces(j)%cell(:)
    !                     select case(sgg%med(i)%conformalSurface(1)%faces(j)%direction)
    !                     case(FACE_Z)
    !                         sgg%med(i)%conformalSurface(1)%faces(j)%region_I_fields%H => Hz(cell(1), cell(2), cell(3))
    !                         k = [cell, FACE_Z]
    !                         ! 1: Ey(cell(1)    , cell(2)    , cell(3))
    !                         ! si existe
    !                         call edge_map%get(key(k), edge)
    !                         if (edge%ratio == 0) then 
    !                             sgg%med(i)%conformalSurface(1)%faces(j)%region_I_fields%E1 = 0
    !                             sgg%med(i)%conformalSurface(1)%faces(j)%region_II_fields%E1 = edge%region_II_fields%E
    !                         else if (edge%ratio == 1) then 
    !                             sgg%med(i)%conformalSurface(1)%faces(j)%region_I_fields%E1 => Ey(cell(1)    , cell(2)    , cell(3))
                                
    !                             allocate(sgg%med(i)%conformalSurface(1)%faces(j)%region_II_fields%E1%owned, 0.0)
    !                             sgg%med(i)%conformalSurface(1)%faces(j)%region_II_fields%E1%p => % 
    !                                 sgg%med(i)%conformalSurface(1)%faces(j)%region_II_fields%E1%owned
    !                             ! sgg%med(i)%conformalSurface(1)%faces(j)%region_I_fields%E1 = Ey(cell(1)    , cell(2)    , cell(3))
    !                             ! sgg%med(i)%conformalSurface(1)%faces(j)%region_II_fields%E1 = 0
    !                         else 
    !                             sgg%med(i)%conformalSurface(1)%faces(j)%region_I_fields%E1 => Ey(cell(1)    , cell(2)    , cell(3))
    !                             ! sgg%med(i)%conformalSurface(1)%faces(j)%region_II_fields%E1 => edge%region_II_fields%E
    !                             allocate(sgg%med(i)%conformalSurface(1)%faces(j)%region_II_fields%E1%owned, 0.0)
    !                             sgg%med(i)%conformalSurface(1)%faces(j)%region_II_fields%E1%p => % 
    !                                 sgg%med(i)%conformalSurface(1)%faces(j)%region_II_fields%E1%owned
    !                         end if
    !                         ! Ex(cell(1)    , cell(2) + 1, cell(3))
    !                         ! Ey(cell(1) + 1, cell(2)    , cell(3))
    !                         ! Ex(cell(1)    , cell(2)    , cell(3))                            
    !                     end select
    !                 end do
    !             end if
    !         end  if
    !     end do

    !     do i = 1, conformal_media%n_faces_media
    !         do j = 1, conformal_media%face_media(i)%size
    !             cell(:) = conformal_media%face_media(i)%faces(j)%cell(:)
    !             select case(conformal_media%face_media(j)%faces(k)%direction)
    !             case(FACE_X)
    !                 conformal_media%face_media(i)%faces(j)%region_I_fields%H  => Hx(cell(1), cell(2)    , cell(3))
    !                 conformal_media%face_media(i)%faces(j)%region_I_fields%E1 => Ez(cell(1), cell(2)    , cell(3))
    !                 conformal_media%face_media(i)%faces(j)%region_I_fields%E2 => Ey(cell(1), cell(2)    , cell(3) + 1)
    !                 conformal_media%face_media(i)%faces(j)%region_I_fields%E3 => Ez(cell(1), cell(2) + 1, cell(3))
    !                 conformal_media%face_media(i)%faces(j)%region_I_fields%E4 => Ey(cell(1), cell(2)    , cell(3))
    !             case(FACE_Y)
    !                 conformal_media%face_media(i)%faces(j)%region_I_fields%H  => Hy(cell(1)    , cell(2), cell(3))
    !                 conformal_media%face_media(i)%faces(j)%region_I_fields%E1 => Ex(cell(1)    , cell(2), cell(3))
    !                 conformal_media%face_media(i)%faces(j)%region_I_fields%E2 => Ez(cell(1) + 1, cell(2), cell(3))
    !                 conformal_media%face_media(i)%faces(j)%region_I_fields%E3 => Ex(cell(1)    , cell(2), cell(3) + 1)
    !                 conformal_media%face_media(i)%faces(j)%region_I_fields%E4 => Ez(cell(1)    , cell(2), cell(3))
    !             case(FACE_Z)
    !                 conformal_media%face_media(i)%faces(j)%region_I_fields%H  => Hz(cell(1)    , cell(2)    , cell(3))
    !                 conformal_media%face_media(i)%faces(j)%region_I_fields%E1 => Ey(cell(1)    , cell(2)    , cell(3))
    !                 conformal_media%face_media(i)%faces(j)%region_I_fields%E2 => Ex(cell(1)    , cell(2) + 1, cell(3))
    !                 conformal_media%face_media(i)%faces(j)%region_I_fields%E3 => Ey(cell(1) + 1, cell(2)    , cell(3))
    !                 conformal_media%face_media(i)%faces(j)%region_I_fields%E4 => Ex(cell(1)    , cell(2)    , cell(3))
    !            end select
    !         end do
    !     end do
    !     ! and fields from region II?

    ! end subroutine

    subroutine advanceConformalE(sgg)
        ! type(SGGFDTDINFO), intent(in) :: sgg
        ! type(media_maps_t), intent(in) :: media_maps
        real (kind = rkind), pointer :: E, H1, H2, H3, H4
        real (kind = rkind) :: g1, g2, gm1, gm2
        real (kind = rkind) :: id1, id2, ratio
        ! real (kind = rkind), pointer :: rIIH, rIIE1, rIIE2, rIIE3, rIIE4
        integer :: i, j
        real (kind = rkind) :: id1, id2

        do i = 1, conformal%nConformalMedia
            if (allocated(conformal%medium(i)%edged)) then 
                do j = 1, size(conformal%medium(i)%edges)
                    ! E update in region 1
                    E  => conformal%medium(i)%edges(j)%r1fields%E%p
                    H1 => conformal%medium(i)%edges(j)%r1fields%H1%p
                    H2 => conformal%medium(i)%edges(j)%r1fields%H2%p
                    H3 => conformal%medium(i)%edges(j)%r1fields%H3%p
                    H4 => conformal%medium(i)%edges(j)%r1fields%H4%p

                    g1 = conformal%g%g1(i)
                    g2 = conformal%g%g2(i)
                    id1 = conformal%id1(i)
                    id2 = conformal%id2(i)

                    E = g1*E + g2*((H2 - H4)*id1 - (H3-H1)*id2)

                    ! E update in region 2
                    E  => conformal%medium(i)%edges(j)%r2fields%E%p
                    H1 => conformal%medium(i)%edges(j)%r2fields%H1%p
                    H2 => conformal%medium(i)%edges(j)%r2fields%H2%p
                    H3 => conformal%medium(i)%edges(j)%r2fields%H3%p
                    H4 => conformal%medium(i)%edges(j)%r2fields%H4%p

                    ratio = conformal%ratio(i)
                    g2 = (1-ratio)*g2/ratio
                    E = g1*E + g2*((H2 - H4)*id1 - (H3-H1)*id2)
                end do
            end if
        end do
    end subroutine

    subroutine advanceConformalH()
        ! type(SGGFDTDINFO), intent(in) :: sgg
        ! type(media_maps_t), intent(in) :: media_maps
        real (kind = rkind), pointer :: H, E1, E2, E3, E4
        real (kind = rkind) :: g1, g2, gm1, gm2
        real (kind = rkind) :: id1, id2, ratio
        ! real (kind = rkind), pointer :: rIIH, rIIE1, rIIE2, rIIE3, rIIE4
        integer :: i, j

        do i = 1, conformal%nConformalMedia
            if (allocated(conformal%medium(i)%faces)) then 
                do j = 1, size(conformal%medium(i)%faces)
                    ! H update in region 1
                    H  => conformal%medium(i)%faces(j)%r1fields%H%p
                    E1 => conformal%medium(i)%faces(j)%r1fields%E1%p
                    E2 => conformal%medium(i)%faces(j)%r1fields%E2%p
                    E3 => conformal%medium(i)%faces(j)%r1fields%E3%p
                    E4 => conformal%medium(i)%faces(j)%r1fields%E4%p

                    gm1 = conformal%g%gm1(i)
                    gm2 = conformal%g%gm2(i)
                    id1 = conformal%id1(i)
                    id2 = conformal%id2(i)

                    H = gm1*H + gm2*((E2 - E4)*id1 - (E3-E1)*id2)

                    ! H update in region 2
                    H  => conformal%medium(i)%faces(j)%r2fields%H%p
                    E1 => conformal%medium(i)%faces(j)%r2fields%E1%p
                    E2 => conformal%medium(i)%faces(j)%r2fields%E2%p
                    E3 => conformal%medium(i)%faces(j)%r2fields%E3%p
                    E4 => conformal%medium(i)%faces(j)%r2fields%E4%p

                    ratio = conformal%ratio(i)
                    gm2 = gm2*ratio/(1-ratio)
                    H = gm1*H + gm2*((E2 - E4)*id1 - (E3-E1)*id2)
                end do
            end if
        end do
    end subroutine

end module