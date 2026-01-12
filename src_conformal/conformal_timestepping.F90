module conformal_mod

    ! use conformal_mod
    use fhash, only: fhash_tbl_t, key=>fhash_key

    integer, parameter :: CONFORMAL_FACE  =  0
    integer, parameter :: CONFORMAL_EDGE  =  1

    type, public :: media_maps_t
        type(fhash_tbl_t) :: Ex, Ey, Ez, Hz, Hy, Hz
    contains
        procedure :: getMedium => media_maps_getMedium
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

    function buildMediaMaps(conformal_surface, edge_ratios, face_ratios, num_media) result(res)
        type(ConformalMedia_t), intent(in) :: conformal_surface
        real(kind=rkind), dimension(:), allocatable, intent(in) :: edge_ratios, face_ratios
        type(media_maps_t) :: res
        integer (kind=4) :: med_number
        integer (kind=4) :: cell(3)
        integer :: i, j

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

    subroutine initConformal(sgg, Ex, Ey, Ez, Hx, Hy, Hz, conformal_media)
        type (sggfdtdinfo), intent(inout) :: sgg
        real (kind=rkind), intent(in) , target :: &
        Ex(sgg%alloc(iEx)%XI : sgg%alloc(iEx)%XE,sgg%alloc(iEx)%YI : sgg%alloc(iEx)%YE,sgg%alloc(iEx)%ZI : sgg%alloc(iEx)%ZE),&
        Ey(sgg%alloc(iEy)%XI : sgg%alloc(iEy)%XE,sgg%alloc(iEy)%YI : sgg%alloc(iEy)%YE,sgg%alloc(iEy)%ZI : sgg%alloc(iEy)%ZE),&
        Ez(sgg%alloc(iEz)%XI : sgg%alloc(iEz)%XE,sgg%alloc(iEz)%YI : sgg%alloc(iEz)%YE,sgg%alloc(iEz)%ZI : sgg%alloc(iEz)%ZE),&
        Hx(sgg%alloc(iHx)%XI : sgg%alloc(iHx)%XE,sgg%alloc(iHx)%YI : sgg%alloc(iHx)%YE,sgg%alloc(iHx)%ZI : sgg%alloc(iHx)%ZE),&
        Hy(sgg%alloc(iHy)%XI : sgg%alloc(iHy)%XE,sgg%alloc(iHy)%YI : sgg%alloc(iHy)%YE,sgg%alloc(iHy)%ZI : sgg%alloc(iHy)%ZE),&
        Hz(sgg%alloc(iHz)%XI : sgg%alloc(iHz)%XE,sgg%alloc(iHz)%YI : sgg%alloc(iHz)%YE,sgg%alloc(iHz)%ZI : sgg%alloc(iHz)%ZE)
        type(conformal_media_t), intent(in) :: conformal_media
        type (edge_t) :: edge
        integer (kind=4) :: cell(3), k(4), ratioII
        ! make fields in faces and edges point to the corresponding region I or II fields
        integer :: i, j

        do i = 1, sgg%NumMedia    
            if (sgg%med(i)%Is%conformal .and. sgg%med(i)%Is%surface) then 
                if (featureIs(sgg%med(i)%conformalSurface(1)), CONFORMAL_FACE) then 
                    do j = 1, sgg%med(i)%conformalSurface(1)%size 
                        cell(:) = sgg%med(i)%conformalSurface(1)%faces(j)%cell(:)
                        select case(sgg%med(i)%conformalSurface(1)%faces(j)%direction)
                        case(FACE_Z)
                            sgg%med(i)%conformalSurface(1)%faces(j)%region_I_fields%H => Hz(cell(1), cell(2), cell(3))
                            k = [cell, FACE_Z]
                            ! 1: Ey(cell(1)    , cell(2)    , cell(3))
                            ! si existe
                            call edge_map%get(key(k), edge)
                            if (edge%ratio == 0) then 
                                sgg%med(i)%conformalSurface(1)%faces(j)%region_I_fields%E1 = 0
                                sgg%med(i)%conformalSurface(1)%faces(j)%region_II_fields%E1 = edge%region_II_fields%E
                            else if (edge%ratio == 1) then 
                                sgg%med(i)%conformalSurface(1)%faces(j)%region_I_fields%E1 = Ey(cell(1)    , cell(2)    , cell(3))
                                sgg%med(i)%conformalSurface(1)%faces(j)%region_II_fields%E1 = 0
                            else 
                                sgg%med(i)%conformalSurface(1)%faces(j)%region_I_fields%E1 => Ey(cell(1)    , cell(2)    , cell(3))
                                sgg%med(i)%conformalSurface(1)%faces(j)%region_II_fields%E1 => edge%region_II_fields%E
                            end if
                            ! Ex(cell(1)    , cell(2) + 1, cell(3))
                            ! Ey(cell(1) + 1, cell(2)    , cell(3))
                            ! Ex(cell(1)    , cell(2)    , cell(3))                            
                        end select
                    end do
                end if
            end  if
        end do

        do i = 1, conformal_media%n_faces_media
            do j = 1, conformal_media%face_media(i)%size
                cell(:) = conformal_media%face_media(i)%faces(j)%cell(:)
                select case(conformal_media%face_media(j)%faces(k)%direction)
                case(FACE_X)
                    conformal_media%face_media(i)%faces(j)%region_I_fields%H  => Hx(cell(1), cell(2)    , cell(3))
                    conformal_media%face_media(i)%faces(j)%region_I_fields%E1 => Ez(cell(1), cell(2)    , cell(3))
                    conformal_media%face_media(i)%faces(j)%region_I_fields%E2 => Ey(cell(1), cell(2)    , cell(3) + 1)
                    conformal_media%face_media(i)%faces(j)%region_I_fields%E3 => Ez(cell(1), cell(2) + 1, cell(3))
                    conformal_media%face_media(i)%faces(j)%region_I_fields%E4 => Ey(cell(1), cell(2)    , cell(3))
                case(FACE_Y)
                    conformal_media%face_media(i)%faces(j)%region_I_fields%H  => Hy(cell(1)    , cell(2), cell(3))
                    conformal_media%face_media(i)%faces(j)%region_I_fields%E1 => Ex(cell(1)    , cell(2), cell(3))
                    conformal_media%face_media(i)%faces(j)%region_I_fields%E2 => Ez(cell(1) + 1, cell(2), cell(3))
                    conformal_media%face_media(i)%faces(j)%region_I_fields%E3 => Ex(cell(1)    , cell(2), cell(3) + 1)
                    conformal_media%face_media(i)%faces(j)%region_I_fields%E4 => Ez(cell(1)    , cell(2), cell(3))
                case(FACE_Z)
                    conformal_media%face_media(i)%faces(j)%region_I_fields%H  => Hz(cell(1)    , cell(2)    , cell(3))
                    conformal_media%face_media(i)%faces(j)%region_I_fields%E1 => Ey(cell(1)    , cell(2)    , cell(3))
                    conformal_media%face_media(i)%faces(j)%region_I_fields%E2 => Ex(cell(1)    , cell(2) + 1, cell(3))
                    conformal_media%face_media(i)%faces(j)%region_I_fields%E3 => Ey(cell(1) + 1, cell(2)    , cell(3))
                    conformal_media%face_media(i)%faces(j)%region_I_fields%E4 => Ex(cell(1)    , cell(2)    , cell(3))
               end select
            end do
        end do
        ! and fields from region II?

    end subroutine

    subroutine advanceHz(media_maps)
        type(media_maps_t), intent(in) :: media_maps
        real (kind = rkind), pointer :: H, E1, E2, E3, E4
        real (kind = rkind), pointer :: rIIH, rIIE1, rIIE2, rIIE3, rIIE4
        integer :: i, j, medium, direction
        real (kind = rkind) :: id1, id2
        integer (kind=4) :: cell(3)
        logical :: found

        !check if volume/surface
        do i = 1, conformal_media%n_faces_media
            do j = 1, conformal_media%face_media(i)%size
                H  => conformal_media%face_media(i)%faces(j)%region_I_fields%H
                E1 => conformal_media%face_media(i)%faces(j)%region_I_fields%E1
                E2 => conformal_media%face_media(i)%faces(j)%region_I_fields%E2
                E3 => conformal_media%face_media(i)%faces(j)%region_I_fields%E3
                E4 => conformal_media%face_media(i)%faces(j)%region_I_fields%E4        
                cell(:) = conformal_media%face_media(i)%faces(j)%cell(:)
                direction = conformal_media%face_media(j)%faces(k)%direction
                ! asignacion de medios: esta en otro sitio?  en este esquema, los medios this%g%gm1(medium)
                ! son = 0, porque se hace todo aquí. El numero de sggMiHx se mantiene, pero ahora hay un confg%g1 ?
                ! hay que crear medios para las complementarios de los medios conf?
                medium = media_maps%getMedium(cell, direction, found)
                if (.not. found) write(*,*) 'error'
                ! case(FACE_X)
                !     medium = media_maps%Hx%get(key(cell))
                !     ! medium = sggMiHx(cell(1), cell(2), cell(3))
                !     ! id1 = 
                !     ! id2 = 
                ! case(FACE_Y)
                !     medium =sggMiHy(cell(1), cell(2), cell(3))
                !     ! id1 = 
                !     ! id2 = 
                ! case(FACE_Z)
                !     medium =sggMiHz(cell(1), cell(2), cell(3))
                !     ! id1 = Idye(cell(2))
                !     ! id2 = Idxe(cell(1))
                ! end select
                ! i.e Hz, E2, E4 => Ex
                ! if Ex in region II, g1 = 0, g2 = 0, and field = 0
                H = g%gm1(medium)*H + g%gm2(medium)*((E2 - E4)*id1 - (E3-E1)*id2)
                ! H = this%g%gm1(medium)*H + this%g%gm2(medium)*((E2 - E4)*id1 - (E3-E1)*id2)
                ! update region II 
                if (isSurface) then 
                    rIIH  => conformal_media%face_media(i)%faces(j)%region_II_fields%H
                    rIIE1 => conformal_media%face_media(i)%faces(j)%region_II_fields%E1
                    rIIE2 => conformal_media%face_media(i)%faces(j)%region_II_fields%E2
                    rIIE3 => conformal_media%face_media(i)%faces(j)%region_II_fields%E3
                    rIIE4 => conformal_media%face_media(i)%faces(j)%region_II_fields%E4        

                    rIIH = this%g%gm1(medium)*rIIH + this%g%gm2(medium)*((rIIE2 - rIIE4)*id1 - (rIIE3-rIIE1)*id2)

                end if

            end do
        end do
    end subroutine

end module