module conformal_advance_m

    use conformal_fields_m
    use FDETypes_m, only: face_t, edge_t, FACE_X, FACE_Y, FACE_Z, EDGE_X, EDGE_Y, EDGE_Z, EDGE_REGION_CONFORMAL, EDGE_REGION_I, EDGE_REGION_II, constants_t, INTEGERSIZEOFMEDIAMATRICES

    implicit none

contains

    subroutine advanceConformalEdgesE(conformal_fields, Id_dir1, Id_dir2, sggMed, E_grid, g, direction)
        type(conformal_fields_t), intent(inout) :: conformal_fields
        real(kind=rkind), pointer, dimension(:) :: Id_dir1, Id_dir2
        integer(kind=INTEGERSIZEOFMEDIAMATRICES) , allocatable , dimension(:,:,:) :: sggMed
        real(kind=rkind) , pointer , dimension(:,:,:) :: E_grid
        type(constants_t), intent(in) :: g
        integer(kind=4), intent(in) :: direction
        
        real (kind = rkind), pointer :: E, H1, H2, H3, H4
        real (kind = rkind) :: E_eff, Id1, Id2

        real(kind=rkind) :: corr
        integer(kind=integersizeofmediamatrices) :: med

        integer(kind=4), dimension(3) :: c
        integer :: i, j, k, n, id
        do n = 1, size(conformal_fields%edges_id_list, 2)
            id = conformal_fields%edges_id_list(direction, n)
            if (id == 0) exit            
            corr = 0.0
            E_eff = 0.0
            c = conformal_fields%edges(id)%cell
            i = c(1); j = c(2); k = c(3)
            med = sggMed(i, j, k); 
            select case (direction)
                case (EDGE_X); Id1 = Id_dir1(k); Id2 = Id_dir2(j)
                case (EDGE_Y); Id1 = Id_dir1(i); Id2 = Id_dir2(k)
                case (EDGE_Z); Id1 = Id_dir1(j); Id2 = Id_dir2(i)
            end select
            !region 1
            if (conformal_fields%edges(id)%edge_region == EDGE_REGION_I .or. &
                conformal_fields%edges(id)%edge_region == EDGE_REGION_CONFORMAL) then 
                E  => conformal_fields%edges(id)%region_I_fields%E
                H1 => conformal_fields%edges(id)%region_I_fields%H1
                H2 => conformal_fields%edges(id)%region_I_fields%H2
                H3 => conformal_fields%edges(id)%region_I_fields%H3
                H4 => conformal_fields%edges(id)%region_I_fields%H4
                E = g%g1(med)*E + g%g2(med)*((H1 - H3)*Id2 - (H4-H2)*Id1)
                E_eff = E
                
            end if
            !region 2
            if (conformal_fields%edges(id)%edge_region == EDGE_REGION_II .or. &
                conformal_fields%edges(id)%edge_region == EDGE_REGION_CONFORMAL) then 
                E  => conformal_fields%edges(id)%region_II_fields%E
                H1 => conformal_fields%edges(id)%region_II_fields%H1
                H2 => conformal_fields%edges(id)%region_II_fields%H2
                H3 => conformal_fields%edges(id)%region_II_fields%H3
                H4 => conformal_fields%edges(id)%region_II_fields%H4
                    
                corr = computeEdgeRegionIICorrection(conformal_fields%edges(id))
                E = g%g1(med)*E + g%g2(med)*corr*((H1 - H3)*Id2 - (H4-H2)*Id1)
                E_eff = E_eff + E
            end if
            E_grid(i,j,k) = E_eff
        
        end do

    end subroutine

    subroutine advanceConformalFacesH(conformal_fields, Id_dir1, Id_dir2, sggMed, g, direction)
        type(conformal_fields_t), intent(inout) :: conformal_fields
        real(kind=rkind), pointer, dimension(:) :: Id_dir1, Id_dir2
        integer(kind=INTEGERSIZEOFMEDIAMATRICES), allocatable, dimension(:,:,:) :: sggMed
        type(constants_t), intent(in) :: g
        integer(kind=4), intent(in) :: direction

        real(kind=rkind), pointer :: H, E1, E2, E3, E4, Etan
        real(kind=rkind) :: Id1, Id2, corr
        integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: med
        integer(kind=4), dimension(3) :: c
        integer :: i, j, k, n, id

        do n = 1, size(conformal_fields%faces_id_list, 2)
            id = conformal_fields%faces_id_list(direction, n)
            if (id == 0) exit
            c = conformal_fields%faces(id)%cell
            i = c(1); j = c(2); k = c(3)
            med = sggMed(i, j, k)
            ! contour_length = conformal_fields%faces(id)%intersect%length
            select case (direction)
                case (FACE_X); Id1 = Id_dir1(k); Id2 = Id_dir2(j)
                case (FACE_Y); Id1 = Id_dir1(i); Id2 = Id_dir2(k)
                case (FACE_Z); Id1 = Id_dir1(j); Id2 = Id_dir2(i)
            end select

            ! region 1
            H  => conformal_fields%faces(id)%region_I_fields%H
            ! if (H/=0) then 
            !     call sleep(1)
            ! end if
            E1 => conformal_fields%faces(id)%region_I_fields%E1
            E2 => conformal_fields%faces(id)%region_I_fields%E2
            E3 => conformal_fields%faces(id)%region_I_fields%E3
            E4 => conformal_fields%faces(id)%region_I_fields%E4
            Etan => conformal_fields%faces(id)%region_I_fields%Etan
            H = g%gm1(med)*H + g%gm2(med)*((E2 - E4)*Id1 - (E3 - E1)*Id2 - Etan)

            ! region 2
            H  => conformal_fields%faces(id)%region_II_fields%H
            ! if (H/=0) then 
            !     call sleep(1)
            ! end if
            E1 => conformal_fields%faces(id)%region_II_fields%E1
            E2 => conformal_fields%faces(id)%region_II_fields%E2
            E3 => conformal_fields%faces(id)%region_II_fields%E3
            E4 => conformal_fields%faces(id)%region_II_fields%E4
            Etan => conformal_fields%faces(id)%region_II_fields%Etan
            corr = computeFaceRegionIICorrection(conformal_fields%faces(id))
            H = g%gm1(med)*H + g%gm2(med)*corr*((E2 - E4)*Id1 - (E3 - E1)*Id2  - Etan)
        end do
    end subroutine

    subroutine advanceConformalSGBCE(conformal_fields)
        type(conformal_fields_t), intent(inout) :: conformal_fields
        integer :: i
        do i = 1, size(conformal_fields%sgbcs)
            call conformal_fields%sgbcs(i)%state%advance()
        end do
    end subroutine

end module