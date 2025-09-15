module evolution_operator

    use Resuming
    use Solver_mod
    use fdetypes
    use Report
    use SEMBA_FDTD_mod
    use smbjson_testingTools

    use fhash,  key => fhash_key

    implicit none

    type :: field_array_t
        real(RKIND), allocatable, dimension(:,:,:) :: data
        character(len=2) :: field_type  ! 'Ex', 'Ey', 'Ez', 'Hx', etc.
    end type

    type :: int_array
        integer, allocatable :: data(:)
    end type

    private 

    public :: GenerateElectricalInputBasis,  GenerateMagneticalInputBasis, AddElectricFieldIndices, AddMagneticFieldIndices, fhash_get_int_array, int_array, GenerateRowIndexMap, get_field_bounds_from_json, GenerateOutputFields, field_array_t
    public :: GenerateColumnIndexMap, GenerateEvolutionOperator

contains

    subroutine GenerateElectricalInputBasis(M, dim1, dim2, M_ee, M_eo, M_oe, M_oo)

        real(RKIND), intent(in)  :: M(:,:,:)
        integer,     intent(in)  :: dim1, dim2  

        real(RKIND), intent(out) :: M_ee(size(M,1), size(M,2), size(M,3))
        real(RKIND), intent(out) :: M_eo(size(M,1), size(M,2), size(M,3))
        real(RKIND), intent(out) :: M_oe(size(M,1), size(M,2), size(M,3))
        real(RKIND), intent(out) :: M_oo(size(M,1), size(M,2), size(M,3))

        integer :: i, j, k
        integer :: ijk(3)
        integer :: v1, v2
        integer :: sz(3)
        sz = shape(M)

        M_ee = 0.0_RKIND
        M_eo = 0.0_RKIND
        M_oe = 0.0_RKIND
        M_oo = 0.0_RKIND

        do i = 0, sz(1)-1
            do j = 0, sz(2)-1
            do k = 0, sz(3)-1
                ijk = [i, j, k]
                v1 = ijk(dim1)
                v2 = ijk(dim2)

                select case (2*mod(v1,2) + mod(v2,2))
                case (0); M_ee(i+1,j+1,k+1) = 1.0_RKIND
                case (1); M_eo(i+1,j+1,k+1) = 1.0_RKIND
                case (2); M_oe(i+1,j+1,k+1) = 1.0_RKIND
                case (3); M_oo(i+1,j+1,k+1) = 1.0_RKIND
                end select
            end do
            end do
        end do
        end subroutine

    subroutine GenerateMagneticalInputBasis(A, M1, M2, M3, M)

        real(RKIND), intent(in)  :: A(:,:,:)
        integer, intent(in)      :: M1, M2, M3  

        ! Output
        real(RKIND), allocatable, intent(out) :: M(:,:,:,:,:,:)
        ! M(m1, m2, m3, i, j, k)

        integer :: i, j, k
        integer :: sz(3)
        integer :: idx1, idx2, idx3

        sz = shape(A)
        allocate(M(M1, M2, M3, sz(1), sz(2), sz(3)))
        M = 0.0_RKIND

        do i = 0, sz(1)-1
            do j = 0, sz(2)-1
            do k = 0, sz(3)-1
                idx1 = mod(i, M1)
                idx2 = mod(j, M2)
                idx3 = mod(k, M3)
                M(idx1+1, idx2+1, idx3+1, i+1, j+1, k+1) = 1.0_RKIND
            end do
            end do
        end do
        end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Creation of the basis for the input fields
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine GenerateInputFieldsBasis(b, FieldList)

        type (bounds_t), intent( IN)  ::  b
        type(field_array_t), allocatable, intent(OUT) :: FieldList(:)
        
        ! Generating the basis for the electical fields
        real (kind = RKIND), dimension    ( b%Ex%NX, b%Ex%NY, b%Ex%NZ) ::  Ex
        real (kind = RKIND), dimension    ( b%Ey%NX, b%Ey%NY, b%Ey%NZ) ::  Ey
        real (kind = RKIND), dimension    ( b%Ez%NX, b%Ez%NY, b%Ez%NZ) ::  Ez
        
        ! Generating the basis for the magnetical fields
        real (kind = RKIND), dimension    ( b%HX%NX, b%HX%NY, b%HX%NZ) ::  Hx
        real (kind = RKIND), dimension    ( b%Hy%NX, b%Hy%NY, b%Hy%NZ) ::  Hy
        real (kind = RKIND), dimension    ( b%Hz%NX, b%Hz%NY, b%Hz%NZ) ::  Hz

        ! Allocating the basis for the electrical fields
        real (kind = RKIND), dimension    (b%Ex%NX,b%Ex%NY, b%Ex%NZ ) :: Ex_ee
        real (kind = RKIND), dimension    (b%Ex%NX,b%Ex%NY, b%Ex%NZ ) :: Ex_eo
        real (kind = RKIND), dimension    (b%Ex%NX,b%Ex%NY, b%Ex%NZ ) :: Ex_oe
        real (kind = RKIND), dimension    (b%Ex%NX,b%Ex%NY, b%Ex%NZ ) :: Ex_oo
        real (kind = RKIND), dimension    (b%Ey%NX,b%Ey%NY, b%Ey%NZ ) :: Ey_ee
        real (kind = RKIND), dimension    (b%Ey%NX,b%Ey%NY, b%Ey%NZ ) :: Ey_eo
        real (kind = RKIND), dimension    (b%Ey%NX,b%Ey%NY, b%Ey%NZ ) :: Ey_oe
        real (kind = RKIND), dimension    (b%Ey%NX,b%Ey%NY, b%Ey%NZ ) :: Ey_oo
        real (kind = RKIND), dimension    (b%Ez%NX,b%Ez%NY, b%Ez%NZ ) :: Ez_ee
        real (kind = RKIND), dimension    (b%Ez%NX,b%Ez%NY, b%Ez%NZ ) :: Ez_eo
        real (kind = RKIND), dimension    (b%Ez%NX,b%Ez%NY, b%Ez%NZ ) :: Ez_oe
        real (kind = RKIND), dimension    (b%Ez%NX,b%Ez%NY, b%Ez%NZ ) :: Ez_oo

        ! Allocating the basis for the magnetical fields
        real (kind = RKIND), allocatable, dimension(:,:,:,:,:,:) :: Hx_m
        real (kind = RKIND), allocatable, dimension(:,:,:,:,:,:) :: Hy_m
        real (kind = RKIND), allocatable, dimension(:,:,:,:,:,:) :: Hz_m
        
        integer :: idx, i1, i2, i3

        allocate(FieldList(66))
        
        Ex = 0.0_RKIND
        Ey = 0.0_RKIND
        Ez = 0.0_RKIND

        Hx = 0.0_RKIND
        Hy = 0.0_RKIND
        Hz = 0.0_RKIND

        call GenerateElectricalInputBasis(Ex, 2, 3, Ex_ee, Ex_eo, Ex_oe, Ex_oo)
        call GenerateElectricalInputBasis(Ey, 1, 3, Ey_ee, Ey_eo, Ey_oe, Ey_oo)
        call GenerateElectricalInputBasis(Ez, 1, 2, Ez_ee, Ez_eo, Ez_oe, Ez_oo)

        ! Storing the electrical fields in the FieldList
        FieldList(1)%data = Ex_ee
        FieldList(2)%data = Ex_eo
        FieldList(3)%data = Ex_oe
        FieldList(4)%data = Ex_oo

        FieldList(5)%data = Ey_ee
        FieldList(6)%data = Ey_eo
        FieldList(7)%data = Ey_oe
        FieldList(8)%data = Ey_oo

        FieldList(9)%data  = Ez_ee
        FieldList(10)%data = Ez_eo
        FieldList(11)%data = Ez_oe
        FieldList(12)%data = Ez_oo


        call GenerateMagneticalInputBasis(Hx, 2, 3, 3, Hx_m)
        call GenerateMagneticalInputBasis(Hy, 3, 2, 3, Hy_m)
        call GenerateMagneticalInputBasis(Hz, 3, 3, 2, Hz_m)

        ! Storing the magnetical fields in the FieldList

        idx = 0

        do idx = 1, 12
            if (idx <= 4) then
                FieldList(idx)%field_type = 'Ex'
            else if (idx <= 8) then
                FieldList(idx)%field_type = 'Ey'
            else if (idx <= 12) then
                FieldList(idx)%field_type = 'Ez'
            end if
        end do

        do i1 = 1, 2
            do i2 = 1, 3
                do i3 = 1, 3
                    FieldList(idx)%data = Hx_m(i1, i2, i3, :, :, :)
                    FieldList(idx)%field_type = 'Hx'
                    idx = idx + 1 
                end do
            end do
        end do

        do i1 = 1, 3
            do i2 = 1, 2
                do i3 = 1, 3
                    FieldList(idx)%data = Hy_m(i1, i2, i3, :, :, :)
                    FieldList(idx)%field_type = 'Hy'
                    idx = idx + 1
                end do
            end do
        end do

        do i1 = 1, 3
            do i2 = 1, 3
                do i3 = 1, 2
                    FieldList(idx)%data = Hz_m(i1, i2, i3, :, :, :)
                    FieldList(idx)%field_type = 'Hz'
                    idx = idx + 1
                end do
            end do
        end do

    end subroutine 

    subroutine GenerateOutputFields(input_flags_no_json, input_file_json, fieldInput, fieldOutput)
        
        character (len=*), intent(in) :: input_file_json
        character (len=*), optional, intent(in) :: input_flags_no_json

        type(field_array_t), intent(in) :: fieldInput
        type(field_array_t), intent(out) :: fieldOutput

        type (bounds_t) ::  bounds
        type(semba_fdtd_t) :: semba
        type(solver_t) :: solver

        integer :: i, j, k
        integer, dimension(3) :: dims
        

        call semba%init(input_flags_no_json // ' ' // input_file_json)
        call solver%init_control(semba%l, semba%maxSourceValue, semba%time_desdelanzamiento)
        call solver%init(semba%sgg,semba%eps0, semba%mu0, semba%sggMiNo,& 
                            semba%sggMiEx,semba%sggMiEy,semba%sggMiEz,& 
                            semba%sggMiHx,semba%sggMiHy,semba%sggMiHz, & 
                            semba%sggMtag, semba%SINPML_fullsize, semba%fullsize, semba%tag_numbers)


        call get_field_bounds_from_json(bounds, semba%fullsize)

        
        dims = shape(fieldInput%data)
        
        allocate(fieldOutput%data( &
            size(fieldInput%data, 1), &
            size(fieldInput%data, 2), &
            size(fieldInput%data, 3)))

        fieldOutput%data = 0.0_RKIND
        
        do i = 1, dims(1)
            do j = 1, dims(2)
                do k = 1, dims(3)
                    select case (fieldInput%field_type)
                    case ('Ex')
                        call solver%set_field_value(iEx, [i-1,i-1], [j-1,j-1], [k-1,k-1], fieldInput%data(i,j,k))
                    case ('Ey')
                        call solver%set_field_value(iEy, [i-1,i-1], [j-1,j-1], [k-1,k-1], fieldInput%data(i,j,k))
                    case ('Ez')
                        call solver%set_field_value(iEz, [i-1,i-1], [j-1,j-1], [k-1,k-1], fieldInput%data(i,j,k))
                    case ('Hx')
                        call solver%set_field_value(iHx, [i-1,i-1], [j-1,j-1], [k-1,k-1], fieldInput%data(i,j,k))
                    case ('Hy')
                        call solver%set_field_value(iHy, [i-1,i-1], [j-1,j-1], [k-1,k-1], fieldInput%data(i,j,k))
                    case ('Hz')
                        call solver%set_field_value(iHz, [i-1,i-1], [j-1,j-1], [k-1,k-1], fieldInput%data(i,j,k))
                    end select
                end do
            end do
        end do

        call solver%step(semba%sgg, semba%eps0, semba%mu0, semba%SINPML_FULLSIZE, semba%tag_numbers)

        ! This is wrong, for example with Ex as imput, we should have Ex, Hy and Hz as output
        do i = 1, dims(1)
            do j = 1, dims(2)
                do k = 1, dims(3)
                    select case (fieldInput%field_type)
                    case ('Ex')
                        fieldOutput%data(i, j, k) = solver%get_field_value(iEx, i-1, j-1, k-1)
                    case ('Ey')
                        fieldOutput%data(i, j, k) = solver%get_field_value(iEy, i-1, j-1, k-1)
                    case ('Ez')
                        fieldOutput%data(i, j, k) = solver%get_field_value(iEz, i-1, j-1, k-1)
                    case ('Hx')
                        fieldOutput%data(i, j, k) = solver%get_field_value(iHx, i-1, j-1, k-1)
                    case ('Hy')
                        fieldOutput%data(i, j, k) = solver%get_field_value(iHy, i-1, j-1, k-1)
                    case ('Hz')
                        fieldOutput%data(i, j, k) = solver%get_field_value(iHz, i-1, j-1, k-1)
                    end select
                end do
            end do
        end do

        call solver%set_field_value(iEx, [solver%bounds%Ex%xi, solver%bounds%Ex%xe], [solver%bounds%Ex%yi, solver%bounds%Ex%ye], [solver%bounds%Ex%zi, solver%bounds%Ex%ze], 0.0)
        call solver%set_field_value(iEy, [solver%bounds%Ey%xi, solver%bounds%Ey%xe], [solver%bounds%Ey%yi, solver%bounds%Ey%ye], [solver%bounds%Ey%zi, solver%bounds%Ey%ze], 0.0)
        call solver%set_field_value(iEz, [solver%bounds%Ez%xi, solver%bounds%Ez%xe], [solver%bounds%Ez%yi, solver%bounds%Ez%ye], [solver%bounds%Ez%zi, solver%bounds%Ez%ze], 0.0)
        call solver%set_field_value(iHx, [solver%bounds%Hx%xi, solver%bounds%Hx%xe], [solver%bounds%Hx%yi, solver%bounds%Hx%ye], [solver%bounds%Hx%zi, solver%bounds%Hx%ze], 0.0)
        call solver%set_field_value(iHy, [solver%bounds%Hy%xi, solver%bounds%Hy%xe], [solver%bounds%Hy%yi, solver%bounds%Hy%ye], [solver%bounds%Hy%zi, solver%bounds%Hy%ze], 0.0)
        call solver%set_field_value(iHz, [solver%bounds%Hz%xi, solver%bounds%Hz%xe], [solver%bounds%Hz%yi, solver%bounds%Hz%ye], [solver%bounds%Hz%zi, solver%bounds%Hz%ze], 0.0)

    end subroutine

    subroutine AddElectricFieldIndices(RowIndexMap, Efield, H1field, H2field, startingIndex_Efield, startingIndex_H1field, startingIndex_H2field, shiftDirection_H1, shiftDirection_H2)
        type(fhash_tbl_t), intent(inout) :: RowIndexMap
        type(limit_t), intent(in) :: Efield, H1field, H2field
        integer, intent(in) :: startingIndex_Efield, startingIndex_H1field, startingIndex_H2field
        character(len=1), intent(in) :: shiftDirection_H1, shiftDirection_H2  

        integer :: i, j, k
        integer :: combinedIndex_E, combinedIndex_H1, combinedIndex_H2
        integer :: indexShift_H1, indexShift_H2, auxiliarIndexShift_H1, auxiliarIndexShift_H2

        type(int_array) :: wrapper 
        integer, allocatable :: indexList(:)
        integer :: countIndex, positionList


        do i = 1, Efield%Nx
            do j = 1, Efield%Ny
                do k = 1, Efield%Nz
                    combinedIndex_E  = ((i - 1)*Efield%Ny  + (j - 1))*Efield%Nz   + k
                    combinedIndex_H1 = ((i - 1)*H1field%Ny + (j - 1))*H1field%Nz  + k
                    combinedIndex_H2 = ((i - 1)*H2field%Ny + (j - 1))*H2field%Nz  + k

                    countIndex = 1
                    positionList = 1


                    select case (shiftDirection_H1)
                        case ('i')
                            if (i > 1 .and. i < Efield%Nx) then
                                indexShift_H1 = ((i - 2)*H1field%Ny + (j - 1))*H1field%Nz  + k
                                countIndex = countIndex + 2
                            else if (i == 1) then
                                indexShift_H1 = -1
                                countIndex = countIndex + 1
                            else
                                indexShift_H1 = -2
                                auxiliarIndexShift_H1 = ((i - 2)*H1field%Ny + (j - 1))*H1field%Nz  + k
                                countIndex = countIndex + 1
                            end if
                        case ('j')
                            if (j > 1 .and. j < Efield%Ny) then
                                indexShift_H1 = ((i - 1)*H1field%Ny + (j - 2))*H1field%Nz  + k
                                countIndex = countIndex + 2
                            else if (j == 1) then
                                indexShift_H1 = -1
                                countIndex = countIndex + 1
                            else
                                indexShift_H1 = -2
                                auxiliarIndexShift_H1 = ((i - 1)*H1field%Ny + (j - 2))*H1field%Nz  + k
                                countIndex = countIndex + 1
                            end if
                        case ('k')
                            if (k > 1 .and. k < Efield%Nz) then
                                indexShift_H1 = ((i - 1)*H1field%Ny + (j - 1))*H1field%Nz  + (k - 1)
                                countIndex = countIndex + 2
                            else if (k == 1) then
                                indexShift_H1 = -1
                                countIndex = countIndex + 1
                            else
                                indexShift_H1 = -2
                                auxiliarIndexShift_H1 = ((i - 1)*H1field%Ny + (j - 1))*H1field%Nz  + (k - 1)
                                countIndex = countIndex + 1
                            end if
                    end select

                    select case (shiftDirection_H2)
                        case ('i')
                            if (i > 1 .and. i < Efield%Nx) then
                                indexShift_H2 = ((i - 2)*H2field%Ny + (j - 1))*H2field%Nz  + k
                                countIndex = countIndex + 2
                            else if (i == 1) then
                                indexShift_H2 = -1
                                countIndex = countIndex + 1
                            else
                                indexShift_H2 = -2
                                auxiliarIndexShift_H2 = ((i - 2)*H2field%Ny + (j - 1))*H2field%Nz  + k
                                countIndex = countIndex + 1
                            end if
                        case ('j')
                            if (j > 1 .and. j < Efield%Ny) then
                                indexShift_H2 = ((i - 1)*H2field%Ny + (j - 2))*H2field%Nz  + k
                                countIndex = countIndex + 2
                            else if (j == 1) then
                                indexShift_H2 = -1
                                countIndex = countIndex + 1
                            else
                                indexShift_H2 = -2
                                auxiliarIndexShift_H2 = ((i - 1)*H2field%Ny + (j - 2))*H2field%Nz  + k
                                countIndex = countIndex + 1
                            end if
                        case ('k')
                            if (k > 1 .and. k < Efield%Nz) then
                                indexShift_H2 = ((i - 1)*H2field%Ny + (j - 1))*H2field%Nz  + (k - 1)
                                countIndex = countIndex + 2
                            else if (k == 1) then
                                indexShift_H2 = -1
                                countIndex = countIndex + 1
                            else
                                indexShift_H2 = -2
                                auxiliarIndexShift_H2 = ((i - 1)*H2field%Ny + (j - 1))*H2field%Nz  + (k - 1)
                                countIndex = countIndex + 1
                            end if
                    end select


                    allocate(indexList(countIndex))
                    indexList(positionList) = startingIndex_Efield + combinedIndex_E
                    positionList = positionList + 1


                    if (indexShift_H2 /= -1 .and. indexShift_H2 /= -2) then
                        indexList(positionList)     = startingIndex_H2field + combinedIndex_H2
                        indexList(positionList + 1) = startingIndex_H2field + indexShift_H2
                        positionList = positionList + 2
                    else if (indexShift_H2 == -1) then   ! Border at the beginning
                        indexList(positionList) = startingIndex_H2field + combinedIndex_H2
                        positionList = positionList + 1
                    else                            ! Border at the end
                        indexList(positionList) = startingIndex_H2field + auxiliarIndexShift_H2
                        positionList = positionList + 1
                    end if

                    if (indexShift_H1 /= -1 .and. indexShift_H1 /= -2) then
                        indexList(positionList)     = startingIndex_H1field + combinedIndex_H1
                        indexList(positionList + 1) = startingIndex_H1field + indexShift_H1
                        positionList = positionList + 2
                    else if (indexShift_H1 == -1) then   ! Border at the beginning
                        indexList(positionList) = startingIndex_H1field + combinedIndex_H1
                        positionList = positionList + 1
                    else                            ! Border at the end
                        indexList(positionList) = startingIndex_H1field + auxiliarIndexShift_H1
                        positionList = positionList + 1
                    end if


                    wrapper%data = indexList
                    call RowIndexMap%set(key(startingIndex_Efield + combinedIndex_E), value=wrapper)  
                    deallocate(indexList)
                end do
            end do
        end do
    end subroutine 

    subroutine AddMagneticFieldIndices(RowIndexMap, Hfield, E1field, E2field, startingIndex_Hfield, startingIndex_E1field, startingIndex_E2field, shiftDirection_E1, shiftDirection_E2)
        type(fhash_tbl_t), intent(inout) :: RowIndexMap
        type(limit_t), intent(in) :: Hfield, E1field, E2field
        integer, intent(in) :: startingIndex_Hfield, startingIndex_E1field, startingIndex_E2field
        character(len=1), intent(in) :: shiftDirection_E1, shiftDirection_E2

        integer :: i, j, k
        integer :: combinedIndex_H, combinedIndex_E1, combinedIndex_E2
        integer :: indexShift_E1, indexShift_E2

        integer, allocatable :: fullIndexList(:), uniqueIndexList(:)
        type(int_array) :: relatedIndices_E1, relatedIndices_E1_shift, relatedIndices_E2, relatedIndices_E2_shift, wrapper
        integer :: relatedIndices_maximumSize


        do i = 1, Hfield%Nx
            do j = 1, Hfield%Ny
                do k = 1, Hfield%Nz
                    combinedIndex_H  = ((i - 1)*Hfield%Ny  + (j - 1))*Hfield%Nz   + k
                    combinedIndex_E1 = ((i - 1)*E1field%Ny + (j - 1))*E1field%Nz  + k
                    combinedIndex_E2 = ((i - 1)*E2field%Ny + (j - 1))*E2field%Nz  + k

                    select case (shiftDirection_E1)
                        case ('i')
                            indexShift_E1 = ( i     *E1field%Ny + (j - 1))*E1field%Nz  + k
                        case ('j')
                            indexShift_E1 = ((i - 1)*E1field%Ny + j      )*E1field%Nz  + k
                        case ('k')
                            indexShift_E1 = ((i - 1)*E1field%Ny + (j - 1))*E1field%Nz  + (k + 1)
                    end select

                    select case (shiftDirection_E2)
                        case ('i')
                            indexShift_E2 = ( i     *E2field%Ny + (j - 1))*E2field%Nz  + k
                        case ('j')
                            indexShift_E2 = ((i - 1)*E2field%Ny + j      )*E2field%Nz  + k
                        case ('k')
                            indexShift_E2 = ((i - 1)*E2field%Ny + (j - 1))*E2field%Nz  + (k + 1)
                    end select

                    call fhash_get_int_array(RowIndexMap, key(startingIndex_E1field + combinedIndex_E1), relatedIndices_E1)
                    call fhash_get_int_array(RowIndexMap, key(startingIndex_E1field + indexShift_E1   ), relatedIndices_E1_shift)
                    call fhash_get_int_array(RowIndexMap, key(startingIndex_E2field + combinedIndex_E2), relatedIndices_E2)
                    call fhash_get_int_array(RowIndexMap, key(startingIndex_E2field + indexShift_E2   ), relatedIndices_E2_shift)


                    relatedIndices_maximumSize = size(relatedIndices_E1%data) + size(relatedIndices_E1_shift%data) + size(relatedIndices_E2%data) + size(relatedIndices_E2_shift%data)
                    allocate(fullIndexList(relatedIndices_maximumSize))

                    fullIndexList(1:size(relatedIndices_E1%data)) = relatedIndices_E1%data
                    fullIndexList(size(relatedIndices_E1%data)+1:size(relatedIndices_E1%data)+size(relatedIndices_E1_shift%data)) = relatedIndices_E1_shift%data
                    fullIndexList(size(relatedIndices_E1%data)+size(relatedIndices_E1_shift%data)+1:size(relatedIndices_E1%data)+size(relatedIndices_E1_shift%data)+size(relatedIndices_E2%data)) = relatedIndices_E2%data
                    fullIndexList(size(relatedIndices_E1%data)+size(relatedIndices_E1_shift%data)+size(relatedIndices_E2%data)+1:) = relatedIndices_E2_shift%data

                    call RemoveDuplicates(fullIndexList, uniqueIndexList)

                    wrapper%data = uniqueIndexList

                    call RowIndexMap%set(key(startingIndex_Hfield + combinedIndex_H), value=wrapper)

                    deallocate(fullIndexList, uniqueIndexList)
                end do
            end do
        end do
    end subroutine
        
    subroutine RemoveDuplicates(inputArray, outputArray)
        integer, intent(in) :: inputArray(:)
        integer, allocatable, intent(out) :: outputArray(:)
        integer :: i, j, n
        logical :: found
        integer, allocatable :: fullIndexList(:)

        allocate(fullIndexList(size(inputArray)))
        n = 0

        do i = 1, size(inputArray)
            found = .false.
            do j = 1, n
                if (fullIndexList(j) == inputArray(i)) then
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                n = n + 1
                fullIndexList(n) = inputArray(i)
            end if
        end do

        allocate(outputArray(n))
        outputArray = fullIndexList(1:n)
        deallocate(fullIndexList)
    end subroutine 

    subroutine GenerateRowIndexMap(bounds, RowIndexMap)

        type(bounds_t), intent(IN) :: bounds
        type(fhash_tbl_t), intent(OUT) :: RowIndexMap
        integer :: shiftEx, shiftEy, shiftEz, shiftHx, shiftHy, shiftHz

        shiftEx = 0
        shiftEy = shiftEx + bounds%Ex%Nx * bounds%Ex%Ny * bounds%Ex%Nz
        shiftEz = shiftEy + bounds%Ey%Nx * bounds%Ey%Ny * bounds%Ey%Nz
        shiftHx = shiftEz + bounds%Ez%Nx * bounds%Ez%Ny * bounds%Ez%Nz
        shiftHy = shiftHx + bounds%Hx%Nx * bounds%Hx%Ny * bounds%Hx%Nz
        shiftHz = shiftHy + bounds%Hy%Nx * bounds%Hy%Ny * bounds%Hy%Nz


        call AddElectricFieldIndices(RowIndexMap, bounds%Ex, bounds%Hy, bounds%Hz, shiftEx, shiftHy, shiftHz, 'k', 'j')
        call AddElectricFieldIndices(RowIndexMap, bounds%Ey, bounds%Hx, bounds%Hz, shiftEy, shiftHx, shiftHz, 'k', 'i')
        call AddElectricFieldIndices(RowIndexMap, bounds%Ez, bounds%Hx, bounds%Hy, shiftEz, shiftHx, shiftHy, 'j', 'i')

        call AddMagneticFieldIndices(RowIndexMap, bounds%Hx, bounds%Ez, bounds%Ey, shiftHx, shiftEz, shiftEy, 'j', 'k')
        call AddMagneticFieldIndices(RowIndexMap, bounds%Hy, bounds%Ex, bounds%Ez, shiftHy, shiftEx, shiftEz, 'k', 'i')
        call AddMagneticFieldIndices(RowIndexMap, bounds%Hz, bounds%Ex, bounds%Ey, shiftHz, shiftEx, shiftEy, 'j', 'i')


    end subroutine GenerateRowIndexMap

    subroutine GenerateColumnIndexMap(bounds, ColIndexMap)

        type(bounds_t), intent(IN) :: bounds
        type(fhash_tbl_t), intent(OUT) :: ColIndexMap
        type(fhash_tbl_t) :: RowIndexMap
        integer :: m1, m2, dataIdx, listPosition, countSize, totalElements
        type(int_array) :: wrapper1, wrapper2, wrapperColumn
        integer, allocatable :: tempData(:)


        totalElements = bounds%Ex%NX * bounds%Ex%NY * bounds%Ex%NZ + &
                        bounds%Ey%NX * bounds%Ey%NY * bounds%Ey%NZ + &
                        bounds%Ez%NX * bounds%Ez%NY * bounds%Ez%NZ + &
                        bounds%Hx%NX * bounds%Hx%NY * bounds%Hx%NZ + &
                        bounds%Hy%NX * bounds%Hy%NY * bounds%Hy%NZ + &
                        bounds%Hz%NX * bounds%Hz%NY * bounds%Hz%NZ

        call GenerateRowIndexMap(bounds, RowIndexMap)

        do m1 = 1, totalElements
            call fhash_get_int_array(RowIndexMap, key(m1), wrapper1)

            do dataIdx = 1, size(wrapper1%data)
                listPosition = 1
                countSize = 0

                do m2 = 1, totalElements
                    call fhash_get_int_array(RowIndexMap, key(m2), wrapper2)

                    if (any(wrapper2%data == wrapper1%data(dataIdx))) then
                        countSize = countSize + 1
                    end if
                end do

                allocate(tempData(countSize))

                do m2 = 1, totalElements
                    call fhash_get_int_array(RowIndexMap, key(m2), wrapper2)

                    if (any(wrapper2%data == wrapper1%data(dataIdx))) then
                        tempData(listPosition) = m2
                        listPosition = listPosition + 1
                    end if
                end do

                wrapperColumn%data = tempData
                call ColIndexMap%set(key(wrapper1%data(dataIdx)), value=wrapperColumn)
                deallocate(tempData)
            end do
        end do

    end subroutine GenerateColumnIndexMap

    subroutine fhash_get_int_array(tbl, k, val)
        type(fhash_tbl_t), intent(in) :: tbl
        class(fhash_key_t), intent(in) :: k
        type(int_array), intent(out) :: val

        integer :: stat
        class(*), allocatable :: raw

        call tbl%get_raw(k, raw, stat)

        if (stat /= 0) then
            allocate(val%data(0))
            return
        end if

        select type(d => raw)
        type is (int_array)
            val = d
        class default
            allocate(val%data(0))
        end select
    end subroutine

    subroutine get_field_bounds_from_json(field_bounds, fullsize)
        type(bounds_t), intent(out) :: field_bounds
        TYPE (limit_t), DIMENSION (1:6) :: fullsize
        type(integer) :: Nx, Ny, Nz

        Nx = fullsize(1)%xe - fullsize(1)%xi + 1
        Ny = fullsize(2)%ye - fullsize(2)%yi + 1
        Nz = fullsize(3)%ze - fullsize(3)%zi + 1

        field_bounds%Ex%NX = nX
        field_bounds%Ex%NY = nY + 1
        field_bounds%Ex%NZ = nZ + 1

        field_bounds%Ey%NX = nX + 1
        field_bounds%Ey%NY = nY
        field_bounds%Ey%NZ = nZ + 1

        field_bounds%Ez%NX = nX + 1
        field_bounds%Ez%NY = nY + 1
        field_bounds%Ez%NZ = nZ

        field_bounds%Hx%NX = field_bounds%Ex%Nx + 1
        field_bounds%Hx%NY = field_bounds%Ex%Ny - 1
        field_bounds%Hx%NZ = field_bounds%Ex%Nz - 1

        field_bounds%Hy%NX = field_bounds%Ey%Nx - 1
        field_bounds%Hy%NY = field_bounds%Ey%Ny + 1
        field_bounds%Hy%NZ = field_bounds%Ey%Nz - 1

        field_bounds%Hz%NX = field_bounds%Ez%Nx - 1
        field_bounds%Hz%NY = field_bounds%Ez%Ny - 1
        field_bounds%Hz%NZ = field_bounds%Ez%Nz + 1
    end subroutine

    subroutine GenerateEvolutionOperator(input_flags_no_json, input_file_json, evolutionOperator)
        
        character (len=*), intent(in) :: input_file_json
        character (len=*), optional, intent(in) :: input_flags_no_json
        real (kind = RKIND), allocatable, dimension(:, :), intent(out) ::  evolutionOperator

        type(field_array_t), allocatable :: fieldInputList(:)
        type(int_array) :: wrapper

        type (bounds_t) ::  bounds
        type(fhash_tbl_t) :: ColIndexMap
        type(semba_fdtd_t) :: semba
        type(solver_t) :: solver

        integer :: shiftEx, shiftEy, shiftEz, shiftHx, shiftHy, shiftHz
        integer :: i, j, k, m, totalElements, fieldIdx
        integer :: i_rel, j_rel, k_rel, m_rel, wrapperIdx
        real(kind = RKIND) :: fieldValue
        integer, dimension(3) :: dims
        

        call semba%init(input_flags_no_json // ' ' // input_file_json)
        call solver%init_control(semba%l, semba%maxSourceValue, semba%time_desdelanzamiento)
        call solver%init(semba%sgg,semba%eps0, semba%mu0, semba%sggMiNo,& 
                            semba%sggMiEx,semba%sggMiEy,semba%sggMiEz,& 
                            semba%sggMiHx,semba%sggMiHy,semba%sggMiHz, & 
                            semba%sggMtag, semba%SINPML_fullsize, semba%fullsize, semba%tag_numbers)


        call get_field_bounds_from_json(bounds, semba%fullsize)

        shiftEx = 0
        shiftEy = 0       + bounds%Ex%Nx * bounds%Ex%Ny * bounds%Ex%Nz
        shiftEz = shiftEy + bounds%Ey%Nx * bounds%Ey%Ny * bounds%Ey%Nz
        shiftHx = shiftEz + bounds%Ez%Nx * bounds%Ez%Ny * bounds%Ez%Nz
        shiftHy = shiftHx + bounds%Hx%Nx * bounds%Hx%Ny * bounds%Hx%Nz
        shiftHz = shiftHy + bounds%Hy%Nx * bounds%Hy%Ny * bounds%Hy%Nz

        totalElements = bounds%Ex%NX * bounds%Ex%NY * bounds%Ex%NZ + &
                        bounds%Ey%NX * bounds%Ey%NY * bounds%Ey%NZ + &
                        bounds%Ez%NX * bounds%Ez%NY * bounds%Ez%NZ + &
                        bounds%Hx%NX * bounds%Hx%NY * bounds%Hx%NZ + &
                        bounds%Hy%NX * bounds%Hy%NY * bounds%Hy%NZ + &
                        bounds%Hz%NX * bounds%Hz%NY * bounds%Hz%NZ

        allocate(evolutionOperator(totalElements, totalElements))
        evolutionOperator = 0.0_RKIND

        call GenerateColumnIndexMap(bounds, ColIndexMap)
        call GenerateInputFieldsBasis(bounds, fieldInputList)

        do fieldIdx = 1, size(fieldInputList)
            dims = shape(fieldInputList(fieldIdx)%data)
            
            do i = 1, dims(1)
                do j = 1, dims(2)
                    do k = 1, dims(3)
                        select case (fieldInputList(fieldIdx)%field_type)
                        case ('Ex')
                            call solver%set_field_value(iEx, [i-1,i-1], [j-1,j-1], [k-1,k-1], fieldInputList(fieldIdx)%data(i,j,k))
                        case ('Ey')
                            call solver%set_field_value(iEy, [i-1,i-1], [j-1,j-1], [k-1,k-1], fieldInputList(fieldIdx)%data(i,j,k))
                        case ('Ez')
                            call solver%set_field_value(iEz, [i-1,i-1], [j-1,j-1], [k-1,k-1], fieldInputList(fieldIdx)%data(i,j,k))
                        case ('Hx')
                            call solver%set_field_value(iHx, [i-1,i-1], [j-1,j-1], [k-1,k-1], fieldInputList(fieldIdx)%data(i,j,k))
                        case ('Hy')
                            call solver%set_field_value(iHy, [i-1,i-1], [j-1,j-1], [k-1,k-1], fieldInputList(fieldIdx)%data(i,j,k))
                        case ('Hz')
                            call solver%set_field_value(iHz, [i-1,i-1], [j-1,j-1], [k-1,k-1], fieldInputList(fieldIdx)%data(i,j,k))
                        end select
                    end do
                end do
            end do
            

            call solver%step(semba%sgg, semba%eps0, semba%mu0, semba%SINPML_FULLSIZE, semba%tag_numbers)
            
            
            do i = 1, dims(1)
                do j = 1, dims(2)
                    do k = 1, dims(3)
                        m = ((i - 1)*dims(2) + (j - 1))*dims(3) + k

                        select case (fieldInputList(fieldIdx)%field_type)
                        case ('Ex')
                            call fhash_get_int_array(ColIndexMap, key(shiftEx + m), wrapper)
                        case ('Ey')
                            call fhash_get_int_array(ColIndexMap, key(shiftEy + m), wrapper)
                        case ('Ez')
                            call fhash_get_int_array(ColIndexMap, key(shiftEz + m), wrapper)
                        case ('Hx')
                            call fhash_get_int_array(ColIndexMap, key(shiftHx + m), wrapper)
                        case ('Hy')
                            call fhash_get_int_array(ColIndexMap, key(shiftHy + m), wrapper)
                        case ('Hz')
                            call fhash_get_int_array(ColIndexMap, key(shiftHz + m), wrapper)
                        end select

                        do wrapperIdx = 1, size(wrapper%data)
                            m_rel = wrapper%data(wrapperIdx)

                            select case (fieldInputList(fieldIdx)%field_type)
                            case ('Ex')
                                k_rel = mod(m_rel - shiftEx - 1, bounds%Ex%Nz) + 1
                                j_rel = mod((m_rel - shiftEx - 1) / bounds%Ex%Nz, bounds%Ex%Ny) + 1
                                i_rel = (m_rel - shiftEx - 1) / (bounds%Ex%Nz * bounds%Ex%Ny) + 1

                                fieldValue = solver%get_field_value(iEx, i_rel - 1, j_rel - 1, k_rel - 1)
                            case ('Ey')
                                k_rel = mod(m_rel - shiftEy - 1, bounds%Ey%Nz) + 1
                                j_rel = mod((m_rel - shiftEy - 1) / bounds%Ey%Nz, bounds%Ey%Ny) + 1
                                i_rel = (m_rel - shiftEy - 1) / (bounds%Ey%Nz * bounds%Ey%Ny) + 1

                                fieldValue = solver%get_field_value(iEy, i_rel - 1, j_rel - 1, k_rel - 1)
                            case ('Ez')
                                k_rel = mod(m_rel - shiftEz - 1, bounds%Ez%Nz) + 1
                                j_rel = mod((m_rel - shiftEz - 1) / bounds%Ez%Nz, bounds%Ez%Ny) + 1
                                i_rel = (m_rel - shiftEz - 1) / (bounds%Ez%Nz * bounds%Ez%Ny) + 1

                                fieldValue = solver%get_field_value(iEz, i_rel - 1, j_rel - 1, k_rel - 1)
                            case ('Hx')
                                k_rel = mod(m_rel - shiftHx - 1, bounds%Hx%Nz) + 1
                                j_rel = mod((m_rel - shiftHx - 1) / bounds%Hx%Nz, bounds%Hx%Ny) + 1
                                i_rel = (m_rel - shiftHx - 1) / (bounds%Hx%Nz * bounds%Hx%Ny) + 1

                                fieldValue = solver%get_field_value(iHx, i_rel - 1, j_rel - 1, k_rel - 1)
                            case ('Hy')
                                k_rel = mod(m_rel - shiftHy - 1, bounds%Hy%Nz) + 1
                                j_rel = mod((m_rel - shiftHy - 1) / bounds%Hy%Nz, bounds%Hy%Ny) + 1
                                i_rel = (m_rel - shiftHy - 1) / (bounds%Hy%Nz * bounds%Hy%Ny) + 1

                                fieldValue = solver%get_field_value(iHy, i_rel - 1, j_rel - 1, k_rel - 1)
                            case ('Hz')
                                k_rel = mod(m_rel - shiftHz - 1, bounds%Hz%Nz) + 1
                                j_rel = mod((m_rel - shiftHz - 1) / bounds%Hz%Nz, bounds%Hz%Ny) + 1
                                i_rel = (m_rel - shiftHz - 1) / (bounds%Hz%Nz * bounds%Hz%Ny) + 1

                                fieldValue = solver%get_field_value(iHz, i_rel - 1, j_rel - 1, k_rel - 1)
                            end select

                            evolutionOperator(m_rel, m) = fieldValue
                        end do
                        
                    end do
                end do
            end do

            call solver%set_field_value(iEx, [solver%bounds%Ex%xi, solver%bounds%Ex%xe], [solver%bounds%Ex%yi, solver%bounds%Ex%ye], [solver%bounds%Ex%zi, solver%bounds%Ex%ze], 0.0)
            call solver%set_field_value(iEy, [solver%bounds%Ey%xi, solver%bounds%Ey%xe], [solver%bounds%Ey%yi, solver%bounds%Ey%ye], [solver%bounds%Ey%zi, solver%bounds%Ey%ze], 0.0)
            call solver%set_field_value(iEz, [solver%bounds%Ez%xi, solver%bounds%Ez%xe], [solver%bounds%Ez%yi, solver%bounds%Ez%ye], [solver%bounds%Ez%zi, solver%bounds%Ez%ze], 0.0)
            call solver%set_field_value(iHx, [solver%bounds%Hx%xi, solver%bounds%Hx%xe], [solver%bounds%Hx%yi, solver%bounds%Hx%ye], [solver%bounds%Hx%zi, solver%bounds%Hx%ze], 0.0)
            call solver%set_field_value(iHy, [solver%bounds%Hy%xi, solver%bounds%Hy%xe], [solver%bounds%Hy%yi, solver%bounds%Hy%ye], [solver%bounds%Hy%zi, solver%bounds%Hy%ze], 0.0)
            call solver%set_field_value(iHz, [solver%bounds%Hz%xi, solver%bounds%Hz%xe], [solver%bounds%Hz%yi, solver%bounds%Hz%ye], [solver%bounds%Hz%zi, solver%bounds%Hz%ze], 0.0)
            
        end do
    end subroutine

end module