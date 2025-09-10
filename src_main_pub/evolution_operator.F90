module evolution_operator

    use Resuming
    use Solver_mod
    use fdetypes
    use Report

    use fhash,  key => fhash_key

    implicit none

    type :: field_array_t
        real(RKIND), pointer, dimension(:,:,:) :: data
        character(len=2) :: field_type  ! 'Ex', 'Ey', 'Ez', 'Hx', etc.
    end type

    type :: int_array
        integer, allocatable :: data(:)
    end type

    private 

    public :: GenerateElectricalInputBasis,  GenerateMagneticalInputBasis, AddElectricFieldIndices, AddMagneticFieldIndices, fhash_get_int_array, int_array, GenerateRowIndexMap

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
        real (kind = RKIND), dimension    ( 0 :      b%Ex%NX-1 , 0 :      b%Ex%NY-1 , 0 :      b%Ex%NZ-1 ) ::  Ex
        real (kind = RKIND), dimension    ( 0 :      b%Ey%NX-1 , 0 :      b%Ey%NY-1 , 0 :      b%Ey%NZ-1 ) ::  Ey
        real (kind = RKIND), dimension    ( 0 :      b%Ez%NX-1 , 0 :      b%Ez%NY-1 , 0 :      b%Ez%NZ-1 ) ::  Ez
        
        ! Generating the basis for the magnetical fields
        real (kind = RKIND), dimension    ( 0 :      b%HX%NX-1 , 0 :      b%HX%NY-1 , 0 :      b%HX%NZ-1 ) ::  Hx
        real (kind = RKIND), dimension    ( 0 :      b%Hy%NX-1 , 0 :      b%Hy%NY-1 , 0 :      b%Hy%NZ-1 ) ::  Hy
        real (kind = RKIND), dimension    ( 0 :      b%Hz%NX-1 , 0 :      b%Hz%NY-1 , 0 :      b%Hz%NZ-1 ) ::  Hz

        ! Allocating the basis for the electrical fields
        real (kind = RKIND), allocatable, dimension(:,:,:) :: Ex_ee
        real (kind = RKIND), allocatable, dimension(:,:,:) :: Ex_eo
        real (kind = RKIND), allocatable, dimension(:,:,:) :: Ex_oe
        real (kind = RKIND), allocatable, dimension(:,:,:) :: Ex_oo
        real (kind = RKIND), allocatable, dimension(:,:,:) :: Ey_ee
        real (kind = RKIND), allocatable, dimension(:,:,:) :: Ey_eo
        real (kind = RKIND), allocatable, dimension(:,:,:) :: Ey_oe
        real (kind = RKIND), allocatable, dimension(:,:,:) :: Ey_oo
        real (kind = RKIND), allocatable, dimension(:,:,:) :: Ez_ee
        real (kind = RKIND), allocatable, dimension(:,:,:) :: Ez_eo
        real (kind = RKIND), allocatable, dimension(:,:,:) :: Ez_oe
        real (kind = RKIND), allocatable, dimension(:,:,:) :: Ez_oo

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
                    idx = idx + 1
                    FieldList(idx)%data = Hx_m(i1, i2, i3, :, :, :)
                    FieldList(idx)%field_type = 'Hx'
                end do
            end do
        end do

        do i1 = 1, 3
            do i2 = 1, 2
                do i3 = 1, 3
                    idx = idx + 1
                    FieldList(idx)%data = Hy_m(i1, i2, i3, :, :, :)
                    FieldList(idx)%field_type = 'Hy'
                end do
            end do
        end do

        do i1 = 1, 3
            do i2 = 1, 3
                do i3 = 1, 2
                    idx = idx + 1
                    FieldList(idx)%data = Hz_m(i1, i2, i3, :, :, :)
                    
                end do
            end do
        end do

    end subroutine 

    ! subroutine GenerateOutputFields(b, FieldList)
        
    !     type (bounds_t), intent( IN)  ::  b
    !     type(field_array_t), allocatable, intent(OUT) :: FieldListOutput(:) ! Aquí necesito cambiar el tipo de variable de los outputs para tener en cuenta que con el input i, se generan varios outputs
    !     allocate(FieldListOutput(66))

    !     call GenerateInputFieldsBasis(b, FieldListInput)

    !     integer :: i

    !     do i = 1, size(FieldListInput)
    !     ! Acá es necesario realizar el paso temporal y extraer los campos usando el timestepping/resuming, en todo caso si la función es general, se llama fuera del case y se almacena dependiendo
    !     ! del caso.
    !         select case (trim(FieldListInput(i)%field_type))
    !             case ("Ex")
    !                 call Advance_Ex()
    !                 call Advance_Hy()
    !                 call Advance_Hz()
    !             case ("Ey")
    !                 call Advance_Ey()
    !                 call Advance_Hx()
    !                 call Advance_Hz()
    !             case ("Ez")
    !                 call Advance_Ez()
    !                 call Advance_Hx()
    !                 call Advance_Hy()
    !             case ("Hx")
    !                 call Advance_Ex()
    !                 call Advance_Ey()
    !                 Call Advance_Ez()
    !                 call Advance_Hx()
    !             case ("Hy")
    !                 call Advance_Hy()
    !             case ("Hz")
    !                 call Advance_Hz()
    !         end select
    !     end do

    ! end subroutine

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


    end subroutine

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

end module