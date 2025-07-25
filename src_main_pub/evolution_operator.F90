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

    public :: GenerateElectricalInputBasis,  GenerateMagneticalInputBasis, AddElectricFieldIndices, fhash_get_int_array, int_array

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

    subroutine AddElectricFieldIndices(RowIndexMap, field, shiftE, shiftM1, shiftM2, dirM1, dirM2)
        type(fhash_tbl_t), intent(inout) :: RowIndexMap
        type(limit_t), intent(in) :: field
        integer, intent(in) :: shiftE, shiftM1, shiftM2
        character(len=1), intent(in) :: dirM1, dirM2  

        integer :: i, j, k, m, m_shift1, m_shift2
        integer :: Nx, Ny, Nz

        type(int_array) :: wrapper 
        integer, allocatable :: indexList(:)
        integer :: countIndex, positionList

        Nx = field%Nx
        Ny = field%Ny
        Nz = field%Nz

        do i = 1, Nx
            do j = 1, Ny
                do k = 1, Nz
                    m = ((i - 1)*Ny + (j - 1))*Nz  + k
                    countIndex = 1
                    positionList = 1

                    select case (dirM1)
                        case ('i')
                            if (i > 1 .and. i < Nx) then
                                m_shift1 = ((i - 2)*Ny + (j - 1))*Nz  + k
                                countIndex = countIndex + 2
                            else if (i == 1) then
                                m_shift1 = -1
                                countIndex = countIndex + 1
                            else
                                m_shift1 = -2
                                countIndex = countIndex + 1
                            end if
                        case ('j')
                            if (j > 1 .and. j < Ny) then
                                m_shift1 = ((i - 1)*Ny + (j - 2))*Nz  + k
                                countIndex = countIndex + 2
                            else if (j == 1) then
                                m_shift1 = -1
                                countIndex = countIndex + 1
                            else
                                m_shift1 = -2
                                countIndex = countIndex + 1
                            end if
                        case ('k')
                            if (k > 1 .and. k < Nz) then
                                m_shift1 = ((i - 1)*Ny + (j - 1))*Nz  + (k - 1)
                                countIndex = countIndex + 2
                            else if (k == 1) then
                                m_shift1 = -1
                                countIndex = countIndex + 1
                            else
                                m_shift1 = -2
                                countIndex = countIndex + 1
                            end if
                    end select

                    select case (dirM2)
                        case ('i')
                            if (i > 1 .and. i < Nx) then
                                m_shift2 = ((i - 2)*Ny + (j - 1))*Nz  + k
                                countIndex = countIndex + 2
                            else if (i == 1) then
                                m_shift2 = -1
                                countIndex = countIndex + 1
                            else
                                m_shift2 = -2
                                countIndex = countIndex + 1
                            end if
                        case ('j')
                            if (j > 1 .and. j < Ny) then
                                m_shift2 = ((i - 1)*Ny + (j - 2))*Nz  + k
                                countIndex = countIndex + 2
                            else if (j == 1) then
                                m_shift2 = -1
                                countIndex = countIndex + 1
                            else
                                m_shift2 = -2
                                countIndex = countIndex + 1
                            end if
                        case ('k')
                            if (k > 1 .and. k < Nz) then
                                m_shift2 = ((i - 1)*Ny + (j - 1))*Nz  + (k - 1)
                                countIndex = countIndex + 2
                            else if (k == 1) then
                                m_shift2 = -1
                                countIndex = countIndex + 1
                            else
                                m_shift2 = -2
                                countIndex = countIndex + 1
                            end if
                    end select

                    ! Allocate the indexList with the size of countIndex
                    allocate(indexList(countIndex))

                    indexList(positionList) = shiftE + m
                    positionList = positionList + 1

                    if (m_shift2 /= -1 .and. m_shift2 /= -2) then
                        indexList(positionList) = shiftM1 + m
                        indexList(positionList + 1) = shiftM1 + m_shift2
                        positionList = positionList + 2
                    else if (m_shift2 == -1) then   ! Border at the beginning
                        indexList(positionList) = shiftM1 + m
                        positionList = positionList + 1
                    else                            ! Border at the end
                        indexList(positionList) = shiftM1 + m_shift2
                        positionList = positionList + 1
                    end if


                    if (m_shift1 /= -1 .and. m_shift1 /= -2) then
                        indexList(positionList) = shiftM2 + m
                        indexList(positionList + 1) = shiftM2 + m_shift1
                        positionList = positionList + 2
                    else if (m_shift1 == -1) then   ! Border at the beginning
                        indexList(positionList) = shiftM2 + m
                        positionList = positionList + 1
                    else                            ! Border at the end
                        indexList(positionList) = shiftM2 + m_shift1
                        positionList = positionList + 1
                    end if

                    wrapper%data = indexList
                    call RowIndexMap%set(key(shiftE + m), value=wrapper)  


                    deallocate(indexList)
                end do
            end do
        end do
    end subroutine 

    ! subroutine AddMagneticFieldIndices(RowIndexMap, field, shiftH, shiftE1, shiftE2, dir1, dir2)
    !     type(fhash_tbl_t), intent(inout) :: RowIndexMap
    !     type(limit_t), intent(in) :: field
    !     integer, intent(in) :: shiftH, shiftE1, shiftE2
    !     character(len=1), intent(in) :: dir1, dir2

    !     integer :: i, j, k, m, m_shift1, m_shift2
    !     integer :: Nx, Ny, Nz
    !     integer, allocatable :: temp(:), indexList(:)
    !     integer, allocatable :: aux1(:), aux2(:), aux3(:), aux4(:)
    !     integer :: totalSize

    !     Nx = field%Nx + 2
    !     Ny = field%Ny + 2
    !     Nz = field%Nz + 2

    !     do i = 1, Nx - 2
    !         do j = 1, Ny - 2
    !             do k = 1, Nz - 2
    !                 m = (i*Ny + j)*Nz + k

    !                 select case (dir1)
    !                     case ('i')
    !                         m_shift1 = ((i + 1)*Ny + j)*Nz + k
    !                     case ('j')
    !                         m_shift1 = (i*Ny + (j + 1))*Nz + k
    !                     case ('k')
    !                         m_shift1 = (i*Ny + j)*Nz + (k + 1)
    !                 end select

    !                 select case (dir2)
    !                     case ('i')
    !                         m_shift2 = ((i + 1)*Ny + j)*Nz + k
    !                     case ('j')
    !                         m_shift2 = (i*Ny + (j + 1))*Nz + k
    !                     case ('k')
    !                         m_shift2 = (i*Ny + j)*Nz + (k + 1)
    !                 end select

    !                 call RowIndexMap%get(key(shiftE1 + m), aux1)
    !                 call RowIndexMap%get(key(shiftE1 + m_shift1), aux2)
    !                 call RowIndexMap%get(key(shiftE2 + m), aux3)
    !                 call RowIndexMap%get(key(shiftE2 + m_shift2), aux4)

    !                 totalSize = size(aux1) + size(aux2) + size(aux3) + size(aux4)
    !                 allocate(temp(totalSize))
    !                 temp(1:size(aux1)) = aux1
    !                 temp(size(aux1)+1:size(aux1)+size(aux2)) = aux2
    !                 temp(size(aux1)+size(aux2)+1:size(aux1)+size(aux2)+size(aux3)) = aux3
    !                 temp(size(aux1)+size(aux2)+size(aux3)+1:) = aux4

    !                 call RemoveDuplicates(temp, indexList)


    !                 call RowIndexMap%set(key(shiftH + m), value=indexList)

    !                 deallocate(temp, indexList, aux1, aux2, aux3, aux4)
    !             end do
    !         end do
    !     end do
    ! end subroutine

    ! subroutine AddBoundaryIndices(RowIndexMap, sggBorder, field, shiftField, dir):
    !     type(fhash_tbl_t), intent(inout) :: RowIndexMap
    !     type(Border_t), intent(in) :: sggBorder
    !     type(limit_t), intent(in) :: field

    !     integer, intent(in) :: shiftField
    !     character(len=1), intent(in) :: dir

    !     !Hx Down
    !     if (sggBorder%IsDownPMC) then
    !         if (layoutnumber == 0)      Hx( : , : ,C(iHx)%ZI-1)=-Hx( : , : ,C(iHx)%ZI)
    !     endif
    !     !Hx Up
    !     if (sggBorder%IsUpPMC) then
    !         if (layoutnumber == size-1) Hx( : , : ,C(iHx)%ZE+1)=-Hx( : , : ,C(iHx)%ZE)
    !     endif
    !     !Hx Left
    !     if (sggBorder%IsLeftPMC) then
    !         Hx( : ,C(iHx)%YI-1, : )=-Hx( : ,C(iHx)%YI, : )
    !     endif
    !     !Hx Right
    !     if (sggBorder%IsRightPMC) then
    !         Hx( : ,C(iHx)%YE+1, : )=-Hx( : ,C(iHx)%YE, : )
    !     endif
    !     !Hy Back
    !     if (sggBorder%IsBackPMC) then
    !         Hy(C(iHy)%XI-1, : , : )=-Hy(C(iHy)%XI, : , : )
    !     endif
    !     !Hy Front
    !     if (sggBorder%IsFrontPMC) then
    !         Hy(C(iHy)%XE+1, : , : )=-Hy(C(iHy)%XE, : , : )
    !     endif
    !     !Hy Down
    !     if (sggBorder%IsDownPMC) then
    !         if (layoutnumber == 0)      Hy( : , : ,C(iHy)%ZI-1)=-Hy( : , : ,C(iHy)%ZI)
    !     endif
    !     !Hy Up
    !     if (sggBorder%IsUpPMC) then
    !         if (layoutnumber == size-1) Hy( : , : ,C(iHy)%ZE+1)=-Hy( : , : ,C(iHy)%ZE)
    !     endif
    !     !
    !     !Hz Back
    !     if (sggBorder%IsBackPMC) then
    !         Hz(C(iHz)%XI-1, : , : )=-Hz(C(iHz)%XI, : , : )
    !     endif
    !     !Hz Front
    !     if (sggBorder%IsFrontPMC) then
    !         Hz(C(iHz)%XE+1, : , : )=-Hz(C(iHz)%XE, : , : )
    !     endif
    !     !Hz Left
    !     if (sggBorder%IsLeftPMC) then
    !         Hz( : ,C(iHz)%YI-1, : )=-Hz( : ,C(iHz)%YI, : )
    !     endif
    !     !Hz Right
    !     if (sggBorder%IsRightPMC) then
    !         Hz( : ,C(iHz)%YE+1, : )=-Hz( : ,C(iHz)%YE, : )
    !     endif

    ! end subroutine
        
!     subroutine RemoveDuplicates(inputArray, outputArray)
!         integer, intent(in) :: inputArray(:)
!         integer, allocatable, intent(out) :: outputArray(:)
!         integer :: i, j, n
!         logical :: found
!         integer, allocatable :: temp(:)

!         allocate(temp(size(inputArray)))
!         n = 0

!         do i = 1, size(inputArray)
!             found = .false.
!             do j = 1, n
!                 if (temp(j) == inputArray(i)) then
!                     found = .true.
!                     exit
!                 end if
!             end do
!             if (.not. found) then
!                 n = n + 1
!                 temp(n) = inputArray(i)
!             end if
!         end do

!         allocate(outputArray(n))
!         outputArray = temp(1:n)
!         deallocate(temp)
!     end subroutine 

!     subroutine GenerateRowIndexMap(b, RowIndexMap)

!         type(bounds_t), intent(IN) :: b
!         type(fhash_tbl_t), intent(OUT) :: RowIndexMap
!         integer :: shiftEx, shiftEy, shiftEz, shiftHx, shiftHy, shiftHz

!         shiftEx = 0
!         shiftEy = b%Ex%Nx * b%Ex%Ny * b%Ex%Nz
!         shiftEz = shiftEy + b%Ey%Nx * b%Ey%Ny * b%Ey%Nz
!         shiftHx = shiftEz + b%Ez%Nx * b%Ez%Ny * b%Ez%Nz
!         shiftHy = shiftHx + b%Hx%Nx * b%Hx%Ny * b%Hx%Nz
!         shiftHz = shiftHy + b%Hy%Nx * b%Hy%Ny * b%Hy%Nz

!         call AddElectricFieldIndices(RowIndexMap, b%sweepEx, shiftEx, shiftHy, shiftHz, 'k', 'j')
!         call AddElectricFieldIndices(RowIndexMap, b%sweepEy, shiftEy, shiftHx, shiftHz, 'k', 'i')
!         call AddElectricFieldIndices(RowIndexMap, b%sweepEz, shiftEz, shiftHx, shiftHy, 'j', 'i')

!         ! Before the magnetic fields, it is necessary to create the map of indices related to the boundary conditions

!         call AddMagneticFieldIndices(RowIndexMap, b%sweepHx, shiftHx, shiftEy, shiftEz, 'k', 'j')
!         call AddMagneticFieldIndices(RowIndexMap, b%sweepHy, shiftHy, shiftEx, shiftEz, 'k', 'i')
!         call AddMagneticFieldIndices(RowIndexMap, b%sweepHz, shiftHz, shiftEx, shiftEy, 'j', 'i')

!         ! And also, it seems to be boundary conditions for the magnetic fields, so we need to add them as well


!     end subroutine

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