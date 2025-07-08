module evolution_operator

    use Resuming
    use Solver_mod
    use fdetypes
    use Report

    use fhash, only: fhash_tbl, key => fhash_key

    type :: field_array_t
        real(RKIND), pointer, dimension(:,:,:) :: data
        character(len=2) :: field_type  ! 'Ex', 'Ey', 'Ez', 'Hx', etc.
    end type

    implicit none
    private 

    public :: evolution_operator

contains

    subroutine GenerateElectricalInputBasis(M, dim1, dim2, M_ee, M_eo, M_oe, M_oo)

        real(RKIND), intent(in)  :: M(:,:,:)
        integer,     intent(in)  :: dim1, dim2  

        real(RKIND), intent(out) :: M_ee(size(M,1), size(M,2), size(M,3))
        real(RKIND), intent(out) :: M_eo(size(M,1), size(M,2), size(M,3))
        real(RKIND), intent(out) :: M_oe(size(M,1), size(M,2), size(M,3))
        real(RKIND), intent(out) :: M_oo(size(M,1), size(M,2), size(M,3))

        integer :: i, j, k
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
                v1 = [i, j, k](dim1)
                v2 = [i, j, k](dim2)

                select case (2*mod(v1,2) + mod(v2,2))
                case (0); M_ee(i,j,k) = 1.0_RKIND
                case (1); M_eo(i,j,k) = 1.0_RKIND
                case (2); M_oe(i,j,k) = 1.0_RKIND
                case (3); M_oo(i,j,k) = 1.0_RKIND
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
                M(idx1+1, idx2+1, idx3+1, i, j, k) = 1.0_RKIND
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
        allocate(FieldList(66))

        ! Generating the basis for the electical fields
        real (kind = RKIND), dimension    ( 0 :      b%Ex%NX-1 , 0 :      b%Ex%NY-1 , 0 :      b%Ex%NZ-1 ) ::  Ex
        real (kind = RKIND), dimension    ( 0 :      b%Ey%NX-1 , 0 :      b%Ey%NY-1 , 0 :      b%Ey%NZ-1 ) ::  Ey
        real (kind = RKIND), dimension    ( 0 :      b%Ez%NX-1 , 0 :      b%Ez%NY-1 , 0 :      b%Ez%NZ-1 ) ::  Ez

        Ex = 0.0_RKIND
        Ey = 0.0_RKIND
        Ez = 0.0_RKIND

        call GenerateElectricalInputBasis(Ex, 2, 3, Ex_ee, Ex_eo, Ex_oe, Ex_oo)
        call GenerateElectricalInputBasis(Ey, 1, 3, Ey_ee, Ey_eo, Ey_oe, Ey_oo)
        call GenerateElectricalInputBasis(Ez, 1, 2, Ez_ee, Ez_eo, Ez_oe, Ez_oo)

        ! Storing the electrical fields in the FieldList
        FieldList(1)%data => Ex_ee
        FieldList(2)%data => Ex_eo
        FieldList(3)%data => Ex_oe
        FieldList(4)%data => Ex_oo

        FieldList(5)%data => Ey_ee
        FieldList(6)%data => Ey_eo
        FieldList(7)%data => Ey_oe
        FieldList(8)%data => Ey_oo

        FieldList(9)%data  => Ez_ee
        FieldList(10)%data => Ez_eo
        FieldList(11)%data => Ez_oe
        FieldList(12)%data => Ez_oo

        ! Generating the basis for the magnetical fields
        real (kind = RKIND), dimension    ( 0 :      b%HX%NX-1 , 0 :      b%HX%NY-1 , 0 :      b%HX%NZ-1 ) ::  Hx
        real (kind = RKIND), dimension    ( 0 :      b%Hy%NX-1 , 0 :      b%Hy%NY-1 , 0 :      b%Hy%NZ-1 ) ::  Hy
        real (kind = RKIND), dimension    ( 0 :      b%Hz%NX-1 , 0 :      b%Hz%NY-1 , 0 :      b%Hz%NZ-1 ) ::  Hz

        Hx = 0.0_RKIND
        Hy = 0.0_RKIND
        Hz = 0.0_RKIND

        call GenerateMagneticalInputBasis(Hx, 2, 3, 3, Hx_m)
        call GenerateMagneticalInputBasis(Hy, 3, 2, 3, Hy_m)
        call GenerateMagneticalInputBasis(Hz, 3, 3, 2, Hz_m)

        ! Storing the magnetical fields in the FieldList
        integer :: idx, i1, i2, i3

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
                    FieldList(idx)%data => Hx_m(i1, i2, i3, :, :, :)
                    FieldList(idx)%field_type = 'Hx'
                end do
            end do
        end do

        do i1 = 1, 3
            do i2 = 1, 2
                do i3 = 1, 3
                    idx = idx + 1
                    FieldList(idx)%data => Hy_m(i1, i2, i3, :, :, :)
                    FieldList(idx)%field_type = 'Hy'
                end do
            end do
        end do

        do i1 = 1, 3
            do i2 = 1, 3
                do i3 = 1, 2
                    idx = idx + 1
                    FieldList(idx)%data => Hz_m(i1, i2, i3, :, :, :)
                    
                end do
            end do
        end do

    end subroutine 

    subroutine GenerateOutputFields(b, FieldList)
        
        type (bounds_t), intent( IN)  ::  b
        type(field_array_t), allocatable, intent(OUT) :: FieldListOutput(:) ! Aquí necesito cambiar el tipo de variable de los outputs para tener en cuenta que con el input i, se generan varios outputs
        allocate(FieldListOutput(66))

        call GenerateInputFieldsBasis(b, FieldListInput)

        integer :: i

        do i = 1, size(FieldListInput)
        ! Acá es necesario realizar el paso temporal y extraer los campos usando el timestepping/resuming, en todo caso si la función es general, se llama fuera del case y se almacena dependiendo
        ! del caso.
            select case (trim(FieldListInput(i)%field_type))
                case ("Ex")
                    call Advance_Ex()
                    call Advance_Hy()
                    call Advance_Hz()
                case ("Ey")
                    call Advance_Ey()
                    call Advance_Hx()
                    call Advance_Hz()
                case ("Ez")
                    call Advance_Ez()
                    call Advance_Hx()
                    call Advance_Hy()
                case ("Hx")
                    call Advance_Ex()
                    call Advance_Ey()
                    Call Advance_Ez()
                    call Advance_Hx()
                case ("Hy")
                    call Advance_Hy()
                case ("Hz")
                    call Advance_Hz()
            end select
        end do

    end subroutine

    subroutine AddElectricFieldIndices(RowIndexMap, field, shiftE, shiftM1, shiftM2, dirM1, dirM2)
        type(fhash_tbl), intent(inout) :: RowIndexMap
        type(limit_t), intent(in) :: field
        integer, intent(in) :: shiftE, shiftM1, shiftM2
        character(len=1), intent(in) :: dirM1, dirM2  

        integer :: i, j, k, m, m_shift1, m_shift2
        integer, allocatable :: indexList(:)
        integer :: Nx, Ny, Nz

        Nx = field%Nx
        Ny = field%Ny
        Nz = field%Nz

        do i = 1, Nx - 2
            do j = 1, Ny - 2
                do k = 1, Nz - 2
                    m = (i * (Ny - 1) + j) * (Nz - 1) + k

                    select case (dirM1)
                        case ('i')
                            m_shift1 = ((i - 1)*(Ny - 1) + j)*(Nz - 1) + k
                        case ('j')
                            m_shift1 = (i*(Ny - 1) + (j - 1))*(Nz - 1) + k
                        case ('k')
                            m_shift1 = (i*(Ny - 1) + j)*(Nz - 1) + (k - 1)
                    end select

                    select case (dirM2)
                        case ('i')
                            m_shift2 = ((i - 1)*(Ny - 1) + j)*(Nz - 1) + k
                        case ('j')
                            m_shift2 = (i*(Ny - 1) + (j - 1))*(Nz - 1) + k
                        case ('k')
                            m_shift2 = (i*(Ny - 1) + j)*(Nz - 1) + (k - 1)
                    end select

                    allocate(indexList(5))
                    indexList(1) = shiftE + m
                    indexList(2) = shiftM1 + m
                    indexList(3) = shiftM1 + m_shift2
                    indexList(4) = shiftM2 + m
                    indexList(5) = shiftM2 + m_shift1

                    call RowIndexMap%set(key(shiftE + m), value=indexList)

                    deallocate(indexList)
                end do
            end do
        end do
    end subroutine 

    subroutine AddMagneticFieldIndices(RowIndexMap, field, shiftH, shiftE1, shiftE2, dir1, dir2)
        type(fhash_tbl), intent(inout) :: RowIndexMap
        type(limit_t), intent(in) :: field
        integer, intent(in) :: shiftH, shiftE1, shiftE2
        character(len=1), intent(in) :: dir1, dir2

        integer :: i, j, k, m, m_shift1, m_shift2
        integer :: Nx, Ny, Nz
        integer, allocatable :: temp(:), indexList(:)
        integer, allocatable :: aux1(:), aux2(:), aux3(:), aux4(:)

        Nx = field%Nx
        Ny = field%Ny
        Nz = field%Nz

        do i = 0, Nx - 1
            do j = 0, Ny - 1
                do k = 0, Nz - 1
                    m = (i*(Ny - 1) + j)*(Nz - 1) + k

                    select case (dir1)
                        case ('i')
                            m_shift1 = ((i + 1)*(Ny - 1) + j)*(Nz - 1) + k
                        case ('j')
                            m_shift1 = (i*(Ny - 1) + (j + 1))*(Nz - 1) + k
                        case ('k')
                            m_shift1 = (i*(Ny - 1) + j)*(Nz - 1) + (k + 1)
                    end select

                    select case (dir2)
                        case ('i')
                            m_shift2 = ((i + 1)*(Ny - 1) + j)*(Nz - 1) + k
                        case ('j')
                            m_shift2 = (i*(Ny - 1) + (j + 1))*(Nz - 1) + k
                        case ('k')
                            m_shift2 = (i*(Ny - 1) + j)*(Nz - 1) + (k + 1)
                    end select

                    call RowIndexMap%get(key(shiftE1 + m), aux1)
                    call RowIndexMap%get(key(shiftE1 + m_shift1), aux2)
                    call RowIndexMap%get(key(shiftE2 + m), aux3)
                    call RowIndexMap%get(key(shiftE2 + m_shift2), aux4)

                    integer :: totalSize
                    totalSize = size(aux1) + size(aux2) + size(aux3) + size(aux4)
                    allocate(temp(totalSize))
                    temp(1:size(aux1)) = aux1
                    temp(size(aux1)+1:size(aux1)+size(aux2)) = aux2
                    temp(size(aux1)+size(aux2)+1:size(aux1)+size(aux2)+size(aux3)) = aux3
                    temp(size(aux1)+size(aux2)+size(aux3)+1:) = aux4

                    call RemoveDuplicates(temp, indexList)


                    call RowIndexMap%set(key(shiftH + m), value=indexList)

                    deallocate(temp, indexList, aux1, aux2, aux3, aux4)
                end do
            end do
        end do
    end subroutine
        
    subroutine RemoveDuplicates(inputArray, outputArray)
        integer, intent(in) :: inputArray(:)
        integer, allocatable, intent(out) :: outputArray(:)
        integer :: i, j, n
        logical :: found
        integer, allocatable :: temp(:)

        allocate(temp(size(inputArray)))
        n = 0

        do i = 1, size(inputArray)
            found = .false.
            do j = 1, n
                if (temp(j) == inputArray(i)) then
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                n = n + 1
                temp(n) = inputArray(i)
            end if
        end do

        allocate(outputArray(n))
        outputArray = temp(1:n)
        deallocate(temp)
    end subroutine 

    subroutine GenerateRowIndexMap(b, RowIndexMap)

        type(bounds_t), intent(IN) :: b
        type(fhash_tbl), intent(OUT) :: RowIndexMap
        integer :: shiftEx, shiftEy, shiftEz, shiftHx, shiftHy, shiftHz

        shiftEx = 0
        shiftEy = b%Ex%Nx * b%Ex%Ny * b%Ex%Nz
        shiftEz = shiftEy + b%Ey%Nx * b%Ey%Ny * b%Ey%Nz
        shiftHx = shiftEz + b%Ez%Nx * b%Ez%Ny * b%Ez%Nz
        shiftHy = shiftHx + b%Hx%Nx * b%Hx%Ny * b%Hx%Nz
        shiftHz = shiftHy + b%Hy%Nx * b%Hy%Ny * b%Hy%Nz

        call AddElectricFieldIndices(RowIndexMap, b%Ex, shiftEx, shiftHy, shiftHz, 'k', 'j')
        call AddElectricFieldIndices(RowIndexMap, b%Ey, shiftEy, shiftHx, shiftHz, 'k', 'i')
        call AddElectricFieldIndices(RowIndexMap, b%Ez, shiftEz, shiftHx, shiftHy, 'j', 'i')

        ! Before the magnetic fields, it is necessary to create the map of indices related to the boundary conditions

        call AddMagneticFieldIndices(RowIndexMap, b%Hx, shiftHx, shiftEy, shiftEz, 'k', 'j')
        call AddMagneticFieldIndices(RowIndexMap, b%Hy, shiftHy, shiftEx, shiftEz, 'k', 'i')
        call AddMagneticFieldIndices(RowIndexMap, b%Hz, shiftHz, shiftEx, shiftEy, 'j', 'i')


    end subroutine