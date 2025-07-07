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

    subroutine GenerateRowIndexMap(b, RowIndexMap)

        type(bounds_t), intent(IN) :: b
        type(fhash_tbl), intent(OUT) :: RowIndexMap
        integer :: shiftEx, shiftEy, shiftEz, shiftHx, shiftHy, shiftHz
        integer :: i, j, k, m
        integer, allocatable :: indexList(:)

        shiftEx = 0
        shiftEy = b%Ex%Nx * b%Ex%Ny * b%Ex%Nz
        shiftEz = shiftEy + b%Ey%Nx * b%Ey%Ny * b%Ey%Nz
        shiftHx = shiftEz + b%Ez%Nx * b%Ez%Ny * b%Ez%Nz
        shiftHy = shiftHx + b%Hx%Nx * b%Hx%Ny * b%Hx%Nz
        shiftHz = shiftHy + b%Hy%Nx * b%Hy%Ny * b%Hy%Nz

        do i = 1, b%Ex%Nx-2
            do j = 1, b%Ex%Ny-2
                do k = 1, b%Ex%Nz-2
                    m = (i*(b%Ex%Ny - 1) + j)*(b%Ex%Nz - 1) + k
                    m_shift_j = (i*(b%Ex%Ny - 1) + (j - 1))*(b%Ex%Nz - 1) + k
                    m_shift_k = (i*(b%Ex%Ny - 1) + j)*(b%Ex%Nz - 1) + (k - 1)

                    allocate(indexList(5))
                    indexList(1) = m
                    indexList(2) = shiftHy + m
                    indexList(3) = shiftHy + m_shift_k
                    indexList(4) = shiftHz + m
                    indexList(5) = shiftHz + m_shift_j

                    call fhash_insert(RowIndexMap, m, indexList)

                    deallocate(indexList)
                    
                end do
            end do
        end do

        do i = 1, b%Ey%Nx-2
            do j = 1, b%Ey%Ny-2
                do k = 1, b%Ey%Nz-2
                    m = (i*(b%Ey%Ny - 1) + j)*(b%Ey%Nz - 1) + k
                    m_shift_i = ((i - 1)*(b%Ey%Ny - 1) + j)*(b%Ey%Nz - 1) + k
                    m_shift_k = (i*(b%Ey%Ny - 1) + j)*(b%Ey%Nz - 1) + (k - 1)

                    allocate(indexList(5))
                    indexList(1) = shiftEy + m
                    indexList(2) = shiftHx + m
                    indexList(3) = shiftHx + m_shift_k
                    indexList(4) = shiftHz + m
                    indexList(5) = shiftHz + m_shift_i

                    call fhash_insert(RowIndexMap, shiftEy + m, indexList)

                    deallocate(indexList)
                    
                end do
            end do
        end do

        do i = 1, b%Ez%Nx-2
            do j = 1, b%Ez%Ny-2
                do k = 1, b%Ez%Nz-2
                    m = (i*(b%Ez%Ny - 1) + j)*(b%Ez%Nz - 1) + k
                    m_shift_i = ((i - 1)*(b%Ez%Ny - 1) + j)*(b%Ez%Nz - 1) + k
                    m_shift_j = (i*(b%Ez%Ny - 1) + (j - 1))*(b%Ez%Nz - 1) + k

                    allocate(indexList(5))
                    indexList(1) = shiftEz + m
                    indexList(2) = shiftHx + m
                    indexList(3) = shiftHx + m_shift_j
                    indexList(4) = shiftHy + m
                    indexList(5) = shiftHy + m_shift_i

                    call fhash_insert(RowIndexMap, shiftEz + m, indexList)

                    deallocate(indexList)
                    
                end do
            end do
        end do

    end subroutine
            
