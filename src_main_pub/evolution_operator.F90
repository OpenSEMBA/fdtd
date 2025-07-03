!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Module to handle the creation of the evolution operator
!  Date :  July, 3, 2025
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module evolution_operator

    use Resuming
    use Solver_mod
    use fdetypes
    use Report

    type :: field_array_t
        real(RKIND), pointer, dimension(:,:,:) :: data
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

        idx = 12
        do i1 = 1, 2
            do i2 = 1, 3
                do i3 = 1, 3
                    idx = idx + 1
                    FieldList(idx)%data => Hx_m(i1, i2, i3, :, :, :)
                end do
            end do
        end do

        do i1 = 1, 3
            do i2 = 1, 2
                do i3 = 1, 3
                    idx = idx + 1
                    FieldList(idx)%data => Hy_m(i1, i2, i3, :, :, :)
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