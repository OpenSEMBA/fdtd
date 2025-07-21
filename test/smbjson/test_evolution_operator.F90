integer function test_evolution_operator_dimension_Field_basis() bind (C, name="test_evolution_operator_dimension_Field_basis") result(err)
    use smbjson
    use smbjson_testingTools
    use evolution_operator

    implicit none

    real(RKIND), dimension(2,2,2) :: M_E
    real(RKIND), dimension(2,2,2) :: M_ee, M_eo, M_oe, M_oo
    real(RKIND), dimension(3,3,3) :: A
    real(RKIND), allocatable :: M_H(:,:,:,:,:,:)  
    integer, parameter :: M1 = 2, M2 = 3, M3 = 3

    err = 0

    M_E = 0.0_RKIND
    A = 0.0_RKIND
    call GenerateElectricalInputBasis(M_E, 1, 2, M_ee, M_eo, M_oe, M_oo)
    call GenerateMagneticalInputBasis(A, M1, M2, M3, M_H)

    if (size(M_ee, 1) /= 2 .or. size(M_ee, 2) /= 2 .or. size(M_ee, 3) /= 2) err = err + 1
    if (size(M_eo, 1) /= 2 .or. size(M_eo, 2) /= 2 .or. size(M_eo, 3) /= 2) err = err + 1
    if (size(M_oe, 1) /= 2 .or. size(M_oe, 2) /= 2 .or. size(M_oe, 3) /= 2) err = err + 1
    if (size(M_oo, 1) /= 2 .or. size(M_oo, 2) /= 2 .or. size(M_oo, 3) /= 2) err = err + 1

    if (size(M_H,1) /= M1 .or. size(M_H,2) /= M2 .or. size(M_H,3) /= M3) err = err + 1
    if (size(M_H,4) /= 3 .or. size(M_H,5) /= 3 .or. size(M_H,6) /= 3) err = err + 1

end function test_evolution_operator_dimension_Field_basis

integer function test_evolution_operator_poisition_E_basis() bind(C, name="test_evolution_operator_poisition_E_basis") result(err)
    use smbjson
    use smbjson_testingTools
    use evolution_operator

    implicit none

    integer, parameter :: D1 = 2, D2 = 2, D3 = 2
    real(RKIND), dimension(D1,D2,D3) :: M_E
    real(RKIND), dimension(D1,D2,D3) :: M_ee, M_eo, M_oe, M_oo, M_sum
    integer, parameter :: dim1 = 1, dim2 = 2
    integer :: i, j, k

    err = 0

    call GenerateElectricalInputBasis(M_E, dim1, dim2, M_ee, M_eo, M_oe, M_oo)

    do i = 1, D1
        do j = 1, D2
            do k = 1, D3
                if (abs(M_ee(i,j,k) + M_eo(i,j,k) + M_oe(i,j,k) + M_oo(i,j,k) - 1.0_RKIND) > 1.0e-12_RKIND) then
                    err = err + 1
                end if

                if (count([M_ee(i,j,k), M_eo(i,j,k), M_oe(i,j,k), M_oo(i,j,k)] == 1.0_RKIND) /= 1) then
                    err = err + 1
                end if
            end do
        end do
    end do

    M_sum = M_ee + M_eo + M_oe + M_oo
    if (any(abs(M_sum - 1.0_RKIND) > 1.0e-12_RKIND)) then
        err = err + 1
    end if

end function test_evolution_operator_poisition_E_basis

integer function test_evolution_operator_position_H_basis() bind(C, name="test_evolution_operator_position_H_basis") result(err)
    use smbjson
    use smbjson_testingTools
    use evolution_operator

    implicit none

    integer, parameter :: D1 = 2, D2 = 2, D3 = 2
    real(RKIND), dimension(D1,D2,D3) :: A
    real(RKIND), allocatable :: M_H(:,:,:,:,:,:)  
    integer, parameter :: M1 = 2, M2 = 3, M3 = 3
    integer :: i, j, k
    integer :: i_m1, i_m2, i_m3
    integer :: count_ones

    err = 0

    A = 0.0_RKIND
    call GenerateMagneticalInputBasis(A, M1, M2, M3, M_H)

    ! Loop over the dimensions of the H basis
    do i = 1, D1
        do j = 1, D2
            do k = 1, D3
                count_ones = 0

                ! Loop over all the elements in the H basis
                do i_m1 = 1, M1
                    do i_m2 = 1, M2
                        do i_m3 = 1, M3
                            if (abs(M_H(i_m1,i_m2,i_m3,i,j,k) - 1.0_RKIND) < 1.0e-12_RKIND) then
                                count_ones = count_ones + 1
                            else if (abs(M_H(i_m1,i_m2,i_m3,i,j,k)) > 1.0e-12_RKIND) then
                                err = err + 1
                            end if
                        end do
                    end do
                end do

                if (count_ones /= 1) then
                    err = err + 1
                end if

            end do
        end do
    end do

end function test_evolution_operator_position_H_basis

! One final test related to the basis of the inputs for the construction of the evolution operator could be check if the separation corresponds with the one detected in the map

! integer function test_evolution_operator_oneStep() bind (C) result(err)
!     use smbjson
!     use smbjson_testingTools
!     use evolution_operator

!     implicit none

!     type(evolution_operator) :: evolOp

!     err = 0

!     call evolOp%GenerateOperator()

!     ExternalField_t :: field

!     expected_field = smbjson%step(field)
!     result_field = evolOp%step(field, 1)

!     if (any(expected_field%field /= result_field%field)) then
!         err = err + 1
!     end if

! end function test_evolution_operator_oneStep
