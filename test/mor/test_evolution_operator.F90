integer function test_evolution_operator_dimension_Field_basis() bind (C, name="test_evolution_operator_dimension_Field_basis") result(err)
    use smbjson
    use smbjson_testingTools
    use evolution_operator 
    use fhash, only: fhash_tbl_t, key => fhash_key

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

integer function test_evolution_operator_E_indices_map() bind(C, name="test_evolution_operator_E_indices_map") result(err)
    use smbjson
    use smbjson_testingTools
    use evolution_operator
    use fhash, key => fhash_key

    implicit none

    integer :: i, j, k
    integer :: m
    type(bounds_t) :: bounds
    type(fhash_tbl_t) :: RowIndexMap
    type(int_array) :: wrapper
    integer :: ElementsInMap

    bounds%Ex%NX = 2
    bounds%Ex%NY = 3
    bounds%Ex%NZ = 4

    bounds%Hy%NX = 1
    bounds%Hy%NY = 4
    bounds%Hy%NZ = 3

    bounds%Hz%NX = 1
    bounds%Hz%NY = 3
    bounds%Hz%NZ = 4

    err = 0
    ElementsInMap = 0

    call AddElectricFieldIndices(RowIndexMap, bounds%Ex, bounds%Hy, bounds%Hz, 0, 0, 0, 'k', 'j')

    do i = 1, bounds%Ex%NX
        do j = 1, bounds%Ex%NY
            do k = 1, bounds%Ex%NZ
                m = ((i - 1)*bounds%Ex%NY + (j - 1))*bounds%Ex%NZ  + k

                call fhash_get_int_array(RowIndexMap, key(m), wrapper)

                ! Check if the map has been created correctly for each i, j, k
                if (size(wrapper%data) == 0) then
                    err = err + 1
                    cycle
                else 
                    ElementsInMap = ElementsInMap + 1
                end if

                ! First we check the number of neighbours in the frontier of the first direction
                if (j == 1 .or. j == bounds%Ex%Ny) then
                    if (k == 1 .or. k == bounds%Ex%NZ) then
                        if (size(wrapper%data) /= 3) then
                            err = err + 1
                        end if
                    else
                        if (size(wrapper%data) /= 4) then
                            err = err + 1
                        end if
                    end if
                ! Then we check the number of neighbours in the frontier of the second direction that are not neighbours of the first direction               
                else if (k == 1 .or. k == bounds%Ex%NZ) then
                    if (size(wrapper%data) /= 4) then
                        err = err + 1
                    end if
                ! Finally we check the number of neighbours in the interior of the grid
                else
                    if (size(wrapper%data) /= 5) then
                        err = err + 1
                    end if
                end if

            end do
        end do
    end do

    if (ElementsInMap /= bounds%Ex%NX * bounds%Ex%NY * bounds%Ex%NZ) then
        err = err + 1
    end if

end function test_evolution_operator_E_indices_map

integer function test_evolution_operator_H_indices_map() bind(C, name="test_evolution_operator_H_indices_map") result(err)
    use smbjson
    use smbjson_testingTools
    use evolution_operator
    use fhash, key => fhash_key

    implicit none

    integer :: i, j, k
    integer :: m
    type(bounds_t) :: bounds
    type(fhash_tbl_t) :: RowIndexMap
    type(int_array) :: wrapper
    integer :: ElementsInMap
    integer :: shiftEx, shiftEy, shiftEz, shiftHx, shiftHy, shiftHz

    bounds%Ex%NX = 1
    bounds%Ex%NY = 4
    bounds%Ex%NZ = 4

    bounds%Ey%NX = 2
    bounds%Ey%NY = 3
    bounds%Ey%NZ = 4

    bounds%Ez%NX = 2
    bounds%Ez%NY = 4
    bounds%Ez%NZ = 3

    bounds%Hx%NX = 2
    bounds%Hx%NY = 3
    bounds%Hx%NZ = 3

    bounds%Hy%NX = 1
    bounds%Hy%NY = 4
    bounds%Hy%NZ = 3

    bounds%Hz%NX = 1
    bounds%Hz%NY = 3
    bounds%Hz%NZ = 4

    shiftEx = 0
    shiftEy = 0       + bounds%Ex%Nx * bounds%Ex%Ny * bounds%Ex%Nz
    shiftEz = shiftEy + bounds%Ey%Nx * bounds%Ey%Ny * bounds%Ey%Nz
    shiftHx = shiftEz + bounds%Ez%Nx * bounds%Ez%Ny * bounds%Ez%Nz
    shiftHy = shiftHx + bounds%Hx%Nx * bounds%Hx%Ny * bounds%Hx%Nz
    shiftHz = shiftHy + bounds%Hy%Nx * bounds%Hy%Ny * bounds%Hy%Nz

    err = 0
    ElementsInMap = 0

    ! To verify the H indices map, first I need to create the map of all the Electic fields

    call AddElectricFieldIndices(RowIndexMap, bounds%Ey, bounds%Hx, bounds%Hz, shiftEy, shiftHx, shiftHz, 'k', 'i')
    call AddElectricFieldIndices(RowIndexMap, bounds%Ez, bounds%Hx, bounds%Hy, shiftEz, shiftHx, shiftHy, 'j', 'i')
    call AddMagneticFieldIndices(RowIndexMap, bounds%Hx, bounds%Ez, bounds%Ey, shiftHx, shiftEz, shiftEy, 'j', 'k')

    do i = 1, bounds%Hx%NX
        do j = 1, bounds%Hx%NY
            do k = 1, bounds%Hx%NZ
                m = ((i - 1)*bounds%Hx%NY + (j - 1))*bounds%Hx%NZ  + k

                call fhash_get_int_array(RowIndexMap, key(shiftHx + m), wrapper)

                ! Check if the map has been created correctly for each i, j, k
                if (size(wrapper%data) == 0) then
                    err = err + 1
                    cycle
                else 
                    ElementsInMap = ElementsInMap + 1
                end if

                ! In theory, the indices related to the H field must be maximum 17 and minimum 9, due to the form of the evolution operator, 
                ! it is difficult to check the exact number of neighbours, but we can check the limits
                if (size(wrapper%data) < 9 .or. size(wrapper%data) > 17) then
                    err = err + 1
                end if

                ! In this particular test, the frontiers must have those specific number of neighbours due to the geometry
                ! First we check the number of neighbours in the frontier of the first direction
                if (j == 1 .or. j == bounds%Hx%Ny) then
                    if (k == 1 .or. k == bounds%Hx%NZ) then
                        if (size(wrapper%data) /= 11) then
                            err = err + 1
                        end if
                    else
                        if (size(wrapper%data) /= 12) then
                            err = err + 1
                        end if
                    end if
                ! Then we check the number of neighbours in the frontier of the second direction that are not neighbours of the first direction               
                else if (k == 1 .or. k == bounds%Hx%NZ) then
                    if (size(wrapper%data) /= 12) then
                        err = err + 1
                    end if
                ! Finally we check the number of neighbours in the interior of the grid
                else
                    if (size(wrapper%data) /= 13) then
                        err = err + 1
                    end if
                end if

            end do
        end do
    end do

    if (ElementsInMap /= bounds%Hx%NX * bounds%Hx%NY * bounds%Hx%NZ) then
        err = err + 1
    end if
end function test_evolution_operator_H_indices_map


integer function test_evolution_operator_indices_map_all_fields() bind(C, name="test_evolution_operator_indices_map_all_fields") result(err)
    use smbjson
    use smbjson_testingTools
    use evolution_operator
    use fhash, key => fhash_key

    implicit none

    integer :: m
    type(bounds_t) :: bounds
    type(fhash_tbl_t) :: RowIndexMap
    type(int_array) :: wrapper
    integer :: totalElements

    bounds%Ex%NX = 1
    bounds%Ex%NY = 4
    bounds%Ex%NZ = 4

    bounds%Ey%NX = 2
    bounds%Ey%NY = 3
    bounds%Ey%NZ = 4

    bounds%Ez%NX = 2
    bounds%Ez%NY = 4
    bounds%Ez%NZ = 3

    bounds%Hx%NX = 2
    bounds%Hx%NY = 3
    bounds%Hx%NZ = 3

    bounds%Hy%NX = 1
    bounds%Hy%NY = 4
    bounds%Hy%NZ = 3

    bounds%Hz%NX = 1
    bounds%Hz%NY = 3
    bounds%Hz%NZ = 4

    err = 0

    totalElements = bounds%Ex%NX * bounds%Ex%NY * bounds%Ex%NZ + &
                    bounds%Ey%NX * bounds%Ey%NY * bounds%Ey%NZ + &
                    bounds%Ez%NX * bounds%Ez%NY * bounds%Ez%NZ + &
                    bounds%Hx%NX * bounds%Hx%NY * bounds%Hx%NZ + &
                    bounds%Hy%NX * bounds%Hy%NY * bounds%Hy%NZ + &
                    bounds%Hz%NX * bounds%Hz%NY * bounds%Hz%NZ

    call GenerateRowIndexMap(bounds, RowIndexMap)

    do m = 1, totalElements
        call fhash_get_int_array(RowIndexMap, key(m), wrapper)

        if (size(wrapper%data) == 0) then
            err = err + 1
        end if
    end do

    end function test_evolution_operator_indices_map_all_fields

    integer function test_evolution_operator_column_map_creation() bind(C, name="test_evolution_operator_column_map_creation") result(err)
    use smbjson
    use smbjson_testingTools
    use evolution_operator
    use fhash, key => fhash_key

    implicit none

    integer :: m
    type(bounds_t) :: bounds
    type(fhash_tbl_t) :: ColIndexMap
    type(int_array) :: wrapper
    integer :: totalElements

    bounds%Ex%NX = 1
    bounds%Ex%NY = 4
    bounds%Ex%NZ = 4

    bounds%Ey%NX = 2
    bounds%Ey%NY = 3
    bounds%Ey%NZ = 4

    bounds%Ez%NX = 2
    bounds%Ez%NY = 4
    bounds%Ez%NZ = 3

    bounds%Hx%NX = 2
    bounds%Hx%NY = 3
    bounds%Hx%NZ = 3

    bounds%Hy%NX = 1
    bounds%Hy%NY = 4
    bounds%Hy%NZ = 3

    bounds%Hz%NX = 1
    bounds%Hz%NY = 3
    bounds%Hz%NZ = 4

    err = 0

    totalElements = bounds%Ex%NX * bounds%Ex%NY * bounds%Ex%NZ + &
                    bounds%Ey%NX * bounds%Ey%NY * bounds%Ey%NZ + &
                    bounds%Ez%NX * bounds%Ez%NY * bounds%Ez%NZ + &
                    bounds%Hx%NX * bounds%Hx%NY * bounds%Hx%NZ + &
                    bounds%Hy%NX * bounds%Hy%NY * bounds%Hy%NZ + &
                    bounds%Hz%NX * bounds%Hz%NY * bounds%Hz%NZ

    call GenerateColumnIndexMap(bounds, ColIndexMap)

    do m = 1, totalElements
        call fhash_get_int_array(ColIndexMap, key(m), wrapper)

        if (size(wrapper%data) == 0) then
            err = err + 1
        end if
    end do

    end function test_evolution_operator_column_map_creation

    integer function test_evolution_operator_chaeck_map_consistency() bind(C, name="test_evolution_operator_chaeck_map_consistency") result(err)
    use smbjson
    use smbjson_testingTools
    use evolution_operator
    use fhash, key => fhash_key

    implicit none

    integer :: m, i
    type(bounds_t) :: bounds
    type(fhash_tbl_t) :: ColIndexMap, RowIndexMap
    type(int_array) :: wrapperCol, wrapperRow
    integer :: totalElements

    bounds%Ex%NX = 1
    bounds%Ex%NY = 4
    bounds%Ex%NZ = 4

    bounds%Ey%NX = 2
    bounds%Ey%NY = 3
    bounds%Ey%NZ = 4

    bounds%Ez%NX = 2
    bounds%Ez%NY = 4
    bounds%Ez%NZ = 3

    bounds%Hx%NX = 2
    bounds%Hx%NY = 3
    bounds%Hx%NZ = 3

    bounds%Hy%NX = 1
    bounds%Hy%NY = 4
    bounds%Hy%NZ = 3

    bounds%Hz%NX = 1
    bounds%Hz%NY = 3
    bounds%Hz%NZ = 4

    err = 0

    totalElements = bounds%Ex%NX * bounds%Ex%NY * bounds%Ex%NZ + &
                    bounds%Ey%NX * bounds%Ey%NY * bounds%Ey%NZ + &
                    bounds%Ez%NX * bounds%Ez%NY * bounds%Ez%NZ + &
                    bounds%Hx%NX * bounds%Hx%NY * bounds%Hx%NZ + &
                    bounds%Hy%NX * bounds%Hy%NY * bounds%Hy%NZ + &
                    bounds%Hz%NX * bounds%Hz%NY * bounds%Hz%NZ

    call GenerateRowIndexMap(bounds, RowIndexMap)
    call GenerateColumnIndexMap(bounds, ColIndexMap)

    do m = 1, totalElements
        call fhash_get_int_array(ColIndexMap, key(m), wrapperCol)

        do i = 1, size(wrapperCol%data)
            call fhash_get_int_array(RowIndexMap, key(wrapperCol%data(i)), wrapperRow)

            if (all(wrapperRow%data /= m)) then
                err = err + 1
            end if
        end do
    end do

    do m = 1, totalElements
        call fhash_get_int_array(RowIndexMap, key(m), wrapperRow)

        do i = 1, size(wrapperRow%data)
            call fhash_get_int_array(ColIndexMap, key(wrapperRow%data(i)), wrapperCol)

            if (all(wrapperCol%data /= m)) then
                err = err + 1
            end if
        end do
    end do

    end function test_evolution_operator_chaeck_map_consistency

    integer function test_evolution_operator_read_bounds_from_json() bind(C, name="test_evolution_operator_read_bounds_from_json") result(err)
    use smbjson
    use smbjson_testingTools
    use evolution_operator
    use SEMBA_FDTD_mod
    

    implicit none

    type(bounds_t) :: bounds
    type(semba_fdtd_t) :: semba
    character(len=*),parameter :: filename = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'grid_50x3x3.fdtd.json'

    err = 0

    call semba%init(trim('-i '//filename))
    call get_field_bounds_from_json(bounds, semba%fullsize)

    if (bounds%Ex%Nx /= 50 .or. bounds%Ex%Ny /= 4 .or. bounds%Ex%Nz /= 4) err = err + 1
    if (bounds%Ey%Nx /= 51 .or. bounds%Ey%Ny /= 3 .or. bounds%Ey%Nz /= 4) err = err + 1
    if (bounds%Ez%Nx /= 51 .or. bounds%Ez%Ny /= 4 .or. bounds%Ez%Nz /= 3) err = err + 1
    if (bounds%Hx%Nx /= 51 .or. bounds%Hx%Ny /= 3 .or. bounds%Hx%Nz /= 3) err = err + 1
    if (bounds%Hy%Nx /= 50 .or. bounds%Hy%Ny /= 4 .or. bounds%Hy%Nz /= 3) err = err + 1
    if (bounds%Hz%Nx /= 50 .or. bounds%Hz%Ny /= 3 .or. bounds%Hz%Nz /= 4) err = err + 1


    end function test_evolution_operator_read_bounds_from_json

    integer function test_evolution_operator_get_field_outputs() bind(C, name="test_evolution_operator_get_field_outputs") result(err)
    use smbjson
    use smbjson_testingTools
    use evolution_operator
    use SEMBA_FDTD_mod

    implicit none

    type(field_array_t) :: fieldInput, fieldOutput
    real(RKIND), dimension(3,4,4) :: M_Ex, M_ee, M_eo, M_oe, M_oo

    character(len=*),parameter :: filename = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'grid_3x3x3.fdtd.json'

    call GenerateElectricalInputBasis(M_Ex, 2, 3, M_ee, M_eo, M_oe, M_oo)
    fieldInput%data = M_ee
    fieldInput%field_type = 'Ex'
    
    call GenerateOutputFields('-i', filename, fieldInput, fieldOutput)

    end function test_evolution_operator_get_field_outputs

    integer function test_evolution_operator_comparison_with_solver() bind(C, name="test_evolution_operator_comparison_with_solver") result(err)
    use smbjson
    use smbjson_testingTools
    use evolution_operator
    use fhash, key => fhash_key

    implicit none

    integer :: m
    type(bounds_t) :: bounds
    type(fhash_tbl_t) :: RowIndexMap

    bounds%Ex%NX = 50
    bounds%Ex%NY = 4
    bounds%Ex%NZ = 4

    bounds%Ey%NX = 51
    bounds%Ey%NY = 3
    bounds%Ey%NZ = 4

    bounds%Ez%NX = 51
    bounds%Ez%NY = 4
    bounds%Ez%NZ = 3

    bounds%Hx%NX = bounds%Ex%Nx + 1
    bounds%Hx%NY = bounds%Ex%Ny - 1
    bounds%Hx%NZ = bounds%Ex%Nz - 1

    bounds%Hy%NX = bounds%Ey%Nx - 1
    bounds%Hy%NY = bounds%Ey%Ny + 1
    bounds%Hy%NZ = bounds%Ey%Nz - 1

    bounds%Hz%NX = bounds%Ez%Nx - 1
    bounds%Hz%NY = bounds%Ez%Ny - 1
    bounds%Hz%NZ = bounds%Ez%Nz + 1

    err = 0

    ! Generate the evolution operator with the basis, the map and one step with the solver
    ! I can make a function with these three steps inside the evolution operator module
    call GenerateRowIndexMap(bounds, RowIndexMap)

    ! With the evolution operator and the initial excitation, I can generate for example five steps and then compare with the full solver

    end function test_evolution_operator_comparison_with_solver