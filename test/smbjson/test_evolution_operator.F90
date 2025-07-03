integer function test_evolution_operator_numberOfField_basis() bind (C) result(err)
    use smbjson
    use smbjson_testingTools
    use evolution_operator

    implicit none

    type(evolution_operator) :: evolOp

    err = 0

    call evolOp%GenerateInputFieldsBasis()
    if (evolOp%numberOfField_basis /= 66) err = err + 1


integer function test_evolution_operator_oneStep() bind (C) result(err)
    use smbjson
    use smbjson_testingTools
    use evolution_operator

    implicit none

    type(evolution_operator) :: evolOp

    err = 0

    call evolOp%GenerateOperator()

    ExternalField_t :: field

    expected_field = smbjson%step(field)
    result_field = evolOp%step(field, 1)

    if (any(expected_field%field /= result_field%field)) then
        err = err + 1
    end if

end function test_evolution_operator_oneStep
