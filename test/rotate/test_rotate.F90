integer function test_simple_rotate() bind(C) result(err)
    use smbjson
    use rotate_testingTools
    use nfde_rotate_m
    character(len=*),parameter :: filename = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'asimetricPlanewave.fdtd.json'
    type(parser_t) :: parser
    type(Parseador) :: parsed

    parser = parser_t(filename)
    parsed = parser%readProblemDescription()
    if (parser%isInitialized) then
        err = 0
    else 
        err = 1
    end if



    call nfde_rotate(parsed, 1)
    call expect_eq_int(err, 20, parsed%matriz%totalZ)

end function
