integer function test_create_vtk_file() bind(C) result(err)
    use vtk

    integer :: test_err = 0
    err = test_err
end function test_create_vtk_file

integer function test_write_vtk_file() bind(c) result(err)
end function test_write_vtk_file

integer function test_vtk_file_for_probe() bind(c) result(err)
end function test_vtk_file_for_probe