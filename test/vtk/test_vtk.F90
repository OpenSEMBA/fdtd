module test_vtk_m
    use VTK
    implicit none
    contains
        integer function test_init_vtk() bind(C) result(err)
        end function test_init_vtk
end module