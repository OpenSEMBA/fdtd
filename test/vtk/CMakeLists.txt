message(STATUS "Creating build system for test/vtk")

add_library(
    vtk_test_fortran
    "rotate_testingTools.F90"

)

target_link_libraries(vtk_test_fortran
    semba-outputs
)

add_library(vtk_tests "vtk_tests.cpp")

target_link_libraries(vtk_tests
    vtk_test_fortran
    GTest::gtest
)