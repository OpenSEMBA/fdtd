message(STATUS "Creating build system for test/output")

add_library(
    output_test_fortran
    "test_output.F90"
    "test_output_utils.F90"
    "test_volumic_utils.F90"
)

add_library(vtkAPI_test_fortran
    "test_vtkAPI.F90"
)

add_library(xdmfAPI_test_fortran
    "test_xdmfAPI.F90"
)

add_library(output_tests "output_tests.cpp")
add_library(vtkAPI_tests "vtkAPI_tests.cpp")
add_library(xdmfAPI_tests "xdmfAPI_tests.cpp")

target_link_libraries(output_test_fortran
    semba-outputs
    fdtd-output
    test_utils_fortran
)
target_link_libraries(output_tests
    output_test_fortran
    GTest::gtest
)
target_link_libraries(vtkAPI_test_fortran
    vtkAPI
    test_utils_fortran
)
target_link_libraries(vtkAPI_tests
    vtkAPI_test_fortran
    GTest::gtest
)

target_link_libraries(xdmfAPI_test_fortran
    xdmfAPI
    test_utils_fortran
    ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES}
)

target_link_libraries(xdmfAPI_tests
    xdmfAPI_test_fortran
    GTest::gtest
)

target_include_directories(xdmfAPI_test_fortran PRIVATE ${HDF5_INCLUDE_DIRS})

