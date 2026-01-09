message(STATUS "Creating build system for test/output")

add_library(
    output_test_fortran
    "test_output.F90"
    "test_output_utils.F90"
)

target_link_libraries(output_test_fortran
    semba-outputs
    fdtd-output
    test_utils_fortran
)

add_library(output_tests "output_tests.cpp")

target_link_libraries(output_tests
    output_test_fortran
    GTest::gtest
)