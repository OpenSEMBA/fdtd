message(STATUS "Creating build system for border tests")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library(border_test_fortran
    "test_periodic.F90"
)

target_link_libraries(border_test_fortran
    semba-components
    semba-types
)

add_library(border_tests "border_tests.cpp")

target_link_libraries(border_tests
    border_test_fortran
    GTest::gtest
)
