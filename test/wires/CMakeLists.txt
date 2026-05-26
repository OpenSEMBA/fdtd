message(STATUS "Creating build system for test/wires")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library(
    wires_test_fortran
    "test_wires.F90"
)

target_link_libraries(wires_test_fortran
    semba-components
)

add_library(wires_tests "wires_tests.cpp")

target_link_libraries(wires_tests
    wires_test_fortran
    GTest::gtest
)

