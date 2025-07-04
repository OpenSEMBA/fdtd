
message(STATUS "Creating build system for test/system")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library(
    system_test_fortran
    "system_testingTools.F90"
    "test_init_solver.F90"
)

target_link_libraries(system_test_fortran
    fhash
    mtlnsolver
    ngspice_interface
    smbjson
    ${JSONFORTRAN_LIB}
)

add_library(system_tests "system_tests.cpp")

target_link_libraries(system_tests
    semba-main
    system_test_fortran
    GTest::gtest
)