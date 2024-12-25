
message(STATUS "Creating build system for test/system")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library(
    system_test_fortran
    "system_testingTools.F90"
    "test_map_wires_to_bundles.F90"
    "test_json_to_solver_input_shielded_pair.F90"
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
    system_test_fortran
    GTest::gtest
)