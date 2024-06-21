cmake_minimum_required(VERSION 3.5)
enable_language (CXX Fortran)
enable_testing ()

message(STATUS "Creating build system for test/mtln")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library(
    mtln_solver_tests_fortran
    "mtln_testingTools.F90"
    "test_mtln_types.F90"
    "test_preprocess.F90"
    "test_spice.F90"
    "test_math.F90"
    "test_mtl.F90"
    "test_mtl_bundle.F90"
    "test_Paul.F90"
    "test_fhash.F90"
    "test_dispersive.F90"
    "test_wire_panel.F90"
)
target_link_libraries(mtln_solver_tests_fortran
    mtlnsolver 
    fhash
    ngspice_interface
    ${NGSPICE_LIB}
)

add_library(mtln_tests "mtln_tests.cpp")

target_link_libraries(mtln_tests
    mtln_solver_tests_fortran
    GTest::gtest
)