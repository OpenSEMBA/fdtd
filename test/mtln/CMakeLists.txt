message(STATUS "Creating build system for test/mtln")

add_library(
    mtln_solver_tests_fortran
    "mtln_testingTools.F90"
    "test_mtln_types.F90"
    "test_preprocess.F90"
    "test_spice.F90"
    "test_math.F90"
    "test_mtl.F90"
    "test_mtl_bundle.F90"
    # "test_Paul.F90"
    "test_fhash.F90"
    "test_dispersive.F90"
    "test_multipolar_expansion.F90"
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
# The codemodels test loads the code models built from the ngspice submodule.
configure_file(ngspice_codemodels_dir.inc.in
    "${CMAKE_CURRENT_BINARY_DIR}/ngspice_codemodels_dir.inc" @ONLY)
target_include_directories(mtln_solver_tests_fortran PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
if (TARGET ngspice-codemodels)
    add_dependencies(mtln_solver_tests_fortran ngspice-codemodels)
endif()
