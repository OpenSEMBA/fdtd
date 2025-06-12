message(STATUS "Creating build system for test/rotate")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library(
    rotate_test_fortran
    "rotate_testingTools.F90"
    "test_rotate_generateSpaceSteps.F90"
    "test_rotate_generateCurrentFieldSources.F90"
    "test_rotate_generatePlaneWaves.F90"
    "test_rotate_generateBoxSources.F90"
    "test_rotate_generateFronteras.F90"
    "test_rotate_generatePECs.F90"
    "test_rotate_generateNONMetals.F90"
    "test_rotate_generateAnisotropic.F90"
    "test_rotate_generateThinWires.F90"
    "test_rotate_generateSlantedWires.F90"
    "test_rotate_generateThinSlots.F90"
    "test_rotate_generateLossyThinSurface.F90"
    "test_rotate_generateFDMs.F90"
    "test_rotate_generateSONDAs.F90"
    "test_rotate_generateMasSondas.F90"
    "test_rotate_generateBloqueProbes.F90"
    "test_rotate_generateVolumicProbes.F90"
)

target_link_libraries(rotate_test_fortran
    fhash
    smbjson
    semba-main
    ${JSONFORTRAN_LIB}
)

add_library(rotate_tests "rotate_tests.cpp")

target_link_libraries(rotate_tests
    rotate_test_fortran
    GTest::gtest
)