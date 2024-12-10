message(STATUS "Creating build system for test/smbjson")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library (smbjson_test_fortran 
    "smbjson_testingTools.F90" 
    "test_idchildtable.F90"
    "test_cells.F90"
    "test_mesh.F90"
    "test_parser.F90"    
    "test_read_planewave.F90"
    "test_read_dielectricslab.F90"
    "test_read_sgbc.F90"
    "test_read_holland1981.F90"
    "test_read_towelHanger.F90"
    "test_read_connectedWires.F90"
    "test_read_currentinjection.F90"
    "test_read_shieldedPair.F90"
    "test_read_mtln.F90"
    "test_read_sphere.F90"
    "test_read_airplane.F90"
    "test_read_large_airplane_mtln.F90"
)
target_link_libraries(smbjson_test_fortran
    smbjson 
    fhash 
    ${JSONFORTRAN_LIB}
)

add_library(smbjson_tests "smbjson_tests.cpp")

target_link_libraries(smbjson_tests
    smbjson 
    smbjson_test_fortran
    GTest::gtest
)
