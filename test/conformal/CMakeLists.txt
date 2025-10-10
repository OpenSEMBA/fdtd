message(STATUS "Creating build system for test/conformal")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library (conformal_test_fortran 
    "test_geometry.F90"
    "test_fhash.F90"
    "test_filling.F90"
    "test_pec_volume.F90"
)


target_link_libraries(conformal_test_fortran
    conformal 
    fhash 
)

add_library(conformal_tests "conformal_tests.cpp")

target_link_libraries(conformal_tests
    conformal
    conformal_test_fortran
    GTest::gtest
)
