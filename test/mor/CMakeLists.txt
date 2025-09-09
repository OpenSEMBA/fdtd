message(STATUS "Creating build system for test/mor")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library (mor_test_fortran 
    "test_evolution_operator.F90"
)

target_link_libraries(mor_test_fortran
    fhash 
    ${JSONFORTRAN_LIB}
    semba-main
)

add_library(mor_tests "mor_tests.cpp")

target_link_libraries(mor_tests
    mor_test_fortran
    GTest::gtest
)
