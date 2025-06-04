message(STATUS "Creating build system for test/rotate")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library(
    rotate_test_fortran
    "rotate_testingTools.F90"
    "test_rotate.F90"
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