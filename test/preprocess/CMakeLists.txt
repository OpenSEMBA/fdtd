message(STATUS "Creating build system for test/preprocess")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library(preprocess_test_fortran
    "test_searchtag.F90"
    "test_check_dielectric_tags.F90"
    "test_check_animated_tags.F90"
    "test_check_lossy_tags.F90"
)

target_link_libraries(preprocess_test_fortran
    semba-main
)

add_library(preprocess_tests "preprocess_tests.cpp")

target_link_libraries(preprocess_tests
    preprocess_test_fortran
    GTest::gtest
)
