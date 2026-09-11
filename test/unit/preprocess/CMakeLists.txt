message(STATUS "Creating build system for test/unit/preprocess")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library (preprocess_test_fortran 
    "test_preprocess_geom.F90"
)

target_link_libraries(preprocess_test_fortran
    semba-main
    semba-reports
    fhash
)

add_library(preprocess_tests "preprocess_tests.cpp")

target_link_libraries(preprocess_tests
    semba-main
    preprocess_test_fortran
    GTest::gtest
)
