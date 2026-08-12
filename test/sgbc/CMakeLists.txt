message(STATUS "Creating build system for test/sgbc")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library(sgbc_test_fortran
    "test_sgbc.F90"
)

target_link_libraries(sgbc_test_fortran
    semba-components
)

add_library(sgbc_tests "sgbc_tests.cpp")

target_link_libraries(sgbc_tests
    sgbc_test_fortran
    GTest::gtest
)
